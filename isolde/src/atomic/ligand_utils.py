def recluster_ligands(session, model):
    from chimerax.atomic import Residue, Residues, concatenate
    m = model
    polymer = m.residues[m.residues.polymer_types!=Residue.PT_NONE]
    non_polymer = m.residues[m.residues.polymer_types==Residue.PT_NONE]
    # Ligands bound to polymeric residues should keep their existing chain IDs
    unbound, bound_map = sort_bound_and_unbound_ligands(non_polymer)
    chain_map, unclassified = cluster_unbound_ligands(m, unbound)
    non_polymer.chain_ids = 'Xtmp'
    m.renumber_residues(non_polymer, 0)
    assign_bound_ligand_chain_ids_and_residue_numbers(m, bound_map)
    assign_unbound_residues(m, chain_map, unclassified)
    
def assign_unbound_residues(m, chain_map, unclassified):
    from chimerax.atomic import Residue, concatenate
    for cid, residues in chain_map.items():
        first_ligand_number = next_available_ligand_number(m, cid)
        residues.chain_ids = cid+'tmp'
        residues = concatenate([residues[residues.names!='HOH'], residues[residues.names=='HOH']])
        m.renumber_residues(residues, first_ligand_number)
        residues.chain_ids=cid
    if len(unclassified):
        # Try reclassifying with a more permissive cutoff, now that all the other ligands 
        # are assigned to chains
        chain_map, unclassified = cluster_unbound_ligands(m, unclassified, cutoff=8)
        assign_unbound_residues(m, chain_map, [])
        if len(unclassified):
            session = m.session
            warn_str = ('{} residues could not be automatically assigned to chains. '
                'These have been given the chain ID UNK.').format(len(unclassified))
            session.logger.warning(warn_str)
            unclassified.chain_ids='UNK'
            m.renumber_residues(unclassified, 1)
    
def find_duplicate_residues(m):
    seen=set()
    duplicates=[]
    for r in reversed(m.residues):
        sig = (r.chain_id, r.number)
        if sig in seen:
            duplicates.append(r)
        seen.add(sig)
    return duplicates



_metal_weight=5
_carbon_weight=0.5

def cluster_unbound_ligands(model, unbound, cutoff=5):
    from chimerax.geometry import find_close_points
    from chimerax.atomic import Residue, Residues
    from collections import defaultdict
    import numpy
    m = model
    chain_ids = m.residues.unique_chain_ids
    other_residues = m.residues.subtract(unbound)
    #polymeric = m.residues[m.residues.polymer_types!=Residue.PT_NONE]
    ligand_atoms = unbound.atoms[unbound.atoms.element_names != 'H']
    chain_map = {}
    for cid in chain_ids:
        cres = other_residues[other_residues.chain_ids==cid]
        catoms = cres.atoms[cres.atoms.element_names!='H']
        ci, li = find_close_points(catoms.coords, ligand_atoms.coords, cutoff)
        close_ligand_atoms = ligand_atoms[li]
        weights = numpy.ones(len(close_ligand_atoms), numpy.double)
        weights[close_ligand_atoms.element_names == 'C'] = _carbon_weight
        weights[close_ligand_atoms.elements.is_metal] = _metal_weight
        chain_map[cid] = Weighted_Counter([a.residue for a in close_ligand_atoms], weights)
    unclassified = []
    closest_chain_map = defaultdict(list)
    for r in unbound:
        max_atoms = 0
        closest = None
        for cid in chain_ids:
            close = chain_map[cid].get(r, None)
            if close is not None:
                if close > max_atoms:
                    closest = cid
                    max_atoms = close
        if closest is not None:
            closest_chain_map[closest].append(r)
        else:
            unclassified.append(r)
    return {cid: Residues(residues) for cid, residues in closest_chain_map.items()}, Residues(unclassified)

def create_merged_residue(residues, base_residue=None, new_name=None):
    if base_residue is None:
        base_residue = residues[0]
    if base_residue not in residues:
        raise TypeError('base_residue must be a member of residues!')
    if new_name is None:
        new_name = base_residue.name
    from chimerax.atomic import Residues
    from chimerax.atomic.struct_edit import add_atom, add_bond
    residues = Residues(residues)
    seen = set()
    for r in residues:
        for n in r.neighbors:
            seen.add(n)
    if not all([r in seen for r in residues]):
        raise TypeError('All residues must be connected through covalent bonds!')
    atoms = residues.atoms
    other_residues = residues.subtract(Residues([base_residue]))
    bonds = atoms.intra_bonds
    external_bonds = []
    for r in seen:
        if r not in other_residues:
            for n in r.neighbors:
                if n in residues:
                    external_bonds.extend(r.bonds_between(n))
    next_atom_number = highest_atom_number(base_residue.atoms)+1

    atom_map = {}
    other_atoms = other_residues.atoms
    # Do heavy atoms first, so we can assign hydrogen names based on them
    other_heavy_atoms = other_atoms[other_atoms.element_names!='H']
    for a in other_heavy_atoms:
        name = a.element.name + str(next_atom_number)
        next_atom_number += 1
        atom_map[a] = add_atom(name, a.element, base_residue, a.coord, occupancy=a.occupancy, bfactor=a.bfactor)
    # Add the hydrogens
    for a, na in atom_map.items():
        base_number = na.name[1:]
        i=1
        for n in a.neighbors:
            if n.element.name=='H':
                name = 'H'+base_number+str(i)
                i += 1
                add_atom(name, n.element, base_residue, n.coord, bonded_to=na, occupancy=n.occupancy, bfactor=n.bfactor)
    # Add all internal bonds
    for a, na in atom_map.items():
        for n in a.neighbors:
            nn = atom_map.get(n, None)
            if nn is not None and nn not in na.neighbors:
                add_bond(na, nn)
    # Now we just need to add the bonds between base_residue and the merged portion, and any 
    # other bonds between the merged portion and the wider model. To avoid valence errors we'll need
    # to first delete the existing bonds.
    for b in external_bonds:
        atoms = b.atoms
        if atoms[1] in atom_map.keys():
            atoms = atoms[::-1]
        na = atom_map[atoms[0]]
        a = atoms[1]
        b.delete()
        add_bond(na, a)
    other_residues.delete()
    base_residue.name=new_name
        




def highest_atom_number(atoms):
    import re
    atoms = atoms[atoms.element_names!='H']
    highest=0
    for a in atoms:
        name = a.name
        number = re.findall(r'\d+', name)
        if len(number):
            num = int(number[0])
            if num > highest:
                highest=num
    return highest





    


class Weighted_Counter(dict):
    def __init__(self, elements, weights):
        super().__init__(self)
        for e, w in zip(elements, weights):
            v = self.get(e, 0)
            v += w
            self[e] = v                

def next_available_ligand_number(m, chain_id):
    from chimerax.atomic import Residue
    cres = m.residues[m.residues.chain_ids==chain_id]
    lres = cres[cres.polymer_types==Residue.PT_NONE]
    if len(lres):
        return max(lres.numbers)+1
    pres = cres[cres.polymer_types!=Residue.PT_NONE]
    if len(pres):
        max_polymeric_residue_number = max(pres.numbers)
    else:
        return 0
    first_ligand_number = round(max_polymeric_residue_number+1000,-3)
    if first_ligand_number - max_polymeric_residue_number < 100:
        first_ligand_number += 1000
    return first_ligand_number

known_sugars = ["NGC","NGA","RM4","FCB","GLC","GCS","GTR","GAL","BDR","RIB","BDF","BGC","BXX","XYZ","FUL","FUB","Z6H","PSV","A2G","LDY","RIP","SHD","NDG","ARB","SOE","SLB","BM3","LFR","Z6J","GIV","PA1","ABE","AHR","XXR","ARA","AFD","ADA","IDR","MAN","RAM","32O","NAG","GUP","T6T","G6D","FUC","Z9N","BMA","XYP","FCA","BDP","TYV","BM7","GCU","GXL","XYS","HSY","LXC","FRU","WOO","ALL","0MK","SIA","BXY","64K","GL0","GLA"]
def assign_bound_ligand_chain_ids_and_residue_numbers(m, bound_map):
    # The wwPDB steering committee has dictated that for protein-linked glycans,
    # the following rules apply:
    #   - if the modelled portion of the glycan is a single residue, it should have
    #     the same chain ID as the protein.
    #   - if more than one residue, it should have a unique chain ID.
    from chimerax.atomic import Residues, Residue, concatenate
    import numpy
    for cid, bound_residues in bound_map.items():
        first_ligand_number = next_available_ligand_number(m, cid)
        new_glycan_cid = []
        assign_to_chain = []
        groups = independent_groupings(bound_residues)
        for g in groups:
            if len(g)==1:
                assign_to_chain.append(g)
            elif any (r.name in known_sugars for r in g):
                new_glycan_cid.append(g)
            else:
                assign_to_chain.append(g)
        new_glycan_cid = list(sorted(new_glycan_cid, key=lambda g:linked_polymer_residue(g).number))
        assign_to_chain = list(sorted(assign_to_chain, key=lambda g:linked_polymer_residue(g).number))
        for i,g in enumerate(new_glycan_cid):
            gid = cid+'gl'+str(i)
            residues = sort_glycan_residues(g)
            residues.chain_ids=gid
            m.renumber_residues(residues, 1)
        if len(assign_to_chain):
            assign_to_chain = concatenate([Residues(g) for g in assign_to_chain])
            assign_to_chain.chain_ids = 'XXtmp'
            m.renumber_residues(assign_to_chain, first_ligand_number)
            assign_to_chain.chain_ids = cid


def sort_glycan_residues(residues):
    from chimerax.atomic import Residues
    polymer_stem_res = linked_polymer_residue(residues)
    for r in polymer_stem_res.neighbors:
        if r in residues:
            break
    order = [r]
    def _sort_walk(r, residues, order):
        bonds = [r.bonds_between(n)[0] for n in r.neighbors if n in residues and n not in order]
        ordered_bonds = []
        for b in bonds:
            atoms = b.atoms
            if atoms[0].residue == r:
                ordered_bonds.append(atoms)
            else:
                ordered_bonds.append(atoms[::-1])
        ordered_bonds = list(sorted(ordered_bonds, key=lambda b:
            (int(b[0].name[-1]))))
        for b in ordered_bonds:
            next_r = b[1].residue
            order.append(next_r)
            _sort_walk(next_r, residues, order)
    _sort_walk(r, residues, order)
    return Residues(order)



        
def linked_polymer_residue(residue_group):
    from chimerax.atomic import Residue
    for r in residue_group:
        for n in r.neighbors:
            if n.polymer_type != Residue.PT_NONE:
                return n
    return None

def independent_groupings(residues):
    residues = list(residues)
    groups = []
    while len(residues) > 0:
        r = residues[0]
        connected = set([r])
        walk_tree(r, residues, connected)
        groups.append(connected)
        for r in connected:
            residues.remove(r)

    return groups        

def walk_tree(r, residues, connected):
    for n in r.neighbors:
        if n in connected or n not in residues:
            continue
        connected.add(n)
        walk_tree(n, residues, connected)

def sort_bound_and_unbound_ligands(residues):
    unbound = []
    from collections import defaultdict
    bound_map = defaultdict(set)
    for r in residues:
        if not len(r.neighbors):
            unbound.append(r)
        cid = bound_to_polymer(r)
        if cid is not None:
            bound_map[cid].add(r)
        else:
            unbound.append(r)
    from chimerax.atomic import Residues
    return (Residues(unbound), {cid: Residues(rset) for cid, rset in bound_map.items()})

def bound_to_polymer(residue, seen=None):
    from chimerax.atomic import Residue
    if seen is None:
        seen=set()
    seen.add(residue)
    for n in residue.neighbors:
        if n in seen:
            continue
        if n.polymer_type!=Residue.PT_NONE:
            return n.chain_id
        seen.add(n)
        return bound_to_polymer(n, seen)
    return None


