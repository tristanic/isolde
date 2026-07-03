# @Author: Tristan Croll
# @Date:   02-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 02-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

from .param_provider import LigandBackedParameterisationProvider

_OPENFF_FORCEFIELD = 'openff-2.2.0.offxml'


def _cx_residue_to_offmol(residue):
    '''
    Convert a ChimeraX residue to an OpenFF Molecule.

    Builds an RDKit RWMol from atom elements and 3D coordinates, uses
    rdDetermineBonds to infer bond orders from geometry, then converts to
    an OpenFF Molecule via ``from_rdkit``.

    Requires the residue to have explicit hydrogen atoms for reliable bond-order
    inference.  ``charge=0`` is assumed; charged ligands may fail and fall back
    to AMBER14.
    '''
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDetermineBonds
    from openff.toolkit import Molecule as OFFMolecule

    cx_atoms = list(residue.atoms)
    atom_index = {a: i for i, a in enumerate(cx_atoms)}
    n = len(cx_atoms)

    mol = Chem.RWMol()
    conf = Chem.Conformer(n)
    for i, a in enumerate(cx_atoms):
        mol.AddAtom(Chem.Atom(a.element.number))
        coord = a.coord  # Angstrom
        conf.SetAtomPosition(i, (float(coord[0]), float(coord[1]), float(coord[2])))
    mol.AddConformer(conf, assignId=True)

    bonds_obj = residue.atoms.intra_bonds
    a1_col, a2_col = bonds_obj.atoms
    for a1, a2 in zip(a1_col, a2_col):
        mol.AddBond(atom_index[a1], atom_index[a2], Chem.BondType.SINGLE)

    rdDetermineBonds.DetermineBondOrders(mol, charge=0)
    mol_ro = mol.GetMol()
    AllChem.AssignStereochemistryFrom3D(mol_ro)

    return OFFMolecule.from_rdkit(mol_ro, allow_undefined_stereo=True)


def _assign_partial_charges(off_mol):
    '''
    Try AM1-BCC (needs ambertools or openeye), fall back to Gasteiger (RDKit).
    Gasteiger is lower quality but available with a pure pip install.
    '''
    for method in ('am1bcc', 'gasteiger'):
        try:
            off_mol.assign_partial_charges(method)
            return
        except Exception:
            continue
    raise RuntimeError(
        'No partial charge method available — install ambertools or '
        'openeye-toolkits for AM1-BCC charges, or openff-nagl for a '
        'neural-network alternative.')


def _openff_system_to_ffxml(residue, system) -> str:
    '''
    Convert an OpenFF-produced OpenMM System to a ForceField XML template
    string suitable for
    ``ForceField.loadFile(io.StringIO(xml), resname_prefix='USER_')``.

    OpenFF (like Espaloma) stores charges, sigma, and epsilon all in
    NonbondedForce — no CustomNonbondedForce.
    '''
    import xml.etree.ElementTree as ET
    from openmm import unit as omm_unit
    from openmm.openmm import (HarmonicBondForce, HarmonicAngleForce,
                               PeriodicTorsionForce, NonbondedForce)

    cx_atoms = list(residue.atoms)
    n_atoms = system.getNumParticles()
    rname = residue.name
    type_names = [f'OPN_{rname}_{i}' for i in range(n_atoms)]

    charges = [0.0] * n_atoms
    sigma   = [0.35] * n_atoms
    epsilon = [0.40] * n_atoms

    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            for ni in range(force.getNumParticles()):
                q, sig, eps = force.getParticleParameters(ni)
                charges[ni] = q.value_in_unit(omm_unit.elementary_charge)
                sigma[ni]   = sig.value_in_unit(omm_unit.nanometer)
                epsilon[ni] = eps.value_in_unit(omm_unit.kilojoule_per_mole)
            break

    root = ET.Element('ForceField')

    at_el = ET.SubElement(root, 'AtomTypes')
    for i, cx_a in enumerate(cx_atoms):
        ET.SubElement(at_el, 'Type', {
            'name':    type_names[i],
            'class':   type_names[i],
            'element': cx_a.element.name.capitalize(),
            'mass':    str(cx_a.element.mass),
        })

    res_el = ET.SubElement(ET.SubElement(root, 'Residues'), 'Residue', name=rname)
    for i, cx_a in enumerate(cx_atoms):
        ET.SubElement(res_el, 'Atom',
            name=cx_a.name,
            type=type_names[i],
            charge=str(charges[i]))
    bonds_obj = residue.atoms.intra_bonds
    a1_col, a2_col = bonds_obj.atoms
    for a1_name, a2_name in zip(a1_col.names, a2_col.names):
        ET.SubElement(res_el, 'Bond', atomName1=a1_name, atomName2=a2_name)

    nb_el = ET.SubElement(root, 'NonbondedForce',
        coulomb14scale='0.8333333333333333', lj14scale='0.5')
    ET.SubElement(nb_el, 'UseAttributeFromResidue', name='charge')
    for i in range(n_atoms):
        ET.SubElement(nb_el, 'Atom',
            type=type_names[i],
            sigma=str(sigma[i]),
            epsilon=str(epsilon[i]))

    for force in system.getForces():
        if isinstance(force, HarmonicBondForce):
            bf_el = ET.SubElement(root, 'HarmonicBondForce')
            for bi in range(force.getNumBonds()):
                i1, i2, r0, k = force.getBondParameters(bi)
                ET.SubElement(bf_el, 'Bond',
                    type1=type_names[i1], type2=type_names[i2],
                    length=str(r0.value_in_unit(omm_unit.nanometer)),
                    k=str(k.value_in_unit(
                        omm_unit.kilojoule_per_mole / omm_unit.nanometer**2)))

        elif isinstance(force, HarmonicAngleForce):
            af_el = ET.SubElement(root, 'HarmonicAngleForce')
            for ai in range(force.getNumAngles()):
                i1, i2, i3, theta0, k = force.getAngleParameters(ai)
                ET.SubElement(af_el, 'Angle',
                    type1=type_names[i1], type2=type_names[i2],
                    type3=type_names[i3],
                    angle=str(theta0.value_in_unit(omm_unit.radian)),
                    k=str(k.value_in_unit(
                        omm_unit.kilojoule_per_mole / omm_unit.radian**2)))

        elif isinstance(force, PeriodicTorsionForce):
            tf_el = ET.SubElement(root, 'PeriodicTorsionForce')
            for ti in range(force.getNumTorsions()):
                i1, i2, i3, i4, per, phase, k = force.getTorsionParameters(ti)
                ET.SubElement(tf_el, 'Proper',
                    type1=type_names[i1], type2=type_names[i2],
                    type3=type_names[i3], type4=type_names[i4],
                    periodicity1=str(per),
                    phase1=str(phase.value_in_unit(omm_unit.radian)),
                    k1=str(k.value_in_unit(omm_unit.kilojoule_per_mole)))

    return ET.tostring(root, encoding='unicode')


def _openff_parameterise_residue(residue, logger) -> str:
    '''Orchestrate the full OpenFF parameterisation pipeline for one residue.'''
    from openff.toolkit import ForceField

    off_mol = _cx_residue_to_offmol(residue)
    _assign_partial_charges(off_mol)
    ff = ForceField(_OPENFF_FORCEFIELD)
    off_top = off_mol.to_topology()
    system = ff.create_openmm_system(off_top,
                                     charge_from_molecules=[off_mol],
                                     constraints=None)
    return _openff_system_to_ffxml(residue, system)


class OpenFFParameterisationProvider(LigandBackedParameterisationProvider):
    '''
    Parameterisation provider using AMBER14 for standard residues and
    OpenFF Sage (openff-2.2.0) for ligands.

    Partial charges are assigned via AM1-BCC if ambertools or openeye-toolkits
    are available, falling back to Gasteiger (RDKit) for pure pip installs.
    '''

    @property
    def backend_forcefield_key(self):
        return 'amber14+openff'

    def _parameterise_one(self, residue, logger):
        xml = _openff_parameterise_residue(residue, logger)
        return xml, {}
