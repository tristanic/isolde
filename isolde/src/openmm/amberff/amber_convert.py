# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 10-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll



'''
Conversion of AMBER mol2+frcmod files to OpenMM XML, doing its best to weed out
any dud entries. Requires parmed.
'''

# Incorrectly-parameterised residues
_blacklist = set((
'HOH', 'DIS', 'MTO', 'DOD', # Water
# 'NH3', # Incorrectly protonated as NH4
'PI', '2HP', # Should be HPO4(2-), H2PO4(-), both PO4(3-)
'SOH', # Should be HSO4-, modelled as SO42-
'ACE', # N-terminal acetylation - not a free ligand
))

_obsolete = set((
"CL2", "PS0", "FS4", "NMO", "MEP", "BFP", "2LU", "GTT", "NA6", "PHC", "HSZ",
"RIC", "XL1", "MH3", "MPS", "OME", "BA4", "GS4", "PHM", "ILG", "PLY", "OC1",
"OCM", "N2B", "HAP", "PAS", "1NI", "CYP", "THL", "CRC", "4RR", "CD3", "MT7",
"LIN", "GSA", "PYH", "T", "IP4", "DHO", "PNL", "577", "BDN", "OC6", "7MQ",
"148", "CS0", "QEH", "T36", "MDG", "LT4", "MO6", "ATW", "NAC", "ACL", "HQ",
"MW2", "OTB", "BGG", "1CU", "GLB", "ASQ", "DCG", "CPP", "PVL", "MO3", "TFH",
"IDT", "0AT", "T37", "SPY", "CD5", "MGY", "INT", "SBU", "ETH", "EOX", "BTA",
"NAN", "FPC", "MEO", "PLE", "395", "GTE", "BOP", "MBA", "OCO", "FLH", "NC",
"NH", "2EP", "FIB", "5AX", "543", "GMY", "0EZ", "CLO", "LAU", "MTY", "CPG",
"MFA", "0SN", "HMP", "6MC", "CRW", "PR5", "PGC", "ZN2", "PPL", "ATI", "4EX",
"FAT", "RAW", "PKR", "CM", "HP3", "KO4", "NED", "IU6", "BRO", "NAH", "PSS",
"9ZT", "HYD", "6CT", "9HE", "TTH", "NA5", "6JZ", "ASK", "13H", "DML", "F22",
"KGR", "LA4", "M4C", "AA3", "OC4", "ABK", "LCH", "0AU", "MGA", "TBM", "FDC",
"NIV", "GCL", "HDP", "MAM", "OCN", "1NA", "TIB", "EYS", "HGC", "QB4", "QND",
"CRY", "OXO", "R51", "MB0", "IOX", "TRB", "2OF", "B1P", "BUT", "I5C", "NOZ",
"MP3", "DGH", "NFB", "NEV", "SEO", "FLO", "G32", "DG0", "PA3", "OLE", "CAY",
"DSP", "GTN", "MM4", "SAP", "ZN3", "CO5", "0SP", "AFL", "SI2", "PTP", "TY3",
"FRF", "TZO", "191", "VAF", "OIP", "CB2", "RT", "2PL", "7BO", "MCL", "MIP",
"TS", "BTC", "TOS", "G33", "BAC", "RDP", "CD1", "CEC", "SEU", "MJQ", "PR0",
"ZH3", "INY", "QTR", "CSW", "IP3", "NYC", "PY2", "6MT", "ISB", "AYG", "H",
"EHN", "OF2", "PGQ", "SOH", "EHP", "GCM", "Q72", "DDB", "GLW", "O4M", "OBS",
"GSD", "DHN", "T31", "TMN", "CL1", "BOX", "OC8", "CMN", "HSR", "SEA", "CGL",
"OWL", "S", "FA", "ADG", "IT1", "BRM", "CEA", "CBM", "GLQ", "OCL", "S4U", "CSE",
"LG", "TNP", "GUR", "GSW", "COE", "GTO", "GLI", "NI3", "A1P", "DIS", "LC",
"EGL", "T0M", "OF3", "TML", "PDL", "LLA", "HII", "ITS", "OTE", "OC3", "MCB",
"NAM", "DFG", "ASI", "CYO", "LAA", "OLI", "4SR", "PCC", "SPG", "KOL", "WRS",
"CH2", "BH4", "0ZO", "HAA", "ATG", "CQR", "TGU", "OMB", "AMA", "TSO", "DXN",
"PO0", "DGC", "BP", "B51", "DHU", "HHP", "3IH", "DRT", "CHA", "NGL", "MAI",
"39N", "YAP", "RAA", "638", "D1P", "1ZD", "OC7", "DLA", "UBR", "SUL", "0A5",
"STY", "FRO", "LOL", "UIC", "NNS", "U18", "ACU", "OXA", "PF", "0AV", "EOM",
"ICI", "XPP", "B1F", "HNP", "DOH", "2AB", "CSY", "C25", "A34", "0AM", "STO",
"PA4", "GRP", "CYA", "A39", "1ZT", "HPB", "MHM", "TS3", "OET", "AGC", "GLM",
"A35", "ANE", "A9D", "DGM", "EGG", "OF1", "MAD", "SAA", "SFN", "MMR", "2SI",
"DMR", "LYM", "NCP", "PTC", "IOH", "TMB", "MNV", "BUG", "2PI", "PC2", "CBZ",
"VAS", "CLE", "HAC", "DOX", "GUD", "450", "DDM", "U2N", "3MU", "WRR", "NGN",
"119", "PCZ", "0AA", "BH1", "BPC", "MC1", "MO4", "MOL", "LCX", "LTR", "2OG",
"EKE", "DUM", "POC", "SCP", "THQ", "LA5", "PGY", "ETD", "CH3", "L0H", "UPR",
"BCC", "3MD", "ABG", "BAR", "PPM", "TPH", "CY7", "LYW", "OR5", "LP2", "ITM",
"XBO", "C32", "FS3", "TRG", "0AL", "PI1", "THB", "F2O", "KTH", "E4N", "NI2",
"NEM", "SRI", "MN6", "WSQ", "BRI", "PEI", "PTF", "3B7", "ALM", "YJC", "1PY",
"O2", "NWB", "3OF", "45H", "FMR", "P2K", "B4M", "NC2", "345", "RNY", "HBL",
"GLH", "MH1", "0P0", "IPS", "0AZ", "OX", "MO5", "MW3", "AC4", "OC5", "LEM",
"CN", "DNJ", "IDO", "CBX", "AKL", "NA2", "PBK", "SOM", "NAW", "NEW", "ZO3",
"3DB", "MN5", "NI1", "SOY", "2AO", "TF6", "CYL", "MPN", "ROB", "HPG", "OC2",
"0AC", "6AB", "CFS", "CYM", "G42", "HV5", "MO2", "MW1", "ETO", "DFC", "K00",
"NIK", "SCC", "BOT", "RON", "2OA", "ZNO", "PAG", "CUC", "HSU", "102", "PPP",
"5IT", "AB7", "OPH", "PIG", "FCY", "URY", "NTY", "LHU", "PG7", "53P", "2MB",
"NAO", "I42", "MO1", "YF1", "0UH", "LOF", "QX", "U25", "GDB", "267", "DPH",
"2SP", "3KV", "B7D", "5SA", "FU4", "MTO", "G25", "OHE", "DSX", "SEG", "BZO",
"0A3", "5HP", "MCE", "HPC", "HAD", "CNM", "YH", "AHA", "HEQ", "AB5", "LEP"
))

_iron_sulfur = set(("SF4", "CYF"))

_leap_files = {
'FF14SB':   'leaprc.ff14sb.redq',
# 'WATER':    'leaprc.water.tip3p',
'GAFF2':    'leaprc.gaff2',
# 'DNA':      'leaprc.DNA.bsc1',
# 'RNA':      'leaprc.RNA.OL3',
# 'MODRNA08': 'leaprc.modrna08',
'PHOSAA10': 'leaprc.phosaa10',
'GLYCAM06': 'leaprc.GLYCAM_06j-1'
}

def _find_mol2_frcmod_pairs(input_dir, blacklist=set(), search_subdirs=True):
    import os
    from glob import glob

    file_dict = {}
    if search_subdirs:
        mol2files = glob(os.path.join(input_dir, '**/*.mol2'), recursive=True)
    else:
        mol2files = glob(os.path.join(input_dir, '*.mol2'), recursive=False)
    for m in mol2files:
        name = os.path.splitext(os.path.basename(m))[0].upper()
        if name in blacklist:
            continue
        f = os.path.splitext(m)[0]+'.frcmod'
        if os.path.isfile(f):
            file_dict[name] = (m,f)

    return file_dict

def amber_to_ffxml(frcmod_file, mol2_file, output_name=None, atom_names_from=None):
    import os
    import parmed as pmd
    mol2 = pmd.load_file(mol2_file)
    ff = pmd.openmm.OpenMMParameterSet.from_parameterset(
        pmd.amber.AmberParameterSet(frcmod_file)
    )
    ff.residues[mol2.name] = mol2
    if atom_names_from is not None:
        from chimerax.core.errors import UserError
        from chimerax.atomic import Residue
        if not isinstance(atom_names_from, Residue):
            raise UserError('atom_names_from must be a single residue')
        from chimerax.isolde.graph import make_graph_from_residue, make_graph_from_parmed_residue
        rg = make_graph_from_residue(atom_names_from)
        pr = ff.residues[mol2.name]
        pg = make_graph_from_parmed_residue(pr)
        ri, pi, _ = rg.maximum_common_subgraph(pg)
        if len(ri) != len(atom_names_from.atoms):
            raise UserError('Failed to map all atoms in the residue to the mol2 file. This may be because the mol2 file contains multiple residues, or because of a mismatch between the residue and mol2 file. If your mol2 file contains multiple residues, try splitting it into individual files for each residue and converting them separately.')
        for a, p in zip(ri, pi):
            pr.atoms[p].name = atom_names_from.atoms[a].name

    if output_name is None:
        output_name = os.path.splitext(os.path.basename(mol2_file))[0]+'.xml'
    ff.write(output_name)
    return output_name




def amber_to_xml_individual(input_dir, output_dir, resname_prefix=None,
        dict_filename = None, compress_to = None, blacklist=False, search_subdirs=True):
    '''
    Convert all .mol2/.frcmod pairs found under input_dir into OpenMM XML files
    in output_dir. If resname_prefix is given, it will be prepended to all
    residue names in the xml files.
    '''
    import os
    import parmed as pmd
    from glob import glob

    if blacklist:
        blacklist = _blacklist.union(_obsolete)
    else:
        blacklist=set()

    file_dict = _find_mol2_frcmod_pairs(input_dir, blacklist=blacklist, search_subdirs=search_subdirs)

    resname_to_ff = {}
    xmlfiles = []
    fails = []
    with open(os.path.join(output_dir, 'failures.log'), 'wt') as logfile:

        for key, (m2f, frcmod) in file_dict.items():
            try:
                mol2 = pmd.load_file(m2f)
                if resname_prefix is not None:
                    mol2.name = resname_prefix+key
                else:
                    mol2.name = key
            except:
                print('Failed to open {}.'.format(m2f))
                logfile.write(m2f+'\n')
                fails.append(key)
                continue
            ff = pmd.openmm.OpenMMParameterSet.from_parameterset(
                pmd.amber.AmberParameterSet(frcmod)
            )
            resname_to_ff[key] = mol2.name
            ff.residues[mol2.name] = mol2
            xml_name = key+'.xml'
            xmlfiles.append(xml_name)
            ff.write(os.path.join(output_dir, xml_name))

    for f in fails:
        del file_dict[f]

    if dict_filename is not None:
        import json
        if os.path.splitext(dict_filename)[1] != 'json':
            dict_filename+='.json'
        with open(dict_filename, 'wt') as jfile:
            json.dump(resname_to_ff, jfile)

    if compress_to is not None:
        import zipfile
        if os.path.splitext(compress_to)[1].lower() != '.zip':
            compress_to += '.zip'
        with zipfile.ZipFile(compress_to, 'w', zipfile.ZIP_DEFLATED) as z:
            for xf in xmlfiles:
                z.write(xf, os.path.basename(xf))


def make_combined_forcefield(dirname, output_xml, resname_prefix=None):
    '''
    Combine all parameters and residue topologies into a single file.
    '''
    import os
    import parmed as pmd
    from glob import glob

    file_dict = _find_mol2_frcmod_pairs(dirname, blacklist=_blacklist.union(_obsolete))

    rdict = {}

    fails = []
    for key, (m2f, _) in file_dict.items():
        try:
            mol2 = pmd.load_file(m2f)
            if resname_prefix is not None:
                mol2.name = resname_prefix+key
            else:
                mol2.name = key
            rdict[mol2.name] = mol2
        except:
            fails.append(key)

    for name in fails:
        del file_dict[name]

    aff = pmd.amber.AmberParameterSet()
    for fp in file_dict.values():
        frcmod = fp[1]
        name = os.path.splitext(os.path.basename(frcmod))[0].upper()
        try:
            aff.load_parameters(frcmod)
        except pmd.exceptions.ParameterError:
            print('WARNING: Could not load parameters for {}. Removing this entry'
                ' from the residue dictionary'.format(name))
            del rdict[resname_prefix+name]
    ff = pmd.openmm.OpenMMParameterSet.from_parameterset(aff)


    # ff = pmd.openmm.OpenMMParameterSet.from_parameterset(
    #     pmd.amber.AmberParameterSet([fp[1] for fp in file_dict.values()])
    # )
    ff.residues = rdict
    ff.condense()

    ff.write(output_xml, skip_duplicates=False)

def convert_bryce_parms(dirname, output_xml, resname_prefix='BRYCE_'):
    '''
    Convert the AMBER parameterisations found at
    http://research.bmh.manchester.ac.uk/bryce/amber/
    '''
    files = {
    ('ATP', 'ADP', 'GTP', 'GDP'):   ('phos.lib', 'frcmod.phos'), # Nucleotide di/triphosphates
    ('FAD',):                       ('FADH-.lib', 'FADH-.frcfld'), # FAD-
    ('FMN',):                       ('fmn.off', 'fmn.frcfld'), # Flavin mononucleotide
    ('NAD', 'NAH', 'NP2', 'NP3',
        'NPD', 'NPH'):              ('nad_variants.lib', 'nad.frcmod'),
    ('H1D',):                       ('H1D.off', 'frcmod_h1d'), # Phosphohistidine (ND1, 1-)
    ('H1E',):                       ('H1E.off', 'frcmod_h1d'), # Phosphohistidine (NE2, 1-)
    ('H2D',):                       ('H2D.off', 'frcmod_h2d'), # Phosphohistidine (ND1, 2-)
    ('H2E',):                       ('H2E.off', 'frcmod_h2e'), # Phosphohistidine (NE2, 2-)
    ('HEM',):                       ('hem.off', 'frcmod.hemall'), # Heme
    }
    import parmed as pmd
    import os
    aff = pmd.amber.AmberParameterSet()
    from collections import OrderedDict
    residues = OrderedDict()
    for resnames, (template_file, param_file) in files.items():
        aff.load_parameters(os.path.join(dirname, param_file))
        inres = pmd.load_file(os.path.join(dirname, template_file))
        for r, res in inres.items():
            res.name = resname_prefix+res.name
            residues[res.name] = res
    off = pmd.openmm.OpenMMParameterSet.from_parameterset(aff)
    off.residues = residues
    off.condense()
    off.write(output_xml, skip_duplicates=False)

def find_duplicate_definitions(dirname):
    import parmed as pmd
    from glob import glob
    import os
    from collections import defaultdict
    known_angles = defaultdict(lambda: list())
    known_bonds = defaultdict(lambda: list())
    known_propers = defaultdict(lambda: list())
    known_impropers = defaultdict(lambda: list())
    import numpy

    file_dict = _find_mol2_frcmod_pairs(dirname, blacklist=_blacklist.union(_obsolete))
    for key, (m2f, frcmod) in file_dict.items():
        ff = pmd.openmm.OpenMMParameterSet.from_parameterset(
            pmd.amber.AmberParameterSet(frcmod)
        )
        for batoms, bdef in ff.bond_types.items():
            sba = frozenset(batoms)
            params = numpy.array((bdef.k, bdef.req))
            kbl = known_bonds[sba]
            if not len(kbl):
                kbl.append(params)
                continue
            for kb in kbl:
                if all(numpy.isclose(params, kb)):
                    break
            else:
                kbl.append(params)

        for aatoms, adef in ff.angle_types.items():
            saa = frozenset((aatoms, aatoms[::-1]))
            params = numpy.array((adef.k, adef.theteq))
            kal = known_angles[saa]
            if not len(kal):
                kal.append(params)
                continue
            for ka in kal:
                if all(numpy.isclose(params, ka)):
                    break
            else:
                kal.append(params)

        for datoms, ddeflist in ff.dihedral_types.items():
            sda = frozenset((datoms, datoms[::-1]))
            kpl = known_propers[sda]
            ddeflist = [numpy.array((ddef.phi_k, ddef.per, ddef.phase, ddef.scee, ddef.scnb)) for ddef in ddeflist]
            if not len(kpl):
                kpl.append(ddeflist)
                continue
            any_match = False
            for kpd in kpl:
                if len(kpd) == len(ddeflist):
                    match = True
                    for kp, ddef in zip(kpd, ddeflist):
                        if not all(numpy.isclose(kp, ddef)):
                            match = False
                            break
                    if match:
                        any_match = True
                        break
            if not any_match:
                kpl.append(ddeflist)


            # for ddef in ddeflist:
            #     params = numpy.array((ddef.phi_k, ddef.per, ddef.phase, ddef.scee, ddef.scnb))
            #     kpl = known_propers[sda]
            #     if not len(kpl):
            #         kpl.append(params)
            #         continue
            #     for kp in kpl:
            #         if all(numpy.isclose(params, kp)):
            #             break
            #     else:
            #         kpl.append(params)

        for iatoms, idef in ff.improper_periodic_types.items():
            params = numpy.array((idef.phi_k, idef.per, idef.phase, idef.scee, idef.scnb))
            kil = known_impropers[iatoms]
            if not len(kil):
                kil.append(params)
                continue
            for ki in kil:
                if all(numpy.isclose(params, ki)):
                    break
            else:
                kil.append(params)
    return known_bonds, known_angles, known_propers, known_impropers



#: AMBER -> OpenMM unit conversions (frcmod uses kcal/mol, Angstrom, degrees;
#: OpenMM ffXML uses kJ/mol, nm, radians). Harmonic force constants also gain a
#: factor of 2 because AMBER's energy is k(x-x0)^2 while OpenMM's is (1/2)k(x-x0)^2.
_KCAL = 4.184                       # kcal/mol -> kJ/mol
_BOND_K = _KCAL * 100.0 * 2.0       # kcal/mol/A^2 -> kJ/mol/nm^2, x2
_ANGLE_K = _KCAL * 2.0              # kcal/mol/rad^2 -> kJ/mol/rad^2, x2


def _frcmod_openmm_params(frcmod_file):
    '''Parse an AMBER frcmod directly into OpenMM-unit parameter tables (kJ/mol,
    nm, rad), keyed by GAFF atom-type tuples. Parsing the (small, fixed-column)
    text ourselves avoids a ParmEd/NumPy dependency that has proven fragile.

    Returns ``{'bonds': {(t1,t2): (length, k)},
               'angles': {(t1,t2,t3): (angle, k)},
               'propers': [ (t1,t2,t3,t4, [(per, phase, k), ...]) ],
               'impropers': [ (t1,t2,t3,t4, [(per, phase, k), ...]) ] }``.

    frcmod atom types occupy 2 columns each, dash-separated (``c3-ss``,
    ``c -ss``), so they are read by fixed offsets. Multi-term propers (signalled by
    a negative periodicity) are accumulated into one entry.
    '''
    import math
    out = {'bonds': {}, 'angles': {}, 'propers': [], 'impropers': []}
    _HEADERS = {'MASS', 'BOND', 'ANGLE', 'ANGL', 'DIHE', 'DIHEDRAL',
                'IMPROPER', 'IMPROP', 'NONBON', 'NONB', 'LJEDIT', 'CMAP'}
    _CANON = {'ANGL': 'ANGLE', 'DIHEDRAL': 'DIHE', 'IMPROP': 'IMPROPER',
              'NONB': 'NONBON'}
    section = None
    prev_proper_open = False        # last proper line had negative periodicity
    with open(frcmod_file) as f:
        lines = f.readlines()
    for raw in lines[1:]:           # first line is a free-text title
        line = raw.rstrip('\n')
        s = line.strip()
        up = s.upper()
        if up in _HEADERS:
            section = _CANON.get(up, up)
            prev_proper_open = False
            continue
        if not s:
            section = None if section in ('MASS', 'NONBON') else section
            prev_proper_open = False
            continue
        try:
            if section == 'BOND':
                t1, t2 = line[0:2].strip(), line[3:5].strip()
                p = line[5:].split()
                out['bonds'][(t1, t2)] = (float(p[1]) * 0.1, float(p[0]) * _BOND_K)
            elif section == 'ANGLE':
                t1, t2, t3 = line[0:2].strip(), line[3:5].strip(), line[6:8].strip()
                p = line[8:].split()
                out['angles'][(t1, t2, t3)] = (math.radians(float(p[1])),
                                               float(p[0]) * _ANGLE_K)
            elif section == 'DIHE':
                t = (line[0:2].strip(), line[3:5].strip(),
                     line[6:8].strip(), line[9:11].strip())
                p = line[11:].split()
                idivf, pk, phase, pn = (float(p[0]), float(p[1]),
                                        float(p[2]), float(p[3]))
                term = (abs(int(round(pn))), math.radians(phase),
                        pk / idivf * _KCAL)
                if prev_proper_open and out['propers'] and out['propers'][-1][:4] == t:
                    out['propers'][-1][4].append(term)
                else:
                    out['propers'].append((t[0], t[1], t[2], t[3], [term]))
                prev_proper_open = pn < 0
            elif section == 'IMPROPER':
                t = (line[0:2].strip(), line[3:5].strip(),
                     line[6:8].strip(), line[9:11].strip())
                p = line[11:].split()
                pk, phase, pn = float(p[0]), float(p[1]), float(p[2])
                out['impropers'].append(
                    (t[0], t[1], t[2], t[3],
                     [(abs(int(round(pn))), math.radians(phase), pk * _KCAL)]))
        except (IndexError, ValueError):
            continue                # tolerate stray/short lines
    return out


def _lookup_bond(params, g1, g2):
    return params['bonds'].get((g1, g2)) or params['bonds'].get((g2, g1))


def _lookup_angle(params, g1, g2, g3):
    return params['angles'].get((g1, g2, g3)) or params['angles'].get((g3, g2, g1))


def _lookup_proper(params, g1, g2, g3, g4):
    '''Match a proper by exact GAFF quad, its reverse, or with wildcard ('' /'X')
    end atoms -- the frcmod's usual generic form.'''
    def _match(pt):
        a, b, c, d = pt[0], pt[1], pt[2], pt[3]
        def wc(x):
            return x in ('', 'X')
        for (w, x, y, z) in ((g1, g2, g3, g4), (g4, g3, g2, g1)):
            if (wc(a) or a == w) and b == x and c == y and (wc(d) or d == z):
                return True
        return False
    for pt in params['propers']:
        if _match(pt):
            return pt[4]
    return None


def _lookup_improper(params, central, peripherals):
    '''Match a frcmod improper to a centre of GAFF type ``central`` with the three
    ``peripherals`` (GAFF types). frcmod impropers put the central atom THIRD
    (a1-a2-a3-a4, central=a3), with wildcards ('' / 'X') common on a1/a2. Returns
    ``(entry_tuple, terms)`` (entry as stored, central still 3rd) or ``None``.'''
    peri = list(peripherals)
    for entry in params['impropers']:
        a1, a2, a3, a4, terms = entry
        if a3 not in ('', 'X') and a3 != central:
            continue
        remaining = list(peri)
        ok = True
        for s in (a1, a2, a4):
            if s in ('', 'X'):
                continue
            if s in remaining:
                remaining.remove(s)
            else:
                ok = False
                break
        wilds = sum(1 for s in (a1, a2, a4) if s in ('', 'X'))
        if ok and wilds >= len(remaining):
            return (a1, a2, a3, a4), terms
    return None


def _parse_gaff2(gaff2_xml):
    '''Read GAFF2 per-type element/mass (from ``<AtomTypes>``) and LJ
    sigma/epsilon (from ``<NonbondedForce>``), so the emitted self-contained ffXML
    can define its prefixed types fully. Returns
    ``{gaff_type: {'element','mass','sigma','epsilon'}}``.'''
    import xml.etree.ElementTree as ET
    info = {}
    root = ET.parse(gaff2_xml).getroot()
    at = root.find('AtomTypes')
    if at is not None:
        for t in at.findall('Type'):
            info[t.get('name')] = {'element': t.get('element'), 'mass': t.get('mass'),
                                   'sigma': None, 'epsilon': None}
    nb = root.find('NonbondedForce')
    if nb is not None:
        for a in nb.findall('Atom'):
            ty = a.get('type') or a.get('class')
            if ty in info:
                info[ty]['sigma'] = a.get('sigma')
                info[ty]['epsilon'] = a.get('epsilon')
    return info


def covalent_to_ffxml(per_res, template_names, mol, frcmod_file, output_path,
                      gaff2_xml, gaff_equivalent=None, lookup_seam=None):
    '''Emit a SELF-CONTAINED OpenMM ffXML for a covalent unit.

    Every atom that carries a GAFF type is given a UNIQUE, per-unit **prefixed**
    atom type (e.g. ``MCOV_UNL_ce``), and the file defines those types in full --
    ``<AtomTypes>``, ``<NonbondedForce>`` LJ, and EVERY bonded parameter touching
    them -- so it stands alone and cannot collide with another ligand parameterised
    in the same session (the openmmforcefields trick; OpenMM has no atom-type-name
    length limit). Atoms of the frozen ff14SB core keep their standard ff14SB types
    (shared with the surrounding protein), so the external peptide bonds stay
    parameterised by the base force field.

    A bonded term is emitted iff it touches at least one prefixed (GAFF) atom;
    terms wholly among standard ff14SB atoms are left to the base force field.
    Values come from the ``-a Y`` frcmod (the molecule's complete GAFF parameter
    set), looked up by the atoms' actual GAFF types and re-keyed onto the emitted
    (prefixed / ff14SB) types. Bonds + angles are complete (so ``createSystem``
    succeeds); propers and impropers are included where the frcmod defines them.
    Impropers are emitted OpenMM-style with the central atom FIRST (the frcmod
    stores it third), matching the shipped amberff14SB.xml convention.

    Args:
        gaff2_xml: path to gaff2.xml (source of prefixed types' element/mass/LJ).
        frcmod_file: -a Y frcmod from run_antechamber; None emits templates only.
        lookup_seam: reserved for curated seam overrides (currently unused with the
            prefixed-type scheme).
    '''
    import xml.etree.ElementTree as ET
    if gaff_equivalent is None:
        from .boundary_params import gaff_equivalent as gaff_equivalent

    gaff_info = _parse_gaff2(gaff2_xml)
    params = _frcmod_openmm_params(frcmod_file) if frcmod_file else \
        {'bonds': {}, 'angles': {}, 'propers': [], 'impropers': []}

    # Decide each atom's EMITTED type: a unique prefixed type for GAFF atoms
    # (ligand / retyped shell), the standard ff14SB type otherwise. An atom is
    # GAFF-typed exactly when its final type equals the type ANTECHAMBER assigned.
    rd2spec = {}
    for r, d in per_res.items():
        pfx = template_names[r] + '_'
        for spec in d['atoms']:
            spec['is_gaff'] = (spec['type'] == spec['gaff_type'])
            spec['etype'] = (pfx + spec['gaff_type']) if spec['is_gaff'] else spec['type']
            rd2spec[spec['rd']] = spec

    prefixed = {}       # emitted type -> source GAFF type
    for d in per_res.values():
        for spec in d['atoms']:
            if spec['is_gaff']:
                prefixed[spec['etype']] = spec['gaff_type']

    root = ET.Element('ForceField')

    if prefixed:
        at_el = ET.SubElement(root, 'AtomTypes')
        for etype, g in sorted(prefixed.items()):
            gi = gaff_info.get(g, {})
            ET.SubElement(at_el, 'Type', {'name': etype, 'class': etype,
                'element': gi.get('element') or 'C', 'mass': gi.get('mass') or '0.0'})

    residues_el = ET.SubElement(root, 'Residues')
    for r, d in per_res.items():
        rel = ET.SubElement(residues_el, 'Residue', {'name': template_names[r]})
        for spec in d['atoms']:
            ET.SubElement(rel, 'Atom', {'name': spec['name'], 'type': spec['etype'],
                                        'charge': '%.5f' % spec['charge']})
        for (n1, n2) in d['bonds']:
            ET.SubElement(rel, 'Bond', {'atomName1': n1, 'atomName2': n2})
        for nm in sorted(d['external']):
            ET.SubElement(rel, 'ExternalBond', {'atomName': nm})

    if prefixed:
        # Charge is supplied per-atom in the templates via the base force field's
        # UseAttributeFromResidue; here we only add the LJ for the prefixed types.
        nb_el = ET.SubElement(root, 'NonbondedForce',
            {'coulomb14scale': '0.8333333333333334', 'lj14scale': '0.5'})
        # Charge is taken per-atom from the residue templates (like gaff2.xml);
        # without this, OpenMM's standalone parser rejects charge-less <Atom>s.
        ET.SubElement(nb_el, 'UseAttributeFromResidue', {'name': 'charge'})
        for etype, g in sorted(prefixed.items()):
            gi = gaff_info.get(g, {})
            if gi.get('sigma') is not None:
                ET.SubElement(nb_el, 'Atom', {'type': etype, 'sigma': gi['sigma'],
                                              'epsilon': gi['epsilon']})

    def _emit(specs):
        # Emit a term iff it touches a prefixed (GAFF) atom; all-ff14SB terms are
        # covered by the base force field.
        return any(s['is_gaff'] for s in specs)

    def _gt(spec):
        return spec.get('gaff_type') or gaff_equivalent(spec['type'])

    # A bond or angle whose force constant came back zero (or absent) means
    # parmchk2/GAFF could not parameterise it (an "ATTN, need revision" placeholder)
    # -- leaving that internal coordinate unrestrained. We fail loudly rather than
    # emit junk. (Zero-barrier *torsions* are legitimate GAFF and are NOT checked.)
    # Skipped when there is no frcmod (template-only emission, e.g. in tests).
    check_zero = frcmod_file is not None
    unparam = set()

    # --- bonds ---
    bond_out = {}
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        si, sj = rd2spec.get(i), rd2spec.get(j)
        if si is None or sj is None or not _emit((si, sj)):
            continue
        key = frozenset((si['etype'], sj['etype']))
        if key in bond_out:
            continue
        p = _lookup_bond(params, _gt(si), _gt(sj))
        if check_zero and (p is None or p[1] == 0.0):
            unparam.add('bond %s-%s (GAFF %s-%s)'
                        % (si['name'], sj['name'], _gt(si), _gt(sj)))
        if p is not None:
            bond_out[key] = (si['etype'], sj['etype'], p[0], p[1])

    # --- angles (a-b-c around each central atom) ---
    angle_out = {}
    for b_atom in mol.GetAtoms():
        bi = b_atom.GetIdx()
        if rd2spec.get(bi) is None:
            continue
        nbrs = [n.GetIdx() for n in b_atom.GetNeighbors() if rd2spec.get(n.GetIdx())]
        for x in range(len(nbrs)):
            for y in range(x + 1, len(nbrs)):
                ai, ci = nbrs[x], nbrs[y]
                specs = (rd2spec[ai], rd2spec[bi], rd2spec[ci])
                if not _emit(specs):
                    continue
                e1, e2, e3 = (s['etype'] for s in specs)
                key = (e1, e2, e3) if e1 <= e3 else (e3, e2, e1)
                if key in angle_out:
                    continue
                p = _lookup_angle(params, _gt(specs[0]), _gt(specs[1]), _gt(specs[2]))
                if check_zero and (p is None or p[1] == 0.0):
                    unparam.add('angle %s-%s-%s (GAFF %s-%s-%s)'
                                % (specs[0]['name'], specs[1]['name'], specs[2]['name'],
                                   _gt(specs[0]), _gt(specs[1]), _gt(specs[2])))
                if p is not None:
                    angle_out[key] = (e1, e2, e3, p[0], p[1])

    if unparam:
        from chimerax.core.errors import UserError
        raise UserError(
            'GAFF/parmchk2 could not assign force constants for the following '
            'internal coordinate(s) of %s -- they came back zeroed, which would '
            'leave them unrestrained:\n  %s\n\n'
            'This chemistry is outside GAFF2\'s coverage. Parameterisation was '
            'aborted rather than produce a template with unrestrained bonds/angles; '
            'the affected atoms need parameters supplied by hand.'
            % (', '.join(sorted(set(template_names.values()))),
               '\n  '.join(sorted(unparam))))

    # --- propers ---
    proper_out = {}
    for b in mol.GetBonds():
        j, k = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if rd2spec.get(j) is None or rd2spec.get(k) is None:
            continue
        bj = mol.GetAtomWithIdx(j)
        bk = mol.GetAtomWithIdx(k)
        for na in bj.GetNeighbors():
            ai = na.GetIdx()
            if ai == k or rd2spec.get(ai) is None:
                continue
            for nd in bk.GetNeighbors():
                di = nd.GetIdx()
                if di == j or rd2spec.get(di) is None:
                    continue
                specs = (rd2spec[ai], rd2spec[j], rd2spec[k], rd2spec[di])
                if not _emit(specs):
                    continue
                e = tuple(s['etype'] for s in specs)
                key = min(e, e[::-1])
                if key in proper_out:
                    continue
                terms = _lookup_proper(params, *[_gt(s) for s in specs])
                if terms:
                    proper_out[key] = (e, terms)

    # --- impropers (centre with 3 neighbours; emitted OpenMM central-FIRST) ---
    improper_out = {}
    for c_atom in mol.GetAtoms():
        cs = rd2spec.get(c_atom.GetIdx())
        if cs is None or not cs['is_gaff']:
            continue
        nspecs = [rd2spec[n.GetIdx()] for n in c_atom.GetNeighbors()
                  if rd2spec.get(n.GetIdx()) is not None]
        if len(nspecs) != 3:
            continue
        res = _lookup_improper(params, _gt(cs), [_gt(n) for n in nspecs])
        if not res:
            continue
        (a1, a2, a3, a4), terms = res
        # Map the frcmod peripheral slots (a1,a2,a4) to neighbour etypes, keeping
        # wildcards as '' -- then promote the centre to first (OpenMM convention).
        remaining = list(nspecs)
        slots = []
        for sgaff in (a1, a2, a4):
            if sgaff in ('', 'X'):
                slots.append('')
                continue
            match = next((n for n in remaining if _gt(n) == sgaff), None)
            if match is not None:
                slots.append(match['etype'])
                remaining.remove(match)
            else:
                slots.append('')
        key = (cs['etype'], tuple(sorted(slots)))
        if key not in improper_out:
            improper_out[key] = ((cs['etype'], slots[0], slots[1], slots[2]), terms)

    if bond_out:
        hb = ET.SubElement(root, 'HarmonicBondForce')
        for (e1, e2, length, k) in bond_out.values():
            ET.SubElement(hb, 'Bond', {'type1': e1, 'type2': e2,
                                       'length': '%.6f' % length, 'k': '%.4f' % k})
    if angle_out:
        ha = ET.SubElement(root, 'HarmonicAngleForce')
        for (e1, e2, e3, ang, k) in angle_out.values():
            ET.SubElement(ha, 'Angle', {'type1': e1, 'type2': e2, 'type3': e3,
                                        'angle': '%.6f' % ang, 'k': '%.4f' % k})
    if proper_out or improper_out:
        pt = ET.SubElement(root, 'PeriodicTorsionForce')
        for (e, terms) in proper_out.values():
            attrs = {'type1': e[0], 'type2': e[1], 'type3': e[2], 'type4': e[3]}
            for n, (per, phase, k) in enumerate(terms, start=1):
                attrs['periodicity%d' % n] = str(per)
                attrs['phase%d' % n] = '%.6f' % phase
                attrs['k%d' % n] = '%.4f' % k
            ET.SubElement(pt, 'Proper', attrs)
        for (e, terms) in improper_out.values():
            attrs = {'type1': e[0], 'type2': e[1], 'type3': e[2], 'type4': e[3]}
            for n, (per, phase, k) in enumerate(terms, start=1):
                attrs['periodicity%d' % n] = str(per)
                attrs['phase%d' % n] = '%.6f' % phase
                attrs['k%d' % n] = '%.4f' % k
            ET.SubElement(pt, 'Improper', attrs)

    try:
        ET.indent(root, space='  ')   # human-readable multi-line layout (Py>=3.9)
    except AttributeError:
        pass
    ET.ElementTree(root).write(output_path, encoding='unicode', xml_declaration=False)
    return output_path


if __name__ == "__main__":
    import sys
    amber_to_ffxml(sys.argv[1], sys.argv[2])

# if __name__ == "__main__":
#     import sys
#     make_combined_forcefield(*sys.argv[1:])
