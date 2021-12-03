# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tic20@cam.ac.uk
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

def amber_to_ffxml(frcmod_file, mol2_file, output_name=None):
    import os
    import parmed as pmd
    mol2 = pmd.load_file(mol2_file)
    ff = pmd.openmm.OpenMMParameterSet.from_parameterset(
        pmd.amber.AmberParameterSet(frcmod_file)
    )
    ff.residues[mol2.name] = mol2
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



if __name__ == "__main__":
    import sys
    amber_to_ffxml(sys.argv[1], sys.argv[2])

# if __name__ == "__main__":
#     import sys
#     make_combined_forcefield(*sys.argv[1:])
