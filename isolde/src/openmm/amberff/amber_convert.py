'''
Conversion of AMBER mol2+frcmod files to OpenMM XML, doing its best to weed out
any dud entries. Requires parmed.
'''

# Incorrectly-parameterised residues
_blacklist = set((
'HOH', 'DIS', 'MTO', 'DOD', # Water
'NH3', # Incorrectly protonated as NH4
'PI', '2HP', # Should be HPO4(2-), H2PO4(-), both PO4(3-)
'SOH', # Should be HSO4-, modelled as SO42-
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


_leap_files = {
'FF14SB':   'leaprc.ff14sb.redq',
# 'WATER':    'leaprc.water.tip3p',
'GAFF2':    'leaprc.gaff2',
# 'DNA':      'leaprc.DNA.bsc1',
# 'RNA':      'leaprc.RNA.OL3',
'MODRNA08': 'leaprc.modrna08',
'PHOSAA10': 'leaprc.phosaa10',
'GLYCAM06': 'leaprc.GLYCAM_06j-1'
}

def make_combined_forcefield(dirname, output_xml):
    '''
    Combine all parameters and residue topologies into a single file.
    '''
    import os
    import parmed as pmd
    from glob import glob

    # BASE_DIR = os.path.abspath(__file__)
    #
    # _moriarty_lib_dir = os.path.join(BASE_DIR, 'src', 'moriarty_lib')
    #
    # AMBERHOME = os.environ['AMBERHOME']
    #
    #
    # LEAP_DIR = os.path.join(AMBERHOME, 'dat', 'leap', 'cmd')


    cwd = os.path.abspath(os.curdir)
    os.chdir(dirname)

    file_dict = {}

    mol2files = glob('**/*.mol2', recursive = True)
    for m in mol2files:
        name = os.path.splitext(os.path.basename(m))[0].upper()
        if name in _blacklist or name in _obsolete:
            continue
        name = 'CCD_'+name
        f = os.path.splitext(m)[0]+'.frcmod'
        if os.path.isfile(f):
            file_dict[name] = (m,f)

    rdict = {}

    fails = []
    for key, (m2f, _) in file_dict.items():
        try:
            mol2 = pmd.load_file(m2f)
            mol2.name = key
            rdict[key] = mol2
        except:
            fails.append(key)

    for name in fails:
        del file_dict[name]


    ff = pmd.openmm.OpenMMParameterSet.from_parameterset(
        pmd.amber.AmberParameterSet([fp[1] for fp in file_dict.values()])
    )
    ff.residues = rdict

    ff.write(output_xml)
    os.chdir(cwd)

if __name__ == "__main__":
    import sys
    make_combined_forcefield(*sys.argv[1:])
