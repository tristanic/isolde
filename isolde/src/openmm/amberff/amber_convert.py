'''
Conversion of AMBER mol2+frcmod files to OpenMM XML, doing its best to weed out
any dud entries. Requires parmed.
'''

# Incorrectly-parameterised residues
_blacklist = (
'HOH', 'DIS', 'MTO', 'DOD', # Water
'NH3', # Incorrectly protonated as NH4
'PI', '2HP', # Should be HPO4(2-), H2PO4(-), both PO4(3-)
'SOH', # Should be HSO4-, modelled as SO42-
'DTU', 'DTV', # enantiomers of DTT that somehow get missed by ParmED
)

_obsolete = (
# Annotated as obsolete in components.cif
'0A3', '0A5', '0AZ', '0EZ', '0P0', '0SN', '0UH', '0ZO', '102',
'191', '2AB', '2AO', '2EP', '2MB', '39N', '3B7', '3IH', '3MU', '45H', '4EX',
'5AX', '6JZ', '7BO', '9HE', '9ZT', 'A1P', 'A9D', 'AB5', 'ABG', 'ACL', 'AKL',
'ALM', 'ASK', 'ATG', 'ATI', 'ATW', 'AYG', 'B4M', 'B51', 'B7D', 'BA4', 'BAC',
'BCC', 'BGG', 'BH1', 'BOP', 'BUG', 'CB2', 'CFS', 'CGL', 'CH2', 'CHA', 'CL1',
'CLE', 'CPP', 'CRW', 'CSE', 'CSW', 'CSY', 'CY7', 'DGC', 'DMR', 'DPH', 'EGG',
'EHN', 'EKE', 'EOM', 'EYS', 'FDC', 'FMR', 'FPC', 'FRF', 'FRO', 'FU4', 'G25',
'GLH', 'GLM', 'GLQ', 'GMY', 'GRP', 'GSW', 'GTT', 'GUD', 'HEQ', 'HII', 'I42',
'ILG', 'IP3', 'IP4', 'IT1', 'ITM', 'K00', 'KOL', 'KTH', 'L0H', 'L0H', 'LA4',
'LA5', 'LAA', 'LEM', 'LEP', 'LOL', 'LT4', 'LYM', 'LYW', 'M4C', 'MAD', 'MAI',
'MAI', 'MC1', 'MCL', 'MJQ', 'MM4', 'MMR', 'MNV', 'MPN', 'MT7', 'NAM', 'NYC',
'OPH', 'OWL', 'PA3', 'PBK', 'PCZ', 'PHC', 'PHM', 'PI1', 'PKR', 'PLE', 'PO0',
'PPL', 'PPP', 'PR0', 'PS0', 'PSS', 'PTC', 'PTF', 'PY2', 'PYH', 'QB4', 'QND',
'RNY', 'SAP', 'SCC', 'SOH', 'SOY', 'SPY', 'SRI', 'SRI', 'T0M', 'TBM', 'TF6',
'TGU', 'THL', 'TNP', 'TOS', 'TZO', 'U2N', 'UPR', 'URY', 'VAF', 'VAS', 'WRR',
'WRS', 'WSQ', 'XBO', 'YAP', 'YJC',
# Annotated as obsolete on ligand expo
'267', 'AC4', 'IPS', 'ANE', 'BDN', 'CM', 'CBM', 'CYL', 'EGL', 'CBX', 'GUR',
'FIB', 'ICI', 'KGR', 'LAU', 'BUT', 'OLI', '119', 'HPG', 'PGC', 'PGQ', 'RAW',
'SCP', 'NAN', 'SPG', 'SUL', 'TFH', 'TMN', 'SOM', 'CRY',
)

def make_combined_forcefield(dirname, output_xml):
    import os
    import parmed as pmd
    from glob import glob

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
