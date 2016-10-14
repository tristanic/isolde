
class available_forcefields():
    # Main force field files
    main_files = [
        'amber99sbildn.xml',
        'amber99sbnmr.xml',
        'amber10.xml'
        ]
    
    main_file_descriptions = [
        'AMBER99 with improved backbone & sidechain torsions',
        'AMBER99 with modifications to fit NMR data',
        'AMBER10'
        ]
    
    # Implicit solvent force field files. The main file and the implicit
    # solvent file must match.
    implicit_solvent_files = [
        'amber99_obc.xml',
        'amber99_obc.xml',
        'amber10_obc.xml'
        ]
    
    # Explicit water models
    explicit_water_files = [
        'tip3pfb.xml',
        'tip4pfb.xml',
        'tip3p.xml'
        ]
    
    explicit_water_descriptions = [
        'TIP3P-FB (DOI: 10.1021/jz500737m)',
        'TIP4P-FB (DOI: 10.1021/jz500737m)',
        'Original TIP3P water (not recommended)'
        ]
    
    
    
