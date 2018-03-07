import os

cwd = os.path.dirname(os.path.abspath(__file__))
forcefields = {
    'amber14':  [os.path.join(cwd, 'amberff', f) for f in
        ['amberff14SB.xml','tip3p_standard.xml',
            'tip3p_HFE_multivalent.xml', 'tip3p_IOD_multivalent.xml']],
}
