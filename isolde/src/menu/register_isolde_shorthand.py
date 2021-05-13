
tooltip = ('Initialise and/or list a series of shorthand ISOLDE commands')

def run_script(session):
    from chimerax.core.commands import run
    run(session, 'isolde shorthand')