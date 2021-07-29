
aliases = (
    ('st', 'isolde step $*'),
    ('aw', 'isolde add water $*'),
    ('awsf', 'isolde add water $* sim false'),
    ('al', 'isolde add ligand $*'),
    ('aa', 'isolde add aa $1 sel $*'),
    ('ht', 'isolde mod his sel $*'),
    ('so', 'setattr sel atoms occupancy $*'),
    ('ab', 'isolde adjust bfactors $*'),
    ('ss', 'isolde sim start sel'),
    ('rt', 'isolde release torsions sel $*'),
    ('rd', 'isolde release distances sel $*'),
    ('ra', 'rd; rt'),
    ('pf', 'isolde pepflip sel'),
    ('cf', 'isolde cisflip sel'),
    ('cbb', 'color bfactor $*'),
    ('cbo', 'color byattr occupancy $*'),
    ('cbc', 'color bychain; color byhet'),
    ('cs', 'clipper set contourSensitivity $*')

)

global _shorthand_initialized
_shorthand_initialized = False

def register_isolde_shorthand_commands(session):
    global _shorthand_initialized
    from chimerax.core.commands import create_alias
    if not _shorthand_initialized:
        log_string = ('<Pre>Initialising ISOLDE-specific command aliases:\n')
    else:
        log_string = ('<Pre>Registered ISOLDE-specific command aliases:\n')
    log_string += ('Alias\tEquivalent full command\n'
                   '-------------------------------------------------\n'
    )
    for (alias, command) in aliases:
        if not _shorthand_initialized:
            create_alias(alias, command, logger=session.logger)
        log_string += f'{alias}\t{command}\n'.replace('$*', f'{{arguments}}')
    log_string += '</Pre>'
    session.logger.info(log_string, is_html=True)
    _shorthand_initialized = True