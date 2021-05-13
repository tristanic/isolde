
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
    ('pf', 'isolde pepflip sel'),
    ('cf', 'isolde cisflip sel'),
    ('cbb', 'color bfactor $*'),
    ('cbo', 'color byattr occupancy $*'),
    ('cbc', 'color bychain; color byhet'),
    ('cs', 'clipper set contourSensitivity $*')

)

def register_isolde_shorthand_commands(session):
    from chimerax.core.commands import create_alias
    log_string = ('<Pre>Initialising ISOLDE-specific command aliases:\n'
                  'Alias\tEquivalent full command\n'
                  '-------------------------------------------------\n'
    )
    for (alias, command) in aliases:
        create_alias(alias, command, logger=session.logger)
        log_string += f'{alias}\t{command}\n'.replace('$*', f'{{arguments}}')
    log_string += '</Pre>'
    session.logger.info(log_string, is_html=True)