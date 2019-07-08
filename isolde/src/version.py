
def version(session):
    installed_bundle_info = session.toolshed.bundle_info(session.logger)
    for b in installed_bundle_info:
        if b.name == 'ChimeraX-ISOLDE':
            return b.version
 
