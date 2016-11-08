# vim: set expandtab shiftwidth=4 softtabstop=4:

def get_singleton(session, create=True):
    if not session.ui.is_gui:
        return None
    from chimerax.core import tools
    from .tool import ISOLDE_ToolUI
    return tools.get_singleton(session, ISOLDE_ToolUI, 'ISOLDE', create=create)
    
def isolde(session, sim = False, atoms = None, map = None, mask_radius = 3.0, range = 5, softbuffer = 5.0, 
            hardbuffer = 8.0, T = 100.0, k = 1.0):
    ''' Run an interactive MD simulation
    
    Parameters
    ----------
    sim : Bool
        Start an interactive simulation
    atoms : Atomspec
        Core atom selection to base the simulation around. Note: this will
        be automatically extended to whole residues
    map : Model
        Map to couple the simulation to (optional)
    mask_radius : Float (default: 3.0)
        Distance (in Angstroms) from mobile atoms to mask the map down to
    range : Integer (default: 5)
        Number of residues before and after the selection to include in
        the simulation. Note that very small simulations surrounded by
        fixed atoms tend to be quite unstable, so in most cases it's best
        to have at least a little padding.
    softbuffer : Float (default: 5.0)
        Residues with atoms approaching within this distance of those
        defined by the "atoms" and "range" terms will be automatically
        included in the simulation. It is recommended to provide at least
        the default 5 Angstroms buffer.
    hardbuffer : Float (default: 8.0)
        Residues with atoms approaching within this distance of those
        defined by the "atoms", "range" and "softbuffer" terms will be
        included in the simulation, but rigidly fixed in space. Note that
        if this term is reduced to below the length of a bond, then the
        simulation will entirely ignore surrounding atoms (even those with
        direct bonds to mobile atoms). This is undesirable in almost all
        cases.
    T : Float (default: 100.0)
        Temperature of the simulation in Kelvin
    k : Float (default: 1.0)
        Arbitrary scaling constant defining how strongly the map pulls on
        atoms. It is recommended to experiment with this to find a value that
        guides the structure into place without introducing severe distortions.
        
    '''            
                
    log = session.logger
    if not session.ui.is_gui:
        log.warning("Sorry, ISOLDE currently requires ChimeraX to be in GUI mode")
        return
    if sim:
        if atoms is None:
            log.warning("You need to define an atom selection in order to start a simulation!")
            return
    
    ISOLDE = get_singleton(session)
    iobj = ISOLDE.isolde
    
    # Reset ISOLDE to defaults.
    #iobj.__init__(ISOLDE)
    
    iobj.b_and_a_padding = range
    iobj.soft_shell_cutoff = softbuffer
    iobj.hard_shell_cutoff = hardbuffer
    iobj.simulation_temperature = T
    
    
    if sim:
        if iobj._simulation_running:
            log.warning("You already have a simulation running!")
            return
        iobj.set_sim_selection_mode('from_picked_atoms')
        us = atoms.unique_structures
        if len(us) != 1:
            e = "Selection text must define atoms from exactly one model."
            raise Exception(e)
        sel_model = us[0]
        sel_model.selected = False
        atoms.selected = True
        
        if map is not None:
            iobj.set_sim_mode('em')
            iobj.master_map_list = {}
            iobj.add_map('map0', map, mask_radius, k, mask=True)
        else:
            iobj.set_sim_mode('free')
        
        iobj.start_sim()
        
        
            
    
                 
    
    
def register_isolde():
    from chimerax.core.commands import register, CmdDesc, AtomsArg, FloatArg, ModelArg, IntArg, BoolArg, NoArg
    desc = CmdDesc(
        required  = [],
        optional  = [('sim', NoArg),
                    ('atoms', AtomsArg),
                    ('map', ModelArg),
                    ('mask_radius', FloatArg),
                    ('range', IntArg),
                    ('softbuffer', FloatArg),
                    ('hardbuffer', FloatArg),
                    ('T', FloatArg),
                    ('k', FloatArg)]
                )
    register('isolde', desc, isolde)

_fps_tracker = None
