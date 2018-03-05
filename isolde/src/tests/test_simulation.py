from time import time
class TestSimulation:
    def __init__(self, structure):
        self.structure = structure
        self.session = structure.session
        self.atoms = structure.atoms
        initialize_openmm()

        self._pause = True
        self._stop = False
        self._sim_running = False

        self._topology = None
        self._system = None
        self._platform = None
        self._simulation = None
        self._sim_thread_handler = None

        self._sim_steps_per_gui_update = 50
        from simtk import unit
        self._temperature = 100*unit.kelvin
        from simtk import openmm
        self._integrator_type = openmm.VariableLangevinIntegrator
        self._integrator_tol = 1e-4
        self._constraint_tol = 1e-4
        self._friction = 5.0/unit.picoseconds
        self._platform_name = 'OpenCL'
        self._create_openmm_system()
        self._loop_start_time = 0;
        self._last_loop_period = 0;

    @property
    def sim_rate(self):
        return 1/self._last_loop_period

    def prepare_sim(self):
        integrator = self._integrator_type(self._temperature, self._friction,
            self._integrator_tol)
        integrator.setConstraintTolerance(self._constraint_tol)

        from simtk.openmm import app
        s = self._simulation = app.Simulation(self._topology, self._system,
            integrator, self._platform)
        c = s.context
        c.setPositions(0.1*self.atoms.coords)
        c.setVelocitiesToTemperature(self._temperature)
        from ..openmm import openmm_interface
        self._sim_thread_handler = openmm_interface.OpenMM_Thread_Handler(s.context)

    def start_sim(self):
        if self._sim_running:
            raise RuntimeError('Simulation is already running!')
        if self._simulation is None:
            self.prepare_sim()
        self._pause = False
        self._stop = False
        self._sim_running = True
        self._minimize_and_go()

    def _minimize_and_go(self):
        th = self._sim_thread_handler
        from ..delayed_reaction import delayed_reaction
        delayed_reaction(self.session.triggers, 'new frame', th.minimize, [], th.thread_finished,
            self._update_coordinates_and_repeat, [True])

    def _repeat_step(self):
        self._last_loop_period = time()-self._loop_start_time
        self._loop_start_time = time()
        th = self._sim_thread_handler
        from ..delayed_reaction import delayed_reaction
        if th.unstable():
            f = th.minimize
            f_args = []
        else:
            f = th.step
            f_args = (self._sim_steps_per_gui_update,)
        delayed_reaction(self.session.triggers, 'new frame', f, f_args,
            th.thread_finished, self._update_coordinates_and_repeat, [])

    def _update_coordinates_and_repeat(self, reinit_vels = False):
        th = self._sim_thread_handler
        self.atoms.coords = th.coords
        if reinit_vels:
            th.reinitialize_velocities()
        if self._stop:
            self._sim_thread_handler.delete()
            self._sim_thread_handler = None
            self._simulation = None
            self._sim_running = False
            return
        if not self._pause:
            self._repeat_step()

    @property
    def thread_handler(self):
        return self._sim_thread_handler

    @property
    def pause(self):
        return self._pause

    @pause.setter
    def pause(self, flag):
        if not self._sim_running:
            raise RuntimeError('No simulation running!')
        if flag != self._pause:
            self._pause = flag
            if not flag:
                self._repeat_step()

    def stop(self):
        self._stop = True



    def _create_openmm_system(self):
        from simtk.openmm import app
        from simtk import openmm as mm
        from simtk import unit

        self._topology, self._particle_positions = openmm_topology_and_coordinates(self.structure)

        forcefield = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')

        system = forcefield.createSystem(self._topology,
                                        nonbondedMethod = app.CutoffNonPeriodic,
                                        nonbondedCutoff=1.0*unit.nanometers,
                                        constraints=app.HBonds,
                                        rigidWater=True)
        self._system = system
        platform = self._platform = mm.Platform.getPlatformByName(self._platform_name)

def run_test(session):
    import os
    base_dir = os.path.dirname(__file__)
    from chimerax.core.commands import open
    m = open.open(session, os.path.join(base_dir,'1pmx_1.pdb'))
    test_sim = TestSimulation(m[0])
    test_sim.prepare_sim()
    test_sim.start_sim()
    return test_sim


def openmm_topology_and_coordinates(mol):
    '''Make OpenMM topology and positions from ChimeraX AtomicStructure.'''
    a = mol.atoms
    n = len(a)
    r = a.residues
    aname = a.names
    ename = a.element_names
    rname = r.names
    rnum = r.numbers
    cids = r.chain_ids
    from simtk.openmm.app import Topology, Element
    top = Topology()
    cmap = {}
    rmap = {}
    atoms = {}
    for i in range(n):
        cid = cids[i]
        if not cid in cmap:
            cmap[cid] = top.addChain()	# OpenMM chains have no name.
        rid = (rname[i], rnum[i], cid)
        if not rid in rmap:
            rmap[rid] = top.addResidue(rname[i], cmap[cid])
        element = Element.getBySymbol(ename[i])
        atoms[i] = top.addAtom(aname[i], element, rmap[rid])
    a1, a2 = mol.bonds.atoms
    for i1, i2 in zip(a.indices(a1), a.indices(a2)):
        top.addBond(atoms[i1], atoms[i2])
    from simtk.openmm import Vec3
    pos = a.coords
    return top, pos


_openmm_initialized = False
def initialize_openmm():
    # On linux need to set environment variable to find plugins.
    # Without this it gives an error saying there is no "CPU" platform.
    global _openmm_initialized
    if not _openmm_initialized:
        _openmm_initialized = True
        from sys import platform
        if platform == 'linux' or platform == 'darwin':
            from os import environ, path
            from chimerax import app_lib_dir
            environ['OPENMM_PLUGIN_DIR'] = path.join(app_lib_dir, 'plugins')
