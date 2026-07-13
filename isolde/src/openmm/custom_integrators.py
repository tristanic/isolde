from openmm.openmm import CustomIntegrator

class VelocityChecker(CustomIntegrator):
    '''
    Performs an in-line sanity check on atom velocities, useful to catch "explosions" before they get too far.
    Does not affect the simulation.
    '''
    MAX_SPEED = 50.0 # nm ps-1 (50 km/s)
    def __init__(self):
        super().__init__(0)
        self.addGlobalVariable('fast_count', 0.0)
        # The (1-step(-m)) factor is a mass-based mask: it is 1 for real
        # (mass>0) particles and 0 for any zero-mass particle. This makes it
        # impossible *by construction* to count a virtual site (e.g. a
        # crystallographic symmetry copy) or a fixed atom as "fast" - important
        # because a virtual site slaved to a fast-moving parent inherits its
        # speed and would otherwise trip the instability check spuriously. It
        # keys on the intrinsic mass rather than any particle-ordering
        # assumption, so it remains correct for any future virtual particles.
        self.addComputeSum('fast_count',
            f'(1-step(-m))*step(sqrt(_x(v)^2+_y(v)^2+_z(v)^2)/({self.MAX_SPEED}*3)-1)')

class Smoother(CustomIntegrator):
    '''
    Maintains an exponential moving average of atomic coordinates. Does not affect the simulation.
    '''
    def __init__(self, smoothing_alpha, enabled=True):
        super().__init__(0)
        self.addGlobalVariable('reset_smooth', 1.0)
        self.addGlobalVariable('enabled', int(enabled))
        self.addGlobalVariable('smoothing_alpha', smoothing_alpha)
        self.addPerDofVariable('smoothed', 0)
    
        self.beginIfBlock('enabled > 0')

        self.beginIfBlock('reset_smooth>0')
        self.addComputePerDof('smoothed', 'x')
        self.addComputeGlobal('reset_smooth', '0')
        self.endBlock() # reset_smooth>0

        self.beginIfBlock('reset_smooth=0')
        self.addComputePerDof('smoothed', 'x*smoothing_alpha + smoothed*(1-smoothing_alpha)')
        self.endBlock() # reset_smooth=0

        self.endBlock() # enabled > 0
    
        self.beginIfBlock('enabled = 0')
        self.addComputePerDof('smoothed', 'x')
        self.endBlock()
    
    @property
    def enabled(self):
        return bool(self.getGlobalVariableByName('enabled'))
    
    @enabled.setter
    def enabled(self, flag):
        self.setGlobalVariableByName('enabled', int(flag))
    
    @property
    def smoothing_alpha(self):
        return self.getGlobalVariableByName('smoothing_alpha')
    
    @smoothing_alpha.setter
    def smoothing_alpha(self, val):
        self.setGlobalVariableByName('smoothing_alpha', val)

