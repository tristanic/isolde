'''
Utilities for shifting chains in register.
'''
import numpy
from scipy import interpolate
from .util import is_continuous_protein_chain, find_polymer
from chimerax.core.atomic import Atoms

class ProteinRegisterShifter:
    '''
    Shifts a protein chain in register by dragging CA and CB (where 
    present) atom targets along splines. 
    '''
    def __init__(self, session, isolde, atoms):
        if not is_continuous_protein_chain(atoms):
            raise TypeError(
                'Selection must be atoms from a continuous peptide chain!')
        
        self.spring_constant = 100 # kJ/mol/A2
        
        # Number of GUI update steps between each residue along the spline
        self.spline_steps_per_residue = 10
        self._spline_step = 1/self.spline_steps_per_residue
        # Current value of the spline parameter
        self._current_position_on_spline = 0
        
        self.finished = False
        
        # Trigger handler to update position along spline
        self._handler = None
         
        self.session = session
        self.isolde = isolde
        self.polymer = find_polymer(atoms)
        residues = self.residues = atoms.unique_residues
        self._extended_atoms = None
        
        nres = len(residues)
        # We need to build up the array of atoms procedurally here, since
        # we want to know where there's a gap in the CB positions.
        atoms = self._make_atom_arrays(residues)
        coords = []
        n_atoms = self._n_atoms = Atoms(atoms[:,0])
        ncoords = n_atoms.coords
        
        ca_atoms = self._ca_atoms = Atoms(atoms[:,1])
        cacoords = ca_atoms.coords
        
        c_atoms = self._c_atoms = Atoms(atoms[:,2])
        ccoords = c_atoms.coords
        
        nocb = numpy.equal(atoms[:,3], None)
        nocb_indices = numpy.argwhere(nocb).ravel()
        cb_indices = self._cb_indices = numpy.argwhere(numpy.invert(nocb)).ravel()
        cbcoords = numpy.empty([nres, 3])
        cb_atoms = self._cb_atoms = Atoms(atoms[cb_indices,3])
        cbcoords[cb_indices] = cb_atoms.coords
        
        # Fill in the missing CB positions. If CB is missing we'll assume
        # we're dealing with a glycine and fudge it by using the HA3 
        # position
        
        glyres = residues[nocb_indices]
        glyha3 = glyres.atoms.filter(glyres.atoms.names == 'HA3')
        cbcoords[nocb_indices] = glyha3.coords
        
        u_vals = numpy.arange(0, nres)
        
        # Prepare the splines
        nspl = self.n_spline = interpolate.splprep(ncoords.transpose(), u=u_vals)
        caspl = self._ca_spline = interpolate.splprep(cacoords.transpose(), u=u_vals)
        cspl = self._c_spline = interpolate.splprep(ccoords.transpose(), u=u_vals)
        cbspl = self._cb_spline = interpolate.splprep(cbcoords.transpose(), u=u_vals)

    def _make_atom_arrays(self, residues):
        nres = len(residues)
        atoms = numpy.empty([nres, 4], object)
        for i, r in enumerate(residues):
            ratoms = r.atoms
            atoms[i][0] = ratoms.filter(ratoms.names=='N')[0]
            atoms[i][1] = ratoms.filter(ratoms.names=='CA')[0]
            atoms[i][2] = ratoms.filter(ratoms.names=='C')[0]
            try:
                atoms[i][3] = ratoms.filter(ratoms.names=='CB')[0]
            except IndexError:
                atoms[i][3] = None
        return atoms

        
    def shift_register(self, nres):
        '''
        Shift the atoms in register by the given number of residues. A 
        positive nres yields a shift towards the C-terminus, while a 
        negative nres gives a shift towards the N-terminus. Residues that
        fall off the end of the target region will become unrestrained,
        while residues that enter the start will become restrained.
        '''
        isolde = self.isolde
        p = self.polymer
        residues = self.residues
        polymer_length = len(p)
        self._spline_length = len(residues)
        self._shift_length = nres
        if nres < 0:
            self._spline_step = -self._spline_step
        indices_in_polymer = p.indices(residues)
        spline_start_index = first_index_to_consider = indices_in_polymer[0]
        spline_end_index = last_index_to_consider = indices_in_polymer[-1]
        if nres > 0:
            first_index_to_consider -= nres
            if first_index_to_consider <0:
                first_index_to_consider = 0
        else:
            last_index_to_consider -= nres
            if last_index_to_consider >= polymer_length:
                last_index_to_consider = polymer_length-1
        
        xr = self._extended_residues = p[first_index_to_consider:last_index_to_consider+1]
        
        xa = self._extended_atoms = self._make_atom_arrays(xr)
        
        # clear any rotamer, backbone and distance restraints on 
        # extended selection
        isolde.release_rotamers(xr)
        isolde.clear_secondary_structure_restraints_for_selection(residues = xr)
        isolde.release_xyz_restraints_on_selected_atoms(sel = xr.atoms)
        
        self._positions_along_spline = (p.indices(xr) - spline_start_index).astype(float)
        
        self._handler = isolde.triggers.add_handler('completed simulation step', self._step_forward)
    
    def release_all(self):
        if self._extended_residues is not None:
            self.isolde.release_xyz_restraints_on_selected_atoms(sel = self._extended_residues.atoms)
    
    def _step_forward(self, *_):
        '''
        Move one step forward along the spline, and adjust the target 
        positions accordingly.
        '''
        isolde = self.isolde
        xr = self._extended_residues
        xa = self._extended_atoms
        if abs(self._current_position_on_spline) >= abs(self._shift_length):
            # Our target has reached the end of the spline. Switch to a 
            # final "polishing" routine
            isolde.triggers.remove_handler(self._handler)
            self._handler = None
            self.finished = True
            #self._handler = self.isolde.triggers.add_handler('completed simulation step', self._final_polish)
        self._current_position_on_spline += self._spline_step
        sp = self._positions_along_spline 
        sp += self._spline_step
        outside = numpy.logical_or(sp < 0, sp > self._spline_length-1)
        inside = numpy.argwhere(numpy.invert(outside)).ravel()
        outside = numpy.argwhere(outside).ravel()
        # Release restraints on atoms outside of the spline
        out_a = xr[outside].atoms
        if len(out_a):
            isolde.release_xyz_restraints_on_selected_atoms(sel = out_a)
        
        # Update the restraints on atoms inside of the spline
        in_a = xa[inside]
        in_pos = sp[inside]
        
        nspl = self.n_spline
        caspl = self._ca_spline
        cspl = self._c_spline
        cbspl = self._cb_spline
        
        splev = interpolate.splev
        restraints = isolde._sim_pos_restr
        k = self.spring_constant
        
        n_targets = numpy.column_stack(splev(in_pos, nspl[0]))
        ca_targets = numpy.column_stack(splev(in_pos, caspl[0]))
        c_targets = numpy.column_stack(splev(in_pos, cspl[0]))
        cb_targets = numpy.column_stack(splev(in_pos, cbspl[0]))
        #all_targets = numpy.column_stack((n_targets, ca_targets, c_targets, cb_targets))
        
        for n_t, ca_t, c_t, cb_t, atoms in zip(n_targets, ca_targets, c_targets, cb_targets, in_a):
            n = atoms[0]
            rest = restraints[n]
            rest.target = n_t
            rest.spring_constant = k
            ca = atoms[1]
            rest = restraints[ca]
            rest.target = ca_t
            rest.spring_constant = k
            c = atoms[2]
            rest = restraints[c]
            rest.target = c_t
            rest.spring_constant = k
            cb = atoms[3]
            if cb is not None:
                rest = restraints[cb]
                rest.target = cb_t
                rest.spring_constant = k
            
        
        
        
        
            
        
        
        
        
       
        
            
        
        
        