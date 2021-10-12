# @Author: Tristan Croll <tic20>
# @Date:   24-Oct-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 05-May-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



class Clash_Table_Mgr:

    def __init__(self, isolde):
        self.isolde = isolde
        self.session = isolde.session
        iw = isolde.iw
        mf = self._main_frame = iw._validate_clashes_main_frame
        mf.hide()
        self._stub_frame = iw._validate_clashes_stub_frame
        self._show_button = iw._validate_clashes_show_button
        self._hide_button = iw._validate_clashes_hide_button
        self._update_button = iw._validate_clashes_update_button
        self._geom_clash_table = iw._validate_clashes_geom_table
        self._energy_clash_table = iw._validate_clashes_energy_table

        self._clash_detection_handler = None

        self._sim_start_handler = isolde.triggers.add_handler(
            'simulation started', self._sim_start_cb
        )
        self._sim_end_handler = isolde.triggers.add_handler(
            'simulation terminated', self._sim_end_cb
        )
        self._connect_functions()

    def _connect_functions(self):
        self._show_button.clicked.connect(self._show_clash_frame)
        self._hide_button.clicked.connect(self._hide_clash_frame)
        self._update_button.clicked.connect(self.populate_clash_table)
        self._geom_clash_table.itemClicked.connect(self._show_selected_clash)
        self._energy_clash_table.itemClicked.connect(self._show_selected_clash)

    def _show_clash_frame(self):
        self._stub_frame.hide()
        self._main_frame.show()
        self.populate_clash_table()

    def _hide_clash_frame(self):
        self._stub_frame.show()
        self._main_frame.hide()

    def _set_table_visibility(self):
        sr = self.isolde.simulation_running
        self._geom_clash_table.setVisible(not sr)
        self._energy_clash_table.setVisible(sr)

    def _clear_tables(self):
        self._geom_clash_table.setRowCount(0)
        self._energy_clash_table.setRowCount(0)

    def populate_clash_table(self):
        self._set_table_visibility()
        isolde = self.isolde
        self._clear_tables()
        if isolde.selected_model is None:
            return
        sm = isolde.sim_manager
        sr = isolde.simulation_running
        if sr:
            sc = sm.sim_construct
            atoms = sc.all_atoms
            self.populate_clash_table_from_energies(atoms, *sm.sim_handler.find_clashing_atoms())
        else:
            self.populate_clash_table_from_geometry(isolde.selected_model)

    def _sim_start_cb(self, *_):
        sh = self.isolde.sim_handler
        if self._main_frame.isVisible():
            # Need to wait until the first coordinate update to get clashes
            sh.triggers.add_handler('coord update', self._populate_on_first_coord_update)
        # No need to track this handler. The triggerset will be deleted when the simulation ends
        sh.triggers.add_handler('clash detected', self._severe_clash_cb
        )

    def _sim_end_cb(self, *_):
        if self.isolde.gui_mode:
            if self._main_frame.isVisible():
                # re-running _show_clash_frame() will select and populate the correct
                # table
                self._show_clash_frame()

    def _severe_clash_cb(self, *_):
        self.isolde.sim_manager.pause=True
        from ..dialog import generic_warning
        msg_string = ('ISOLDE has detected severe clashes in the model that the '
            'minimiser is unable to reconcile on its own. The simulation '
            'is still initialised, but cannot continue until these are '
            'corrected. Clicking OK will open the clash validation panel with '
            'a list of atoms currently experiencing extreme forces. In most '
            'cases these can be dealt with by choosing more appropriate rotamers. '
            'For more extreme issues you may need to stop the simulation and '
            'tell ISOLDE to temporarily ignore some residues using the command '
            '"isolde ignore {selection}". Once you have moved the adjacent '
            'atoms into non-clashing positions, you can stop the simulation and '
            'stop ignoring the residues using "isolde ~ignore". \n'
            'NOTE: the true culprit may not necessarily be at the top of the '
            'list. Look for clashes that have no direct path away from each '
            'other (e.g. bonds threaded through rings). After rearranging atoms '
            'you may check to see if this has solved the problem by pressing '
            'the play button.'
            )
        generic_warning(msg_string)
        iw = self.isolde.iw
        st = iw._sim_tab_widget
        st.setCurrentIndex(2)
        self._show_clash_frame()

    def _populate_on_first_coord_update(self, *_):
        if self._main_frame.isVisible():
            self._show_clash_frame()
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def populate_clash_table_from_geometry(self, model):
        '''
        If no simulation is running, use ChimeraX's built-in clash toolset to find
        and list geometric clashes in decreasing order of overlap severity.
        '''
        atoms = model.atoms
        from chimerax.clashes import clashes
        clash_dict = clashes.find_clashes(self.session, atoms, inter_model=False)
        t = self._geom_clash_table
        t.setRowCount(0)
        if not clash_dict:
            return
        from functools import reduce
        # Clash dict is bi-directional, so number of clashes is half the total count
        clash_count = 1/2 * reduce(lambda x,y: x+y, (len(d.keys()) for d in clash_dict.values()))
        t.setRowCount(clash_count)
        seen = set()
        from chimerax.atomic import Atoms
        # Make list of unique clashing atom pairs and their distances
        clash_list = []
        for a1, clashes in clash_dict.items():
            seen.add(a1)
            for a2, dist in clashes.items():
                if a2 in seen:
                    continue
                clash_list.append((Atoms([a1,a2]), dist))
        # Sort clashes in decreasing order of overlap
        clash_list = sorted(clash_list, key=lambda x: x[1], reverse=True)

        from Qt.QtWidgets import QTableWidgetItem
        for i, (catoms, overlap) in enumerate(clash_list):
            a1, a2 = catoms
            r1, r2 = catoms.residues
            data = (
            "{} {}{}: {}".format(r1.name, r1.chain_id, r1.number, a1.name),
            "{} {}{}: {}".format(r2.name, r2.chain_id, r2.number, a2.name),
            "{:0.2f}".format(overlap)
            )
            for j, d in enumerate(data):
                item = QTableWidgetItem(d)
                item.data = catoms
                t.setItem(i, j, item)
        t.resizeColumnsToContents()

    def populate_clash_table_from_energies(self, atoms, clashing_indices, forces):
        '''
        When running a simulation, list atoms that are experiencing extreme forces.
        This is preferable, since clashes in a static molecule can lead to
        extremely stretched bonds on energy minimization, with no van der Waals
        overlap to detect.

        Args:
            * atoms:
                - All currently-simulated atoms, in the order they were added to
                  the simulation
            * clashing_indices:
                - Indices of atoms experiencing excessive forces, sorted in
                  decreasing order of force magnitude
            * forces:
                - The forces on the clashing atoms, in decreasing order
            * clash_table:
                - The table to populate with details of the clashing atoms
        '''
        catoms =atoms[clashing_indices]
        t = self._energy_clash_table
        t.setRowCount(0)
        if not len(catoms):
            return

        forces /= 10 # kJ mol-1 nm-1 to kJ mol-1 A-1
        residues = catoms.residues
        chain_ids = residues.chain_ids
        resnames = residues.names
        resnums = residues.numbers
        anames = catoms.names
        from Qt.QtWidgets import QTableWidgetItem
        from chimerax.atomic import Atoms
        t.setRowCount(len(catoms))
        for i, (a, cid, rname, rnum, aname, force) in enumerate(zip(
                catoms, chain_ids, resnames, resnums, anames, forces)):
            data = (
                cid,
                "{}-{}".format(rname, rnum),
                aname,
                "{:0.3g}".format(force)
            )
            aa = Atoms([a])
            for j, d in enumerate(data):
                item = QTableWidgetItem(d)
                item.data = aa
                t.setItem(i, j, item)
        t.resizeColumnsToContents()

    def _show_selected_clash(self, item):
        atoms = item.data
        residues=atoms.residues
        self.session.selection.clear()
        residues.atoms.selected=True
        from ..navigate import get_stepper
        get_stepper(residues[0].structure).step_to(residues[0])
        # session.selection.clear()
        # atoms.selected = True
