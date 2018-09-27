from __future__ import division

# TODO: lots of Windows fixes
# TODO: separate window for controling visibility of many models?

# this module is designed to be run only by Coot
# Coot 0.6-pre-1 or greater is required!

from random import random
import subprocess
import traceback
import cPickle
import string
import sys
import os

# =============================================================================
# Copy of code from libtbx.utils for handling unicode
# UnicodeDecodeError is passed instead of Sorry
def to_unicode(text, codec='utf8'):
  '''
  Function for handling text when it is first encountered
  Changes str type to unicode type
  The input is returned unmodified if it is already unicode
  Will convert other types (e.g. int, float) to str
  '''
  if (isinstance(text, unicode)):
    return text
  elif (isinstance(text, str)):
    new_text = text
    try:
      new_text = text.decode(codec)
    except UnicodeDecodeError:
      raise
    finally:
      return new_text
  elif (text is not None):
    return unicode(text)
  else:
    return None

def to_str(text, codec='utf8'):
  '''
  Function for handling text when it is passed to cctbx functions that expect
  the str type
  Changes unicode type to str type
  The input is returned unmodified if it is already str
  Will convert other types (e.g. int, float) to str
  '''
  if (isinstance(text, str)):
    return text
  elif (isinstance(text, unicode)):
    new_text = text
    try:
      new_text = text.encode(codec)
    except UnicodeDecodeError:
      raise
    finally:
      return new_text
  elif (text is not None):
    return str(text)
  else:
    return None
# =============================================================================

#if "COOT_PYTHON_DIR" in os.environ :
#  sys.path.append(os.environ["COOT_PYTHON_DIR"])
phenix_dir = to_unicode(os.environ.get("COOT_PHENIX_PATH", None))
phenix_build = to_unicode(os.environ.get("COOT_PHENIX_BUILD_PATH", None))
if (sys.platform != "win32") :
  home_dir = to_unicode(os.environ.get("HOME"))
  if (phenix_build is not None) :
    os.environ["PATH"]+=":%s:%s:%s" % (
      os.path.join(phenix_build, "probe", "exe"),
      os.path.join(phenix_build, "reduce", "exe"),
      os.path.join(phenix_build, "bin"))
else :
  home_dir = to_unicode(to_unicode(os.environ.get("HOMEPATH")))
  if (phenix_build is not None) :
    os.environ["PATH"]+=";%s;%s;%s" % (
      os.path.join(phenix_build, "probe", "exe"),
      os.path.join(phenix_build, "reduce", "exe"),
      os.path.join(phenix_build, "bin"))

import gtk
import gobject
import coot # XXX: coot_utils broken?
try :
  import coot_python
except Exception, e :
  print "Could not import coot_python module!"
  print "PHENIX GUI extensions will be disabled."
  class empty (object) :
    def main_menubar (self) :
      return None
    def main_toolbar (self) :
      return None
  coot_python = empty()

import socket
#socket.setdefaulttimeout(0.01)
from SimpleXMLRPCServer import SimpleXMLRPCServer
from xmlrpclib import ServerProxy

def coot_log (f) :
  def log_wrapper (self, *args, **kwds) :
    _args = list(args)
    _kwds = dict(kwds)
    str_args = ", ".join([ str(arg) for arg in _args ])
    str_kwds = ", ".join([ "%s=%s" % (kwd, str(val))
                           for kwd, val in _kwds.iteritems() ])
    call_signature = []
    if str_args != "" : call_signature.append(str_args)
    if str_kwds != "" : call_signature.append(str_kwds)
    try :
      print "%s.%s(%s)" % (self.__class__.__name__, f.__name__,
        ", ".join(call_signature))
      sys.stdout.flush()
    except IOError, e :
      pass
    return f(self, *args, **kwds)
  return log_wrapper

def coot_startup () :
  print "Loading PHENIX Coot extensions..."
  gui = coot_phenix_interface()
  set_console_display_commands_hilights(0, 0, 0)
  set_find_hydrogen_torsions(1)
  #set_do_probe_dots_on_rotamers_and_chis(1)
  #set_do_probe_dots_post_refine(1)
  set_rotamer_lowest_probability(0.5)
  try :
    set_nomenclature_errors_on_read("ignore")
  except Exception :
    pass

########################################################################
# PHENIX INTERFACE

#--- XML-RPC server object
class coot_xmlrpc_server (SimpleXMLRPCServer) :
  def __init__ (self, addr, phenix_interface) :
    self.phenix_interface = phenix_interface
    SimpleXMLRPCServer.__init__(self, addr, logRequests=0)

  def _dispatch (self, method, params) :
    if not self.phenix_interface.enable_xmlrpc :
      return -1
    result = -1
    func = None
    if hasattr(self.phenix_interface, method) :
      func = getattr(self.phenix_interface, method)
    elif hasattr(coot, method) :
      func = getattr(coot, method)
    if not hasattr(func, "__call__") :
      print "%s is not a callable object!" % method
    else :
      try :
        result = func(*params)
      except Exception, e :
        traceback_str = "\n".join(traceback.format_tb(sys.exc_info()[2]))
        raise Exception("%s\nOriginal traceback:%s" % (to_str(e), traceback_str))
      else :
        if result is None :
          result = -1
    return result

#--- manager for communication with phenix
class coot_phenix_interface (object) :
  def __init__ (self) :
    self._client_id = 0
    self._guis = {}
    self._btns = {}
    self._menu = None
    self._results = {}
    self._molprobity_gui = None
    self._molprobity_btn = None
    self._molprobity_results = None
    self._project_list = []
    self._current_project = None
    self._update_btn = None
    self._probe_file = None
    self.enable_xmlrpc = True
    self.xmlrpc_server = None
    self._start_model = None
    self._current_model = None
    self._current_maps = None
    self._compared_models_and_maps = []
    self._highlighted_mol = None
    self._tmp_dir = None
    self.phenix_server = None
    self.settings = { "map_colour" : (0.0, 0.5, 1.0),
                      "diff_map_colour" : (0.0, 1.0, 0.0),
                      "n_map_colour" : (0.25, 0.0, 1.0),
                      "n_diff_map_colour" : (0.5, 1.0, 0.0),
                      "anom_map_colour" : (1.0, 1.0, 0.0),
                      "iso_diff_map_colour" : (0.5, 0.0, 1.0), }
    set_tip_of_the_day_flag(0)
    menubar = coot_python.main_menubar()
    toolbar = coot_python.main_toolbar()

    # add GUI objects
    if menubar is not None :
      menu = coot_menubar_menu("PHENIX")
      #add_simple_coot_menu_menuitem(menu, "AutoBuild...",
      #  self.launch_autobuild)
      add_simple_coot_menu_menuitem(menu, "Reload MolProbity to-do list",
        self.reload_molprobity_gui)
      #add_simple_coot_menu_menuitem(menu,
      #  "Reload real-space correlation to-do list", self.reload_rsc_gui)
      self._menu = menu
    if toolbar is not None :
      try :
        toolbar.insert(gtk.SeparatorToolItem(), -1)
      except Exception, e :
        print "Error: %s" % to_str(e)
      self.h_button = gtk.ToggleToolButton()
      self.h_button.set_label("Hydrogens off")
      self.h_button.set_is_important(True)
      toolbar.insert(self.h_button, -1)
      self.h_button.connect("clicked", self.toggle_hydrogens)
      self.h_button.set_active(True)
      self.h_button.show()

    # start XML-RPC server
    if "CCTBX_COOT_PORT" in os.environ :
      port = string.atoi(os.environ["CCTBX_COOT_PORT"])
      try :
        self.xmlrpc_server = coot_xmlrpc_server(("127.0.0.1", port), self)
        self.xmlrpc_server.socket.settimeout(0.01)
      except Exception, e :
        print "Error starting XML-RPC server:"
        print str(e)
      else :
        print "xml-rpc server running on port %d" % port
        # timeout used to be set to whatever the Phenix preferences have it as,
        # but 250ms seems to be a universally good choice, and much shorter
        # intervals screw up the Coot GUI (at least on Windows)
        gobject.timeout_add(250, self.timeout_func)
        if toolbar is not None :
          self._update_btn = gtk.ToggleToolButton()
          self._update_btn.set_label("Connected to PHENIX")
          self._update_btn.set_is_important(True)
          toolbar.insert(self._update_btn, -1)
          self._update_btn.connect("clicked", self.toggle_update)
          self._update_btn.set_active(True)
          self._update_btn.show()

    if False : # FIXME #"CCTBX_XMLRPC_PORT" in os.environ :
      port = string.atoi(os.environ["CCTBX_XMLRPC_PORT"])
      self.phenix_server = ServerProxy(uri="http://127.0.0.1:%d/RPC2" % port)
      if self._menu is not None :
        add_simple_coot_menu_menuitem(self._menu, "Run phenix.refine...",
          self.launch_phenix_refine)
      if toolbar is not None :
        pass
        #phenix_refine_button = gtk.ToolButton(icon_widget=None,
        #                                      label="Refine in PHENIX")
        #phenix_refine_button.set_is_important(True)
        #toolbar.insert(phenix_refine_button, -1)
        #phenix_refine_button.connect("clicked", self.launch_phenix_refine)
        #phenix_refine_button.show()

  #---------------------------------------------------------------------
  # EVENT HANDLERS
  def toggle_hydrogens (self, *args) :
    if self.h_button.get_active() :
      self.h_button.set_label("Hydrogens on")
      for imol in model_molecule_list() :
        set_draw_hydrogens(imol, True)
    else :
      self.h_button.set_label("Hydrogens off")
      for imol in model_molecule_list() :
        set_draw_hydrogens(imol, False)

  def toggle_update (self, *args) :
    if self._update_btn is not None :
      if self._update_btn.get_active() :
        self._update_btn.set_label("Connected to PHENIX")
        self.enable_xmlrpc = True
      else :
        self._update_btn.set_label("Connect to PHENIX")
        self.enable_xmlrpc = False

  def launch_phenix_refine (self, *args) :
    molecule_chooser_gui("Select model for refinement",
      self._launch_phenix_refine)

  def _launch_phenix_refine (self, pdb_mol) :
    if pdb_mol is not None :
      out_file = strip_extension(molecule_name(pdb_mol))
      out_file += "_coot.pdb"
      write_pdb_file(pdb_mol, out_file)
      try :
        print "starting phenix.refine via XML-RPC..."
        self.phenix_server.Refine_refine_new_model(out_file)
      except Exception, e :
        print e
      else :
        show_message(
"""Coot has saved your model as %s and sent it back to PHENIX for \
refinement.  Any previously specified non-PDB input files will be left alone. \
You should now return to the phenix.refine window and adjust the settings, \
then click the "Run" button.""" % out_file)

  def launch_autobuild (self, *args) :
    pass # TODO

  # we have to handle requests to remove them from the queue, but if
  # self.enable_xmlrpc is False, they will simply be ignored.
  def timeout_func (self, *args) :
    if self.xmlrpc_server is not None :
      self.xmlrpc_server.handle_request()
    return True

  def reload_molprobity_gui (self, *args) :
    if self._results.get("molprobity") is not None :
      gui = self._guis.get("molprobity")
      if gui is None or gui.window is None :
        self.load_molprobity_gui(self._results.get("molprobity"))

  def reload_rsc_gui (self, *args) :
    pass

  #---------------------------------------------------------------------
  # XML-RPC methods
  @coot_log
  def is_alive (self) :
    return True

  @coot_log
  def quit (self) :
    gtk.main_quit()

  @coot_log
  def clear (self) :
    pass # TODO: reset everything

  def set_client_id (self, client_id) :
    self._client_id = client_id

  def update_settings (self, settings_pkl) :
    self.settings = cPickle.loads(settings_pkl)

  def set_projects (self, projects) :
    self._project_list = projects.split(";")

  def set_current_project (self, project_id) :
    self._current_project = project_id

  def set_tmp_dir (self, tmp_dir) :
    self._tmp_dir = to_unicode(tmp_dir)

  @coot_log
  def clear_refinement (self) :
    self.close_tmp_model()
    for object_number in range(number_of_generic_objects()) :
      set_display_generic_object(object_number, 0)
    molprobity_gui = self._guis.get("molprobity")
    if molprobity_gui is not None and molprobity_gui.window is not None :
      molprobity_gui.destroy_window()
    try :
      self._guis.pop("molprobity", None)
      self._results.pop("molprobity", None)
    except KeyError : # FIXME why does this happen???
      pass
    return True

  def read_pdb_no_recenter (self, file_name) :
    file_name = to_str(file_name)
    return handle_read_draw_molecule_with_recentre(file_name, 0)

  def read_pdb_recenter (self, file_name) :
    file_name = to_str(file_name)
    return handle_read_draw_molecule_with_recentre(file_name, 1)

  def read_pdb_as_trace (self, file_name) :
    imol = read_pdb(to_str(file_name))
    graphics_to_ca_representation(imol)
    return imol

  def load_alt_orig_ref_model (self, file_name) :
    imol = self.read_pdb_as_trace(file_name)
    set_show_unit_cell(imol, 1)

  @coot_log
  def update_model (self, pdb_out) :
    pdb_out = to_unicode(pdb_out)
    if not os.path.isfile(pdb_out) :
      print "***error: %s not found" % pdb_out
      return False
    if self._current_model is None :
      imol = read_pdb(to_str(pdb_out))
      set_molecule_bonds_colour_map_rotation(imol, 140)
      self._current_model = imol
    else :
      clear_and_update_model_molecule_from_file(
        self._current_model, to_str(pdb_out))
    return True

  def get_newest_model (self) :
    all_mols = molecule_number_list()
    current_mol = None
    for imol in old_mols :
      if is_valid_model_molecule(imol) :
        current_mol = imol
    return current_mol

  @coot_log
  def close_maps (self) :
    old_maps = map_molecule_list()
    for imol in old_maps :
      close_molecule(imol)
    return True

  def close_existing_map (self, file_name, f_label) :
    for imol in map_molecule_list() :
      map_name = molecule_name(imol)
      fields = map_name.split()
      if fields[0] == file_name and fields[1] == f_label :
        close_molecule(imol)

  @coot_log
  def close_all (self, maps=False) :
    old_mols = molecule_number_list()
    for imol in old_mols :
      if (not have_unsaved_changes_p(imol)) :
        close_molecule(imol)
      else :
        print "Molecule %d has been modified, will not close" % imol
        print molecule_name(imol)
    self._current_model = None
    self._start_model = None
    if maps :
      self.close_maps()
    return True

  @coot_log
  def recenter_and_zoom (self, x, y, z) :
    set_rotation_centre(x, y, z)
    set_zoom(30)
    # rotation_center()
    # spin_zoom_trans(axis, nstep, stepsize, zoom_by, x_rel, y_rel, z_rel)
    graphics_draw()

  def highlight_model_selection (self, selection) :
    if self._highlighted_mol is not None :
      close_molecule(self._highlighted_mol)
      self._highlighted_mol = None
    imol = self.get_newest_model()
    imol2 = new_molecule_by_atom_selection(imol, "// /")
    #additional_representation_by_attributes(imol, chain_id, resno_start,
    #  resno_end, ins_code, representation_type, bonds_box_type=0,
    #  bond_width, draw_hydrogens_flag)

  #--- phenix.refine stuff
  @coot_log
  def show_start_model (self, file_name) :
    if self._start_model is None :
      imol = read_pdb(to_str(file_name))
      set_molecule_bonds_colour_map_rotation(imol, 280)
      self._start_model = imol
    else :
      clear_and_update_model_molecule_from_file(
        self._start_model, to_str(file_name))

  @coot_log
  def close_tmp_model (self) :
    if self._current_model is not None :
      close_molecule(self._current_model)
      self._current_model = None
      if (self._current_maps is not None) :
        for imol in self._current_maps :
          close_molecule(imol)
      self._current_maps = None
      return True
    return False

  @coot_log
  def auto_load_maps (self, map_file, f="2FOFCWT", delf="FOFCWT", phi="PH",
      use_filled=True) :
    map_file = to_str(map_file)
    set_colour_map_rotation_for_map(0)
    imol1 = make_and_draw_map(map_file, f, "%s%s" % (phi, f), "", 0, 0)
    imol2 = make_and_draw_map(map_file, delf, "%s%s" % (phi, delf), "", 0, 1)
    set_scrollable_map(imol1)
    return (imol1, imol2)

  @coot_log
  def load_phenix_refine_temp_files (self, tmp_dir, run_name) :
    tmp_dir = to_unicode(tmp_dir)
    run_name = to_unicode(run_name)
    if None in [tmp_dir, run_name] or not os.path.isdir(tmp_dir) :
      return False
    pdb_tmp = os.path.join(tmp_dir, "tmp.refine.pdb")
    map_tmp = os.path.join(tmp_dir, "tmp.refine_maps.mtz")
    if (not os.path.isfile(pdb_tmp)) or (not os.path.isfile(map_tmp)) :
      print "One or more output files not found:"
      print "  %s" % to_str(pdb_tmp)
      print "  %s" % to_str(map_tmp)
      return False
    self.update_model(to_str(pdb_tmp))
    self._current_maps = self.auto_load_maps(to_str(map_tmp), phi="PH")
    set_colour_map_rotation_for_map(0)
    graphics_draw()
    add_status_bar_text("Current model and maps for %s" % to_str(run_name))
    return True

  @coot_log
  def load_phenix_refine_final_pdb (self, output_dir, file_base) :
    output_dir = to_unicode(output_dir)
    file_base = to_unicode(file_base)
    output_base = os.path.join(output_dir, file_base)
    pdb_out = "%s.pdb" % output_base
    if (not os.path.isfile(pdb_out)) :
      print "ERROR: refinement output file not found!"
      print "       expected %s" % pdb_out
      print "       directory contents:"
      subprocess.call("ls -l %s/" % to_str(output_dir))
      return -1
    pdb_mol = read_pdb(to_str(pdb_out))
    set_molecule_bonds_colour_map_rotation(pdb_mol, 30)
    return True

  @coot_log
  def load_phenix_refine_final_files (self, output_dir, file_base) :
    output_dir = to_unicode(output_dir)
    file_base = to_unicode(file_base)
    if (None in [output_dir, file_base]) or (not os.path.isdir(output_dir)) :
      return False
    if (not os.path.isdir(output_dir)) :
      raise RuntimeError("Output directory does not exist!")
    self.load_phenix_refine_final_pdb(output_dir, file_base)
    output_base = os.path.join(output_dir, file_base)
    map_out_old = "%s_map_coeffs.mtz" % output_base
    map_out_new = output_base + ".mtz"
    map_out = None
    if os.path.isfile(map_out_new) :
      map_out = map_out_new
    elif os.path.isfile(map_out_old) :
      map_out = map_out_old
    if (map_out is not None) :
      self.load_refinement_maps(map_out)
    else :
      print "ERROR: can't find map coefficients!"
      print "listing %s:" % output_dir
      for file_name in os.listdir(output_dir) :
        print "  %s" % file_name
    graphics_draw()
    add_status_bar_text("Showing refinement results in %s" % to_str(output_dir))
    return True

  def load_refinement_maps (self, map_file) :
    map_file = to_unicode(map_file)
    (imol1, imol2) = self.auto_load_maps(map_file, phi="PH")
    fcol = self.settings["map_colour"]
    dcol = self.settings["diff_map_colour"]
    set_map_colour(imol1, *fcol)
    set_map_colour(imol2, *dcol)
    set_colour_map_rotation_for_map(0)
    set_scrollable_map(imol1)
    set_imol_refinement_map(imol1)

  @coot_log
  def load_phenix_refine_final_xn_files (self, output_dir, file_base) :
    output_dir = to_unicode(output_dir)
    file_base = to_unicode(file_base)
    self.load_phenix_refine_final_pdb(output_dir, file_base)
    output_base = os.path.join(output_dir, file_base)
    map_out_old = "%s_map_coeffs.mtz" % output_base
    map_out_new = output_base + ".mtz"
    map_out = None
    if os.path.isfile(map_out_new) :
      map_out = map_out_new
    elif os.path.isfile(map_out_old) :
      map_out = map_out_old
    if (map_out is not None) :
      (imol1, imol2) = self.auto_load_maps(map_out, f="2FOFCWT_xray",
        delf="FOFCWT_xray", phi="PH")
      fcol = self.settings["map_colour"]
      dcol = self.settings["diff_map_colour"]
      set_map_colour(imol1, *fcol)
      set_map_colour(imol2, *dcol)
      (imol3, imol4) = self.auto_load_maps(map_out, f="2FOFCWT_neutron",
        delf="FOFCWT_neutron", phi="PH")
      fcol = self.settings["n_map_colour"]
      dcol = self.settings["n_diff_map_colour"]
      set_map_colour(imol3, *fcol) #dcol[0], dcol[1], dcol[2]) # purple
      set_map_colour(imol4, *dcol)
      set_scrollable_map(imol1)
    else :
      print "Sorry, can't find map coefficients!"
    graphics_draw()
    add_status_bar_text("Showing refinement results in %s" % to_str(output_dir))
    return True

  @coot_log
  def load_ensemble_refinement (self, file_base) :
    file_base = to_unicode(file_base)
    pdb_file = file_base + ".pdb"
    mtz_file = file_base + ".mtz"
    imol = read_pdb(to_str(pdb_file))
    graphics_to_b_factor_representation(imol)
    self.load_refinement_maps(mtz_file)

  @coot_log
  def load_molprobity_gui (self, tmp_dir) :
    tmp_dir = to_unicode(tmp_dir)
    molprobity_gui = self._guis.get("molprobity")
    if getattr(molprobity_gui, "window", None) is not None :
      molprobity_gui.destroy_window()
    # data_file = os.path.join(tmp_dir, "molprobity.pkl")
    probe_file = os.path.join(tmp_dir, "probe.txt")
    self.load_probe_results(probe_file, force_reload=True)
    # self._guis["molprobity"] = coot_molprobity_todo_list_gui(self, data_file)
    self._results["molprobity"] = tmp_dir
    return True

  @coot_log
  def load_probe_results (self, probe_file, force_reload=False) :
    probe_file = to_unicode(probe_file)
    if (self._probe_file is None) or force_reload :
      if os.path.isfile(probe_file) :
        self._probe_file = probe_file
        handle_read_draw_probe_dots_unformatted(to_str(probe_file), 0, 0)
        self.show_probe_dots(True, True)

  @coot_log
  def show_probe_dots (self, show_dots, overlaps_only) :
    n_objects = number_of_generic_objects()
    sys.stdout.flush()
    if show_dots :
      for object_number in range(n_objects) :
        obj_name = generic_object_name(object_number)
        if overlaps_only and not obj_name in ["small overlap", "bad overlap"] :
          sys.stdout.flush()
          set_display_generic_object(object_number, 0)
        else :
          set_display_generic_object(object_number, 1)
    else :
      sys.stdout.flush()
      for object_number in range(n_objects) :
        set_display_generic_object(object_number, 0)

  #--- wizard stuff
  @coot_log
  def load_solve_map (self, map_file) :
    map_file = to_unicode(map_file)
    if not os.path.isfile(map_file) :
      return False
    self.close_existing_map(map_file, "/crystal/dataset/FP")
    set_colour_map_rotation_for_map(0)
    imol = make_and_draw_map(to_str(map_file), "/crystal/dataset/FP", "/crystal/dataset/PHIB", "/crystal/dataset/FOM", 1, 0)
    #set_map_colour(imol, 0, 0.5, 1)
    graphics_draw()
    return True

  @coot_log
  def load_resolve_map (self, map_file, difference_map=0, contour_level=1.0,
      prefix="/crystal/dataset") :
    map_file = to_unicode(map_file)
    if not os.path.isfile(map_file) :
      return False
    f_name = "%s/FP" % prefix
    phi_name = "%s/PHIM" % prefix
    fom_name = "%s/FOMM" % prefix
    self.close_existing_map(map_file, f_name)
    set_colour_map_rotation_for_map(0)
    imol = make_and_draw_map(to_str(map_file), f_name, phi_name, fom_name, 1,
      difference_map)
    set_contour_level_in_sigma(imol, contour_level)
    #set_map_colour(imol, 0, 0.5, 1)
    graphics_draw()
    return True

  @coot_log
  def load_ccp4_style_map (self,
                           map_file,
                           difference_map=0,
                           contour_level=1.0,
                           prefix="/crystal/dataset") :
    map_file = to_unicode(map_file)
    if (not os.path.isfile(map_file)) :
      return False
    self.close_existing_map(map_file, "%s/FWT" % prefix)
    set_colour_map_rotation_for_map(0)
    imol = make_and_draw_map(to_str(map_file), "%s/FWT" % prefix, "%s/PHWT" % prefix,
      "", 0, difference_map)
    set_contour_level_in_sigma(imol, contour_level)
    graphics_draw()
    return True

  def load_new_resolve_map (self, *args, **kwds) :
    return self.load_ccp4_style_map(*args, **kwds)

  @coot_log
  def load_autobuild_map (self, map_file) :
    map_file = to_unicode(map_file)
    return self.load_resolve_map(map_file)

  # phenix.phaser (EP), phenix.autosol
  @coot_log
  def load_phaser_map (self, map_file, phi_name="PHIFWT") :
    map_file = to_unicode(map_file)
    if not os.path.isfile(map_file) :
      return False
    set_colour_map_rotation_for_map(0)
    map1 = self.load_any_map(map_file, "FWT", phi_name)
    graphics_draw()
    add_status_bar_text("Showing weighted F(obs) map from Phaser")
    return True

  # phenix.phaser (MR), phenix.automr
  @coot_log
  def load_phaser_mr_maps (self, map_file) :
    map_file = to_unicode(map_file)
    if not os.path.isfile(map_file) :
      return False
    set_colour_map_rotation_for_map(0)
    map1 = self.load_any_map(map_file, "FWT", "PHWT")
    map2 = self.load_any_map(map_file, "DELFWT", "PHDELWT", is_diff_map=1)
#    set_contour_level_in_sigma(map2, 3.0)
    graphics_draw()
    add_status_bar_text("Showing 2mFo-DFC and mFo-DFC maps from Phaser")
    return True

  @coot_log
  def load_autobuild_overall_best (self, output_dir) :
    self.load_current_overall_best(output_dir)

  # phenix.autosol, phenix.autobuild
  @coot_log
  def load_current_overall_best (self, output_dir) :
    output_dir = to_unicode(output_dir)
    run_name = os.path.basename(output_dir)
    pdb_out = os.path.join(output_dir, "overall_best.pdb")
    pdb_out_fallback = os.path.join(output_dir, "working_best.pdb")
    map_out = os.path.join(output_dir, "overall_best_denmod_map_coeffs.mtz")
    self.close_maps()
    if os.path.exists("%s/FINISHED" % output_dir) :
      if self._current_model is not None :
        close_molecule(self._current_model)
        self._current_model = None
      if os.path.isfile(pdb_out) :
        pdb_mol = read_pdb(to_str(pdb_out))
        set_molecule_bonds_colour_map_rotation(pdb_mol, 30)
      elif (os.path.isfile(pdb_out_fallback)) :
        pdb_mol = read_pdb(to_str(pdb_out_fallback))
        set_molecule_bonds_colour_map_rotation(pdb_mol, 30)
    elif os.path.isfile(pdb_out) :
      self.update_model(pdb_out)
    elif (os.path.isfile(pdb_out_fallback)) :
      self.update_model(pdb_out_fallback)
    self.load_ccp4_style_map(map_out)
    if os.path.isfile(pdb_out) :
      try :
        (r_work, r_free) = extract_phenix_refine_r_factors(to_str(pdb_out))
        add_status_bar_text(
          "Showing current best map and RESOLVE model from %s (R-free = %s)" %
          (to_str(run_name), str(r_free)))
      except AssertionError, e :
        pass
    return True

  # phenix.phase_and_build
  def load_current_build (self, output_dir, pdb_file, map_file,
      program_name="phenix.phase_and_build") :
    output_dir = to_unicode(output_dir)
    pdb_file = to_unicode(pdb_file)
    run_name = os.path.basename(output_dir)
    self.close_maps()
    self.update_model(pdb_file)
    self.load_resolve_map(map_file)
    add_status_bar_text("Loaded current best model from %s at %s" %
      (to_str(run_name), time.strftime("%H:%M:%S", time.localtime())))

  @coot_log
  def load_fobs_map (self, map_file, f_col, phi_col, fom_col="", contour=1.0,
      colour=None) :
    map_file = to_unicode(map_file)
    imol = self.load_any_map(map_file, f_col, phi_col, fom_col)
    if imol is not None :
      set_contour_level_in_sigma(imol, contour)
      if (colour is None) :
        colour = self.settings["map_colour"]
      set_map_colour(imol, *colour)
    return imol

  def load_comparison_map (self, *args, **kwds) :
    kwds['colour'] = (0.5, 0., 1.)
    return self.load_fobs_map(*args, **kwds)

  # phenix.development.fem
  def load_fem_map (self, map_file, f_col) :
    map_file = to_unicode(map_file)
    self.load_fobs_map(map_file=map_file, f_col=f_col, phi_col="PHI"+f_col)
    self.load_comparison_map(map_file=map_file, f_col="2mFoDFc",
      phi_col="PHI2mFoDFc")

  def load_difference_map_coeffs (self, map_file, f_col, phi_col) :
    map_file = to_unicode(map_file)
    imol = self.load_any_map(map_file, f_col, phi_col, is_diff_map=1)
    if imol is not None :
      set_contour_level_in_sigma(imol, 3.0)
      set_map_colour(imol, *self.settings["diff_map_colour"])
    return imol

  def load_anom_map_coeffs (self, map_file, f_col, phi_col) :
    map_file = to_unicode(map_file)
    imol = self.load_any_map(map_file, f_col, phi_col)
    if imol is not None :
      set_contour_level_in_sigma(imol, 3.0)
      set_map_colour(imol, *self.settings["anom_map_colour"])
    return imol

  # phenix.fobs_minus_fobs_map
  def load_iso_diff_map_coeffs (self, map_file, f_col, phi_col) :
    map_file = to_unicode(map_file)
    imol = self.load_any_map(map_file, f_col, phi_col, is_diff_map=1)
    if imol is not None :
      set_contour_level_in_sigma(imol, 3.0)
      set_map_colour(imol, *self.settings["iso_diff_map_colour"])
    return imol

  @coot_log
  def load_any_map (self, map_file, f_col, phi_col, fom_col="", is_diff_map=0) :
    map_file = to_unicode(map_file)
    if not os.path.isfile(map_file) :
      return None
    use_weights = 0
    if fom_col != "" :
      use_weights = 1
    imol = make_and_draw_map(to_str(map_file), f_col, phi_col, fom_col,
                             use_weights, is_diff_map)
    return imol

  def load_ccp4_map (self, map_file, is_difference_map=False) :
    map_file = to_unicode(map_file)
    imol = handle_read_ccp4_map(to_str(map_file), 0)
    if (is_difference_map) :
      set_contour_level_in_sigma(imol, 3.0)
      set_map_colour(imol, *self.settings["diff_map_colour"])
    else :
      set_contour_level_in_sigma(imol, 1.5)
      set_map_colour(imol, *self.settings["map_colour"])

  def update_color_setting (self, color_name, *rgb) :
    self.settings[color_name] = tuple(rgb)

  def start_sculptor_gui (self) :
    script_file = os.path.join(phenix_dir,"phaser","phaser","coot_sculptor.py")
    if os.path.isfile(script_file) :
      run_script(script_file)

  def load_superposed_models_and_maps (self, output_dir) :
    output_dir = to_unicode(output_dir)
    pickle_file = os.path.join(output_dir, ".output.pickle")
    assert os.path.isfile(pickle_file)
    output = cPickle.load(open(pickle_file, "rb"))
    (pdb_1, map_1, pdb_2, map_2, is_difference_map) = output
    read_pdb(to_str(os.path.join(output_dir, pdb_1)))
    read_pdb(to_str(os.path.join(output_dir, pdb_2)))
    if is_difference_map :
      for map_file in [map_1, map_2] :
        self.load_ccp4_map(os.path.join(output_dir, map_file), True)
    else :
      for map_file in [map_1, map_2] :
        self.load_ccp4_map(os.path.join(output_dir, map_file), False)

  def recalculate_probe_dots (self, file_name) :
    file_name = to_unicode(file_name)
    for pdb_mol in molecule_number_list() :
      mol_name = molecule_name(pdb_mol)
      if (mol_name == file_name) :
        tmp_dir = self._tmp_dir
        if (tmp_dir is None) :
          if ("PHENIX_TMP" in os.environ) :
            tmp_dir = os.environ["PHENIX_TMP"]
          else :
            tmp_dir = os.path.join(home_dir, ".phenix", "tmp")
          if (not os.path.isdir(tmp_dir)) :
            raise RuntimeError("Temporary directory (%s) does not exist." %
                               to_str(tmp_dir))
        out_file = os.path.join(tmp_dir, "%d.pdb" % int(random() * 1000000))
        safe_delete(out_file)
        write_pdb_file(pdb_mol, out_file)
        probe_file = os.path.join(tmp_dir, "probe.txt")
        subprocess.call(args="""phenix.reduce %s | phenix.probe -u -q -mc -het -once "alta ogt33 not water" "alta ogt33" - > %s""" %
          (to_str(out_file), to_str(probe_file)), shell=True)
        if (not os.path.isfile(probe_file)) :
          print "Missing output file from PROBE:"
          print probe_file
          raise RuntimeError("Missing output file from PROBE (%s)" %
                             to_str(probe_file))
        safe_delete(out_file)
        self.load_probe_results(probe_file, force_reload=True)

  #---------------------------------------------------------------------
  # methods for structure comparison GUI
  def clear_compared_models (self) :
    self._compared_models_and_maps = []

  def load_partial_model_and_map (self, pdb_file, map_file=None,
      diff_map_file=None) :
    pdb_file = to_unicode(pdb_file)
    map_file = to_unicode(map_file)
    diff_map_file = to_unicode(diff_map_file)
    files_and_ids = []
    pdb_mol = self.read_pdb_recenter(to_str(pdb_file))
    map_mol1 = map_mol2 = None
    if (map_file is not None) and (map_file != "") :
      map_mol1 = handle_read_ccp4_map(to_str(map_file), 0)
      set_contour_level_in_sigma(map_mol1, 1.5)
      set_map_colour(map_mol1, *self.settings["map_colour"])
    if (diff_map_file is not None) and (diff_map_file != "") :
      map_mol2 = handle_read_ccp4_map(to_str(diff_map_file), 1)
      set_contour_level_in_sigma(map_mol2, 3.0) # ???
      set_map_colour(map_mol2, *self.settings["diff_map_colour"])
    structure = comparison_structure(
      pdb_file=pdb_file,
      pdb_mol=pdb_mol,
      map_mol1=map_mol1,
      map_mol2=map_mol2)
    self._compared_models_and_maps.append(structure)

  def load_partial_model_and_map_coeffs (self, pdb_file, mtz_file, phi) :
    pdb_file = to_str(pdb_file)
    mtz_file = to_str(mtz_file)
    files_and_ids = []
    pdb_mol = self.read_pdb_recenter(pdb_file)
    imol1 = imol2 = None
    if mtz_file is not None :
      (imol1, imol2) = self.auto_load_maps(mtz_file, phi=phi)
    structure = comparison_structure(
      pdb_file=pdb_file,
      pdb_mol=pdb_mol,
      map_mol1=imol1,
      map_mol2=imol2)
    self._compared_models_and_maps.append(structure)

  def set_active_map_for_comparison (self, pdb_file) :
    pdb_file = to_str(pdb_file)
    for structure in self._compared_models_and_maps :
      if (structure.pdb_file == pdb_file) :
        structure.set_active_map()
        break

  def set_compared_model_visibility (self, file_name, model, maps) :
    file_name = to_str(file_name)
    for structure in self._compared_models_and_maps :
      if (structure.pdb_file == file_name) :
        structure.set_visibility(model, maps)
        break
    else :
      print "Can't find %s" % file_name
    active = []
    for other in self._compared_models_and_maps :
      if (other.maps_visible) :
        active.append(other)
    if (len(active) == 1) :
      active[0].set_active_map()

  def save_compared_models (self, save_dir=None) :
    save_dir = to_unicode(save_dir)
    files = []
    for structure in self._compared_models_and_maps :
      imol = structure.pdb_mol
      file_name = structure.pdb_file
      file_base, ext = os.path.splitext(os.path.basename(file_name))
      if save_dir is None :
        save_dir = os.path.dirname(file_name)
      out_file = os.path.join(save_dir, "%s_coot.pdb" % file_base)
      write_pdb_file(imol, to_str(out_file))
      files.append(to_str(out_file))
    return ";".join(files)

class comparison_structure (object) :
  def __init__ (self,
                pdb_file,
                pdb_mol,
                map_mol1=None,
                map_mol2=None) :
    self.pdb_file = to_unicode(pdb_file)
    self.pdb_mol = pdb_mol
    self.map_mol1 = map_mol1
    self.map_mol2 = map_mol2
    self.maps_visible = True

  def set_visibility (self, model, maps) :
    set_mol_displayed(self.pdb_mol, model)
    for imol in [ self.map_mol1, self.map_mol2 ] :
      if (imol is not None) :
        set_map_displayed(imol, maps)
      self.maps_visible = maps

  def set_active_map (self) :
    if (self.map_mol1 is not None) :
      print "Setting RSR map to 2mFo-DFc associated with %s" % \
        to_str(self.pdb_file)
      set_imol_refinement_map(self.map_mol1)
      set_scrollable_map(self.map_mol1)

#---------------------------------------------------------------------
# phenix.refine launcher
class coot_phenix_molecule_picker (object) :
  def __init__ (self, hint_text, action_label) :
    self.pdb_mol = None
    self.hkl_file = None
    mol_list = molecule_number_list()
    self.model_mol_list = []
    for imol in mol_list :
      if is_valid_model_molecule(imol) :
        self.model_mol_list.append(imol)
    if len(self.model_mol_list) == 0 :
      pass # TODO: error message?
    elif len(self.model_mol_list) == 1 :
      self.pdb_mol = self.model_mol_list[0]
    else :
      self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
      hbox = gtk.HBox(False, 0)
      vbox = gtk.VBox(False, 0)
      h_sep = gtk.HSeparator()
      self.window.add(vbox)
      label = gtk.Label(hint_text)
      molbox = gtk.VBox(False, 0)
      molbox.pack_start(label, False, False, 5)
      self.mol_choices = gtk.combo_box_new_text()
      fill_option_menu_with_coordinates_mol_options(self.mol_choices)
      molbox.pack_start(self.mol_choices, True, True, 5)
      vbox.pack_start(molbox, True, False, 5)
      go_button = gtk.Button("  %s  " % action_label)
      cancel_button = gtk.Button("  Cancel  ")
      go_button.connect("clicked", self.OnRefine)
      cancel_button.connect("clicked", self.OnDelete)
      vbox.pack_start(h_sep, False, False, 5)
      vbox.pack_start(hbox, False, False, 5)
      hbox.pack_start(go_button, False, False, 5)
      hbox.pack_start(cancel_button, False, False, 5)
      self.window.show_all()

  def select_hkl_file (self, *args) :
    self.fs_window = gtk.FileSelection("file selection")
    self.fs_window.ok_button.connect("clicked", self.fs_ok_sel)
    self.fs_window.cancel_button.connect("clicked",
      lambda w: self.fs_window.destroy())
    self.fs_window.show()

  def fs_ok_sel (self, *args) :
    t = self.fs_window.get_filename()
    self.hkl_entry.set_text(t)
    self.fs_window.destroy()

  def OnDelete (self, *args) :
    self.window.destroy()

  def OnRefine (self, *args) :
    if self.pdb_mol is None :
      import operator # import dependency
      hkl_path = None #self.hkl_entry.get_text()
      mol_info = get_option_menu_active_molecule(self.mol_choices,
        self.model_mol_list)
      (imol_str, pdb_path) = mol_info.split()
      try :
        imol = string.atoi(imol_str)
      except Exception :
        pass
      else :
        self.pdb_mol = imol
        self.hkl_file = hkl_path
    if self.window is not None :
      self.window.destroy()

  def get_refinement_inputs (self) :
    return (self.pdb_mol, self.hkl_file)

#---------------------------------------------------------------------
class coot_extension_gui (object) :
  def __init__ (self, phenix_interface, title) :
    self.phenix_interface = phenix_interface
    self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    scrolled_win = gtk.ScrolledWindow()
    self.outside_vbox = gtk.VBox(False, 2)
    self.inside_vbox = gtk.VBox(False, 0)
    self.window.set_title(title)
    self.inside_vbox.set_border_width(0)
    self.window.add(self.outside_vbox)
    self.outside_vbox.pack_start(scrolled_win, True, True, 0)
    scrolled_win.add_with_viewport(self.inside_vbox)
    scrolled_win.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)

  def finish_window (self) :
    self.outside_vbox.set_border_width(2)
    ok_button = gtk.Button("  Close  ")
    self.outside_vbox.pack_end(ok_button, False, False, 0)
    ok_button.connect("clicked", lambda b: self.destroy_window())
    self.window.connect("delete_event", lambda a, b: self.destroy_window())
    self.window.show_all()

  def destroy_window (self, *args) :
    self.window.destroy()
    self.window = None

  def confirm_data (self, data) :
    for data_key in self.data_keys :
      outlier_list = data.get(data_key)
      if outlier_list is not None and len(outlier_list) > 0 :
        return True
    return False

  def create_property_lists (self, data) :
    for data_key in self.data_keys :
      outlier_list = data[data_key]
      if outlier_list is None or len(outlier_list) == 0 :
        continue
      else :
        frame = gtk.Frame(self.data_titles[data_key])
        vbox = gtk.VBox(False, 2)
        frame.set_border_width(6)
        frame.add(vbox)
        self.add_top_widgets(data_key, vbox)
        self.inside_vbox.pack_start(frame, False, False, 5)
        list_obj = residue_properties_list(
          phenix_interface=self.phenix_interface,
          columns=self.data_names[data_key],
          column_types=self.data_types[data_key],
          rows=outlier_list,
          box=vbox)

# Molprobity result viewer
class coot_molprobity_todo_list_gui (coot_extension_gui) :
  data_keys = [ "rama", "rota", "cbeta", "probe" ]
  data_titles = { "rama"  : "Ramachandran outliers",
                  "rota"  : "Rotamer outliers",
                  "cbeta" : "C-beta outliers",
                  "probe" : "Severe clashes" }
  data_names = { "rama"  : ["Chain", "Residue", "Name", "Score"],
                 "rota"  : ["Chain", "Residue", "Name", "Score"],
                 "cbeta" : ["Chain", "Residue", "Name", "Conf.", "Deviation"],
                 "probe" : ["Atom 1", "Atom 2", "Overlap"] }
  data_types = { "rama" : [gobject.TYPE_STRING, gobject.TYPE_STRING,
                           gobject.TYPE_STRING, gobject.TYPE_FLOAT,
                           gobject.TYPE_PYOBJECT],
                 "rota" : [gobject.TYPE_STRING, gobject.TYPE_STRING,
                           gobject.TYPE_STRING, gobject.TYPE_FLOAT,
                           gobject.TYPE_PYOBJECT],
                 "cbeta" : [gobject.TYPE_STRING, gobject.TYPE_STRING,
                            gobject.TYPE_STRING, gobject.TYPE_STRING,
                            gobject.TYPE_FLOAT, gobject.TYPE_PYOBJECT],
                 "probe" : [gobject.TYPE_STRING, gobject.TYPE_STRING,
                            gobject.TYPE_FLOAT, gobject.TYPE_PYOBJECT] }

  def __init__ (self, phenix_interface, data_file) :
    data_file = to_unicode(data_file)
    data = load_pkl(data_file)
    if not self.confirm_data(data) :
      return
    coot_extension_gui.__init__(self,phenix_interface,"MolProbity to-do list")
    self.dots_btn = None
    self.dots2_btn = None
    self._overlaps_only = True
    self.window.set_default_size(420, 600)
    self.create_property_lists(data)
    self.finish_window()

  def add_top_widgets (self, data_key, box) :
    if data_key == "probe" :
      hbox = gtk.HBox(False, 2)
      self.dots_btn = gtk.CheckButton("Show Probe dots")
      hbox.pack_start(self.dots_btn, False, False, 5)
      self.dots_btn.connect("toggled", self.toggle_probe_dots)
      self.dots2_btn = gtk.CheckButton("Overlaps only")
      hbox.pack_start(self.dots2_btn, False, False, 5)
      self.dots2_btn.connect("toggled", self.toggle_all_probe_dots)
      self.dots2_btn.set_active(True)
      self.toggle_probe_dots()
      box.pack_start(hbox, False, False, 0)

  def toggle_probe_dots (self, *args) :
    if self.dots_btn is not None :
      show_dots = self.dots_btn.get_active()
      overlaps_only = self.dots2_btn.get_active()
      if show_dots :
        self.dots2_btn.set_sensitive(True)
      else :
        self.dots2_btn.set_sensitive(False)
      self.phenix_interface.show_probe_dots(show_dots, overlaps_only)

  def toggle_all_probe_dots (self, *args) :
    if self.dots2_btn is not None :
      self._overlaps_only = self.dots2_btn.get_active()
      self.toggle_probe_dots()

class rsc_todo_list_gui (coot_extension_gui) :
  data_keys = ["by_res", "by_atom"]
  data_titles = ["Real-space correlation by residue",
                 "Real-space correlation by atom"]
  data_names = {}
  data_types = {}

class residue_properties_list (object) :
  def __init__ (self, phenix_interface, columns, column_types, rows, box,
      default_size=(380,200)) :
    self.phenix_interface = phenix_interface
    assert len(columns) == (len(column_types) - 1)
    assert len(rows) == 0 or len(rows[0]) == len(column_types)
    self.liststore = gtk.ListStore(*column_types)
    self.listmodel = gtk.TreeModelSort(self.liststore)
    self.listctrl = gtk.TreeView(self.listmodel)
    self.listctrl.column = [None]*len(columns)
    self.listctrl.cell = [None]*len(columns)
    for i, column_label in enumerate(columns) :
      cell = gtk.CellRendererText()
      column = gtk.TreeViewColumn(column_label)
      self.listctrl.append_column(column)
      column.set_sort_column_id(i)
      column.pack_start(cell, True)
      column.set_attributes(cell, text=i)
    self.listctrl.get_selection().set_mode(gtk.SELECTION_SINGLE)
    for row in rows :
      self.listmodel.get_model().append(row)
    self.listctrl.connect("cursor-changed", self.OnChange)

    sw = gtk.ScrolledWindow()
    w, h = default_size
    if len(rows) > 10 :
      sw.set_size_request(w, h)
    else :
      sw.set_size_request(w, 30 + (20 * len(rows)))
    sw.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
    box.pack_start(sw, False, False, 5)
    inside_vbox = gtk.VBox(False, 0)
    sw.add(self.listctrl)

  def OnChange (self, treeview) :
    selection = self.listctrl.get_selection()
    (model, tree_iter) = selection.get_selected()
    if tree_iter is not None :
      row = model[tree_iter]
      xyz = row[-1]
      if isinstance(xyz, tuple) and len(xyz) == 3 :
        self.phenix_interface.recenter_and_zoom(*xyz)

########################################################################
# eLBOW extensions


########################################################################
# Misc. GUI objects and functions
def show_message (message) :
  dialog = gtk.MessageDialog(None, gtk.DIALOG_MODAL, gtk.MESSAGE_INFO,
    gtk.BUTTONS_OK, message)
  dialog.run()
  dialog.destroy()

def load_pkl (file_name) :
  file_name = to_unicode(file_name)
  pkl = open(file_name, "rb")
  data = cPickle.load(pkl)
  pkl.close()
  return data

# copied from Model.py
def extract_phenix_refine_r_factors (file_name) :
  file_name = to_unicode(file_name)
  assert os.path.isfile(file_name)
  (r_work, r_free) = (None, None)
  lines = open(file_name).readlines()
  for line in lines :
    if line.startswith("REMARK Final: r_work =") :
      r_work = line[23:29]
      r_free = line[39:45]
      break
  return (r_work, r_free)

def safe_delete (file_name) :
  file_name = to_unicode(file_name)
  if os.path.exists(file_name) :
    try :
      os.remove(file_name)
    except OSError :
      print "Can't delete %s" % file_name
      return False
    else :
      return True
  else :
    return True

coot_startup()

#---end
