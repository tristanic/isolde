# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll

_openmm_initialized = False
def initialize_openmm():
    # On linux need to set environment variable to find plugins.
    # Without this it gives an error saying there is no "CPU" platform.
    global _openmm_initialized
    if not _openmm_initialized:
        _openmm_initialized = True
        from sys import platform
        if platform in ['linux', 'darwin']:
            from os import environ, path
            from chimerax import app_lib_dir
            environ['OPENMM_PLUGIN_DIR'] = path.join(app_lib_dir, 'plugins')
# initialize_openmm()


def _version():
    import pkg_resources
    return pkg_resources.require('ChimeraX-ISOLDE')[0].version

__version__ = _version()

from chimerax.core.toolshed import BundleAPI
from . import geometry


class _MyAPI(BundleAPI):

    @staticmethod
    def initialize(session, bundle_info):
        initialize_openmm()

    @staticmethod
    def get_class(class_name):
        # 'get_class' is called by session code to get class from bundle that
        # was saved in a session
        if class_name == 'ISOLDE_ToolUI':
            from . import tool
            return tool.ISOLDE_ToolUI
        from .validation import RamaAnnotator, RotamerAnnotator
        from .molobject import (
            _RestraintMgr, MDFFMgr, PositionRestraintMgr,
            TuggableAtomsMgr, _DistanceRestraintMgrBase,
            DistanceRestraintMgr, AdaptiveDistanceRestraintMgr,
            ChiralRestraintMgr, ProperDihedralRestraintMgr,
            _RotamerPreview, RotamerRestraintMgr,
            _Dihedral, ChiralCenter, ProperDihedral,
            Rama, Rotamer, PositionRestraint,
            TuggableAtom, MDFFAtom, DistanceRestraint,
            AdaptiveDistanceRestraint, ChiralRestraint,
            ProperDihedralRestraint, RotamerRestraint,
        )
        from .molarray import (
            ChiralCenters, ProperDihedrals, Ramas, Rotamers, PositionRestraints,
            TuggableAtoms, MDFFAtoms, DistanceRestraints,
            AdaptiveDistanceRestraints, ChiralRestraints,
            ProperDihedralRestraints, RotamerRestraints,
        )
        ct = {
            'RamaAnnotator':                RamaAnnotator,
            'RotamerAnnotator':             RotamerAnnotator,
            '_RestraintMgr':                _RestraintMgr,
            'MDFFMgr':                      MDFFMgr,
            'PositionRestraintMgr':         PositionRestraintMgr,
            'TuggableAtomsMgr':             TuggableAtomsMgr,
            '_DistanceRestraintMgrBase':    _DistanceRestraintMgrBase,
            'DistanceRestraintMgr':         DistanceRestraintMgr,
            'AdaptiveDistanceRestraintMgr': AdaptiveDistanceRestraintMgr,
            'ChiralRestraintMgr':           ChiralRestraintMgr,
            'ProperDihedralRestraintMgr':   ProperDihedralRestraintMgr,
            '_RotamerPreview':              _RotamerPreview,
            'RotamerRestraintMgr':          RotamerRestraintMgr,
            '_Dihedral':                    _Dihedral,
            'ChiralCenter':                 ChiralCenter,
            'ChiralCenters':                ChiralCenters,
            'ProperDihedral':               ProperDihedral,
            'ProperDihedrals':              ProperDihedrals,
            'Rama':                         Rama,
            'Ramas':                        Ramas,
            'Rotamer':                      Rotamer,
            'Rotamers':                     Rotamers,
            'PositionRestraint':            PositionRestraint,
            'PositionRestraints':           PositionRestraints,
            'TuggableAtom':                 TuggableAtom,
            'TuggableAtoms':                TuggableAtoms,
            'MDFFAtom':                     MDFFAtom,
            'MDFFAtoms':                    MDFFAtoms,
            'DistanceRestraint':            DistanceRestraint,
            'DistanceRestraints':           DistanceRestraints,
            'AdaptiveDistanceRestraint':    AdaptiveDistanceRestraint,
            'AdaptiveDistanceRestraints':   AdaptiveDistanceRestraints,
            'ChiralRestraint':              ChiralRestraint,
            'ChiralRestraints':             ChiralRestraints,
            'ProperDihedralRestraint':      ProperDihedralRestraint,
            'ProperDihedralRestraints':     ProperDihedralRestraints,
            'RotamerRestraint':             RotamerRestraint,
            'RotamerRestraints':            RotamerRestraints,
        }
        return ct.get(class_name)

    @staticmethod
    def start_tool(session, tool_name):
        # 'start_tool' is called to start an instance of the tool
        from .tool import ISOLDE_ToolUI
        from chimerax.core import tools
        return tools.get_singleton(session, ISOLDE_ToolUI, 'ISOLDE', create=True)


    @staticmethod
    def register_command(command_name, logger):
        # 'register_command' is lazily called when the command is referenced
        # from chimerax.core.commands import register
        if command_name == 'isolde':
            from . import cmd
            cmd.register_isolde(logger)
        elif command_name in ('rama', '~rama'):
            from .validation import cmd
            cmd.register_rama(logger)
        elif command_name in ('rota', '~rota'):
            from .validation import cmd
            cmd.register_rota(logger)

bundle_api = _MyAPI()
