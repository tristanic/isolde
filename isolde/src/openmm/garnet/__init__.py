# @Author: Tristan Croll
# @Date:   24-Aug-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 24-Aug-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Interim integration of the **garnet-isolde** graph-ML force field into ISOLDE.

This subpackage lets an ISOLDE interactive simulation be parameterised entirely
by garnet-isolde (`garnet_core`), bypassing OpenMM's AMBER template matching.
It is the demonstration-grade first step of the "force field on the model"
project (see ``FORCEFIELD_ON_MODEL_BRIEFING.md``); the parameter store here is
pure Python, deferring the C++ ``FFAtom``/``FFBond`` objects to full Track B.

Design (see the plan the code implements):

* **Parameterise the whole model, then subset.** garnet parameterises from
  topology alone (no coordinates), so we run it once over the complete model and
  cache the result keyed to ChimeraX atoms (:class:`GarnetParameters`). Each
  simulation's OpenMM ``System`` is then built from ISOLDE's usual mobile+fixed
  selection by *reading* that cache — the chemically-incomplete outer edge of the
  fixed buffer needs no capping, exactly as with the AMBER path.
* **AMBER stays the untouched legacy default.** Nothing here runs unless the
  session's force field is ``'garnet'``.

Requires ``garnet_core`` (+ torch, torch_geometric) importable in ChimeraX's
Python; torch is imported lazily so this module loads even when it is absent.
'''

from .params_cache import GarnetParameters, default_checkpoint_path, get_garnet_parameters
from .system_builder import build_garnet_system

__all__ = [
    'GarnetParameters',
    'default_checkpoint_path',
    'get_garnet_parameters',
    'build_garnet_system',
]
