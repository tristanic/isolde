# @Author: Tristan Croll
# @Date:   27-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 27-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Standalone MCP server that exposes a running ISOLDE session as typed MCP tools.

This package is a *separate process* from ChimeraX: it talks to ISOLDE only over
the authenticated localhost REST API (``isolde remote rest start``), so it has no
ChimeraX dependency and works identically against a live-GUI or ``--offscreen``
ISOLDE. See server.py and README.md.
'''
