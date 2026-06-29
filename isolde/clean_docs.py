# Remove the Sphinx HTML output tree (src/docs/user) and its .doctrees cache so
# the next docs build is a full rebuild. Sphinx builds incrementally and never
# deletes orphaned output, so pages whose source was renamed/removed otherwise
# linger -- and would ship in the packaged bundle. This is a clean-out ONLY; run
# a docs build afterwards (make_docs.bat / make docs).
#
# Run with ChimeraX's Python (kept platform-agnostic that way), e.g.:
#   ChimeraX-console --nogui --safemode --exit --script clean_docs.py
import os
import shutil

_here = os.path.dirname(os.path.abspath(__file__))
_out = os.path.join(_here, "src", "docs", "user")
if os.path.isdir(_out):
    shutil.rmtree(_out, ignore_errors=True)
    print("clean_docs: removed %s" % _out)
else:
    print("clean_docs: nothing to remove (%s)" % _out)
