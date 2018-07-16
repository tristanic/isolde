#!/bin/csh -f
./rfltest << eof
1
eof
./rfltest << eof
146
eof
./symtest << eof
1
eof
./symtest << eof
19
eof
./symtest << eof
210
eof
./mtztest
mtzdump hklin junk.mtz < /dev/null
mtzdump hklin junk2.mtz < /dev/null
