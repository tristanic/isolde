import sphinx
import os

from chimerax import app_bin_dir
sphinx.main(argv=[os.path.join(app_bin_dir, 'ChimeraX'),
    os.path.join(os.getcwd(), 'source'),
    os.path.join(os.getcwd(), '..', 'src', 'doc')])
