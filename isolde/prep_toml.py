import openmm, os
lib_dir = openmm.version.openmm_library_path
include_dir = os.path.abspath(os.path.join(lib_dir, '..', 'include'))
with open("pyproject.toml.in", "rt") as f:
    output = f.read()\
        .replace("OPENMM_LIB_DIR", repr(lib_dir).replace("'",""))\
        .replace("OPENMM_INCLUDE_DIR", repr(include_dir).replace("'",""))
with open("pyproject.toml", "wt") as f:
    f.write(output)
