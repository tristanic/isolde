import openmm
lib_dir = openmm.version.openmm_library_path
with open("pyproject.toml.in", "rt") as f:
    output = f.read().replace("OPENMM_LIB_DIR", lib_dir)
with open("pyproject.toml", "wt") as f:
    f.write(output)
