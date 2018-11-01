
def compiled_lib_extension():
    import platform
    pname = platform.system()
    if pname == "Windows":
        return "dll"
    elif pname == "Darwin":
        return "dylib"
    return "so"

def voxel_volume(volume):
    import numpy
    from math import sqrt
    a,b,c = volume.data.step
    angles = numpy.radians(volume.data.cell_angles)
    cos_alpha, cos_beta, cos_gamma = numpy.cos(angles)
    return a*b*c*sqrt(1-cos_alpha**2-cos_beta**2-cos_gamma**2
        + 2*cos_alpha*cos_beta*cos_gamma)
