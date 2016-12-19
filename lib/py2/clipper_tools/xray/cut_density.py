import clipper
from lxml import etree

def cut_density ( mmol=None, mapin=None, radius=1.5, bulk_solvent=True ) :

    log_string = "\n  >> clipper_tools: structure_factors"
    log_string += "\n     bulk_solvent: %s" % bulk_solvent

    xml_root = etree.Element('structure_factors')
    xml_root.attrib['bulk_solvent'] = str(bulk_solvent)

    atoms = mmol.atom_list()
    
    ## continue computing a mask when we have faster access to maps/masks

    mapout = clipper.Xmap_float()
    

    return log_string, xml_root, mapout