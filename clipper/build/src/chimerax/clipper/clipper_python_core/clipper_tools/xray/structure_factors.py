import clipper
from lxml import etree

## Computes structure factors with or without bulk solvent correction
#  @param fsig a clipper.HKL_data_F_sigF_float object
#  @param hklinfo a clipper.HKL_info object
#  @param mmol a clipper.MiniMol object
#  @param bulk_solvent boolean parameter, set True for turning bulk solvent correction on
#  @return a plain text log string, an XML etree and a clipper.HKL_data_F_phi_float object

def structure_factors ( fsigf=None, hklinfo=None, mmol=None, bulk_solvent=True ) :

    log_string = "\n  >> clipper_tools: structure_factors"
    log_string += "\n     bulk_solvent: %s" % bulk_solvent

    xml_root = etree.Element('structure_factors')
    xml_root.attrib['bulk_solvent'] = str ( bulk_solvent )

    crystal = clipper.MTZcrystal()
    atoms = mmol.atom_list()
    
    fc = clipper.HKL_data_F_phi_float ( hklinfo, crystal )
    
    if bulk_solvent :
        sfcb = clipper.SFcalc_obs_bulk_float()
        sfcb ( fc, fsigf, atoms )
        
        bulkfrc = sfcb.bulk_frac();
        bulkscl = sfcb.bulk_scale();
        
        etree.SubElement(xml_root, 'bulk_fraction').text = str ( bulkfrc )
        etree.SubElement(xml_root, 'bulk_scale').text    = str ( bulkscl )
        
        log_string += "\n    bulk_fraction: %f " % bulkfrc
        log_string += "\n    bulk_scale: %f " % bulkscl
    
    else :
        sfc = clipper.SFcalc_obs_base_float()
        sfc ( fc, fsigf, atoms )
    
    log_string += "\n  << structure_factors has finished\n"
    xml_root.attrib['ok']    = 'yes'
    
    return log_string, xml_root, fc
