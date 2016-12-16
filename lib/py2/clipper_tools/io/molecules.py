import clipper
from lxml import etree


## Reads a file in PDB format from the location specified by the input string 
#  @param pdbin a string containing a fully-qualified path to a PDB file
#  @return a log string, an etree (XML) and a clipper.MiniMol object (or None in case of error)

def read_pdb ( pdbin = "undefined" ) :

    log_string = "\n  >> clipper_tools: read_pdb"
    log_string += "\n     pdbin: %s" % pdbin

    xml_root = etree.Element('read_pdb')
    xml_root.attrib['pdbin'] = pdbin
    
    if pdbin is not "undefined" :

        f = clipper.MMDBfile()
        f.read_file ( pdbin )
        mmol = clipper.MiniMol ()
        f.import_minimol ( mmol )
        
        log_string += "\n  << read_pdb has finished\n"
        xml_root.attrib['ok']    = 'yes'
        
        return log_string, xml_root, mmol
    
    else :
        log_string += "\n  ERROR: No input PDB was supplied"
        log_string += "\n  << read_pdb has finished \n"
        xml_root.attrib['ok'] = 'no'
        return log_string, xml_root, None


## Writes molout to a file in PDB format at the location specified by the input string 
#  @param pdbout a string containing a fully-qualified path to a PDB file
#  @param molout the molecule to write
#  @return a log string and an etree (XML)

def write_pdb ( pdbout = "xyzout.pdb", molout = None ) :

    log_string = "\n  >> clipper_tools: write_pdb"
    log_string += "\n     pdbout: %s" % pdbout

    xml_root = etree.Element('write_pdb')
    xml_root.attrib['pdbout'] = pdbout
    
    if molout is None :
        xml_root.attrib['ok']    = 'no'
        log_string += "\n  ERROR: No input MiniMol was supplied"
        log_string += "\n  << write_pdb has finished \n"
        return log_string,xml_root
    else :
        mmdb = clipper.MMDBfile()
        mmdb.export_minimol( molout );
        mmdb.write_file( pdbout, 0 );
        log_string += "\n  << write_pdb has finished \n"
        xml_root.attrib['ok']    = 'yes'

    return log_string,xml_root
