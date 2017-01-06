## @package coordinate_kicks
# This module provides functions for guided/random distortion 
# of models, with a focus on debiasing after Molecular Replacement 

import clipper
from lxml import etree
from clipper_tools.general import is_aminoacid


## Randomly kicks the atoms in mmol, with the specified frequency and rmsd
#  @param mmol input molecule
#  @param amplitude maximum rmsd to kick the atoms (random orientation)
#  @param frequency target frequency with which atoms will be kicked
#  @return a plain text log string and an XML etree
 
def atom_kicks ( mmol=None, amplitude=0.0, frequency=0.0 ) :

    log_string = "\n  >> clipper_tools: atom_kicks"
    log_string += "\n     amplitude: %f" % amplitude
    log_string += "\n     frequency: %f" % frequency

    xml_root = etree.Element('atom_kicks')
    xml_root.attrib['amplitude']    = str(amplitude)
    xml_root.attrib['frequency']    = str(frequency)

    if mmol is None :
        log_string += "\n     ERROR: no valid molecule object supplied\n\n"
        log_string += "\n  << atom_kicks has finished\n"
        xml_root.attrib['ok']    = 'no'
        return log_string, xml_root

    if amplitude == 0.0 or frequency == 0.0 :
        log_string += "\n     ERROR: cannot compute kicks with zero amplitude and/or frequency\n\n"
        log_string += "\n  << atom_kicks has finished\n"
        xml_root.attrib['ok']    = 'no'
        return log_string, xml_root

    if frequency < 0.0 or frequency > 100.0 :
        log_string += "\n     ERROR: frequency is not in the (0,100] range\n\n"
        log_string += "\n  << atom_kicks has finished\n"
        xml_root.attrib['ok']    = 'no'
        return log_string, xml_root

    import random

    model = mmol.model()

    kicked = not_kicked = 0

    for chain in model :
        for residue in chain :
            for atom in residue :
                if random.uniform (0.0, 100.0) < frequency :
                    sign = random.choice ( [-1, 1] )
                    x = atom.coord_orth().x() + sign * random.uniform ( 0.01, amplitude )
                    sign = random.choice ( [-1, 1] )
                    y = atom.coord_orth().y() + sign * random.uniform ( 0.01, amplitude )
                    sign = random.choice ( [-1, 1] )
                    z = atom.coord_orth().z() + sign * random.uniform ( 0.01, amplitude )
                    coords = clipper.Coord_orth(x,y,z)
                    atom.set_coord_orth ( coords )
                    kicked += 1
                else :
                    not_kicked += 1

    log_string += "\n     %i atoms have been kicked, while %i remain in their original positions" % (kicked, not_kicked )
    log_string += "\n     That means %f percent of the atoms have been kicked" % ((kicked / (kicked + not_kicked)) *100.0 )

    log_string += "\n  << atom_kicks has finished\n"
    xml_root.attrib['ok']    = 'yes'

    return log_string, xml_root



## Randomly kicks fragments in mmol, with the specified frequency and rmsd
#  @param mmol input molecule, gets modified
#  @param min_length minimum length for the fragments
#  @param max_length maximum length for the fragments
#  @param amplitude maximum rmsd to kick the fragments (random orientation)
#  @param frequency target frequency with which atoms will be kicked
#  @return a plain text log string and an XML etree

def fragment_kicks ( mmol=None, min_length=0, max_length=0, frequency=0.0 ) :
    
    log_string = "\n  >> clipper_tools: fragment_kicks"
    log_string += "\n     frequency: %f" % frequency

    xml_root = etree.Element('fragment_kicks')
    xml_root.attrib['frequency']    = str(frequency)
    
    if mmol is None :
        log_string += "\n     ERROR: no valid molecule object supplied\n\n"
        log_string += "\n  << fragment_kicks has finished\n"
        xml_root.attrib['ok']    = 'no'
        return log_string, xml_root

    if frequency < 0.0 or frequency > 100.0 :
        log_string += "\n     ERROR: frequency is not in the (0,100] range\n\n"
        log_string += "\n  << fragment_kicks has finished\n"
        xml_root.attrib['ok']    = 'no'
        return log_string, xml_root
    
    if min_length > max_length or max_length < min_length :
        log_string += "\n     ERROR: interval is badly defined\n\n"
        log_string += "\n  << fragment_kicks has finished\n"
        xml_root.attrib['ok']    = 'no'
        return log_string, xml_root
    
    model = mmol.model()
    
    n_fragments_kicked = 0
    import random
    rotation_op = clipper.RTop_orth_identity()
    
    for chain in model :
        ## don't want to propagate kick to another chain
        kick_extent = 0
        for residue in chain :
            if kick_extent > 0:

                ## continue kicking until kick_extent is 0
                residue.transform ( rotation_op )
                kick_extent -= 1
                
                if kick_extent == 0 :
                    n_fragments_kicked += 1 

            elif random.uniform (0.0, 100.0) < frequency and is_aminoacid ( residue.type().trim() ) :
                
                ## decide on the characteristics of the displacement
                kick_extent = random.choice ( range ( min_length, max_length ) )
                kick_attribs = etree.SubElement(xml_root, "kicked_fragment")
                kick_attribs.attrib['start_residue'] = str(residue.id())
                kick_attribs.attrib['kick_extent'] = str ( kick_extent )
                
                ## physically and chemically meaningless - we're not trying to generate
                #  a sensible model, but to explore other conformations quickly
                vertex_c_idx  = residue.lookup("C", clipper.UNIQUE)
                vertex_ca_idx = residue.lookup("CA",clipper.UNIQUE)
                if vertex_c_idx is not -1 and vertex_ca_idx is not -1:
                    vertex_c  = residue[vertex_c_idx]
                    vertex_ca = residue[vertex_ca_idx]
                    
                    vec_rotate = clipper.vec3_float(
                                        vertex_ca.coord_orth().x() - vertex_c.coord_orth().x(),
                                        vertex_ca.coord_orth().y() - vertex_c.coord_orth().y(),
                                        vertex_ca.coord_orth().z() - vertex_c.coord_orth().z())
                    
                    vec_ = vec_rotate.unit()
                    import math
                    
                    rotation = random.uniform ( 0, 6.28 )
                    sin_ = math.sin(rotation)
                    cos_ = math.cos(rotation)
                    rotation_matrix = clipper.mat33_ftype(
                                        cos_ + math.pow(vec_[0],2)*(1-cos_), 
                                        vec_[0]*vec_[1]*(1-cos_)-vec_[2]*sin_, 
                                        vec_[0]*vec_[2]*(1-cos_)+vec_[1]*sin_,                
                                        vec_[0]*vec_[1]*(1-cos_)+vec_[3]*sin_,
                                        cos_ + math.pow(vec_[1],2)*(1-cos_),
                                        vec_[1]*vec_[2]*(1-cos_)-vec_[0]*sin_, 
                                        vec_[0]*vec_[2]*(1-cos_)-vec_[1]*sin_,                
                                        vec_[1]*vec_[2]*(1-cos_)+vec_[0]*sin_,
                                        cos_ + math.pow(vec_[2],2)*(1-cos_)   )
                    
                    rotation_op = clipper.RTop_orth(rotation_matrix)
    log_string += "\n     Number of fragments kicked: %i" % n_fragments_kicked
    log_string += "\n  << fragment_kicks has finished\n"
    kick_attribs.attrib['fragments_kicked'] = str ( n_fragments_kicked )                                
    return log_string, xml_root
    
    
    
## Randomly kicks fragments in mmol, with the specified frequency and rmsd
#  @param matrix_of_vectors a 3D matrix of displacement vectors  
#  @return a plain text log string and an XML etree

def directed_fragment_kicks ( mmol=None, matrix_of_vectors=None):
    log_string = "\n  >> clipper_tools: directed_fragment_kicks"
    xml_root = etree.Element('directed_fragment_kicks')
    
    if matrix_of_vectors is None :
        log_string += "\n     ERROR: a matrix containing displacement vectors must be provided"
        log_string += "\n  << directed_fragment_kicks has finished\n"
        xml_root.attrib['ok']    = 'no'
        return log_string, xml_root
    
    xml_root = etree.Element('directed_fragment_kicks')
    xml_root.attrib['shape'] = matrix_of_vectors.shape()
    
    return log_string, xml_root
