#
#  Copyright 2016 Jon Agirre & The University of York
#  Developed at York Structural Biology Laboratory - Cowtan group
#  Distributed under the terms of the LGPL (www.fsf.org)
#
#  Package containing functions for making it easy to do MR with cryoEM maps
#  To do: account for errors in the scale using control point refinement

import clipper
from clipper_tools import callbacks
from lxml import etree

## Reads EM map, sets origin to 0, pads cell and computes finely-sampled structure factors
#  @param mapin string path to a map that will be read into a clipper.NXmap_float object
#  @param resol estimated resolution (float)
#  @param callback a function that takes care of log string and xml flushing
#  @return a plain text log string, an XML etree and a clipper.HKL_data_F_phi_float object

def prepare_map ( mapin = "",
                  resol = 8.0,
                  callback = callbacks.interactive_flush ) :

    ## Reads numpy array, determines the extent of the electron density
    #  @param numpy_in a numpy array containing grid points
    #  @param tolerance number of points in a plane with value greater than 1 sigma
    #  @return a vector of grid indices: (min_u, max_u, min_v, max_v, min_w, max_w)

    def determine_extent (  numpy_in, tolerance ) :
    
        log_string = ""
        min = clipper.Coord_orth()
        max = clipper.Coord_orth()
        
        map_mean = numpy.mean(map_numpy)
        map_std  = numpy.std (map_numpy)
        
        mask = map_numpy > map_mean + map_std
        
        sum_u = sum(sum(mask))
        sum_w = sum(sum(numpy.transpose(mask)))
        sum_v = sum(numpy.transpose(sum(mask)))
        
        log_string += "\n  >> dumping 1D summaries of the map's content:\n\n  >> U:\n %s\n" % sum_u
        log_string += "\n  >> V:\n %s\n" % sum_v
        log_string += "\n  >> W:\n %s\n" % sum_w
        
        point_list = [ ]
    
        for idx_u, val_u in enumerate(sum_u) :
            if val_u > tolerance :
                point_list.append ( idx_u )

        min_u = point_list[0]
        max_u = point_list[-1]
        
        log_string += "\n  >> First meaningful U: %i ; Last meaningful U: %i" \
              % (min_u, max_u)

        point_list = [ ]
    
        for idx_v, val_v in enumerate(sum_v) :
            if val_v > tolerance :
                point_list.append ( idx_v )

        min_v = point_list[0]
        max_v = point_list[-1]

        log_string += "\n  >> First meaningful V: %i ; Last meaningful V: %i" \
              % (min_v, max_v)

        point_list = [ ]
    
        for idx_w, val_w in enumerate(sum_w) :
            if val_w > tolerance :
                point_list.append ( idx_w )

        min_w = point_list[0]
        max_w = point_list[-1]

        log_string += "\n  >> First meaningful W: %i ; Last meaningful W: %i\n" \
              % (min_w, max_w)

        extent = [ min_u, max_u, min_v, max_v, min_w, max_w ]
        
        return extent, log_string
    
        ################# end determine_extent ################

    ############### main function ################

    # create log string so console-based apps get some feedback
    log_string = "\n  >> clipper_tools: mr_from_em.structure_factors"
    log_string += "\n            mapin: %s" % mapin
    log_string += "\n            resol: %s" % resol

    # create XML tree, to be merged in a global structured results file
    xml_root = etree.Element('structure_factors')
    xml_root.attrib['mapin'] = mapin
    xml_root.attrib['resol'] = str ( resol )
    callback( log_string, xml_root  )

    nxmap = clipper.NXmap_double( )
    xmap  = clipper.Xmap_double ( )
    map_file = clipper.CCP4MAPfile( )
    sg = clipper.Spacegroup.p1()
    resol *= 0.9
    resolution = clipper.Resolution ( resol )

    # nothing in, nothing out
    if mapin == "" :
        return log_string,xml_root,None
    
    # read the cryoEM map into nxmap, get map data irrespective of origin
    map_file.open_read ( mapin )
    map_file.import_nxmap_double ( nxmap )
    map_file.close_read()
    log_string += "\n  >> file %s has been read as nxmap" % mapin
    callback( log_string, xml_root )

    # read the cryoEM map into xmap to get cell dimensions, etc.
    map_file.open_read ( mapin )
    map_file.import_xmap_double ( xmap )
    map_file.close_read()
    log_string += "\n  >> file %s has been read as xmap" % mapin
    callback( log_string, xml_root )
    log_string += "\n  >> cell parameters: %s" % xmap.cell().format()
    log_string += "\n     original translation: %s" % nxmap.operator_orth_grid().trn()

    # put map content in a numpy data structure
    import numpy
    map_numpy = numpy.zeros( (nxmap.grid().nu(), nxmap.grid().nv(), nxmap.grid().nw()), dtype='double')
    log_string += "\n  >> exporting a numpy array of %i x %i x %i grid points" \
               % (nxmap.grid().nu(), nxmap.grid().nv(), nxmap.grid().nw())
    callback( log_string, xml_root  )

    data_points = nxmap.export_numpy ( map_numpy )
    log_string += "\n  >> %i data points have been exported" % data_points
    callback ( log_string, xml_root )
    map_mean = numpy.mean(map_numpy)
    map_stdv = numpy.std(map_numpy)

    log_string += "\n  >> map mean (stdev): %.4f (%.4f)" % (map_mean, map_stdv)

    # compute the extent
    extent, temp_log = determine_extent ( map_numpy, 30 )
    log_string += temp_log
    extent_list = [ extent[1] - extent[0], extent[3] - extent[2], extent[5] - extent[4] ]
    max_extent = max(extent_list)
    
    # create padded xmap and import numpy array
    origin_trans = clipper.vec3_double(extent[0]+((extent[1]-extent[0])/2),
                                       extent[2]+((extent[3]-extent[2])/2),
                                       extent[4]+((extent[5]-extent[4])/2))

    large_a = ( xmap.cell().a() * ( max_extent + xmap.grid_asu().nu())) / xmap.grid_asu().nu()
    large_b = ( xmap.cell().b() * ( max_extent + xmap.grid_asu().nv())) / xmap.grid_asu().nv()
    large_c = ( xmap.cell().c() * ( max_extent + xmap.grid_asu().nw())) / xmap.grid_asu().nw()

    cell_desc = clipper.Cell_descr ( large_a, large_b, large_c, \
                  xmap.cell().alpha(), xmap.cell().beta(), xmap.cell().gamma() )

    large_p1_cell = clipper.Cell ( cell_desc )
    large_grid_sampling = clipper.Grid_sampling ( max_extent + xmap.grid_asu().nu(),
                                                  max_extent + xmap.grid_asu().nv(),
                                                  max_extent + xmap.grid_asu().nw() )

    large_xmap = clipper.Xmap_double ( sg, large_p1_cell, large_grid_sampling )

    log_string += "\n  >> new grid: nu=%i nv=%i nw=%i" % ( large_xmap.grid_asu().nu(),
                                                           large_xmap.grid_asu().nv(),
                                                           large_xmap.grid_asu().nw() )

    log_string += "\n  >> putting map into a large p1 cell..."
    log_string += "\n  >> new cell parameters: %s" % large_p1_cell.format()
    callback( log_string, xml_root )

    large_xmap.import_numpy ( map_numpy )

    # dump map to disk
    map_file = clipper.CCP4MAPfile()
    map_file.open_write ( "mapout_padded.mrc" )
    map_file.export_xmap_double ( large_xmap )
    map_file.close_write()
    log_string += "\n  >> map file mapout_padded.mrc written to disk"
    callback( log_string, xml_root )

    # import it back to nxmap so we can trivially shift the origin
    map_file.open_read ( "mapout_padded.mrc" )
    map_file.import_nxmap_double ( nxmap )
    map_file.close_read()
    log_string += "\n  >> file mapout_padded.mrc has been read back as nxmap"
    callback( log_string, xml_root )

    # now shift the origin
    rtop_zero = clipper.RTop_double ( nxmap.operator_orth_grid().rot(), origin_trans )

    log_string += "\n  >> moving origin..."
    log_string += "\n     original translation: %s  new origin: %s" % ( nxmap.operator_orth_grid().trn(), rtop_zero.trn() )
    callback( log_string, xml_root )

    nxmap_zero = clipper.NXmap_double(nxmap.grid(), rtop_zero )
    nxmap_zero.import_numpy ( map_numpy )

    # dump map to disk
    map_file.open_write ( "mapout_padded_zero.mrc" )
    map_file.export_nxmap_double ( nxmap_zero )
    map_file.close_write()
    log_string += "\n  >> map file mapout_padded_zero.mrc written to disk"
    callback( log_string, xml_root )

    # read it back to an xmap so we can fft-it
    new_xmap = clipper.Xmap_double ()
    map_file.open_read ( "mapout_padded_zero.mrc" )
    map_file.import_xmap_double ( new_xmap )
    map_file.close_read()
    log_string += "\n  >> map file mapout_padded_zero.mrc read back as xmap"
    callback ( log_string, xml_root )

    # create HKL_info using user-supplied resolution parameter
    hkl_info = clipper.HKL_info (sg, large_p1_cell, resolution, True )

    # fft the map
    f_phi = clipper.HKL_data_F_phi_float( hkl_info, large_p1_cell )
    log_string += "\n  >> now computing map coefficients to %0.1f A resolution..." % resol
    callback( log_string, xml_root )
    
    new_xmap.fft_to ( f_phi )
    log_string += "\n  >> writing map coefficients to MTZ file mapout_padded_zero.mtz"
    callback( log_string, xml_root )

    # setup an MTZ file so we can export our map coefficients
    mtzout  = clipper.CCP4MTZfile()
    mtzout.open_write ( "mapout_padded_zero.mtz" )
    mtzout.export_hkl_info ( f_phi.hkl_info() )
    mtzout.export_hkl_data ( f_phi, "*/*/[F, PHI]" )
    mtzout.close_write()
    log_string += "\n  >> all done"
    callback( log_string, xml_root )

    return log_string,xml_root,f_phi

