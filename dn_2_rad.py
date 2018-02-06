#!/usr/bin/env python
"""
SYNOPSIS

    dn_2_rad.py [-h,--help] [-v,--verbose] [-i,--input] [-o,--output]
    

DESCRIPTION

  This program is used to extract the gain parameters and to convert 
  Landsat TM "digital numbers" (DNs) to top-of-atmosphere (TOA) radiances.
  The units of the output are W/(m2*ster*um). We assume that the user
  provides an input band and that the band names follow the standard USGS
  naming convention.

EXAMPLES

            $ ./dn_2_rad.py --input LE71660672004161ASN01_B1.TIF --verbose
            Tue Apr 30 14:26:06 2013
            Start reading LE71660672004161ASN01_B1.TIF
            Using blocksize 7831 x 1
            Creating output LE71660672004161ASN01_B1_TOARAD.tif
            Start reading LE71660672004161ASN01_B2.TIF
            Using blocksize 7831 x 1
            Creating output LE71660672004161ASN01_B2_TOARAD.tif
            Start reading LE71660672004161ASN01_B3.TIF
            Using blocksize 7831 x 1
            Creating output LE71660672004161ASN01_B3_TOARAD.tif
            Start reading LE71660672004161ASN01_B4.TIF
            Using blocksize 7831 x 1
            Creating output LE71660672004161ASN01_B4_TOARAD.tif
            Start reading LE71660672004161ASN01_B5.TIF
            Using blocksize 7831 x 1
            Creating output LE71660672004161ASN01_B5_TOARAD.tif
            Start reading LE71660672004161ASN01_B7.TIF
            Using blocksize 7831 x 1
            Creating output LE71660672004161ASN01_B7_TOARAD.tif
            Tue Apr 30 14:26:31 2013
            TOTAL TIME IN MINUTES: 0.429338383675



EXIT STATUS

    -1 if numpy and/or GDAL aren't present

AUTHOR

    J Gomez-Dans <j.gomez-dans@ucl.ac.uk>

LICENSE

    This script is in the public domain, free from copyrights or restrictions.

"""

import sys
import os
import traceback
import optparse
import time

try:
    import numpy as np
except ImportError:
    print "You need numpy installed"
    sys.exit ( -1 )
try:
  from osgeo import gdal
except ImportError:
    print "You need GDAL installed"
    sys.exit ( -1 )
    
    
GDAL_OPTS = ["COMPRESS=LZW", "INTERLEAVE=PIXEL", "TILED=YES",\
        "SPARSE_OK=TRUE", "BIGTIFF=YES" ]
fname = open('C:/Users/hudjimartsu/Documents/PROJECT/FOREST 2020/TRAINING/PyQgis/DATA/source/LC81220652016230LGN00_MTL.txt', 'r')
def process_metadata ( fname ):
    """A function to extract the relelvant metadata from the
    USGS control file. Returns dicionaries with LMAX, LMIN,
    QCAL_LMIN and QCAL_LMAX for each of the bands of interest."""

    fp = open( fname, 'r') # Open metadata file
    lmax = {} # Dicts to store constants
    lmin = {}
    qc_lmax = {}
    qc_lmin = {}
    gain = {}
    bias = {}

    for line in fp: # 
      # Check for LMAX and LMIN strings
      # Note that parse logic is identical to the first case
      # This version of the code works, but is rather inelegant!
      if ( line.find ("RADIANCE_MULT_BAND") >= 0 ):
          s = line.split("=") # Split by equal sign
          the_band = int(s[0].split("_")[3]) # Band number as integer
          if the_band in [1,2,3,4,5,7]: # Is this one of the bands we want?
              gain[the_band] = float ( s[-1] ) # Get constant as float
      elif ( line.find ("RADIANCE_ADD_BAND") >= 0 ):
          s = line.split("=") # Split by equal sign
          the_band = int(s[0].split("_")[3]) # Band number as integer
          if the_band in [1,2,3,4,5,7]: # Is this one of the bands we want?
              bias[the_band] = float ( s[-1] ) # Get constant as float
      elif ( line.find ("QUANTIZE_CAL_MAX_BAND") >= 0 ):
          s = line.split("=") # Split by equal sign
          the_band = int(s[0].split("_")[4]) # Band number as integer
          if the_band in [1,2,3,4,5,7]: # Is this one of the bands we want?
              qc_lmax[the_band] = float ( s[-1] ) # Get constant as float
      elif ( line.find ("QUANTIZE_CAL_MIN_BAND") >= 0 ):
          s = line.split("=") # Split by equal sign
          the_band = int(s[0].split("_")[4]) # Band number as integer
          if the_band in [1,2,3,4,5,7]: # Is this one of the bands we want?
              qc_lmin[the_band] = float ( s[-1] ) # Get constant as float
      elif ( line.find ("RADIANCE_MAXIMUM_BAND") >= 0 ):
          s = line.split("=") # Split by equal sign
          the_band = int(s[0].split("_")[3]) # Band number as integer
          if the_band in [1,2,3,4,5,7]: # Is this one of the bands we want?
              lmax[the_band] = float ( s[-1] ) # Get constant as float
      elif ( line.find ("RADIANCE_MINIMUM_BAND") >= 0 ):
          s = line.split("=") # Split by equal sign
          the_band = int(s[0].split("_")[3]) # Band number as integer
          if the_band in [1,2,3,4,5,7]: # Is this one of the bands we want?
              lmin[the_band] = float ( s[-1] ) # Get constant as float

    return ( bias, gain, lmax, lmin, qc_lmax, qc_lmin )



def get_metadata ( fname ):
    """This function takes `fname`, a filename (opionally with a path), and
    and works out the associated metadata file"""
    original_fname = os.path.basename ( fname )
    metadata_fname = original_fname.split("_")[0] + "_MTL.txt"
    metadata_fname = os.path.join ( os.path.dirname ( fname ), metadata_fname )
    return metadata_fname
tes1=get_metadata(fname)
def extract_chunk ( the_file, n_blocks=1 ):
    
    """A function that extracts a chunk from a datafile"""
    ds_config = {}
    g = gdal.Open ( the_file )
    block_size = g.GetRasterBand(1).GetBlockSize()
    nx = g.RasterXSize
    ny = g.RasterYSize
    the_bands = np.arange ( g.RasterCount ) + 1
    proj = g.GetProjectionRef()
    geoT = g.GetGeoTransform()
    ds_config['nx'] = nx
    ds_config['ny'] = ny
    ds_config['nb'] = g.RasterCount
    ds_config['geoT'] = geoT
    ds_config['proj'] = proj
    block_size = [ block_size[0]*n_blocks, block_size[1]*n_blocks ]
    # store these numbers in variables that may change later
    nx_valid = block_size[0]
    ny_valid = block_size[1]
    # find total x and y blocks to be read
    nx_blocks = (int)((nx + block_size[0] - 1) / block_size[0]);
    ny_blocks = (int)((ny + block_size[1] - 1) / block_size[1]);
    buf_size = block_size[0]*block_size[1]
    print("Using blocksize %s x %s" %(block_size[0], block_size[1]))
    sys.stdout.flush()
    ################################################################
    # start looping through blocks of data
    ################################################################
    # loop through X-lines
    for X in xrange( nx_blocks ):
        # change the block size of the final piece
        if X == nx_blocks-1:
             nx_valid = nx - X * block_size[0]
             buf_size = nx_valid*ny_valid

        # find X offset
        this_X = X*block_size[0]

        # reset buffer size for start of Y loop
        ny_valid = block_size[1]
        buf_size = nx_valid*ny_valid

        # loop through Y lines
        for Y in xrange( ny_blocks ):
            # change the block size of the final piece
            if Y == ny_blocks-1:
                ny_valid = ny - Y * block_size[1]
                buf_size = nx_valid*ny_valid

            # find Y offset
            this_Y = Y*block_size[1]
            
            buf = g.ReadRaster(this_X, this_Y, nx_valid, ny_valid, \
                buf_xsize=nx_valid, buf_ysize=ny_valid,  \
                band_list= the_bands )
            a = np.frombuffer(buf, dtype=np.uint8 )
            if len(the_bands) > 1:
                data_in = a.reshape(( len(the_bands), ny_valid, nx_valid))
            else:
                data_in = a.reshape(( ny_valid, nx_valid))

            yield ( ds_config, this_X, this_Y, nx_valid, ny_valid, data_in )

def convert_to_radiance ( fname, band, gain, bias, lmax, lmin, qc_lmax, qc_lmin, verbose ):
    """
    This function reads the input file chunk by chunk using ``extract_chunk``,
    and converts from DN to radiance using the methodology given by 
    <http://landsat.usgs.gov/how_is_radiance_calculated.php>. 
    """
    # This variable is used to create the output dataset when the first 
    # chunk of data is read.
    first_time = True
    # `extract_chunk` can be used to efficiently read chunks of data
    # from a GDAL dataset. It uses a "generator", which means that we
    # can iterate over it using e.g. a for loop
    if verbose:
        print "Start reading %s" % fname
    # The output filaname.
    output_fname = fname.replace(".TIF", "_TOARAD.tif" )
    for ( ds_config, this_X, this_Y, nx_valid, ny_valid, data_in ) in \
        extract_chunk ( fname ):
        if first_time:
            if verbose:
                print "Creating output %s" % output_fname
            # Create output dataset if `first_time` is true
            drv = gdal.GetDriverByName ( "GTiff" )
            dst_ds = drv.Create ( output_fname, ds_config['nx'], \
                ds_config['ny'], 1, gdal.GDT_Float32, options=GDAL_OPTS )
            dst_ds.SetGeoTransform ( ds_config['geoT'] )
            dst_ds.SetProjection ( ds_config['proj'] )
            first_time = False
        
        # Create the output array. Not needed, but we can
        # conveniently set the output type here to float32
        
        radiance = np.zeros_like ( data_in, dtype=np.float32 )
        # DN to radiance conversion if we have a sensible DN
        passer = np.logical_and ( qc_lmin < data_in, data_in < qc_lmax )
        radiance = np.where ( passer, \
             (( lmax - lmin )/(qc_lmax-qc_lmin))*\
             ( (data_in*gain + bias) - qc_lmin) + lmin, \
              -999 )
        # Write the output chunk
        dst_ds.GetRasterBand ( 1 ).WriteArray( radiance, \
                    xoff=this_X, yoff=this_Y)
    # Add some useful metadata to the dataset
    dst_ds.GetRasterBand( 1 ).SetNoDataValue ( -999 )
    dst_ds.GetRasterBand( 1 ).SetMetadata ( {"Band": "%d" % band, \
                                    "Units": "W/(m2*ster*um)", \
                                    "Data": "TOA Radiance" } )
    # Need to do this to flush the dataset to disk
    dst_ds = None
    
def main ():
    """The main function""" 
  

    global options
    global args
    
    metadata_file = get_metadata ( options.input_f )
    ( gain, bias, lmax, lmin, qc_lmax, qc_lmin ) = process_metadata ( metadata_file )
    original_fname = os.path.basename ( options.input_f)
    prefix = original_fname.split("_")[0] 
    fname_prefix = os.path.join ( os.path.dirname ( options.input_f ), prefix )
    for the_band in [1,2,3,4,5,7]:
        input_file = fname_prefix + "_B%d.TIF" % the_band
        convert_to_radiance ( input_file, the_band, \
                    gain[the_band], bias[the_band], \
                    lmax[the_band], lmin[the_band], qc_lmax[the_band], \
                    qc_lmin[the_band], options.verbose )
    

if __name__ == '__main__':

    start_time = time.time()
    parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), \
            usage=globals()['__doc__'] )
    parser.add_option ('-v', '--verbose', action='store_true', \
            default=False, help='verbose output')
    parser.add_option ('-i', '--input', action='store', dest="input_f",\
            type=str, help='Input filename')
        
    (options, args) = parser.parse_args()
    if options.verbose: print time.asctime()
    main()
    if options.verbose: print time.asctime()
    if options.verbose: print 'TOTAL TIME IN MINUTES:',
    if options.verbose: print (time.time() - start_time) / 60.0
