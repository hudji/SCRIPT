# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 08:47:51 2017

@author: hudjimartsu
"""

# standard imports
from grab_meta import grab_meta
from dnppy import core
import math
import os
import arcpy
if arcpy.CheckExtension('Spatial')=='Available':
    arcpy.CheckOutExtension('Spatial')
    arcpy.env.overwriteOutput = True

__all__=['toa_reflectance_8',       # complete       
         'toa_reflectance_457']     # complete


def toa_reflectance_8(band_nums, meta_path, outdir = None):
[docs]    """
    Converts Landsat 8 bands to Top-of-Atmosphere reflectance. To be performed
    on raw Landsat 8 level 1 data. See link below for details
    see here [http://landsat.usgs.gov/Landsat8_Using_Product.php]

    :param band_nums:   A list of desired band numbers such as [3,4,5]
    :param meta_path:   The full filepath to the metadata file for those bands
    :param outdir:      Output directory to save converted files. If left False it will save ouput
                        files in the same directory as input files.

    :return output_filelist:    List of files created by this function
    """

    output_filelist = []

    # enforce the list of band numbers and grab metadata from the MTL file
    band_nums = core.enf_list(band_nums)
    band_nums = map(str, band_nums)
    OLI_bands = ['1','2','3','4','5','6','7','8','9']
    meta_path = os.path.abspath(meta_path)
    meta = grab_meta(meta_path)

    # cycle through each band in the list for calculation, ensuring each is in the list of OLI bands
    for band_num in band_nums:
        if band_num in OLI_bands:

            # scrape data from the given file path and attributes in the MTL file
            band_path = meta_path.replace("MTL.txt","B{0}.tif".format(band_num))
            Qcal = arcpy.Raster(band_path)                        
            Mp   = getattr(meta,"REFLECTANCE_MULT_BAND_{0}".format(band_num)) # multiplicative scaling factor
            Ap   = getattr(meta,"REFLECTANCE_ADD_BAND_{0}".format(band_num))  # additive rescaling factor
            SEA  = getattr(meta,"SUN_ELEVATION")*(math.pi/180)       # sun elevation angle theta_se

            # get rid of the zero values that show as the black background to avoid skewing values
            null_raster = arcpy.sa.SetNull(Qcal, Qcal, "VALUE = 0")

            # calculate top-of-atmosphere reflectance
            TOA_ref = (((null_raster * Mp) + Ap)/(math.sin(SEA)))


            # save the data to the automated name if outdir is given or in the parent folder if not
            if outdir is not None:
                outdir = os.path.abspath(outdir)
                outname = core.create_outname(outdir, band_path, "TOA_Ref", "tif")
            else:
                folder = os.path.split(meta_path)[0]
                outname = core.create_outname(folder, band_path, "TOA_Ref", "tif")
                
            TOA_ref.save(outname)
            output_filelist.append(outname)
            print("Saved output at {0}".format(outname))

        # if listed band is not an OLI sensor band, skip it and print message
        else:
            print("Can only perform reflectance conversion on OLI sensor bands")
            print("Skipping band {0}".format(band_num))

    return output_filelist