#! /usr/bin/python
# -*- coding: utf8 -*-

"""
Program purpose

Detailed description
"""

__author__= "Nicolai Holzer"
__author_email__ = "first-name dot last-name @ tu-dresden.de"
__date__ ="2013-05-14"
__version__ = "v0.1.1" #MajorVersion(backward_incompatible).MinorVersion(backward_compatible).Patch(Bug_fixes)


#Changelog
#-------------------------------------------------------------------------------
#2013-05-01: v0.1.0 first version
#2013-05-14: v0.1.1 additions: tm1 and slope threshold, etc

#Imported libraries
#-------------------------------------------------------------------------------
#standard libraries
import time

#related libraries
import numpy as np
#print np.__file__, np.__version__

#GDAL
from osgeo import gdal, gdalnumeric
from osgeo.gdalconst import *

#Arcpy
import arcpy
from arcpy import sa #Spatial Analyst
import arcpy.cartography as ca #Cartography
import arcpy.management as dm #Data Management

#ENVI Tools as Arcpy add-on --> MISSING LICENSE
#import envipy
#envipy.Initialize(arcpy) 

#Tkinter GUI, or wxPython
import Tkinter, tkFileDialog 
tkRoot = Tkinter.Tk()
tkRoot.withdraw() #Get rid of python window of Tkinter

#local applications / library specific import

#===============================================================================

#Module default values / constants
#-------------------------------------------------------------------------------
#Initialisation file
INITDIR = 'D:\Data\Temporary data\imagery\Landsat\Mustag Ata\LE7_Mustag_11092000'
INITFILE_RASTER = 'le7_mustag_11092000_mosaic_sub.tif'
INITFILE_DEM = 'srtm_mosaic_51_05_52_05_egm96_utm_subset_mustagata.tif'
OUTPUT_COORDINATESYSTEM = "Coordinate Systems/Projected Coordinate Systems/UTM/WGS 1984/Northern Hemisphere/WGS 1984 UTM Zone 43N.prj" #via ArcGIS, UTM for metric output

#_______________________________________________________________________________


class ProcessingArcPy:
    """Control class. This class is providing all available functions for reading data"""

    def __init__(self):
        """Constructor"""
        
        self.workspace = arcpy.env.workspace = tkFileDialog.askdirectory(initialdir = INITDIR, parent=tkRoot, title = 'Choose workspace (ArcPy)', mustexist = True)
        arcpy.env.scratchWorkspace = tkFileDialog.askdirectory(initialdir = INITDIR+"\scratchworkspace", parent=tkRoot, title = 'Choose scratch workspace (ArcPy)', mustexist = True)
                
        #Arcpy specific
        arcpy.CheckOutExtension("spatial") #Load sa license
                
        arcpy.env.overwriteOutput = True #Enable overwrite existing data
        arcpy.env.outputCoordinateSystem = OUTPUT_COORDINATESYSTEM #"Coordinate Systems/Projected Coordinate Systems/UTM/WGS 1984/Northern Hemisphere/WGS 1984 UTM Zone 43N.prj"
        
        self.pArcPyTools = ArcPyTools() #Init of own tool class
        
        
    #def __del__(self):
        #"""Desctructor"""
        
        
    def glacierMaskProcessingTest(self):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        inFile = tkFileDialog.askopenfilename(defaultextension='TIFF', filetypes=[('ERDAS IMAGINE','*.img'),('TIFF','*.tif')], initialdir=self.workspace, initialfile=INITFILE_RASTER, multiple=False, parent=tkRoot, title='Chose input raster for processing (ArcPy)')
        
        pNDSI = self.pArcPyTools.NDSI(inFile) #Using Normalized-Difference Snow Index (NDSI) 
        
        initOutFileName = "NDSI_arcpy_"+sa.Raster(inFile).name
        
        outRasterFile = tkFileDialog.asksaveasfilename(defaultextension='TIFF', filetypes=[('ERDAS IMAGINE','*.img'),('TIFF','*.tif')], initialdir=self.workspace, initialfile=initOutFileName, parent=tkRoot, title='Save output raster file (ArcPy)')
        pNDSI.save(outRasterFile)
                
        return
        
        
    def glacierMaskProcessing(self):
        """
        Program description:
        Die Schwellwerte für die Ratiobilder sollten individuell angepasst werden, liegen aber in der Regel um 2.0. 
        TM4/TM5 hat weniger Fehler in den Schattenbereichen aber es wird weniger Schutt erfasst als bei TM3/TM5. 
        Daher empfehlen Paul und Kääb (2005) noch einen zusätzlichen Schwellwert im ersten Kanal. 
        Was auch immer am Ende verwendet wird, eine visuelle Kontrolle und manuelle Verbesserung ist am Ende notwendig.
        
        Bei den Schwellwerten geht es um das mögliche Abbrechen von Eismassen und ist nicht als Schwellwert für das Gletschervorkommen zu sehen.
        Dennoch ist 45° zu testen, da bei SRTM die wahren Neigungen durch die Glättung/Auflösung eher höher sind. 60° ist m.E. mit SRTM zu hoch.
  
        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        inFile = tkFileDialog.askopenfilename(defaultextension='TIFF', filetypes=[('ERDAS IMAGINE','*.img'),('TIFF','*.tif')], initialdir=self.workspace, initialfile=INITFILE_RASTER, multiple=False, parent=tkRoot, title='Chose input multispectral raster file for glacier mask processing (ArcPy)')
        
        #-------------------------------------------------------------------------------
        #Image ration
        
        #ratioMethod = 'ndsi' #best with threshold 0.6
        ratioMethod = 'tm3tm5'#best with threshold 2.6 #--> Best method!
        #ratioMethod = 'tm4tm5' #best with threshold 2.0 
        
        print "Start creating ratio image with method '"+str(ratioMethod)+"' ..."
        
        if ratioMethod == 'ndsi':
            pRatio = self.pArcPyTools.NDSI(inFile) #Using Normalized-Difference Snow Index (NDSI) --> See Function 
            thresholdRatio = 0.6 #Thresholding NDSI > threshold = snow/ice ...  0.5 - 0.6, or 0.7 (rocaviteanuetal2008); #best value: 0.6
            
        elif ratioMethod == 'tm3tm5': #Red/SWIR
            pRatio = sa.Divide(sa.Float(sa.Raster(inFile + "\Band_3")), sa.Float(sa.Raster(inFile + "\Band_5")))
            thresholdRatio = 2.6 #tm3tm5: th=2 (rocaviteanuetal2008); CCI: about 1.8; Paul and Andreassen (2009): 2.6; #best value: 2.6
            
        elif ratioMethod == 'tm4tm5': #VNIR/SWIR
            pRatio = sa.Divide(sa.Float(sa.Raster(inFile + "\Band_4")), sa.Float(sa.Raster(inFile + "\Band_5")))
            thresholdRatio = 2.0 #Tim Golletz threshold = 3 (tm4/tm5) cf. Paul 2001; # best value: 2.0
            
        else:
            raise Exception("Error: No ratio approach chosen")
                
        
        outFileName = str(ratioMethod)+"_arcpy_"+sa.Raster(inFile).name
        pRatio.save(outFileName)
        print "Ratio image '"+str(outFileName)+"' with ratio method '"+str(ratioMethod)+"' created."
        
        #-------------------------------------------------------------------------------
        #Threshold on slope: Criteria: 
        #--> Slope > 60° --> no glacier (ICIMOD); 
        #--> Slope <= 24° --> glacier (Bolch et al. 2005, 2006, 2008)
        #--> Alean (1985), cited by Bolch et al. 2011: threshold of 45° for the slope of the detachment zone for cold glaciers and 25° for warm glaciers (cold glacier threashold of Alean (1985): 45; ICIMOD: 60)
                
        thresholdSlope = float(90) #ignoring value: 90
        print "Start creating slope file with threshold '"+str(thresholdSlope)+"' ..."
        
        #Resample input DEM and derived slope to resolution of ratio image: Data Management Tools --> Raster --> Raster Processing: Resample
        inFileNameDem = tkFileDialog.askopenfilename(defaultextension='TIFF', filetypes=[('ERDAS IMAGINE','*.img'),('TIFF','*.tif')], initialdir=self.workspace, initialfile=INITFILE_DEM, multiple=False, parent=tkRoot, title='Chose input elevation raster file for slope determination (ArcPy)')
        pInFileDem = sa.Raster(inFileNameDem)
        outFileNameDem = "res"+str(pRatio.meanCellWidth).replace(".","pt")+"_"+pInFileDem.name
        arcpy.Resample_management(pInFileDem, outFileNameDem, pRatio.meanCellWidth, "CUBIC") #BILINEAR #CUBIC #resample to Input scene size = Landsat 30m
        
        #Spatial Analyst Tools --> Surface: Slope (Aspect)
        pOutSlope = sa.Slope(outFileNameDem, "DEGREE", 1)
        pOutSlope.save("slope_"+outFileNameDem)
        print "Slope file '"+str("slope_"+outFileNameDem)+"' created out of DEM input file '"+str(inFileNameDem)+"' and resampled to a resolution of '"+str(pRatio.meanCellWidth)+"'."
        
        
        #-------------------------------------------------------------------------------
        #Additional threshold on TM1
        pB1 = sa.Raster(inFile + "\Band_1") #Use of additional threshold in TM1 to improve classification  in cast shadow (CCI, Paul2005 --> rocaviteanuetal2008)
        thresholdB1 = float(0) #Paul and Andreassen (2009): TM1 (DNs >59) #ignoring value: 0
        
        #-------------------------------------------------------------------------------
        #Thresholding glacier yes/no
        
        print "Start thresholding ratio image..."
        
        #Spatial Analyst Tools --> Conditional: Con
        pBinRatio = sa.Con(((pRatio > thresholdRatio) & (pOutSlope < thresholdSlope) &  (pB1 > thresholdB1)), 1, 0) #Threshold on ratio
        
        outFileName = "bin"+str(thresholdRatio).replace(".","pt")+"_slope"+str(thresholdSlope).replace(".","pt")+"_1tm"+str(thresholdB1).replace(".","pt")+"_"+str(outFileName)
        pBinRatio.save(outFileName)
        
        print "Binary ratio image '"+str(outFileName)+"' with ratio method '"+str(ratioMethod)+"', ratio threshold '"+str(thresholdRatio)+"', slope threshold '"+str(thresholdSlope)+"' and TM1 threshold '"+str(thresholdB1)+"' created."
        
        
        #-------------------------------------------------------------------------------
        #Raster to Vector
        
        self.pArcPyTools.rasterToVector(outFileName, "median", 10000) #Eliminate areas smaller 0.01 km^2 = 10000 Square Meters (rocaviteanuetal2008); 0.02 km^2 (ICIMOD)
        
        
        #Detect spatial autocorrelation: Spatial Statistic Tools --> Analyzing Patterns: Spatial Autocorrelation (Morans I)
        #--> Not applicable here
        
        
        return
    
    
    def glacierMaskProcessingEnvi(self):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        #--> Missing license, not finished implemented
        
        inFile = tkFileDialog.askopenfilename(defaultextension='TIFF', filetypes=[('ERDAS IMAGINE','*.img'),('TIFF','*.tif')], initialdir=self.workspace, initialfile=INITFILE_RASTER, multiple=False, parent=tkRoot, title='Chose input multispectral raster file for glacier mask processing (ArcPy)')
        pRatio = self.pArcPyTools.NDSI(inFile) 
        
        #Envi Tools in ArcGIS --> Image Processing: Auto-Threshold Difference Raster
        outFileName = "thrNDSI_arcpy_"+sa.Raster(inFile).name
        arcpy.AutoThresholdDiff_envi(pRatio, outFileName, "Otsu's")
        
        #Envi Tools in ArcGIS: ENVI Filter With Convolution (Median)
        
        #Envi Tools in ArcGIS: Classification Raster To Vector
        
        #Envi Tools in ArcGIS: Cleanup Classification Raster (Smoothing + Aggregation)
        
        return
    
    
    def iceDivideMapping(self):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        #Ice Divide Determination by flowdirection and watershed
        
        print "Start creating basin raster image..."
        
        inFileName = tkFileDialog.askopenfilename(defaultextension='TIFF', filetypes=[('ERDAS IMAGINE','*.img'),('TIFF','*.tif')], initialdir=self.workspace, initialfile=INITFILE_DEM, multiple=False, parent=tkRoot, title='Chose input elevation raster file for ice divide mapping (ArcPy)')
        inFile = sa.Raster(inFileName)
        
        #Spatial Analyst Tools --> Hydrology: Flow Direction
        outFlowDirection = sa.FlowDirection(inFile, "NORMAL")#"FORCE"
        
        #Spatial Analyst Tools --> Hydrology: Basin  
        outBasin = sa.Basin(outFlowDirection)
        
        outFileName = "bas_flow_"+inFile.name
        outBasin.save(outFileName)
                
        print "Basin raster image '"+str(outFileName)+"' created."
        
        self.pArcPyTools.rasterToVector(outFileName, "median", 10000) #Eliminate areas smaller 0.01 km^2 = 10000 Square Meters (rocaviteanuetal2008); 0.02 km^2 (ICIMOD)
        
        return
    
    
    def calculateShadow(self):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        #arcpy.CheckOutExtension("3D")
        
        print "Start creating shadow image..."
        
        inFileName = tkFileDialog.askopenfilename(defaultextension='TIFF', filetypes=[('ERDAS IMAGINE','*.img'),('TIFF','*.tif')], initialdir=self.workspace, initialfile=INITFILE_DEM, multiple=False, parent=tkRoot, title='Chose input elevation raster file for processing shadow (ArcPy)')
        inFile = sa.Raster(inFileName)
        
        #Determine Shadow Area: Spatial Analyst Tools --> Solar Radiation: Area Solar Radiation
        #--> Not applicable here
        
        #3D Anaylist Tools --> Raster Surface: HillShade
        azimuth = float(143.1156527) #SUN_AZIMUTH from Landsat
        altitude = float(51.2588465) #SUN_ELEVATION from Landsat
        outFileName = "out_shadow.tif"# Problem length of filename!
        
        arcpy.HillShade_3d(inFile, outFileName, azimuth, altitude, "SHADOWS", 1)
        
        print "Hillshade shadow image '"+str(outFileName)+"' with azimut '"+str(azimuth)+"' and altitude '"+str(altitude)+"' created."
        
        return


#_______________________________________________________________________________

class ArcPyTools:
    """
    Description
    
    Detailed description

    COMMENTS:
    IMPORTANT:
    """

    #def __init__(self):
        #"""Constructor"""
                
                
    #def __del__ (self):
        #"""Destructor"""
    
    
    def NDSI(self, inFile_):
        """
        Program description: 
        Normalized-Difference Snow Index (NDSI) 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        try: 
            b2 = inFile_ + "\Band_2"
            b5 = inFile_ + "\Band_5"
            
            #Arcpy Calculation
            pNum = sa.Float(sa.Raster(b2)-sa.Raster(b5))
            pDenom = sa.Float(sa.Raster(b2) + sa.Raster(b5))
            pNDSI = sa.Divide(pNum, pDenom)
            
        except Exception:
            raise Exception("Error: Error in calculating NDSI")
            
        return pNDSI 
        
        
        
    def rasterToVector(self, inFileName_, filtering_, delArea_):
        """
        Program description: 
        Normalized-Difference Snow Index (NDSI) 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        #-------------------------------------------------------------------------------
        #Raster Filtering
        
        pBinRatio = sa.Raster(inFileName_)
        print "Start filtering raster..."
                
        if filtering_ == 'median': #3x3 Median Filter #--> Best result
            pFiltBinRatio = sa.FocalStatistics(pBinRatio, sa.NbrRectangle(3, 3, "CELL"), "MEDIAN") #Spatial Analyst Tools --> Neighborhood: Focal Statistics 
            
        elif filtering_ == 'closing': #Erosion / Dilatation (Closing) instead of Median
            pShrinkBinRatio = sa.Shrink(pBinRatio, 1, [0])#Spatial Analyst Tools --> Generalization: Shrink ... Expand: One cell
            pFiltBinRatio = sa.Expand(pShrinkBinRatio, 1, [0])# Closing value '0', so no-glacier area, one pixel size
            
        else:
            raise Exception("Error: No filtering approach chosen")
        
        outFileName = str(filtering_)+"_"+pBinRatio.name
        pFiltBinRatio.save(outFileName)
        print "Binary ratio image '"+str(outFileName)+"' with filtering method '"+str(filtering_)+"' created."
        
        #-------------------------------------------------------------------------------
        #Convert Raster to Polygon
        
        outFileName = outFileName.rsplit(".")[0] #Get filename without ".tif" extension
        outFileName = "vec_"+outFileName.rsplit("arcpy")[0]+".shp" #Split filename at code 'arcpy' to shorten filenmame for following methods (methods fail if output filename is too long)
        
        #Conversion Tools --> From Rater: Raster to Polygon
        print "Start converting raster to polygon..."
        arcpy.RasterToPolygon_conversion(pFiltBinRatio, outFileName, "NO_SIMPLIFY", "VALUE") #No simplify of output polygon
        print "Raster converted to polygon '"+str(outFileName)+"'."
        
        #-------------------------------------------------------------------------------
        #Postprocessing
        
        #Data Management Tools --> Generalization: Eliminate Polygon Part --> Eliminates small islands in the polygon and outside
        inFileName = outFileName
        outFileName = "elim_"+inFileName
        
        print "Start eliminating polygon parts..."
        dm.EliminatePolygonPart(inFileName, outFileName, "AREA", delArea_)#, "", "ANY") # not CONTAINED_ONLY #Eliminate areas smaller 0.01 km^2 = 10000 Square Meters
        print "Polygon '"+str(inFileName)+"' eliminated polygon parts of '"+str(delArea_)+"' to '"+str(outFileName)+"'."
        
        
        #Cartography Tools --> Generalization: Aggregate Polygons
        #ca.AggregatePolygons("buildings.shp", "C:/output/output.gdb/aggregated_buildings", 10)
        # Not suitable in this case!
        
        
        #Cartography Tools --> Generalization: Simplify Polygon with Point Remove
        inFileName = outFileName
        outFileName = "simp_"+inFileName
        tolerance = float(pBinRatio.meanCellWidth) #Cell Size: SRTM: 90m, Landsat: 30m
        
        print "Start simplyfing polygon..."
        ca.SimplifyPolygon(inFileName, outFileName, "POINT_REMOVE", tolerance, delArea_, "RESOLVE_ERRORS","KEEP_COLLAPSED_POINTS") #SQR(15^2+15^2) ~ 22m; Eliminate areas smaller 0.01 km^2 = 10000 Square Meters
        print "Polygon '"+str(inFileName)+"' simplified to polygon '"+str(outFileName)+"' with deleted areas of '"+str(delArea_)+"' and maximal tolerance of '"+str(tolerance)+"'."
        
        
        #Cartography Tools --> Generalization: Smooth Polygon with PEAK
        inFileName = outFileName
        outFileName = "smooth_"+inFileName
        tolerance = 2 * float(pBinRatio.meanCellWidth) #2*Cell Size: SRTM: 180m, Landsat: 60m
        
        print "Start smoothing polygon..."
        ca.SmoothPolygon(inFileName, outFileName, "PAEK", tolerance, "FIXED_ENDPOINT", "FLAG_ERRORS")
        print "Polygon '"+str(inFileName)+"' smoothed to polygon '"+str(outFileName)+"' with maximal tolerance of '"+str(tolerance)+"'."
        
        
        return
        
#_______________________________________________________________________________
#_______________________________________________________________________________


class ProcessingGdal:
    """Control class. This class is providing all available functions for reading data"""

    def __init__(self):
        """Constructor"""
        
        self.workspace = tkFileDialog.askdirectory(initialdir = INITDIR, parent=tkRoot, title = 'Choose workspace (GDAL)', mustexist = True)
        
        #Loading GDAL
        gdal.AllRegister() #Register all drivers
        
        self.GdalTools = GdalTools() #Init of own processing class
        
        
    #def __del__(self):
        #"""Desctructor"""
        
        
    def readRasterFile(self):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        inRasterFile = tkFileDialog.askopenfilename(defaultextension='TIFF', filetypes=[('ERDAS IMAGINE','*.img'),('TIFF','*.tif')], initialdir=self.workspace, initialfile=INITFILE_RASTER, multiple=False, parent=tkRoot, title='Chose input raster for processing (GDAL)')
        pInRasterFile = gdal.Open(inRasterFile, GA_ReadOnly) #Open GDAL file
        
        inRasterFileName = inRasterFile.split("/")[-1] #Get filename out of path
        
        return pInRasterFile, inRasterFileName #Return tuple of above
        
        
    def writeRasterFileFromNumpy(self, pOutNumpy_, outFileName_, pInitRasterFile_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
                
        #Write Gdal file from numpy array
        outRasterFileName = tkFileDialog.asksaveasfilename(defaultextension='TIFF', filetypes=[('TIFF','*.tif')], initialdir=self.workspace, initialfile=outFileName_, parent=tkRoot, title='Save output raster file (GDAL)')
        
        #!!! Only for TIFF, one Numpy band and Float data implemented here
        pDriver = gdal.GetDriverByName("GTiff")
        pOutRasterDataset = pDriver.Create(outRasterFileName, pInitRasterFile_.RasterXSize, pInitRasterFile_.RasterYSize, 1, GDT_Float32) #New file, settings from original file
        gdalnumeric.CopyDatasetInfo(pInitRasterFile_, pOutRasterDataset) #Copy metadata information from original file to new file
        gdalnumeric.BandWriteArray(pOutRasterDataset.GetRasterBand(1), pOutNumpy_) # Copy numpy-data to new file
        
        #Close datasets
        pOutRasterDataset = None
        
        return
        
        
    def glacierMaskProcessingTest(self):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        pInFileTuple = self.readRasterFile()
        inFileName = pInFileTuple[1] #Filename
        pInFile = pInFileTuple[0] #GDAL raster object
        
        pNDSI = self.GdalTools.NDSI(pInFile) #Using Normalized-Difference Snow Index (NDSI)
        
        initOutFileName = "NDSI_gdal_"+inFileName
        self.writeRasterFileFromNumpy(pNDSI, initOutFileName, pInFile) #!!! Only for one Numpy band and Float data implemented
        
        return


#_______________________________________________________________________________


class GdalTools:
    """
    Description

    Detailed description

    COMMENTS:
    IMPORTANT:
    """
    
    #def __init__(self):
        #"""Constructor"""
                
        
    #def __del__ (self):
        #"""Destructor"""
    
    
    def NDSI(self, inFile_):
        """
        Program description: 
        Normalized-Difference Snow Index (NDSI)

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        try: 
            
            pB2 = inFile_.GetRasterBand(2)
            pB5 = inFile_.GetRasterBand(5)
            
            #Read the bands of the dataset into numpy arrays  
            pB2numpy = np.float32(gdalnumeric.BandReadAsArray(pB2))
            pB5numpy = np.float32(gdalnumeric.BandReadAsArray(pB5))
            
            #Calculation in numpy
            pNum = np.subtract(pB2numpy, pB5numpy)
            pDenom = np.add(pB2numpy, pB5numpy)
            pNDSI = np.divide(pNum, pDenom)
            
        except Exception:
            raise Exception("Error: Error in calculating NDSI")
            
        return pNDSI
        
        
#_______________________________________________________________________________

def main():
    """
    Main function.

    Detailed description   
    """
    
    #Initialization
    #-------------------------------------------------------------------------------
    startTime = time.time()
    print("_____________________________________________________________________________________________")
    print("Starting program '" + str(__name__) + "' version '" + str(__version__) + "' from '" + str(__date__) + "':")
    
    #Run program
    #-------------------------------------------------------------------------------
    
    approach = 'ArcPy'
    #approach = 'GDAL'
    
    try:
        if approach == 'ArcPy':
            pProcessingArcPy = ProcessingArcPy() #Initialize
            
            pProcessingArcPy.glacierMaskProcessing()
            pProcessingArcPy.iceDivideMapping()
            pProcessingArcPy.calculateShadow()
            
            #pProcessingArcPy.glacierMaskProcessingTest() #Test for NDVI
            #pProcessingArcPy.glacierMaskProcessingEnvi() # NO LICENSE!
        
        elif approach == 'GDAL':
            pProcessingGdal = ProcessingGdal() #Initialize
            pProcessingGdal.glacierMaskProcessingTest() #Test for NDVI
        
    except Exception: #If Exception occurred in this module or all connected sub-modules
        print('Exception error occurred (see below)! ')
        raise

    finally:
        print("Finished. Total processing time [s]: '" + str(time.time() - startTime) + "'.")
        print("_____________________________________________________________________________________________")
      

if __name__ == "__main__":
    main()
