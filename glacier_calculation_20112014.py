#! /usr/bin/python
# -*- coding: utf8 -*-

"""
Program purpose

Detailed description
"""

__author__= "Nicolai Holzer"
__author_email__ = "first-name dot last-name @ tu-dresden.de"
__date__ ="2014-11-20"
__version__ = "v0.1.7" #MajorVersion(backward_incompatible).MinorVersion(backward_compatible).Patch(Bug_fixes)


#Changelog
#-------------------------------------------------------------------------------
#2013-10-28: v0.1.0 first version
#2013-12-09: v0.1.1 adaption for error calculation and Pleiades
#2014-01-06: v0.1.2 adaption for SRTM-C band penetration and correction of mean dh DTM offset
#2014-01-08: v0.1.3 New NMAD with slope threshold and section wise ablation area nodata mean calculation
#2014-02-28: v0.1.4 Expansion of program and correction of error in uncertainty calculation for ice density assumption
#2014-05-27: v0.1.5 Use of the program for Mustag Ata glacier mass balance, several adaptations and extensions
#2014-08-06: v0.1.6 SRTM-Penetration for accumulation areas, Glacier mass balance output
#2014-09-26: v0.1.7 Adaption after Reviews for Gurla Mandhata


#Imported libraries
#-------------------------------------------------------------------------------
#standard libraries
import time
import csv #Python integrated API for handling CSV data
import datetime as datetime
import math
import gc

#related libraries

#from pylab import *
import numpy as np
print np.__file__, np.__version__
import scipy.stats as stats
import matplotlib.pyplot as plt

import osgeo.ogr as ogr
import osgeo.osr as osr
import osgeo.gdal as gdal
print gdal.__file__, gdal.__version__

try:
    import arcpy
    from arcpy import sa #Spatial Analyst
except ImportError:
    print "No ArcPy package available, not mandatory for most parts..."


#local applications / library specific import


#===============================================================================

#Module default values / constants
#-------------------------------------------------------------------------------
#Labels for entries in SHP attribute tables
LABEL_POINT_X = 'POINT_X'
LABEL_POINT_Y = 'POINT_Y'
LABEL_POINT_Z = 'grid_code' #"GRID_CODE"
LABEL_DH = 'diff_dem' # "RASTERVALU"

#Input data specific constants
NODATA = -9999#-3.40282e+38 #Value to be used for missing or no data
PIXELSIZE = 30 #[m]

#Constants from Huss (2013): Density assumptions for converting geodetic glacier volume change to mass change
ICE_DENSITY = 850 #[kg/m³] #\citep{huss2013}
ICE_DENSITY_ERROR = 60 #[kg/m³] #\citep{huss2013}
WATER_DENSITY = 999.972 #[kg/m³]

#DTM correction and accuracy constants, remaining dh-offset values of DTM to SRTM to correct
DH_THRESHOLD = abs(100.0) #[m] Eventual +-threshold for valid range of dh-Pixels #Gurla Mandhata +-35m
Z_ABL_INTERVAL = 25 #[m] #Section size in ablation area for mean calculation until ELA 



#---
#Area specific constants
#---

#--> MUZTAG ATA

SPATIAL_REF_ESRI = "WGS_1984_UTM_Zone_43N" # Muztag Ata

#DEM Offsets (due to different statistical definition of "stable terrain" and non-excluding of terrain steeper as 10°)
#["ple1_m_alos1", "ple1_m_alos2", "ple1_m_srtm", "alos1_m_srtm", "alos2_m_srtm", "ple1_m_kh91", "alos1_m_kh91", "alos2_m_kh91", "srtm_m_kh91"]
MEAN_DTM_OFFSET_MUSTAG_PLEALOS1 = float(0.5916) #Invert sign!
MEAN_DTM_OFFSET_MUSTAG_PLEALOS2 = float(0.5901) #Invert sign!
MEAN_DTM_OFFSET_MUSTAG_PLESRTM = float(0.6633) #Invert sign!
MEAN_DTM_OFFSET_MUSTAG_ALOS1SRTM = float(0.4409) #Invert sign!
MEAN_DTM_OFFSET_MUSTAG_ALOS2SRTM = float(0.2664) #Invert sign!
MEAN_DTM_OFFSET_MUSTAG_PLEKH9 = float(0.9429) #Invert sign!
MEAN_DTM_OFFSET_MUSTAG_ALOS1KH9 = float(0.1676) #Invert sign!
MEAN_DTM_OFFSET_MUSTAG_ALOS2KH9 = float(-0.1637) #Invert sign!
MEAN_DTM_OFFSET_MUSTAG_SRTMKH9 = float(0.1997) #Invert sign!


#---
#--> GURLA MANDHATA, HALJI

#SPATIAL_REF_ESRI = "WGS_1984_UTM_Zone_44N" # Gurla Mandhata

#DEM Offsets
#["tdx", "ple", "spot"]
MEAN_DTM_OFFSET_GURLA_TDX = float(0.0324) #Invert sign! 
MEAN_DTM_OFFSET_GURLA_PLE = float(1.1498) #Invert sign! 
MEAN_DTM_OFFSET_GURLA_SPOT = float(0.6059) #Invert sign! 

#_______________________________________________________________________________


class MainInterface:
    """Control class for model 'ModelCsvRead'. This class is providing all available functions for reading data"""
    
    
    def __init__(self):
        """Constructor"""
        self.pModel = ModelMassBalance()
     
    #def __del__(self):
        #"""Desctructor"""
    
    
    
    def testMassBalance(self, glacierName_, ela_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        
        pDateNew = datetime.date(2013,9,30) #datetime.date(2013,10,26) #date of test_DEM
        pDateOld = datetime.date(1999,9,30) #datetime.date(2000,2,16) #SRTM 11 – 22 February 2000
         
        #Values derived before from calculateDtmAccuracy()
        dtmError = 4.83 #Nmad
        
        
        #Read DEM
        #-------------------------------------------------------------------------------
        shpFileName = "G081470E30264N_gurla_halji_extendPLE_srtm_ple1_m_srtm.shp"
        pNpGlacier = self.pModel.readShpNumpy("input/"+str(shpFileName))
        
        
        #dh-Adaption for remaining mean offset of stable terrain to SRTM, Nodata and outlier handling
        #-------------------------------------------------------------------------------
        pNpGlacier = self.pModel.adaptElevationValues(pNpGlacier, MEAN_DTM_OFFSET) #Correct remaining offset of input DTM to SRTM, e.g. MEAN_DTM_OFFSET_GURLA_PLE
        
        
        #Detect and replace outliers
        #-------------------------------------------------------------------------------
        #pNpGlacier = self.pModel.deleteNoDataValues(pNpGlacier, "del") #zero_del
        #pNpGlacier = self.pModel.setValidDhRange(pNpGlacier, "q0pt95_nodata") #dhNull_q0pt683_q0pt95_q0pt683_iqr_dhThreshold___nodata_del
        pNpGlacier = self.pModel.setValidDhRangeAccAbl(pNpGlacier, "q0pt95_nodata", "q0pt95", "", ela_) #Accumulation area, Ablation area: dhNull_q0pt683_q0pt95_iqr_dhThreshold___nodata_del
        
        
        #Adapt for SRTM C-band penetration in snow and ice in case of DEM differencing with SRTM
        #-------------------------------------------------------------------------------
        #srtmPenetrationNegativ = map(lambda x: x * -1, srtmPenetration_[0:2])  #--> consider math sign! --> #If DEM differencing with SRTM where SRTM is of NEWER date!
        pNpGlacier = self.pModel.adaptSrtmPenetration(pNpGlacier, ela_, SRTM_PENETRATION, DEBRIS_BOUNDARY)  #dh-Adaption for considering SRTM C-band radar penetration, e.g. SRTM_PENETRATION
        #Seasonal correction here?
        
        
        #Fill nodata areas
        #-------------------------------------------------------------------------------
        pNpGlacier = self.pModel.fillNoDataAreas(pNpGlacier, "accMean", "ablMean", "", ela_) #accMode_accMean_accNull_ablMean_ablNull_manualDhXXX
        
        
        #Save modified pNpGlacier as Shapefile and calculate mass balance
        #-------------------------------------------------------------------------------
        pNpGlacier = self.pModel.calculateMassBalance(pNpGlacier, pDateOld, pDateNew, dtmError, SRTM_PENETRATION_ERROR, glacierName_)
        self.pModel.writeShpNumpy(pNpGlacier, "output/out_"+str(shpFileName))
        self.pModel.rasterizeShpArcPy("output/out_"+str(shpFileName))
        
        
        return pNpGlacier
        
        
        
        
    def testDtmAccuracy(self):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        shpFileName = "Yala_stableterrain_slope45_point_forNMAD1.shp"#"input_difference_image_of_dem_stable_terrain_as_point.shp"
        
        
        #Read DEM
        #-------------------------------------------------------------------------------
        pNpNoGlacier = self.pModel.readShpNumpy("input/"+str(shpFileName))
        
        
        #Nodata and outlier handling, adapt height to set mean offset of DEMs to zero
        #-------------------------------------------------------------------------------
        pNpNoGlacier = self.pModel.deleteNoDataValues(pNpNoGlacier, "del") #zero_del
        #pNpNoGlacier = self.pModel.adaptElevationValues(pNpNoGlacier, MEAN_DTM_OFFSET) #Correct remaining offset of input DTM to SRTM
        pNpNoGlacier = self.pModel.setValidDhRange(pNpNoGlacier, "iqr_del") #dhNull_q0pt683_q0pt95_q0pt683_iqr_dhThreshold___nodata_del
        
        
        #Save modified pNpGlacier as Shapefile and calculate accuracy
        #-------------------------------------------------------------------------------
        dtmAccuracyNmad = self.pModel.calculateDtmAccuracy(pNpNoGlacier, "input_DEM")
        self.pModel.writeShpNumpy(pNpNoGlacier, "output/out_"+str(shpFileName))
        self.pModel.rasterizeShpArcPy("output/out_"+str(shpFileName))
        
        
        return pNpNoGlacier, dtmAccuracyNmad
        
        
        
        
        
    def testPlots(self, pNpGlacier_, demNameNew_, demNameOld_, name_, ela_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        print("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
        print("Creating plots for '" + str(name_) + "' ...")
        
        
        #Plot histogram of glacier elevation changes (dh) of demNameNew_ minus demNameOld_ DTM
        #-------------------------------------------------------------------------------
        #dhHist, binEdges = np.histogram(pNpGlacier_[:,3], bins = 70, range = (-35.0, 35.0))# calculate histogram from numpy
        
        plt.figure()#Necessary to plot multiple figures with plt.show()
        pNpGlacier = np.copy(pNpGlacier_)
        
        
        count, bins, ignored = plt.hist(pNpGlacier[:,3], bins = 70, range = (np.min(pNpGlacier[:,3]), np.max(pNpGlacier[:,3])), normed = True, cumulative = False, histtype = 'bar', log = False, color = 'b')#, range = (-35.0, 35.0), color = "blue")# the histogram of the data
        plt.axvline(np.mean(pNpGlacier[:,3]), color='k', linestyle='dashed', linewidth=1.5) #draw a vertical line from the bottom to the top of the y axis
        plt.axvline(np.median(pNpGlacier[:,3]), color='g', linestyle='dashed', linewidth=1.5) #
        plt.plot(bins, stats.norm.pdf(bins,loc=np.mean(pNpGlacier[:,3]), scale=np.std(pNpGlacier[:,3])), 'k') #loc defines the mean and scale defines the standard deviation
        #plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (bins - mu)**2 / (2 * sigma**2) ), linewidth=2, color='r')
        
        plt.xlabel('Elevation difference $\Delta$h '+str(demNameNew_)+' minus '+str(demNameOld_)+' [m]') #dh
        plt.ylabel('Normed frequency') #Number of grid cells
        plt.title('$\Delta$h-values '+str(demNameNew_)+' minus '+str(demNameOld_)+' \n for ' + str(name_)) #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
        plt.grid(True)
        plt.legend(['Mean = '+str(round(np.mean(pNpGlacier[:,3]), 2))+'m', 'Median = '+str(round(np.median(pNpGlacier[:,3]), 2))+'m', 'Expected PDF'], loc = 'upper right', bbox_to_anchor = (1.1, 1.0)) #round to 2nd decimal
        
        plt.savefig("output/"+str(name_)+'_histogram_'+str(demNameNew_)+'_m_'+str(demNameOld_))
        plt.clf()
        
        
        #Q-Q Plot for both demNameNew_ minus demNameOld_ by using Scipy probplot
        #-------------------------------------------------------------------------------
        
        plt.figure()
        stats.probplot(pNpGlacier[:,3], dist="norm", fit=False, plot=plt) #Scipy Probability plot of the unscaled quantiles 
        plt.axhline(0, color='k')
        #plt.xlabel('Theoretical Quantiles')  
        #plt.ylabel('Sample Quantiles') #dh
        plt.title('Normal Q-Q Plot '+str(demNameNew_)+' minus '+str(demNameOld_)+' for ' + str(name_)) #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
        plt.grid(True)
        plt.legend([str(demNameNew_)+'-minus-'+str(demNameOld_)+'-$\Delta$h quantiles', str(demNameNew_)+'-minus-'+str(demNameOld_)+'-$\Delta$h linear fit'], loc = 'lower right', bbox_to_anchor=(1.1, 0), numpoints = 1)#, borderaxespad = 0) #only show 1 point in legend
        plt.savefig("output/"+str(name_)+'_scipy_qqplot_'+str(demNameNew_)+'_m_'+str(demNameOld_))#.pdf
        plt.clf()
        
        
        #Q-Q Plot for both demNameNew_ by using manual technique
        #-------------------------------------------------------------------------------
        pProb = np.linspace(0, 1, num=101) #probabilities = Quantile values for x
        
        sampleQuants1 = stats.mstats.mquantiles(pNpGlacier[:,3], prob=pProb)#[0.05, 0.5, 0.95]) #Sample Quantiles1
        
        #Theoretical quants = straigth line
        theoreticalPdf1 = stats.norm.rvs(loc=np.mean(pNpGlacier[:,3]), scale=np.std(pNpGlacier[:,3]), size = pNpGlacier[:,3].shape[0]) #Random variates of theoretical probability density function
        theoreticalQuants1 = stats.mstats.mquantiles(theoreticalPdf1, prob=pProb)#[0.05, 0.5, 0.95]) #Quantiles for theoretical probability density function
        
        plt.figure() 
        plt.plot(theoreticalQuants1, sampleQuants1, 'b-')
        plt.axhline(0, color='k')
        plt.xlabel('Theoretical Quantiles')  
        plt.ylabel('Sample Quantiles') #dh
        plt.title('Normal Q-Q Plot '+str(demNameNew_)+' minus '+str(demNameOld_)+' for ' + str(name_)) #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
        plt.grid(True)
        plt.legend([str(demNameNew_)+'-minus-'+str(demNameOld_)+'-$\Delta$h quantiles'], loc = 'lower right', bbox_to_anchor=(1.1, 0), numpoints = 1)#, borderaxespad = 0) #only show 1 point in legend
        plt.savefig("output/"+str(name_)+'_manual_qqplot_'+str(demNameNew_)+'_m_'+str(demNameOld_))#.pdf
        plt.clf()
        
        
        
        #Plot hypsometry height (Z) vs. glacier elevation change (dh) for demNameNew_ minus demNameOld_
        #-------------------------------------------------------------------------------
        
        pX = np.linspace(np.min(pNpGlacier[:,2]), np.max(pNpGlacier[:,2])) #X-Values for polynom x-axes
        
        #Polynomial curve fitting for hypsometry
        polyCoeffs = np.polyfit(pNpGlacier[:,2], pNpGlacier[:,3], deg=3) #Fit the data with a polynomial and get the polynomal coefficients
        pPoly1 = np.poly1d(polyCoeffs, variable='z') # construct the polynomial as Object
        print ("Polynomal fitting for SRTM elevation vs. elevation differences $\Delta$h of "+str(demNameNew_)+' minus '+str(demNameOld_)+' for ' + str(name_) + ":\n")
        print pPoly1
        
        #Scipy.stats.linregress for linear regression
        slope, intercept, rValue1, pValue, stdErr = stats.linregress(pNpGlacier[:,2], pNpGlacier[:,3])
        y1 = slope * pX + intercept
        print "Linear regression for SRTM elevation vs. elevation differences $\Delta$h of "+str(demNameNew_)+' minus '+str(demNameOld_)+' for ' + str(name_) + ":\n y = " + str(round(slope,4)) + " * x +" + str(round(intercept,4)) + "; r_value, p_value, std_err:", rValue1, pValue, stdErr
        
        #Plot
        plt.figure()
        plt.plot(ela_, 0, 'ko', pX) 
        plt.plot(pNpGlacier[:,2], pNpGlacier[:,3], 'b.', ela_, 0, 'ko', pX, pPoly1(pX), 'k-', pX, y1, 'k--', linewidth=1.5) #x1, y1, x2, y2 #plt.plot(t, t, 'r--', t, t**2, 'bs', t, t**3, 'g^')
        plt.axhline(0, color='k', linestyle='--', linewidth=1) #Zero line
        plt.axvline(ela_, color='k', linestyle='--', linewidth=1) #ELA line
        
        plt.xlabel('SRTM elevation (m a.s.l.)') #Z 
        plt.ylabel('Elevation difference $\Delta$h '+str(demNameNew_)+' minus '+str(demNameOld_)+' [m]') #dh
        plt.xlim(np.min(pNpGlacier[:,2]), np.max(pNpGlacier[:,2])) #plt.axis([5400, 6600, -35, 35]) #plt.ylim(-35,35)
        plt.ylim(np.min(pNpGlacier[:,3]), np.max(pNpGlacier[:,3])) #plt.ylim(-DH_THRESHOLD, DH_THRESHOLD)   #plt.axis([5400, 6400, -35, 35]) 
        plt.title('SRTM elevation vs. elevation differences $\Delta$h \n of '+str(demNameNew_)+' minus '+str(demNameOld_)+' for ' + str(name_)) #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
        plt.grid(True)
        #plt.legend([str(demNameNew_)+'-minus-'+str(demNameOld_)+'-$\Delta$h', 'ELA = '+str(ela_)+'m', 'Polynomal fit (deg=3)', 'Linear regression (r='+str(round(rValue1,2))+')'], loc = 'lower right', bbox_to_anchor=(1.1, 0), numpoints = 1) #only show 1 point in legend
        plt.legend(['ELA = '+str(ela_)+'m'], loc = 'lower right', bbox_to_anchor=(1.1, 0), numpoints = 1)#, borderaxespad = 0) #only show 1 point in legend
        plt.savefig(str("output/"+name_)+'_hypsometry_'+str(demNameNew_)+'_m_'+str(demNameOld_))
        plt.clf()
        
        
        del pNpGlacier
        return
        
        
    #-------
    
    
    def haljiMassBalance(self, glacierName_, ela_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        
        pDateNew = datetime.date(2013,9,30) #datetime.date(2013,10,26) #date of PLE_DEM at Halji
        pDateOld = datetime.date(1999,9,30) #datetime.date(2000,2,16) #SRTM 11 – 22 February 2000
        
        #Values derived before from calculateDtmAccuracy()
        dtmError = 4.83 #Nmad
        
        #Read DEM 
        #-------------------------------------------------------------------------------
        shpFileName = str(glacierName_)+"_gurla_halji_extendPLE_srtm_ple1_m_srtm.shp"
        pNpGlacier = self.pModel.readShpNumpy("input/"+str(shpFileName))
        
        
        #dh-Adaption for remaining mean offset of stable terrain to SRTM, Nodata and outlier handling
        #-------------------------------------------------------------------------------
        pNpGlacier = self.pModel.adaptElevationValues(pNpGlacier, MEAN_DTM_OFFSET_GURLA_PLE) #Correct remaining offset of input DTM to SRTM
        
        
        #Detect and replace outliers
        #-------------------------------------------------------------------------------
        #pNpGlacier = self.pModel.deleteNoDataValues(pNpGlacier, "del") #zero_del
        #pNpGlacier = self.pModel.setValidDhRange(pNpGlacier, "q0pt95_nodata") #dhNull_q0pt683_q0pt95_q0pt683_iqr_dhThreshold___nodata_del
        pNpGlacier = self.pModel.setValidDhRangeAccAbl(pNpGlacier, "q0pt95_nodata", "q0pt95", "", ela_) #Accumulation area, Ablation area: dhNull_q0pt683_q0pt95_iqr_dhThreshold___nodata_del
        
        
        #Adapt for SRTM C-band penetration in snow and ice in case of DEM differencing with SRTM
        #-------------------------------------------------------------------------------
        #srtmPenetrationNegativ = map(lambda x: x * -1, srtmPenetration_[0:2])  #--> consider math sign! --> #If DEM differencing with SRTM where SRTM is of NEWER date!
        pNpGlacier = self.pModel.adaptSrtmPenetration(pNpGlacier, ela_, [-2.3, -1.7, 0.6], 0)  #dh-Adaption for considering SRTM C-band radar penetration, e.g. SRTM_PENETRATION
        
        
        #Fill nodata areas
        #-------------------------------------------------------------------------------
        pNpGlacier = self.pModel.fillNoDataAreas(pNpGlacier, "accMean", "ablMean", "", ela_) #accMode_accMean_accNull_ablMean_ablNull_manualDhXXX
        
        
        #Save modified pNpGlacier as Shapefile and calculate mass balance
        #-------------------------------------------------------------------------------
        pNpGlacier = self.pModel.calculateMassBalance(pNpGlacier, pDateOld, pDateNew, dtmError, 0.6, glacierName_) #SRTM_PENETRATION_ERROR
        self.pModel.writeShpNumpy(pNpGlacier, "output/out_"+str(shpFileName))
        self.pModel.rasterizeShpArcPy("output/out_"+str(shpFileName))
        
        
        return pNpGlacier
    
    
    
    
    #-------
        
        
    def mustagMassBalance(self, glacierName_, demName_, dtmAccuracyNmad_, ela_, srtmPenetration_, debrisBoundary_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        
        #Read DEM
        #-------------------------------------------------------------------------------
        shpFileName = str(glacierName_)+"_mustag_srtm_xyzdh_" + str(demName_) + ".shp"
        pNpGlacier = self.pModel.readShpNumpy("input/"+str(shpFileName))
        
        
        #dh-Adaption for remaining mean offset of stable terrain to SRTM, Nodata and outlier handling
        #-------------------------------------------------------------------------------
        #Values derived before from calculateDtmAccuracy()
        # ["ple1_m_alos1", "ple1_m_alos2", "ple1_m_srtm", "alos1_m_srtm", "alos2_m_srtm", "ple1_m_kh91", "alos1_m_kh91", "alos2_m_kh91", "srtm_m_kh91"]
        
        if demName_ in ['ple1_m_alos1']: #Eliminate outliers 
            pDateOld = datetime.date(2009,9,30) #datetime.date(2009,9,10) #ALOS Mustag Ata
            pDateNew = datetime.date(2013,9,30) #datetime.date(2013,6,20) #PLE Mustag Ata
            manualDh = "manualDhMustagPLEALOS"
            pNpGlacier = self.pModel.adaptElevationValues(pNpGlacier, MEAN_DTM_OFFSET_MUSTAG_PLEALOS1) #Correct remaining offset of input DTM to SRTM
            
        elif demName_ in ['ple1_m_alos2']: #Eliminate outliers 
            pDateOld = datetime.date(2009,9,30) #datetime.date(2009,9,10) #ALOS Mustag Ata
            pDateNew = datetime.date(2013,9,30) #datetime.date(2013,6,20) #PLE Mustag Ata
            manualDh = "manualDhMustagPLEALOS"
            pNpGlacier = self.pModel.adaptElevationValues(pNpGlacier, MEAN_DTM_OFFSET_MUSTAG_PLEALOS2) #Correct remaining offset of input DTM to SRTM
            
        elif demName_ in ['ple1_m_srtm']: #Eliminate outliers 
            pDateOld = datetime.date(1999,9,30) #datetime.date(2000,2,16) #SRTM 11 – 22 February 2000
            pDateNew = datetime.date(2013,9,30) #datetime.date(2013,6,20) #PLE Mustag Ata
            manualDh = "manualDhMustagPLESRTM"
            pNpGlacier = self.pModel.adaptElevationValues(pNpGlacier, MEAN_DTM_OFFSET_MUSTAG_PLESRTM) #Correct remaining offset of input DTM to SRTM
            
        elif demName_ in ['alos1_m_srtm']: #Eliminate outliers 
            pDateOld = datetime.date(1999,9,30) #datetime.date(2000,2,16) #SRTM 11 – 22 February 2000
            pDateNew = datetime.date(2009,9,30) #datetime.date(2009,9,10) #ALOS Mustag Ata
            manualDh = "manualDhMustagALOSSRTM"
            pNpGlacier = self.pModel.adaptElevationValues(pNpGlacier, MEAN_DTM_OFFSET_MUSTAG_ALOS1SRTM) #Correct remaining offset of input DTM to SRTM
            
        elif demName_ in ['alos2_m_srtm']: #Eliminate outliers 
            pDateOld = datetime.date(1999,9,30) #datetime.date(2000,2,16) #SRTM 11 – 22 February 2000
            pDateNew = datetime.date(2009,9,30) #datetime.date(2009,9,10) #ALOS Mustag Ata
            manualDh = "manualDhMustagALOSSRTM"
            pNpGlacier = self.pModel.adaptElevationValues(pNpGlacier, MEAN_DTM_OFFSET_MUSTAG_ALOS2SRTM) #Correct remaining offset of input DTM to SRTM
            
        elif demName_ in ['ple1_m_kh91']: #Eliminate outliers 
            pDateOld = datetime.date(1973,9,30) #datetime.date(1973,8,4) #KH-9 Mustag Ata
            pDateNew = datetime.date(2013,9,30) #datetime.date(2013,6,20) #PLE Mustag Ata
            manualDh = "manualDhMustagPLEKH9"
            pNpGlacier = self.pModel.adaptElevationValues(pNpGlacier, MEAN_DTM_OFFSET_MUSTAG_PLEKH9) #Correct remaining offset of input DTM to SRTM
            
        elif demName_ in ['alos1_m_kh91']: #Eliminate outliers 
            pDateOld = datetime.date(1973,9,30) #datetime.date(1973,8,4) #KH-9 Mustag Ata
            pDateNew = datetime.date(2009,9,30) #datetime.date(2009,9,10) #ALOS Mustag Ata
            manualDh = "manualDhMustagALOSKH9"
            pNpGlacier = self.pModel.adaptElevationValues(pNpGlacier, MEAN_DTM_OFFSET_MUSTAG_ALOS1KH9) #Correct remaining offset of input DTM to SRTM
            
        elif demName_ in ['alos2_m_kh91']: #Eliminate outliers 
            pDateOld = datetime.date(1973,9,30) #datetime.date(1973,8,4) #KH-9 Mustag Ata
            pDateNew = datetime.date(2009,9,30) #datetime.date(2009,9,10) #ALOS Mustag Ata
            manualDh = "manualDhMustagALOSKH9"
            pNpGlacier = self.pModel.adaptElevationValues(pNpGlacier, MEAN_DTM_OFFSET_MUSTAG_ALOS2KH9) #Correct remaining offset of input DTM to SRTM
        
        elif demName_ in ['srtm_m_kh91']: #Eliminate outliers 
            pDateOld = datetime.date(1973,9,30) #datetime.date(1973,8,4) #KH-9 Mustag Ata
            pDateNew = datetime.date(1999,9,30) #datetime.date(2000,2,16) #SRTM 11 – 22 February 2000
            manualDh = "manualDhMustagSRTMKH9"
            pNpGlacier = self.pModel.adaptElevationValues(pNpGlacier, MEAN_DTM_OFFSET_MUSTAG_SRTMKH9) #Correct remaining offset of input DTM to SRTM
            
        else:
            raise Exception ("Error: DEM difference image '" + str(demName_) + "' is not implemented.")
        
        
        #Detect and replace outliers
        #-------------------------------------------------------------------------------
        #pNpGlacier = self.pModel.deleteNoDataValues(pNpGlacier, "del") #zero_del
        #pNpGlacier = self.pModel.setValidDhRange(pNpGlacier, "q0pt95_nodata") #dhNull_q0pt683_q0pt95_q0pt683_iqr_dhThreshold___nodata_del
        
        if glacierName_ in ['allglacier']: #Complete glacier area 
            pNpGlacier = self.pModel.setValidDhRangeAccAbl(pNpGlacier, "dhNull_nodata", "dhThreshold", "", ela_) #Accumulation area, Ablation area: dhNull_q0pt683_q0pt95_q0pt683_iqr_dhThreshold___nodata_del
            #pNpGlacier = self.pModel.replaceGlacierAreas(pNpGlacier, demName_, ["kekesayi"]) #Optional: Replace specific areas by already processed glaciers (read from output folder) - Does not work correctly
        else: #individual glaciers
            if glacierName_ in ['kematulejia', 'kuosikulake']: #Poor quality of difference image in ablation area
                pNpGlacier = self.pModel.setValidDhRangeAccAbl(pNpGlacier, "q0pt683_nodata", "q0pt683", "manualZMustag", ela_) #Accumulation area, Ablation area: dhNull_q0pt683_q0pt95_iqr_dhThreshold___nodata_del
            elif glacierName_ in ['muztagata']: #High quality of difference image in accumulation and ablation area
                pNpGlacier = self.pModel.setValidDhRangeAccAbl(pNpGlacier, "q0pt95_nodata", "q0pt95", "manualZMustag", ela_) #Accumulation area, Ablation area: dhNull_q0pt683_q0pt95_iqr_dhThreshold___nodata_del
            else: #High quality of difference image in ablation area --> default
                pNpGlacier = self.pModel.setValidDhRangeAccAbl(pNpGlacier, "q0pt683_nodata", "q0pt95", "manualZMustag", ela_) #Accumulation area, Ablation area: dhNull_q0pt683_q0pt95_iqr_dhThreshold___nodata_del
        
        
        #Adapt for SRTM C-band penetration in snow and ice in case of DEM differencing with SRTM
        #-------------------------------------------------------------------------------
        if demName_ in ['ple1_m_srtm', 'alos1_m_srtm', 'alos2_m_srtm']: #If DEM differencing with SRTM where SRTM is of OLDER date!
            pNpGlacier = self.pModel.adaptSrtmPenetration(pNpGlacier, ela_, srtmPenetration_, debrisBoundary_[glacierName_])  #dh-Adaption for considering SRTM C-band radar penetration, e.g. SRTM_PENETRATION
        elif demName_ in ['srtm_m_kh91']: #If DEM differencing with SRTM where SRTM is of NEWER date!
            srtmPenetrationNegativ = map(lambda x: x * -1, srtmPenetration_[0:2])  #--> consider math sign!
            pNpGlacier = self.pModel.adaptSrtmPenetration(pNpGlacier, ela_ ,srtmPenetrationNegativ, debrisBoundary_[glacierName_])  #dh-Adaption for considering SRTM C-band radar penetration, e.g. SRTM_PENETRATION
        
        
        #Fill nodata areas
        #-------------------------------------------------------------------------------
        if glacierName_ in ['muztagata']: #In case of high quality accumulation dh pixel values for statistical replacement value
            pNpGlacier = self.pModel.fillNoDataAreas(pNpGlacier, "accMean", "ablMean", manualDh, ela_) #accMode_accMean_accNull_ablMean_ablNull_manualDhXXX
        else: #default
            pNpGlacier = self.pModel.fillNoDataAreas(pNpGlacier, "accNull", "ablMean", manualDh, ela_) #accMode_accMean_accNull_ablMean_ablNull_manualDhXXX
        
        
        #Save modified pNpGlacier as Shapefile and calculate mass balance
        #-------------------------------------------------------------------------------
        if demName_ in ['ple1_m_srtm', 'alos1_m_srtm', 'alos2_m_srtm', 'srtm_m_kh91']: #If DEM differencing with SRTM consider srtmPenetrationError_
            pNpGlacier = self.pModel.calculateMassBalance(pNpGlacier, pDateOld, pDateNew, dtmAccuracyNmad_, abs(srtmPenetration_[2]), glacierName_)
        else: #No radar penetration correction --> Set uncertainty to zero
            pNpGlacier = self.pModel.calculateMassBalance(pNpGlacier, pDateOld, pDateNew, dtmAccuracyNmad_, 0, glacierName_)
        self.pModel.writeShpNumpy(pNpGlacier, "output/out_"+str(shpFileName))
        self.pModel.rasterizeShpArcPy("output/out_"+str(shpFileName)) #self.pModel.rasterizeShpGdal("output/out_"+str(shpFileName))
        
        
        return pNpGlacier
        
        
        
        
        
        
        
    def mustagDtmAccuracy(self, demName_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        shpFileName = "stableterrain_mustag_srtm_xyzdh_" + str(demName_) + ".shp"
        
        
        #Read DEM
        #-------------------------------------------------------------------------------
        pNpNoGlacier = self.pModel.readShpNumpy("input/"+str(shpFileName))
        
        
        #Nodata and outlier handling, adapt height to set mean offset of DEMs to zero
        #-------------------------------------------------------------------------------
        pNpNoGlacier = self.pModel.deleteNoDataValues(pNpNoGlacier, "del") #zero_del
        
        
        ["ple1_m_alos1", "ple1_m_alos2", "ple1_m_srtm", "alos1_m_srtm", "alos2_m_srtm", "ple1_m_kh91", "alos1_m_kh91", "alos2_m_kh91", "srtm_m_kh91"]
        if demName_ in ['ple1_m_alos1']: #Eliminate outliers 
            pNpNoGlacier = self.pModel.adaptElevationValues(pNpNoGlacier, MEAN_DTM_OFFSET_MUSTAG_PLEALOS1) #Correct remaining offset of input DTM to SRTM
        elif demName_ in ['ple1_m_alos2']: #Eliminate outliers 
            pNpNoGlacier = self.pModel.adaptElevationValues(pNpNoGlacier, MEAN_DTM_OFFSET_MUSTAG_PLEALOS2) #Correct remaining offset of input DTM to SRTM
        elif demName_ in ['ple1_m_srtm']: #Eliminate outliers 
            pNpNoGlacier = self.pModel.adaptElevationValues(pNpNoGlacier, MEAN_DTM_OFFSET_MUSTAG_PLESRTM) #Correct remaining offset of input DTM to SRTM
        elif demName_ in ['alos1_m_srtm']: #Eliminate outliers 
            pNpNoGlacier = self.pModel.adaptElevationValues(pNpNoGlacier, MEAN_DTM_OFFSET_MUSTAG_ALOS1SRTM) #Correct remaining offset of input DTM to SRTM
        elif demName_ in ['alos2_m_srtm']: #Eliminate outliers 
            pNpNoGlacier = self.pModel.adaptElevationValues(pNpNoGlacier, MEAN_DTM_OFFSET_MUSTAG_ALOS2SRTM) #Correct remaining offset of input DTM to SRTM
        elif demName_ in ['ple1_m_kh91']: #Eliminate outliers 
            pNpNoGlacier = self.pModel.adaptElevationValues(pNpNoGlacier, MEAN_DTM_OFFSET_MUSTAG_PLEKH9) #Correct remaining offset of input DTM to SRTM
        elif demName_ in ['alos1_m_kh91']: #Eliminate outliers 
            pNpNoGlacier = self.pModel.adaptElevationValues(pNpNoGlacier, MEAN_DTM_OFFSET_MUSTAG_ALOS1KH9) #Correct remaining offset of input DTM to SRTM
        elif demName_ in ['alos2_m_kh91']: #Eliminate outliers 
            pNpNoGlacier = self.pModel.adaptElevationValues(pNpNoGlacier, MEAN_DTM_OFFSET_MUSTAG_ALOS2KH9) #Correct remaining offset of input DTM to SRTM
        elif demName_ in ['srtm_m_kh91']: #Eliminate outliers 
            pNpNoGlacier = self.pModel.adaptElevationValues(pNpNoGlacier, MEAN_DTM_OFFSET_MUSTAG_SRTMKH9) #Correct remaining offset of input DTM to SRTM
        
        pNpNoGlacier = self.pModel.setValidDhRange(pNpNoGlacier, "iqr_del") #dhNull_q0pt683_q0pt95_q0pt683_iqr_dhThreshold___nodata_del
        
        
        #Save modified pNpGlacier as Shapefile and calculate accuracy
        #-------------------------------------------------------------------------------
        dtmAccuracyNmad = self.pModel.calculateDtmAccuracy(pNpNoGlacier, demName_)
        self.pModel.writeShpNumpy(pNpNoGlacier, "output/out_"+str(shpFileName))
        self.pModel.rasterizeShpArcPy("output/out_"+str(shpFileName))
        
        
        return pNpNoGlacier, dtmAccuracyNmad
        
        
        
        
    #---------
        
        
    def gurlaMassBalance(self, glacierName_, demName_, dtmAccuracyNmad_, ela_, srtmPenetration_, debrisBoundary_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        
        
        #Read DEM
        #-------------------------------------------------------------------------------
        readMethod = 'shp' #'csv'
        
        if readMethod == 'shp':
            shpFileName = str(glacierName_)+"_xyzdh_gurlaglacier_20001013_09122013_" + str(demName_) + ".shp"
            pNpGlacier = self.pModel.readShpNumpy("input/"+str(shpFileName))
        #elif readMethod == 'csv':
        #    csvFileName = str(glacierName_)+"_xyzdh_gurlaglacier_20001013_09122013_" + str(demName_) + ".csv"
        #    pDocCsvNumpy = self.pModel.readCsvNumpy("input/"++str(csvFileName))
        #    pNpGlacier = self.pModel.adaptCsvNumpy(pDocCsvNumpy)
        #else:
        #    raise Exception ("Error Please select read method 'csv' or 'shp'.")
        
        
        #dh-Adaption for remaining mean offset of stable terrain to SRTM, Nodata and outlier handling
        #-------------------------------------------------------------------------------
        #Values derived before from calculateDtmAccuracy()
        # ["tdx", "ple", "spot"]
        
        
        pDateOld = datetime.date(1999,9,30) #datetime.date(2000,2,16) #SRTM 11 – 22 February 2000
        manualDh = "manualDhGurla" 
        
        if demName_ in ['tdx']: #Eliminate outliers 
            pDateNew = datetime.date(2012,9,30) #datetime.date(2012,04,26) #TDX Gurla Mandhata
            pNpGlacier = self.pModel.adaptElevationValues(pNpGlacier, MEAN_DTM_OFFSET_GURLA_TDX) #Correct remaining offset of input DTM to SRTM
        elif demName_ in ['ple']: #Eliminate outliers 
            pDateNew = datetime.date(2013,9,30) #datetime.date(2013,10,18) #Pleiades Gurla Mandhata: Primarily 18 and partly 26 October 2013 (image of 12 October not used)
            pNpGlacier = self.pModel.adaptElevationValues(pNpGlacier, MEAN_DTM_OFFSET_GURLA_PLE) #Correct remaining offset of input DTM to SRTM
        elif demName_ in ['spot']: #Eliminate outliers 
            pDateNew = datetime.date(2010,9,30) #datetime.date(2010,11,04) #SPOT Gurla Mandhata
            pNpGlacier = self.pModel.adaptElevationValues(pNpGlacier, MEAN_DTM_OFFSET_GURLA_SPOT) #Correct remaining offset of input DTM to SRTM
        else:
            raise Exception ("Error: DEM difference image '" + str(demName_) + "' is not implemented.")
        
        
        
        #Detect and replace outliers
        #-------------------------------------------------------------------------------
        if demName_ in ['spot']:
            pNpGlacier = self.pModel.deleteNoDataValues(pNpGlacier, "del") #zero_del
            pNpGlacier = self.pModel.setValidDhRange(pNpGlacier, "dhThreshold_del") #dhNull_q0pt683_q0pt95_q0pt683_iqr_dhThreshold___nodata_del
        else:
            #pNpGlacier = self.pModel.setValidDhRange(pNpGlacier, "dhThreshold_nodata") #dhNull_q0pt683_q0pt95_q0pt683_iqr_dhThreshold___nodata_del
            pNpGlacier = self.pModel.setValidDhRangeAccAbl(pNpGlacier, "q0pt95_nodata", "q0pt95", manualDh, ela_) #Accumulation area, Ablation area: dhNull_q0pt683_q0pt95_iqr_dhThreshold___nodata_del
        
        
        #Adapt for SRTM C-band penetration in snow and ice in case of DEM differencing with SRTM
        #-------------------------------------------------------------------------------
        #srtmPenetrationNegativ = map(lambda x: x * -1, srtmPenetration_[0:2])  #--> consider math sign! --> If DEM differencing with SRTM where SRTM is of NEWER date!
        pNpGlacier = self.pModel.adaptSrtmPenetration(pNpGlacier, ela_, srtmPenetration_, debrisBoundary_[glacierName_])  #dh-Adaption for considering SRTM C-band radar penetration, e.g. SRTM_PENETRATION
        
        #Seasonal correction for TDX here?
        
        
        #Fill nodata areas
        #-------------------------------------------------------------------------------
        if demName_ in ['tdx', 'ple']:
            #pNpGlacier = self.pModel.fillNoDataAreas(pNpGlacier, "accNull", "ablMean", "manualDh", ela_) #accMode_accMean_accNull_ablMean_ablNull_manualDhXXX
            pNpGlacier = self.pModel.fillNoDataAreas(pNpGlacier, "accMean", "ablMean", manualDh, ela_) #accMode_accMean_accNull_ablMean_ablNull_manualDhXXX
        
        
        #Save modified pNpGlacier as Shapefile and calculate mass balance
        #-------------------------------------------------------------------------------
        pNpGlacier = self.pModel.calculateMassBalance(pNpGlacier, pDateOld, pDateNew, dtmAccuracyNmad_, abs(srtmPenetration_[2]), glacierName_) #If DEM differencing with SRTM consider srtmPenetrationError_
        self.pModel.writeShpNumpy(pNpGlacier, "output/out_"+str(shpFileName))
        self.pModel.rasterizeShpArcPy("output/out_"+str(shpFileName)) #self.pModel.rasterizeShpGdal("output/out_"+str(shpFileName))
        
        
        return pNpGlacier
        
        
        
        
        
    def gurlaDtmAccuracy(self, demName_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        
        shpFileName = "Gurla_extendSmallNew_srtm_xyzdh_09122013_noglacier_slope35_" + str(demName_) + ".shp"
        
        
        #Read DEM
        #-------------------------------------------------------------------------------
        pNpNoGlacier = self.pModel.readShpNumpy("input/"+str(shpFileName))
        
        
        #Nodata and outlier handling, adapt height to set mean offset of DEMs to zero
        #-------------------------------------------------------------------------------
        pNpNoGlacier = self.pModel.deleteNoDataValues(pNpNoGlacier, "del") #zero_del
        
        if demName_ in ['tdx']: #Eliminate outliers 
            pNpNoGlacier = self.pModel.adaptElevationValues(pNpNoGlacier, MEAN_DTM_OFFSET_GURLA_TDX) #Correct remaining offset of input DTM to SRTM
        elif demName_ in ['ple']:
            pNpNoGlacier = self.pModel.adaptElevationValues(pNpNoGlacier, MEAN_DTM_OFFSET_GURLA_PLE) #Correct remaining offset of input DTM to SRTM
        elif demName_ in ['spot']:
            pNpNoGlacier = self.pModel.adaptElevationValues(pNpNoGlacier, MEAN_DTM_OFFSET_GURLA_SPOT) #Correct remaining offset of input DTM to SRTM
        
        pNpNoGlacier = self.pModel.setValidDhRange(pNpNoGlacier, "iqr_del") #dhNull_q0pt683_q0pt95_q0pt683_iqr_dhThreshold___nodata_del
        
        
        #Save modified pNpGlacier as Shapefile and calculate accuracy
        #-------------------------------------------------------------------------------
        dtmAccuracyNmad = self.pModel.calculateDtmAccuracy(pNpNoGlacier, demName_)
        self.pModel.writeShpNumpy(pNpNoGlacier, "output/out_"+str(shpFileName))
        self.pModel.rasterizeShpArcPy("output/out_"+str(shpFileName))
        
        
        return pNpNoGlacier, dtmAccuracyNmad
        
        
        
        
        
        
        
        
    def gurlaPlots(self, pNpGlacier1_, pNpGlacier2_, demName1_, demName2_, name_, ela_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        print("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
        print("Creating plots for '" + str(name_) + "' ...")
        
        
        #Plot histogram of glacier elevation changes (dh) of demName1_ DTM
        #-------------------------------------------------------------------------------
        #dhHist, binEdges = np.histogram(pNpGlacier_[:,3], bins = 70, range = (-35.0, 35.0))# calculate histogram from numpy
        
        plt.figure()#Necessary to plot multiple figures with plt.show()
        pNpGlacier = np.copy(pNpGlacier1_)
        
        count, bins, ignored = plt.hist(pNpGlacier[:,3], bins = 70, range = (-15, 15), normed = True, cumulative = False, histtype = 'bar', log = False, color = 'b')# range = (np.min(pNpGlacier[:,3]), np.max(pNpGlacier[:,3])) # the histogram of the data
        plt.axvline(np.mean(pNpGlacier[:,3]), color='k', linestyle='dashed', linewidth=1.5) #draw a vertical line from the bottom to the top of the y axis
        plt.axvline(np.median(pNpGlacier[:,3]), color='g', linestyle='dashed', linewidth=1.5) #
        plt.plot(bins, stats.norm.pdf(bins,loc=np.mean(pNpGlacier[:,3]), scale=np.std(pNpGlacier[:,3])), 'k') #loc defines the mean and scale defines the standard deviation
        #plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (bins - mu)**2 / (2 * sigma**2) ), linewidth=2, color='r')
        
        plt.xlabel('Elevation difference $\Delta$h '+str(demName1_)+' minus SRTM [m]') #dh
        plt.ylabel('Normed frequency') #Number of grid cells
        plt.title('$\Delta$h-values '+str(demName1_)+' minus SRTM \n for ' + str(name_)) #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
        plt.grid(True)
        plt.legend(['Mean = '+str(round(np.mean(pNpGlacier[:,3]), 2))+'m', 'Median = '+str(round(np.median(pNpGlacier[:,3]), 2))+'m', 'Expected PDF'], loc = 'upper right', bbox_to_anchor = (1.1, 1.0)) #round to 2nd decimal
        
        plt.savefig("output/"+str(name_)+'_histogram_'+str(demName1_))
        plt.clf()
        
        
        #Plot histogram of glacier elevation changes (dh) of demName2_ DTM
        #-------------------------------------------------------------------------------
        plt.figure()#Necessary to plot multiple figures with plt.show()
        pNpGlacier = np.copy(pNpGlacier2_)
        
        count, bins, ignored = plt.hist(pNpGlacier[:,3], bins = 70, range = (-15, 15), normed = True, cumulative = False, histtype = 'bar', log = False, color = 'r')# range = (np.min(pNpGlacier[:,3]), np.max(pNpGlacier[:,3])) # the histogram of the data#the histogram of the data
        plt.axvline(np.mean(pNpGlacier[:,3]), color='k', linestyle='dashed', linewidth=1.5) #draw a vertical line from the bottom to the top of the y axis 
        plt.axvline(np.median(pNpGlacier[:,3]), color='g', linestyle='dashed', linewidth=1.5) #
        plt.plot(bins, stats.norm.pdf(bins,loc=np.mean(pNpGlacier[:,3]), scale=np.std(pNpGlacier[:,3])), 'k') #loc defines the mean and scale defines the standard deviation
        
        plt.xlabel('Elevation differences $\Delta$h '+str(demName2_)+' minus SRTM [m]') #dh
        plt.ylabel('Normed frequency') #Number of grid cells
        plt.title('$\Delta$h-values ' + str(demName2_) + ' minus SRTM \n for ' + str(name_)) 
        plt.grid(True)
        plt.legend(['Mean = '+str(round(np.mean(pNpGlacier[:,3]), 2))+'m', 'Median = '+str(round(np.median(pNpGlacier[:,3]), 2))+'m', 'Expected PDF'], loc = 'upper right', bbox_to_anchor = (1.1, 1.0)) #round to 2nd decimal
        
        plt.savefig(str("output/"+str(name_)+'_histogram_'+str(demName2_)))
        plt.clf()
        
        
        #Q-Q Plot for both demName1_ and demName2_ by using Scipy probplot
        #-------------------------------------------------------------------------------
        pNpGlacier1 = np.copy(pNpGlacier1_)
        pNpGlacier2 = np.copy(pNpGlacier2_)
        
        plt.figure()
        stats.probplot(pNpGlacier1[:,3], dist="norm", fit=False, plot=plt) #Scipy Probability plot of the unscaled quantiles 
        stats.probplot(pNpGlacier2[:,3], dist="norm", fit=False, plot=plt) #Scipy Probability plot of the unscaled quantiles 
        plt.axhline(0, color='k')
        #plt.xlabel('Theoretical Quantiles')  
        #plt.ylabel('Sample Quantiles') #dh
        plt.title('Normal Q-Q Plot for ' + str(name_)) #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
        plt.grid(True)
        plt.legend([str(demName1_)+'-$\Delta$h quantiles', str(demName1_)+'-$\Delta$h linear fit', str(demName2_) + '-$\Delta$h quantiles', str(demName2_) + '-$\Delta$h linear fit'], loc = 'lower right', bbox_to_anchor=(1.1, 0), numpoints = 1)#, borderaxespad = 0) #only show 1 point in legend
        plt.savefig("output/"+str(name_)+'_scipy_qqplot_'+str(demName1_)+'_'+str(demName2_))#.pdf
        plt.clf()
        
        
        #Q-Q Plot for both demName1_ and demName2_ by using manual technique
        #-------------------------------------------------------------------------------
        pProb = np.linspace(0, 1, num=101) #probabilities = Quantile values for x
        
        sampleQuants1 = stats.mstats.mquantiles(pNpGlacier1[:,3], prob=pProb)#[0.05, 0.5, 0.95]) #Sample Quantiles1
        sampleQuants2 = stats.mstats.mquantiles(pNpGlacier2[:,3], prob=pProb)#[0.05, 0.5, 0.95]) #Sample Quantiles2
        
        #Theoretical quants = straigth line
        theoreticalPdf1 = stats.norm.rvs(loc=np.mean(pNpGlacier1[:,3]), scale=np.std(pNpGlacier1[:,3]), size = pNpGlacier1[:,3].shape[0]) #Random variates of theoretical probability density function
        theoreticalQuants1 = stats.mstats.mquantiles(theoreticalPdf1, prob=pProb)#[0.05, 0.5, 0.95]) #Quantiles for theoretical probability density function
        
        theoreticalPdf2 = stats.norm.rvs(loc=np.mean(pNpGlacier2[:,3]), scale=np.std(pNpGlacier1[:,3]), size = pNpGlacier1[:,3].shape[0]) #Random variates of theoretical probability density function
        theoreticalQuants2 = stats.mstats.mquantiles(theoreticalPdf2, prob=pProb)#[0.05, 0.5, 0.95]) #Quantiles for theoretical probability density function
        
        plt.figure() 
        plt.plot(theoreticalQuants1, sampleQuants1, 'b-', theoreticalQuants2, sampleQuants2, 'r-')
        plt.axhline(0, color='k')
        plt.xlabel('Theoretical Quantiles')  
        plt.ylabel('Sample Quantiles') #dh
        plt.title('Normal Q-Q Plot for ' + str(name_)) #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
        plt.grid(True)
        plt.legend([str(demName1_)+'-$\Delta$h quantiles', str(demName2_) +'-$\Delta$h quantiles'], loc = 'lower right', bbox_to_anchor=(1.1, 0), numpoints = 1)#, borderaxespad = 0) #only show 1 point in legend
        plt.savefig("output/"+str(name_)+'_manual_qqplot_'+str(demName1_)+'_'+str(demName2_))#.pdf
        plt.clf()
        
        
        #Plot hypsometry height (Z) vs. glacier elevation change (dh) for demName1_
        #-------------------------------------------------------------------------------
        pNpGlacier = np.copy(pNpGlacier1_)
        #pX = np.linspace(np.min(pNpGlacier[:,2]), np.max(pNpGlacier[:,2])) #X-Values for polynom x-axes
        pX = np.linspace(np.min(pNpGlacier[:,2]), np.max(pNpGlacier[:,2])) #X-Values for polynom x-axes
        
        #Polynomial curve fitting for hypsometry
        polyCoeffs = np.polyfit(pNpGlacier[:,2], pNpGlacier[:,3], deg=3) #Fit the data with a polynomial and get the polynomal coefficients
        pPoly1 = np.poly1d(polyCoeffs, variable='z') # construct the polynomial as Object
        print ("Polynomal fitting for SRTM elevation vs. elevation differences $\Delta$h of "+str(demName1_)+' minus SRTM for ' + str(name_) + "':\n")
        print pPoly1
        
        #Scipy.stats.linregress for linear regression
        slope, intercept, rValue1, pValue, stdErr = stats.linregress(pNpGlacier[:,2], pNpGlacier[:,3])
        y1 = slope * pX + intercept
        print "Linear regression for SRTM elevation vs. elevation differences $\Delta$h of "+str(demName1_)+' minus SRTM for ' + str(name_) + "':\n y = " + str(round(slope,4)) + " * x +" + str(round(intercept,4)) + "; r_value, p_value, std_err:", rValue1, pValue, stdErr
        
        #Plot
        plt.figure()
        #plt.axhline(0, color='k', linestyle='solid', linewidth=1) #Zero line
        #plt.axvline(ela_, color='k', linestyle='solid', linewidth=1) #ELA line
        plt.plot(pNpGlacier[:,2], pNpGlacier[:,3], 'b.', ela_, 0, 'ko', pX, pPoly1(pX), 'k-', pX, y1, 'k--', linewidth = 1.5) #x1, y1, x2, y2 #plt.plot(t, t, 'r--', t, t**2, 'bs', t, t**3, 'g^')
        plt.xlabel('SRTM elevation (m a.s.l.)') #Z 
        plt.ylabel('Elevation difference $\Delta$h '+str(demName1_)+' to SRTM [m]') #dh 
        plt.xlim(np.min(pNpGlacier[:,2]), np.max(pNpGlacier[:,2]))
        #plt.axis([5400, 6600, -35, 35]) #plt.ylim(-35,35)
        plt.title('SRTM elevation vs. elevation differences $\Delta$h \n of '+str(demName1_)+' minus SRTM for ' + str(name_)) #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
        plt.grid(True)
        plt.legend([str(demName1_)+'-$\Delta$h', 'ELA = '+str(ela_)+'m', 'Polynomal fit (deg=3)', 'Linear regression (r='+str(round(rValue1,2))+')'], loc = 'lower right', bbox_to_anchor=(1.1, 0), numpoints = 1) #only show 1 point in legend
        plt.savefig(str("output/"+name_)+'_hypsometry_'+str(demName1_))
        plt.clf()
        
        
        #Plot hypsometry height (Z) vs. glacier elevation change (dh) for demName2_
        #-------------------------------------------------------------------------------
        pNpGlacier = np.copy(pNpGlacier2_)
        #pX = np.linspace(np.min(pNpGlacier[:,2]), np.max(pNpGlacier[:,2])) #X-Values for polynom x-axes
        pX = np.linspace(np.min(pNpGlacier[:,2]), np.max(pNpGlacier[:,2])) #X-Values for polynom x-axes
        
        #Polynomial curve fitting for hypsometry
        polyCoeffs = np.polyfit(pNpGlacier[:,2], pNpGlacier[:,3], deg=3) #Fit the data with a polynomial and get the polynomal coefficients
        pPoly2 = np.poly1d(polyCoeffs, variable='z') # construct the polynomial as Object
        print ("Polynomal fitting for SRTM elevation vs. elevation differences $\Delta$h of "+str(demName2_)+' minus SRTM for ' + str(name_) + "':\n ")
        print pPoly2
        
        #Scipy.stats.linregress for linear regression
        slope, intercept, rValue2, pValue, stdErr = stats.linregress(pNpGlacier[:,2], pNpGlacier[:,3])
        y2 = slope * pX + intercept
        print "Linear regression for SRTM elevation vs. elevation differences $\Delta$h of "+str(demName2_)+' minus SRTM for ' + str(name_) + "': \n y = " + str(round(slope,4)) + " * x +" + str(round(intercept,4)) + "; r_value, p_value, std_err:", rValue2, pValue, stdErr
        
        #Plot
        plt.figure() 
        #plt.axhline(0, color='k', linestyle='solid', linewidth=1) #Zero line
        #plt.axvline(ela_, color='k', linestyle='solid', linewidth=1) #ELA line
        plt.plot(pNpGlacier[:,2], pNpGlacier[:,3], 'r.', ela_, 0, 'ko', pX, pPoly2(pX), 'k-', pX, y2, 'k--', linewidth = 1.5) #x1, y1, x2, y2 #plt.plot(t, t, 'r--', t, t**2, 'bs', t, t**3, 'g^')
        plt.xlabel('SRTM elevation (m a.s.l.)') #Z 
        plt.ylabel('Elevation difference $\Delta$h '+str(demName2_)+' to SRTM [m]') #dh
        plt.xlim(np.min(pNpGlacier[:,2]), np.max(pNpGlacier[:,2]))
        #plt.axis([5400, 6600, -35, 35]) #plt.ylim(-35,35)
        plt.title('SRTM elevation vs. elevation differences $\Delta$h \n of '+str(demName2_)+' minus SRTM for ' + str(name_)) #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
        plt.grid(True)
        plt.legend([str(demName2_)+'-$\Delta$h', 'ELA = '+str(ela_)+'m', 'Polynomal fit (deg=3)', 'Linear regression (r='+str(round(rValue2,2))+')'], loc = 'lower right', bbox_to_anchor=(1.1, 0), numpoints = 1) #only show 1 point in legend
        plt.savefig("output/"+str(name_)+'_hypsometry_'+str(demName2_))
        plt.clf()
        
        
        #Plot hypsometry height (Z) vs. glacier elevation change (dh) for demName1_ and demName2_
        #-------------------------------------------------------------------------------
        pNpGlacier1 = np.copy(pNpGlacier1_)
        pNpGlacier2 = np.copy(pNpGlacier2_)
        
        #Plot
        plt.figure() 
        plt.plot(ela_, 0, 'ko', pX)
        plt.plot(pNpGlacier1[:,2], pNpGlacier1[:,3], 'b.') 
        plt.plot(pNpGlacier2[:,2], pNpGlacier2[:,3], 'ro', mec = 'r', mfc='None', alpha = 0.5)
        plt.plot(ela_, 0, 'ko', pX, pPoly1(pX), 'k-', pX, y1, 'k--', pX, pPoly2(pX), 'g-', pX, y2, 'g--', linewidth = 1.5)
        plt.axhline(0, color='k', linestyle='--', linewidth=1) #Zero line
        plt.axvline(ela_, color='k', linestyle='--', linewidth=1) #ELA line
        
        
        plt.xlabel('SRTM elevation (m a.s.l.)') #Z 
        plt.ylabel('Elevation difference $\Delta$h to SRTM [m]') #dh
        plt.xlim(np.min(pNpGlacier2[:,2]), np.max(pNpGlacier2[:,2])) #plt.axis([5400, 6400, -35, 35]) 
        #plt.xlim(np.min(pNpGlacier2[:,2]), 6400)# 
        plt.ylim(-40, 40) #plt.axis([5400, 6400, -35, 35]) 
        plt.title('SRTM elevation vs. elevation differences $\Delta$h \n of '+str(demName2_)+' and '+str(demName1_)+' minus SRTM for ' + str(name_)) #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
        plt.grid(True)
        #plt.legend([str(demName1_)+'-$\Delta$h', str(demName2_)+'-$\Delta$h', 'ELA = '+str(ela_)+'m', 'Polynomal fit (deg=3) '+str(demName1_), 'Linear regression (r='+str(round(rValue1,2))+') '+str(demName1_), 'Polynomal fit (deg=3) '+str(demName2_), 'Linear regression (r='+str(round(rValue2,2))+') '+str(demName2_)], loc = 'lower right', bbox_to_anchor=(1.1, 0), numpoints = 1, prop={'size':10})#, borderaxespad = 0) #only show 1 point in legend
        plt.legend(['ELA = '+str(ela_)+'m'], loc = 'lower right', bbox_to_anchor=(1.1, 0), numpoints = 1)#, borderaxespad = 0) #only show 1 point in legend
        plt.savefig("output/"+str(name_)+'_hypsometry_'+str(demName1_)+'_'+str(demName2_))#.pdf
        plt.clf()
        
        
        
        #Make scatterplot of elevation differences from demName1_ vs demName2_
        #-------------------------------------------------------------------------------
        pNpGlacier1 = np.copy(pNpGlacier1_)
        pNpGlacier2 = np.copy(pNpGlacier2_)
        
        
        if pNpGlacier1[:,3].shape == pNpGlacier2[:,3].shape: #Arrays must have same shape (number) for dh-values, does only work for glaciers, not stable terrain!
            
            #Scipy.stats.linregress for linear regression
            pX3 = np.linspace(np.min(pNpGlacier1[:,3]), np.max(pNpGlacier1[:,3])) #X-Values for polynom x-axes
            slope, intercept, rValue3, pValue, stdErr = stats.linregress(pNpGlacier1[:,3], pNpGlacier2[:,3])
            y3 = slope * pX3 + intercept
            print "Linear regression for elevation differences $\Delta$h of "+str(demName1_)+' vs. '+str(demName2_)+' for ' + str(name_)+ ":\n y = " + str(round(slope,4)) + " * x +" + str(round(intercept,4)) + "; r_value, p_value, std_err:", rValue3, pValue, stdErr
            
            #Plot
            plt.figure() 
            #matplotlib.pyplot.scatter(x, y, s=20, c=u'b', marker=u'o', cmap=None, norm=None, vmin=None, vmax=None, alpha=None, linewidths=None, verts=None, hold=None, **kwargs)
            plt.scatter(pNpGlacier1[:,3], pNpGlacier2[:,3],  s=20, marker='+', edgecolor='b')# c=u'b', marker = 'o')
            plt.plot(pX3, y3, 'k-', linewidth=1.5) 
            plt.axhline(0, color='k', linestyle='--', linewidth=1) #Zero line
            plt.axvline(0, color='k', linestyle='--', linewidth=1) #Zero line
            plt.xlabel('$\Delta$h of ' + str(demName1_) + ' to SRTM [m]') #Z 
            plt.ylabel('$\Delta$h of ' + str(demName2_) + ' to SRTM [m]') #dh
            plt.title('Elevation differences $\Delta$h of '+str(demName1_)+' vs. '+str(demName2_)+' \n for ' + str(name_)) #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
            plt.grid(True)
            plt.text(np.min(pNpGlacier1[:,3]),np.min(pNpGlacier2[:,3]),'Linear regression (r='+str(round(rValue3,2))+')',bbox={'facecolor':'white', 'pad':10})
            plt.savefig("output/"+str(name_)+'_scatterplot_'+str(demName1_)+'_'+str(demName2_))#.pdf
            plt.clf()
        
       
        
        del pNpGlacier, pNpGlacier1, pNpGlacier2
        
        return
        
        
        
        
#_______________________________________________________________________________

class ModelMassBalance:
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
        
        
    def readCsvNumpy(self, csvFileName_):
        """Creates a empty numpy array with same shape as CSV file and return this numpy array.
        Argument dataType defines the data type of the resulting numpy array."""
                
        
        print("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
        print("Reading input csv-file '" + str(csvFileName_) + "' ...")
        
        #Get shape of CSV-file, so that numpy array can be created
        #-------------------------------------------------------------------------------
        pDocCsv = self.__openCsvFile(csvFileName_)
                
        #Iterate trough all rows and enumerate
        nCsvRows = 0
        for csvRow in pDocCsv:
            nCsvRows = nCsvRows + 1
        nCsvCols = len(csvRow)
        
        #Second step: Create empty numpy file for reading input data
        pDataType = np.float64
        pDocCsvNumpy = np.zeros((nCsvRows, nCsvCols), dtype = pDataType)
        
        
        #Read CSV-file cell by cell and save it to numpy array
        #-------------------------------------------------------------------------------
        #Rows represent time, columns data
        pDocCsv = self.__openCsvFile(csvFileName_)
        isVarName = False #If first row has var_names, determine later
        
        cntRow = 0 #Counter rows, index is cntRow-1
        for row in pDocCsv: #Each row in file
            cnCol = 0 #Counter columns, index is cnCol-1
            cntRow = cntRow + 1
            
            for col in row: #Each cell in row
                cnCol = cnCol + 1
                #print cnCol

                try: #Only read float values = data values 
                    inputValue = np.cast[pDataType](col) #input = e.g. float(col)
                    pDocCsvNumpy[cntRow-1, cnCol-1] = inputValue   #save each cell value to numpy array
                    #print 'Value at row:', cntRow, '; column: ', cnCol, '; value: ' \
                    #, pDocCsvNumpy[cntRow-1, cnCol-1] #Index values beginning with '0', not number of cols, rows

                except: #Values that are not numbers, like alphabetic values or ''
                    inputValue = np.cast[pDataType](NODATA) #input = e.g. float(nodata_)
                    pDocCsvNumpy[cntRow-1, cnCol-1] = inputValue   #Save default no-data value to cell
                    print("Value '" + str(col) + "' at row '" + str(cntRow) + "' and column '" + str(cnCol) + \
                    "' is not of data type '" + str(pDataType) + "'. Use nodata value '" + str(pDocCsvNumpy[cntRow-1, cnCol-1]) + "' instead.")
                    #Index values beginning with '0', not number of cols, rows
                    
                    if cntRow == 1: #If first row has no numbers, it is assumed that these are the column names
                        isVarName = True
                    else:
                        isVarName = False
                        
        if isVarName:   #If first row has var_names, delete first row that consists of nodata_values in numpy array
            pDocCsvNumpy = np.delete(pDocCsvNumpy, 0, 0 )
            print("Detected variable names in first row of file, deleted row in numpy file.")
            
        #Summary
        print("Number of rows found in cvs-file: '" + str(nCsvRows) + "'") #number of rows, e.g. 2493
        print("Number of columns found in cvs-file: '" + str(nCsvCols) + "'") #number of columns, e.g. 12
        print("Number of dimension in input pDocCsvNumpy: '" + str(pDocCsvNumpy.ndim) + "'") 
        print("Shape of input pDocCsvNumpy: '" + str(pDocCsvNumpy.shape) + "'") 
        print("Data Type of input pDocCsvNumpy: '" + str(pDocCsvNumpy.dtype) + "'") 
        print ("Data of registered input in pDocCsvNumpy:")
        print pDocCsvNumpy
        
        return pDocCsvNumpy
        
        
        
    def adaptCsvNumpy(self, pDocCsvNumpy_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
       
        #Create empty numpy file for reading adapted input data
        #-------------------------------------------------------------------------------
        nNpRows = pDocCsvNumpy_.shape[0]
        nNpCols = 4 #X,Y,Z,dh
        pNpGlacier = np.zeros((nNpRows, nNpCols), dtype = pDocCsvNumpy_.dtype)
                
        #Adapt columns for correct order here
        #-------------------------------------------------------------------------------
        #X-Value
        pNpGlacier[:,0] = pDocCsvNumpy_[:,4]
        #Y-Value
        pNpGlacier[:,1] = pDocCsvNumpy_[:,5]
        #Z-Value
        pNpGlacier[:,2] = pDocCsvNumpy_[:,2]
        #dh-Value
        pNpGlacier[:,3] = pDocCsvNumpy_[:,3]
        
        
        #Print Summary of numpy file 
        #-------------------------------------------------------------------------------
        print("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
        print("Number of dimension in input pNpGlacier: '" + str(pNpGlacier.ndim) + "'") 
        print("Shape of input pNpGlacier: '" + str(pNpGlacier.shape) + "'") 
        print("Data Type of input pNpGlacier: '" + str(pNpGlacier.dtype) + "'") 
        print ("Data of input numpy file for dimensions X, Y, Z, dh:")
        print pNpGlacier
        
        return pNpGlacier
        
        
        
               
    def readShpNumpy(self, shpFileName_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        print("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
        print("Reading input shapefile '" + str(shpFileName_) + "' ...")
        
        #Read OGR datasource (Shapefile) and check if exists
        #-------------------------------------------------------------------------------
        pOgrDriver = ogr.GetDriverByName('ESRI Shapefile') #register all the format drivers that are desired
        #Open the input OGR datasource, can be files, RDBMSes, directories full of files, or even remote web services depending on the driver being used. However, the datasource name is always a single string
        pFileShp = pOgrDriver.Open(shpFileName_, 0) # 0 means read-only. 1 means writeable.
        if pFileShp is None:
            raise Exception("Opening of file '" + str(shpFileName_) + "' failed. Check if it exists and if filename suffix is set.")
        
        
        #Get layer from Shapefile and check (only one layer allowed)
        #-------------------------------------------------------------------------------
        pLayer = pFileShp.GetLayer() #An OGRDataSource can potentially have many layers associated with it #lyr = ds.GetLayerByName( "point" )
        nLayers = pFileShp.GetLayerCount() #The number of layers available 
        if nLayers != 1: #Only one layer allowed
            raise Exception("Only one layer is allowed, '" + str(shpFileName_) + "' has '" + str(nLayers) + "' layers.")
        pLayer.ResetReading()
        
        
        #Getting some information from layer and print
        #-------------------------------------------------------------------------------
        pFeatureCount = pLayer.GetFeatureCount() #Number of features in layer
        print ("Number of features in layer: '" + str(pFeatureCount) + "' features.")
        pLayerDefinition = pLayer.GetLayerDefn() #Object, associated with the layer, containing the definitions of all the fields
        for i in range(pLayerDefinition.GetFieldCount()): #Print all the field names
            fieldTypeCode = pLayerDefinition.GetFieldDefn(i).GetType()
            print ("Field number '" + str(i) + "' with field name '" + str(pLayerDefinition.GetFieldDefn(i).GetName()) + "', field type '" + str(pLayerDefinition.GetFieldDefn(i).GetFieldTypeName(fieldTypeCode)) + \
                    "', field width '" + str(pLayerDefinition.GetFieldDefn(i).GetWidth()) + "', and field precision '" + str(pLayerDefinition.GetFieldDefn(i).GetWidth()) + "'.")
        
        
        #Create empty numpy file for reading input data
        #-------------------------------------------------------------------------------
        nNpRows = pFeatureCount
        nNpCols = 4 #X,Y,Z,dh
        pDataType = np.float64
        pNpGlacier = np.zeros((nNpRows, nNpCols), dtype = pDataType)
                
        
        #Loop over all the fields in features, and fetch and report the attributes based on their type: for dH and Z --> Get X,Y value from point Geometry, not from attribute table
        #-------------------------------------------------------------------------------
        cntRow = 0 #counter for rows = features
        for pFeature in pLayer: #Iterate over all Features
            cntRow += 1
            pFeatureDef = pLayer.GetLayerDefn()
            for cnCol in range(pFeatureDef.GetFieldCount()): #column = counter for each field name in features
                pFieldDef = pFeatureDef.GetFieldDefn(cnCol)
                
                if pFeatureDef.GetFieldDefn(cnCol).GetName() == LABEL_POINT_Z: #Z-Value
                    if (pFieldDef.GetType() == ogr.OFTInteger) or (pFieldDef.GetType() == ogr.OFTReal): #Only numbers allowed
                        fieldValue = pFeature.GetFieldAsDouble(cnCol)
                    else:
                        raise Exception("Can not read value in field '" + str(pFeatureDef.GetFieldDefn(cnCol).GetName()) + "' since it has the type '" + str(pFieldDef.GetType()) + "' and only integer or double is allowed.")
                    try: #Only read float values = data values 
                        inputValue = np.cast[pDataType](fieldValue) #input = e.g. float(col)
                        pNpGlacier[cntRow-1, 2] = inputValue   #save each cell value to numpy array
                    except: #Values that are not numbers, like alphabetic values or ''
                        raise Exception("Value '" + str(fieldValue) + "' at row '" + str(cntRow) + "' and column '" + str(cnCol) + \
                        "' is not of data type '" + str(pDataType) + "'. Cannot cast to '" + str(pDataType) + "'.")
                    
                elif pFeatureDef.GetFieldDefn(cnCol).GetName() == LABEL_DH:  #dh-Value
                    if (pFieldDef.GetType() == ogr.OFTInteger) or (pFieldDef.GetType() == ogr.OFTReal): #Only numbers allowed
                        fieldValue = pFeature.GetFieldAsDouble(cnCol)
                    else:
                        raise Exception("Can not read value in field '" + str(pFeatureDef.GetFieldDefn(cnCol).GetName()) + "' since it has the type '" + str(pFieldDef.GetType()) + "' and only integer or double is allowed.")
                    try: #Only read float values = data values 
                        inputValue = np.cast[pDataType](fieldValue) #input = e.g. float(col)
                        pNpGlacier[cntRow-1, 3] = inputValue   #save each cell value to numpy array
                    except: #Values that are not numbers, like alphabetic values or ''
                        raise Exception("Value '" + str(fieldValue) + "' at row '" + str(cntRow) + "' and column '" + str(cnCol) + \
                        "' is not of data type '" + str(pDataType) + "'. Cannot cast to '" + str(pDataType) + "'.")
                
                
                    #elif pFieldDef.GetType() == ogr.OFTString: fieldValue = pFeature.GetFieldAsString(cnCol)
                #elif pFeatureDef.GetFieldDefn(cnCol).GetName() == LABEL_POINT_X: #X-Value
                    #print "Z", fieldValue
                #elif pFeatureDef.GetFieldDefn(cnCol).GetName() == LABEL_POINT_Y: #Y-Value
                    #print "Y", fieldValue
                    
                #else:
                #    print("Ingnoring Field name '" + str(pFeatureDef.GetFieldDefn(cnCol).GetName()) + "', since only allowed is '" + str(LABEL_POINT_X) + "','" + str(LABEL_POINT_Y) + "', '" + str(LABEL_POINT_Z) + "', or '" + str(LABEL_DH) + "'")
            
            
            #Extract the geometry (X and Y) for each feature, and write out the point geometry x and y
            #-------------------------------------------------------------------------------
            geometryRef = pFeature.GetGeometryRef() #geometries are returned as a generic OGRGeometry pointer.
            if geometryRef is not None and (geometryRef.GetGeometryType() == ogr.wkbPoint or geometryRef.GetGeometryType() == ogr.wkbPoint25D): #if it is a point, wkbPoint is used above to convert the type for a wkbPoint25D (a point with a z coordinate) into the base 2D geometry type code (wkbPoint)
                #print "%.3f, %.3f" % ( geometryRef.GetX(), geometryRef.GetY() ) 
                pointX = geometryRef.GetX()
                pointY = geometryRef.GetY()
                
                try: #Only read float values = data values 
                    inputX = np.cast[pDataType](pointX) #input = e.g. float(col)
                    pNpGlacier[cntRow-1, 0] = inputX   #save each cell value to numpy array
                    inputY = np.cast[pDataType](pointY) #input = e.g. float(col)
                    pNpGlacier[cntRow-1, 1] = inputY   #save each cell value to numpy array
                except: #Values that are not numbers, like alphabetic values or ''
                    raise Exception("Coordinate '" + str(pointX) + "' or '" + str(pointY) + "' at row '" + str(cntRow) + \
                    "' is not of data type '" + str(pDataType) + "'. Cannot cast to '" + str(pDataType) + "'.")
            else:
                raise Exception("Only point geometry allowed, '" + str(shpFileName_) + "' has the geometry '" + str(geometryRef.GetGeometryType()) + "'.")
        
        pFileShp = None
        del pFileShp
        
        #Print Summary of numpy file 
        #-------------------------------------------------------------------------------
        print("Number of dimension in input pNpGlacier: '" + str(pNpGlacier.ndim) + "'") 
        print("Shape of input pNpGlacier: '" + str(pNpGlacier.shape) + "'") 
        print("Data Type of input pNpGlacier: '" + str(pNpGlacier.dtype) + "'") 
        print ("Data of input numpy file for dimensions X, Y, Z, dh:")
        print pNpGlacier
        
        return pNpGlacier
        
    
    
    
    def writeShpNumpy(self, pNpData_, shpFileName_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        print("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
        print("Save output data to layer '" + str(shpFileName_) + "'.")
        
        #Get the shapefile driver and create shapefile
        pOgrDriver = ogr.GetDriverByName('ESRI Shapefile')
        pFileShp = pOgrDriver.CreateDataSource(shpFileName_) # 0 means read-only. 1 means writeable.
        if pFileShp is None:
            raise Exception("Creation of output file '" + str(shpFileName_) + "' failed. Check if it exists and if filename suffix is set.")
        
        #Set spatial reference
        pSpatialReference = osr.SpatialReference()
        pSpatialReference.ImportFromESRI(SPATIAL_REF_ESRI) #or .ImportFromEPSG
        
        #Create layer
        pLayer = pFileShp.CreateLayer( "point_out", pSpatialReference, ogr.wkbPoint )
                
        #Create all fields for attribute table
        fieldX = ogr.FieldDefn(LABEL_POINT_X, ogr.OFTReal )
        fieldY = ogr.FieldDefn(LABEL_POINT_Y, ogr.OFTReal )
        fieldZ = ogr.FieldDefn(LABEL_POINT_Z, ogr.OFTReal )
        fieldDh = ogr.FieldDefn(LABEL_DH, ogr.OFTReal )
        #fieldDh.SetWidth( 32 )
        
        #Add the fields to the shapefile
        pLayer.CreateField(fieldX)
        pLayer.CreateField(fieldY)
        pLayer.CreateField(fieldZ)
        pLayer.CreateField(fieldDh)
                
        #Get the FeatureDefn for the shapefile
        pLayerDefinition = pLayer.GetLayerDefn()
        
        
        #Loop through the rows of the numpy array and create each point
        for cntRow in range(pNpData_.shape[0]):
            
            #Create empty point geometry
            pGeometry = ogr.Geometry(ogr.wkbPoint)
            pGeometry.SetPoint(0, pNpData_[cntRow,0], pNpData_[cntRow,1])
            
            #Create a new feature and set its geometry and attribute
            pFeature = ogr.Feature(pLayerDefinition)
            pFeature.SetGeometry(pGeometry)
            #pFeature.SetFID(cntRow)
            pFeature.SetField(LABEL_POINT_X, float(pNpData_[cntRow,0]))
            pFeature.SetField(LABEL_POINT_Y, float(pNpData_[cntRow,1]))
            pFeature.SetField(LABEL_POINT_Z, float(pNpData_[cntRow,2]))
            pFeature.SetField(LABEL_DH, float(pNpData_[cntRow,3]))
            
            #Add the feature to the shapefile
            pLayer.CreateFeature(pFeature)
            
            #Cleanup --> destroy the geometries and feature
            pGeometry.Destroy()
            pFeature.Destroy()
        
        
        #close and flush data
        pFileShp.Destroy()
                
        return
        
        
        
    def rasterizeShpGdal(self, shpFileName_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        print("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
                
        #Raster filename
        rasterFileName = shpFileName_.split(".")[0]+".tif"
        print("Rasterize input layer '" + str(shpFileName_) + "' to '" + str(rasterFileName) + "' by using GDAL (Does not work correctly).")
        
        # Open the data source and read in the extent
        pFileShp = ogr.Open(shpFileName_)
        pFileShpLayer = pFileShp.GetLayer()
        xMin, xMax, yMin, yMax = pFileShpLayer.GetExtent()
        
        # Create the destination data source
        xRes = int((xMax - xMin) / PIXELSIZE)
        yRes = int((yMax - yMin) / PIXELSIZE)
        pFileRaster = gdal.GetDriverByName('GTiff').Create(rasterFileName, xRes, yRes, 1, gdal.GDT_Float64)
        pFileRaster.SetGeoTransform((xMin, float(PIXELSIZE), 0, yMax, 0, -float(PIXELSIZE)))
        pBand = pFileRaster.GetRasterBand(1)
        pBand.SetNoDataValue(float(NODATA))
        
        # Get projection of shapefile and assigned to raster
        srs = osr.SpatialReference()
        srs.ImportFromWkt(pFileShpLayer.GetSpatialRef().__str__())
        pFileRaster.SetProjection(srs.ExportToWkt())
        
        #Rasterize
        #!!! Doesn't work correctly!!!
        gdal.RasterizeLayer(pFileRaster, [1], pFileShpLayer, options=['ALL_TOUCHED = True', 'ATTRIBUTE=%s' %LABEL_DH]) # burn_values=[0]
        
        del pFileShpLayer
        
        return
        
        
        
        
    def rasterizeShpArcPy(self, shpFileName_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        print("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
                
        #Raster filename
        rasterFileName = shpFileName_.split(".")[0]+".tif"
        print("Rasterize input layer '" + str(shpFileName_) + "' to '" + str(rasterFileName) + "' by using Arcpy.")
        
        #Convert Point Shapefile to raster
        arcpy.PointToRaster_conversion(shpFileName_, LABEL_DH, rasterFileName, "MEAN", "", float(PIXELSIZE))# Problem with dh value -1???
        
        #Set values below -9999 to Nodata (necessary due to previous step) 
        arcpy.CheckOutExtension("spatial") #Load sa license
        arcpy.env.overwriteOutput = True #Enable overwrite
        pOutSetNull = sa.SetNull(rasterFileName, rasterFileName, "VALUE < -9999.0") #SetNull(inRaster, inFalseRaster, whereClause)
        pOutSetNull.save(rasterFileName)
        
        gc.collect()
                
        return
        
        
        
        
    def replaceGlacierAreas(self, pNpGlacier_, demName_, glacierNamesListReplace_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        #Optional: Replace specific areas by already processed glaciers (read from output folder) - Does not work correctly
        pNpGlacierOrig = np.copy(pNpGlacier_)
        
        for glacierNameReplace in glacierNamesListReplace_: 
            
            shpFileNameReplace = "output/out_"+str(glacierNameReplace)+"_mustag_srtm_xyzdh_" + str(demName_) + ".shp"
            print ("\nReplace dh-values by glacier shapefile '" + str(shpFileNameReplace) + "' for difference image '" + str(demName_) + "'...")
            pNpGlacierReplace = self.readShpNumpy(str(shpFileNameReplace))
            
            cntManRepl = 0
            
            #for iterPNpGlacierReplace in np.nditer(pNpGlacierReplace):
            for entryReplList in pNpGlacierReplace: #For all entries in numpy array containing replacing glacier values
                
                x = float(entryReplList[0])
                y = float(entryReplList[1])
                
                test = ((pNpGlacierOrig[:,0] == np.array([x])) & (pNpGlacierOrig[:,1] == np.array([y])))
                #print "test", pNpGlacierOrig[test,:]#, pNpGlacierOrig[testY,1]
                #mask = np.in1d(pNpGlacierOrig[:,[0,1]].ravel(), np.array([x,y])).reshape(pNpGlacierOrig.shape[0],2)
                #print "mask", pNpGlacierOrig[mask,[0,1]]#, pNpGlacierOrig[maskY,1]
                
                #indexList = np.where((pNpGlacierOrig[:,0] == float(entryReplList[0])) & (pNpGlacierOrig[:,1] == float(entryReplList[1]))) #Find index where x and y coordinate match with x and y coordinates of replacing glacier
                indexList = np.where(np.logical_and(pNpGlacierOrig[:,0] == float(entryReplList[0]), pNpGlacierOrig[:,1] == float(entryReplList[1]))) #Find index where x and y coordinate match with x and y coordinates of replacing glacier
                
                if pNpGlacierOrig[indexList[0],0].size != 0: #If index in form of list/tuple is not [], meaning that it found a matching value for x  and y  (index has one entry, so size !=0 )
                    
                    cntManRepl += 1
                    print ("Replace dh-value '" + str( pNpGlacierOrig[indexList[0],3]) + "' by '" + str(entryReplList[3]) + "' at position x: '"  + str( pNpGlacierOrig[indexList[0],0]) + "', y: '"  + str( pNpGlacierOrig[indexList[0],1]) + "'.")
                    pNpGlacierOrig[indexList[0],3] = float(entryReplList[3]) #replace dh value
                    
            
                #Does work, but takes too much processing time:
            #    for index in range(np.size(pNpGlacierOrig,axis=0)): #Iterate through all indices until it finds the first matching entry
                
            #        if ((float(pNpGlacierOrig[index,0]) == float(entryReplList[0])) and (float(pNpGlacierOrig[index,1]) == float(entryReplList[1]))):# and (pNpGlacierOrig[index,3] != entryReplList[3])): #If x and y coordinates are the same, but not dh
                        
            #            print ("Replace dh-value '" + str(pNpGlacierOrig[index,3]) + "' by '" + str(entryReplList[3]) + "' at position x: '"  + str(pNpGlacierOrig[index,0]) + "', y: '"  + str(pNpGlacierOrig[index,1]) + "'.")
            #            pNpGlacierOrig[index,3] = float(entryReplList[3]) #replace dh value
            #            cntManRepl += 1
            #            index = np.size(pNpGlacierOrig,axis=0) #Abortion of loop
            
            
            print ("Done. Replaced '" + str(cntManRepl) + "'  dh-values by glacier shapefile '" + str(shpFileNameReplace) + "' for difference image '" + str(demName_) + "'.")
            
            
            
        return pNpGlacierOrig
        
        
        
        
    def setValidDhRangeAccAbl(self, pNpGlacierQuants_, methodAcc_, methodAbl_, methodManual_, ela_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        
        #Manually adapting z-values before outlier detection
        #-------------------------------------------------------------------------------
        if methodManual_ != "":
            pNpGlacierQuants_ = self.__manualZForElaShift(pNpGlacierQuants_, methodManual_, "before")
        
        
        #Sort values in numpy arrays for accumulation and ablation area for statistics
        #-------------------------------------------------------------------------------
        pListGlacierAcc = [] #Python list accumulation
        pListGlacierAbl = [] #Python list ablation
        
        #Extract values depending on ELA 
        cntRow = 0
        for zValue in pNpGlacierQuants_[:,2]:
            cntRow+=1
            if (float(zValue) >= float(ela_)): #Accumulation area
                pListGlacierAcc.append(pNpGlacierQuants_[cntRow-1,:])
            elif (float(zValue) < float(ela_)):  #Ablation area
                pListGlacierAbl.append(pNpGlacierQuants_[cntRow-1,:])
                
        
        
        #Accumulation area
        #-------------------------------------------------------------------------------
        print("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
        if (pListGlacierAcc): #If data for accumulation area exists and outlier method is mentioned
            print("Set valid dh range for accumulation area from '" + str(np.min(np.asarray(pListGlacierAcc)[:,2])) + "'m to '" + str(np.max(np.asarray(pListGlacierAcc)[:,2])) + "'m with method '" + str(methodAcc_) + "':")
            pNpGlacierExclNodataAcc = self.setValidDhRange(np.asarray(pListGlacierAcc), methodAcc_)
            
        else:
            print (" No accumulation area with nodata-values, ELA '" + str(ela_) + "' is above maximum glacier height of '" + str(np.max(pNpGlacierQuants_[:,2])) + "'.")
        
        
        
        #Ablation area for each section
        #-------------------------------------------------------------------------------
        
        if (pListGlacierAbl): #If data for ablation area exists and outlier method is mentioned
            print("\nSet valid dh range for ablation area from '" + str(np.min(np.asarray(pListGlacierAbl)[:,2])) + "'m to '" + str(np.max(np.asarray(pListGlacierAbl)[:,2])) + "'m with method '" + str(methodAbl_) + "':")
            #Sort all entries along the zValues
            pNpGlacierAbl = np.asarray(pListGlacierAbl)
            pNpGlacierAblSortZ = pNpGlacierAbl[pNpGlacierAbl[:,2].argsort()] 
        
            
            #Calculate mean values in sorted zValues array for each Z_ABL_INTERVAL section
            #Initialization
            zLocalMax = pNpGlacierAblSortZ[0,2] + float(Z_ABL_INTERVAL) #First local maximal z Value (deepest z-Value + half offset)
            pListGlacierAblSection = [] #List for first saving  dh values
            cntRow = 0 #Counter all rows
            firstRun = True #First run of loop
            
            
            #Calculate all mean dh-Values per section until final Z-Max = ELA
            
            for zValue in pNpGlacierAblSortZ[:,2]: #For all dh-Values within a section
                
                pListGlacierAblSection.append(pNpGlacierAblSortZ[cntRow,:])
                cntRow += 1
                
                if zValue > zLocalMax:  #(zValue <= zLocalMax and zValue != pNpGlacierAblSortZ[-1,2]): #Z-Value reaches next section
                    
                    if firstRun == True:
                        pNpGlacierExclNodataAbl = self.setValidDhRange(np.asarray(pListGlacierAblSection), methodAbl_)
                        firstRun = False
                    else:
                        pNpGlacierExclNodataAbl = np.concatenate((pNpGlacierExclNodataAbl, self.setValidDhRange(np.asarray(pListGlacierAblSection), methodAbl_)), axis=0)
                        
                    pListGlacierAblSection = [] #Empty Python list ablation
                    zLocalMax = zLocalMax + float(Z_ABL_INTERVAL) #Next local Z-max value
                    
            if pListGlacierAblSection: #If list is not empty, so there are remaining values between zLocalMax and ELA
                
                pNpGlacierExclNodataAbl = np.concatenate((pNpGlacierExclNodataAbl, self.setValidDhRange(np.asarray(pListGlacierAblSection), methodAbl_)), axis=0)
                
        else:
            print (" No ablation area with nodata-values, ELA '" + str(ela_) + "' is below minimum glacier height of '" + str(np.min(pNpGlacierQuants_[:,2])) + "'.")
        
        
        #In case: Manually adapting z-values after outlier detection (countermand)
        #-------------------------------------------------------------------------------
        if (pListGlacierAcc and not pListGlacierAbl): #only accumulation area
            if methodManual_ != "":
                return self.__manualZForElaShift(pNpGlacierExclNodataAcc, methodManual_, "after")
            else:
                return pNpGlacierExclNodataAcc
        elif (not pListGlacierAcc and pListGlacierAbl): #only ablation area
            if methodManual_ != "":
                return self.__manualZForElaShift(pNpGlacierExclNodataAbl, methodManual_, "after")
            else:
                return pNpGlacierExclNodataAbl
        else: #both ablation and accumulation area
            if methodManual_ != "":
                return self.__manualZForElaShift(np.concatenate((pNpGlacierExclNodataAcc, pNpGlacierExclNodataAbl), axis=0), methodManual_, "after")
            else:
                return np.concatenate((pNpGlacierExclNodataAcc, pNpGlacierExclNodataAbl), axis=0)
        
        
        
    
    def setValidDhRange(self, pNpGlacierQuants_, method_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        
        #Delete all nodata values in numpy file for following calculation
        #-------------------------------------------------------------------------------
        pNpGlacierExclNodata = np.copy(pNpGlacierQuants_)
        pNpGlacierExclNodata = self.deleteNoDataValues(pNpGlacierExclNodata, "del") #not use pNpGlacierQuants_ as argument, otherwise will be nans!!!!
        
        
        #print "Set new upper and lower bounds for input data based on quantiles:"
        #print " dH-Threshold found in global constants: ", str(DH_THRESHOLD)
        
        
        #Calculate valide bounds and eliminate outliers 
        #-------------------------------------------------------------------------------
        #Calculate 5% and 95% as well as 31.7% and 68.3% quantiles in dh set as lower and upper bounds to eliminate most extreme outliers
        
        q_0pt317 = stats.mstats.mquantiles(pNpGlacierExclNodata[:,3], alphap =0.05, betap=0.05, prob=[0.317])#), limit=(pNpGlacierExclNodata[0,3],pNpGlacierExclNodata[-1,3]))# # lower bound #alphap =0.05, betap=0.05
        q_0pt683 = stats.mstats.mquantiles(pNpGlacierExclNodata[:,3], alphap =0.05, betap=0.05, prob=[0.683])#, limit=(pNpGlacierExclNodata[0,3],pNpGlacierExclNodata[-1,3])) # # upper bound
        q_0pt05 = stats.mstats.mquantiles(pNpGlacierExclNodata[:,3], alphap =0.05, betap=0.05, prob=[0.05])#), limit=(pNpGlacierExclNodata[0,3],pNpGlacierExclNodata[-1,3]))# # lower bound #alphap =0.05, betap=0.05
        q_0pt95 = stats.mstats.mquantiles(pNpGlacierExclNodata[:,3], alphap =0.05, betap=0.05, prob=[0.95])#, limit=(pNpGlacierExclNodata[0,3],pNpGlacierExclNodata[-1,3])) # # upper bound
        #print " Bounds within 5% and 95% as well as 31.7% and 68.3%  quantiles: ", q_0pt05, q_0pt95, q_0pt317, q_0pt683
        
        #Calculate quantiles for 0.25, 0.5 and 0.75 within bounds of 5% and 95% quantiles
        q1 = stats.mstats.mquantiles(pNpGlacierExclNodata[:,3], alphap =0.05, betap=0.05, prob=[0.25], limit=(q_0pt05,q_0pt95)) #alphap =0.05, betap=0.05
        q2 = stats.mstats.mquantiles(pNpGlacierExclNodata[:,3], alphap =0.05, betap=0.05, prob=[0.5], limit=(q_0pt05,q_0pt95))
        q3 = stats.mstats.mquantiles(pNpGlacierExclNodata[:,3], alphap =0.05, betap=0.05, prob=[0.75], limit=(q_0pt05,q_0pt95))
        #print " Quantiles for 0.25, 0.5 and 0.75 within bounds of 5% and 95% quantiles ", q1, q2, q3
        
        #Calculate 1.5*IQR and calculated new upper and lower bounds lying in between the two-tailed 1,5 time IQR
        iqr = 1.5*abs(q3-q1) # 1.5 interquartile range
        lowerBoundIqr = q2-iqr
        upperBoundIqr = q2+iqr
        #print " 1.5*IQR ('"+str(iqr)+"') within bounds of 5% and 95% quantiles:", lowerBoundIqr, upperBoundIqr
        
        
        #Replace dh-values beyond calculated bounds: dhNull_q0pt683_q0pt95_q0pt683_iqr_dhThreshold___nodata_del 
        #-------------------------------------------------------------------------------
        if "dhNull" in method_: #by only using pixels lying in between the two-tailed 1,5 time IQR
            newLowerBound = float(0) #Set new lower bound for request condition
            newUpperBound = float(0) #Set new upper bound for request condition
        elif "q0pt683" in method_:
            newLowerBound = q_0pt317 #Set new lower bound in for request condition
            newUpperBound = q_0pt683 #Set new upper bound in for request condition
        elif "q0pt95" in method_:
            newLowerBound = q_0pt05 #Set new lower bound in for request condition
            newUpperBound = q_0pt95 #Set new upper bound in for request condition
        elif "iqr" in method_: #by only using pixels lying in between the two-tailed 1,5 time IQR
            newLowerBound = lowerBoundIqr #Set new lower bound for request condition
            newUpperBound = upperBoundIqr #Set new upper bound for request condition
        elif "dhThreshold" in method_: #by only using pixels lying in between the two-tailed 1,5 time IQR
            newLowerBound = (-1) * float(DH_THRESHOLD) #Set new lower bound for request condition
            newUpperBound = float(DH_THRESHOLD) #Set new upper bound for request condition
       
        else: #No replacement, use max and min values as bounds
            newLowerBound = np.min(pNpGlacierExclNodata[:,3])
            newUpperBound = np.max(pNpGlacierExclNodata[:,3])
            #raise Exception("Please set 'method_' correctly, variables 'newLowerBound' and 'newUpperBound' could not be set (method_ = '"+str(method_)+"').")
        
        if "nodata" in method_: #set to value NODATA beyond bounds
            setLowerBound = float(NODATA) #Set new lower bound in dh-Values
            setUpperBound = float(NODATA) #Set new upper bound in dh-Values
        elif "del" in method_: #delete values beyond bounds
            setLowerBound = np.nan #Set new lower bound in dh-Values
            setUpperBound = np.nan #Set new upper bound in dh-Values
        else: #use max and min values beyond bounds
            setLowerBound = newLowerBound #Set new lower bound in dh-Values
            setUpperBound = newUpperBound #Set new upper bound in dh-Values
            #raise Exception("Please set 'method_' correctly, variables 'setLowerBound' and 'setLowerBound' could not be set (method_ = '"+str(method_)+"').")
        
        
        cntRow = 0
        cntLb = 0
        cntUb = 0
        
        for dhValue in pNpGlacierQuants_[:,3]:
            cntRow+=1
            if (dhValue < newLowerBound and dhValue != float(NODATA)) : # dh value lower then lower bound
                cntLb+=1
                pNpGlacierQuants_[cntRow-1,3] = setLowerBound
            
            elif (dhValue > newUpperBound and dhValue != float(NODATA)): #dh value higher then upper bound
                cntUb+=1
                pNpGlacierQuants_[cntRow-1,3] = setUpperBound
                
            #else:
            #    print pNpGlacierQuants_[cntRow-1,3] 
        
        
        pNpGlacierQuants_ = pNpGlacierQuants_[~np.isnan(pNpGlacierQuants_).any(axis=1)] #Delete rows with nan-values in case of method "del"
        
        
        #Summary
        #-------------------------------------------------------------------------------
        #Calculating area per pixel and complete glacier surface area --> A
        areaPerPixel = pow(float(PIXELSIZE),2) #Area per Pixel, depending on pixel size --> (pixel_size * pixel_size) [m²]
        surfaceArea = pNpGlacierQuants_.shape[0] * areaPerPixel  #Complete surface area of glacier --> (number of amount of dh-pixels * area per pixel) [m²]
        
        #if method_ != "":
        #    print("--> '"+str(cntLb)+"' dh-values were lower then lower bound '"+str(newLowerBound)+"' and replaced with '"+str(setLowerBound)+"', '"+str(cntUb)+"' dh-values were higher then upper bound '"+str(newUpperBound)\
        #          +"' and replaced with '"+str(setUpperBound)+"'. New shape: '"+str(pNpGlacierQuants_.shape)+"' corresponds to surface area of '" + str(surfaceArea) + "'m² (~'" + str(round(surfaceArea/pow(10,6),2)) + "'km²).")
        #else:
        #    print ("No outlier replacement employed since method is empty.")
        
        return pNpGlacierQuants_
        
        
        
        
        
        
        
    def fillNoDataAreas(self, pNpGlacierFill_, methodAcc_, methodAbl_, methodManual_, ela_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        
        pNpGlacierCalc = np.copy(pNpGlacierFill_)
        
        #Sort values in numpy arrays for accumulation and ablation area for statistics
        #-------------------------------------------------------------------------------
                
        pListGlacierAcc = [] #Python list accumulation
        pListGlacierAbl = [] #Python list ablation
        
        #Extract values depending on ELA 
        cntRow = 0
        cntNodata = 0
        for zValue in pNpGlacierCalc[:,2]:
            cntRow+=1
            if (float(zValue) > float(ela_) and float(pNpGlacierCalc[cntRow-1,3]) != (float(NODATA))): #Accumulation area
                pListGlacierAcc.append(pNpGlacierCalc[cntRow-1,:])
                #print pNpGlacierCalc[cntRow-1,:]
               
            elif (float(zValue) <= float(ela_) and float(pNpGlacierCalc[cntRow-1,3]) != float(NODATA)):  #Ablation area
                pListGlacierAbl.append(pNpGlacierCalc[cntRow-1,:])
                #print pNpGlacierCalc[cntRow-1,:]
            else:
                cntNodata+=1
            #    print ("Nodata value for dh thickness change entry found at index row '" + str(cntRow-1) + "': '" + str(pNpGlacierCalc[cntRow-1,:]) + "'") 
                
        
        nodataArea = cntNodata * pow(float(PIXELSIZE),2) 
        
        print("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
        print ("Gap filling for '" + str(cntNodata) + "' nodata dh-pixels corresponds to surface area of '" + str(nodataArea) + "'m² (~'" + str(round(nodataArea/pow(10,6),2)) + "'km²):")
        print ("--> Separating ablation and accumulation zone by ELA of '" + str(ela_) + "' meter by ignoring nodata-values:")
        
        
        #Calculate values for accumulation area
        #-------------------------------------------------------------------------------
        if pListGlacierAcc:
            pNpGlacierAcc = np.empty((len(pListGlacierAcc), 4), dtype = pNpGlacierCalc.dtype) #X,Y,Z,dh
            pNpGlacierAcc = np.asarray(pListGlacierAcc)
            
            #Find mean value in accumulation area
            meanAcc = np.mean(pNpGlacierAcc[:,3]) #mean value of entire ablation area 
            
            #Find mode value in accumulation area (most frequent value
            pNpGlacierAccInt = pNpGlacierAcc.astype(int)
            [modeAcc], [cntsModeAcc] = stats.mode(pNpGlacierAccInt[:,3])#dH
            
            accArea = pNpGlacierAcc.shape[0] * pow(float(PIXELSIZE),2)
            
            print (" dh-pixels in accumulation area: '" + str(pNpGlacierAcc.shape[0]) + "' ('" + str(accArea) + "'m² =~'" + str(round(accArea/pow(10,6),2)) + "'km²); Min: '" + str(np.min(pNpGlacierAcc[:,3])) + "'; Max: '" + str(np.max(pNpGlacierAcc[:,3])) \
                + "'; Mean: '" + str(meanAcc) + "'; Standard deviation: '" + str(np.std(pNpGlacierAcc[:,3])) + "'; Median: '" + str(np.median(pNpGlacierAcc[:,3]))\
                 + "'; Integer mode: '" + str(modeAcc) + "' with '" + str(cntsModeAcc) + "' counts.") 
        else:
            print (" Only nodata values in accumulation area or ELA '" + str(ela_) + "' is above maximum glacier height of '" + str(np.max(pNpGlacierCalc[:,2])) + "'.")
            modeAcc = 0
            
            
        #Calculate values for ablation area
        #-------------------------------------------------------------------------------
        if pListGlacierAbl:
            pNpGlacierAbl = np.empty((len(pListGlacierAbl), 4), dtype = pNpGlacierCalc.dtype) #X,Y,Z,dh
            pNpGlacierAbl = np.asarray(pListGlacierAbl)
            
            #Find mean value in ablation area
            meanAbl = np.mean(pNpGlacierAbl[:,3]) #mean value of entire ablation area 
            pMeanAbl = self.__calcMeanDhNodataAbl(pNpGlacierAbl, ela_)  #Calculate mean values in sorted zValues array for each Z_ABL_INTERVAL section
            
            #Find mode value in ablation area (most frequent value)
            pNpGlacierAblInt = pNpGlacierAbl.astype(int)
            [modeAbl], [cntsModeAbl] = stats.mode(pNpGlacierAblInt[:,3])#dH
            
            ablArea = pNpGlacierAbl.shape[0] * pow(float(PIXELSIZE),2)
            
            print (" dh-pixels in ablation area: '" + str(pNpGlacierAbl.shape[0]) + "' ('" + str(ablArea) + "'m² =~'" + str(round(ablArea/pow(10,6),2)) + "'km²); Min: '" + str(np.min(pNpGlacierAbl[:,3])) + "'; Max: '" + str(np.max(pNpGlacierAbl[:,3])) \
                + "'; Mean: '" + str(np.mean(pNpGlacierAbl[:,3])) + "'; Standard deviation: '" + str(np.std(pNpGlacierAbl[:,3])) + "'; Median: '" + str(np.median(pNpGlacierAbl[:,3])) \
                + "' Integer mode: '" + str(modeAbl) + "' with '" + str(cntsModeAbl) + "' counts.") 
        else:
            print ("  Only nodata values in ablation area or ELA '" + str(ela_) + "' is below minimum glacier height of '" + str(np.min(pNpGlacierCalc[:,2])) + "'.")
            modeAbl = 0
        
        #Replacing nodata-values
        #-------------------------------------------------------------------------------
        
        if ("accMode" in methodAcc_) or ("accNull" in methodAcc_) or ("ablNull" in methodAbl_) or ("ablMean" in methodAbl_) or ("manual" in methodManual_):
            #
            
            # Automatically replacing nodata values in accumulation and ablation area with values set above
            cntRow = 0
            cntAccReplc = 0 #counter replaced values for accumulation area
            cntAblReplc = 0 #counter replaced values for ablation area
            cntManReplc = 0 #counter manually replaced values for specific area
            
            
            for dhValue in pNpGlacierCalc[:,3]:
                cntRow+=1
                
                if float(dhValue) == float(NODATA): #for nodata values
                    
                    if float(pNpGlacierCalc[cntRow-1,2]) > float(ela_): # General replacement accumulation area
                        
                        if "accMode" in methodAcc_:
                            nodataReplaceAcc = float(modeAcc) #mode value of accumulation area
                        elif "accMean" in methodAcc_:
                            nodataReplaceAcc = float(meanAcc) #mean value of accumulation area
                        elif "accNull" in methodAcc_:
                            nodataReplaceAcc = 0
                        else:
                            raise Exception ("No replacement method for nodata value in accumulation zone chosen (method '" + str(methodAcc_) + "')")
                        
                        cntAccReplc+=1
                        pNpGlacierCalc[cntRow-1,3] = nodataReplaceAcc
                        
                        
                    elif float(pNpGlacierCalc[cntRow-1,2]) <= float(ela_): #General replacement ablation area
                        
                        if "ablMean" in methodAbl_:
                            #nodataReplaceAbl = float(meanAbl) #mean value of entire ablation area 
                            
                            #Use correspondend mean value calculated for each Z_ABL_INTERVAL section
                            pDiffZ = pMeanAbl[:,0] - pNpGlacierCalc[cntRow-1,2] #allMeanSectorZValues minus Z-Value of Nodata point
                            nodataReplaceAbl = pMeanAbl[np.argmin(abs(pDiffZ)),1] #Get dhMean-value were dhZ-Value to nodata point is minimal
                            
                            
                        elif "ablNull" in methodAbl_:
                            nodataReplaceAbl = 0
                        else:
                            raise Exception ("No replacement method for nodata value in ablation zone chosen (method '" + str(methodAbl_) + "').")
                       
                        cntAblReplc+=1
                        pNpGlacierCalc[cntRow-1,3] = nodataReplaceAbl
                       
                       
                    else:
                        raise Exception ("No replacement for nodata value '" + str(NODATA) + "' at index row '" + str(cntRow-1) + "' was set, please manually do so: " + str(pNpGlacierCalc[cntRow-1,:])+ "'")
                    
                    
                    #!!!!!!!!!!!! Manual replacement per specific area afterwards: !!!!!!!!!!!!!!
                    if methodManual_ != "":
                        pNpGlacierCalc, cntManReplc = self.__manualDhNodata(pNpGlacierCalc, cntRow, cntManReplc, methodManual_, ela_)
                        
                            
            accAreaNodata = cntAccReplc * pow(float(PIXELSIZE),2)
            ablAreaNodata = cntAblReplc * pow(float(PIXELSIZE),2)
            manAreaNodata = cntManReplc * pow(float(PIXELSIZE),2)
            
            
            print ("'" + str(cntAccReplc) + "' nodata values ('" + str(accAreaNodata) + "'m² =~'" + str(round(accAreaNodata/pow(10,6),2)) + "'km²) in accumulation zone replaced with method '" + str(methodAcc_) + "' and '" + str(cntAblReplc) \
                   + "' nodata values ('" + str(ablAreaNodata) + "'m² =~'" + str(round(ablAreaNodata/pow(10,6),2)) + "'km²) in ablation zone replaced with method '" + str(methodAbl_) + "'.'")
            print (" Manually '" + str(cntManReplc) + "' nodata values ('" + str(manAreaNodata) + "'m² =~'" + str(round(manAreaNodata/pow(10,6),2)) + "'km²) replaced with other values for '" + str(methodManual_) + "'.")
            
        else:
            print ("No replacement of nodata values was employed.")
            
            
        #Summary
        #-------------------------------------------------------------------------------
        
        #Find mode value 
        pNpGlacierCalcInt = pNpGlacierCalc.astype(int)
        [modeAll], [cntsModeAll] = stats.mode(pNpGlacierCalcInt[:,3])#dH
        
        print ("All glacier dh pixels after nodata-filling: '" + str(pNpGlacierCalc.shape[0]) + "'; Min: '" + str(np.min(pNpGlacierCalc[:,3])) + "'; Max: '" + str(np.max(pNpGlacierCalc[:,3])) \
                + "'; Mean: '" + str(np.mean(pNpGlacierCalc[:,3])) + "'; Standard deviation: '" + str(np.std(pNpGlacierCalc[:,3])) + "'; Median: '" + str(np.median(pNpGlacierCalc[:,3])) \
                + "' Integer mode: '" + str(modeAll) + "' with '" + str(cntsModeAll) + "' counts.") 
        
        
        
        return pNpGlacierCalc 
        
        
        
        
        
    def deleteNoDataValues(self, pNpGlacierDel_, method_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        
        #Delete complete rows with nodata values in dh in numpy array
        #-------------------------------------------------------------------------------
        if "zero" in method_: #Set nodata value to 0
            nodataValue = 0
        else: #elif "del" in method_: #Set nodata value to NAN and delete row afterwards
            nodataValue = np.NAN
        
        #Set nodata value
        cntRow = 0
        cntNodata = 0 #Counter nodata
        for dhValue in pNpGlacierDel_[:,3]:
            cntRow+=1
            if float(dhValue) == float(NODATA):
                pNpGlacierDel_[cntRow-1,3] = nodataValue
                cntNodata+=1
        
        #Delete all rows that contain np.nan
        #Explanation: np.isnan(a) returns a similar array with True where NaN, False elsewhere. 
        #.any(axis=1) reduces an m*n array to n with an logical or operation on the whole rows,
        # ~ inverts True/False and a[ ] chooses just the rows from the original array, which have True within the brackets.
        pNpGlacierDel_ = pNpGlacierDel_[~np.isnan(pNpGlacierDel_).any(axis=1)] #Delete rows with nan-values
        
        
        #Calculating area per pixel and complete glacier surface area --> A
        areaPerPixel = pow(float(PIXELSIZE),2) #Area per Pixel, depending on pixel size --> (pixel_size * pixel_size) [m²]
        surfaceArea = pNpGlacierDel_.shape[0] * areaPerPixel  #Complete surface area of glacier --> (number of amount of dh-pixels * area per pixel) [m²]
        
        #print("'" + str(cntNodata) + "' nodata values found in numpy file and replaced with '"+str(nodataValue)+"'. New shape: '" + str(pNpGlacierDel_.shape) + "' corresponds to surface area of '" + str(surfaceArea) + "'m² (~'" + str(round(surfaceArea/pow(10,6),2)) + "'km²).")
        
        
        
        return pNpGlacierDel_
        
        
        
        
    def adaptElevationValues(self, pNpGlacier_, offsetValue_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        
        pNpGlacierAdapt = np.copy(pNpGlacier_)
        
        #-------------------------------------------------------------------------------
        #Adapt elevation values on entire area
        cntRow = 0
        cntNodata = 0 #Counter nodata
        
        for dhValue in pNpGlacier_[:,3]:
            if float(dhValue) != float(NODATA):
                pNpGlacierAdapt[cntRow,3] = dhValue + float(offsetValue_)
            else:
                cntNodata+=1
            cntRow+=1
        
        print("\nSet an offset value of '" + str(offsetValue_) + "'m to dH-Values in Numpy array by ignoring '" + str(cntNodata) + "' nodata values.")
        
        return pNpGlacierAdapt
        
        
        
        
    def adaptSrtmPenetration(self, pNpGlacier_, ela_, srtmPenetration_, debrisBoundary_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """    
        
        # penetrationDepthAcc, penetrationDepthAbl, srtmPenetrationError
        penetrationDepthAcc = srtmPenetration_[0]
        penetrationDepthAbl = srtmPenetration_[1]
        
        pNpGlacierAdapt = np.copy(pNpGlacier_)
        
        #-------------------------------------------------------------------------------
        #Adapt elevation values only on accumulation area
        #Correct e.g. SRTM C-Band radar penetration or mean offset of stable terrain
        cntRow = 0
        cntAbl = 0
        cntAcc = 0
        cntDebris = 0
        cntNodata = 0
        
        for zValue in pNpGlacierAdapt[:,2]:
            if float(pNpGlacier_[cntRow,3]) != float(NODATA):
                if (float(zValue) >= float(ela_)): #Accumulation area
                    pNpGlacierAdapt[cntRow,3] += float(penetrationDepthAcc)
                    cntAcc += 1
                elif (float(zValue) >= float(debrisBoundary_)): #Ablation area above eventual debris covered tongue
                    pNpGlacierAdapt[cntRow,3] += float(penetrationDepthAbl)
                    cntAbl += 1
                else: #below debris boundary --> debris covered area
                    cntDebris += 1
            else:
                cntNodata += 1
            
            cntRow+=1
        
        print("\nAdapted SRTM C-band penetration in glacier by '" + str(penetrationDepthAcc) + "'m for '" + str(cntAcc) + "' pixel of accumulation zone and '" + str(penetrationDepthAbl) + \
              "'m for '" + str(cntAbl) + "' pixel of glacier accumulation (ELA altitude of '" + str(ela_) + "'m). Ignoring '" + str(cntNodata) + "' nodata pixel and '" + str(cntDebris) + \
              "' debris covered pixel below altitude of '" + str(debrisBoundary_) + "'m.")

        
        return pNpGlacierAdapt
        
        
        
        
        
        
        
    def calculateMassBalance(self, pNpGlacier_, pDateOld_, pDateNew_, dtmError_, srtmPenetrationError_, glacierName_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
       
        
        #-------------------------------------------------------------------------------
        #Calculate date difference of DTMs from constants
        pDateDiff = pDateNew_ - pDateOld_
        #dateDiffYears = round(abs(float(pDateDiff.days)/365),2) #floating years 
        dateDiffYears = int(abs(float(pDateDiff.days)/365)) #integer years --> years of full mass balance
        
        
        #-------------------------------------------------------------------------------
        #Calculating area per pixel and complete glacier surface area --> A
        areaPerPixel = pow(float(PIXELSIZE),2) #Area per Pixel, depending on pixel size --> (pixel_size * pixel_size) [m²]
        glacierSurfaceArea = pNpGlacier_.shape[0] * areaPerPixel  #Complete surface area of glacier --> (number of amount of dh-pixels * area per pixel) [m²]
        
        
        #-------------------------------------------------------------------------------
        #Calculating volume changes --> dV =  [m³/a]
        
        sumThicknessChange = np.sum(pNpGlacier_[:,3]) #Sum of all glacier elevation difference values [m]
        
        #B_i = Glacier net balance = total net balance = sum of accumulation and ablation over the whole glacier surface = volume change of glacier [m³]
        sumVolumeChange = sumThicknessChange * areaPerPixel #[m]*[m²] = [m³]
        #B_i_dot = Glacier net balance rate = total net balance rate [m³/a] 
        perYearVolumeChange = sumVolumeChange / dateDiffYears # [m³/a] --> Volume change (per year)!!!!
        
        
        #B = Glacier net balance = total net balance = sum of accumulation and ablation over the whole glacier surface = volume change of glacier [kg]
        sumMassBalance = sumVolumeChange * float(ICE_DENSITY) #[m³]*[kg/m³] = [kg] 
        #B_dot = Glacier net balance rate = total net balance rate #[kg/a]
        perYearSumMassBalance = sumMassBalance / dateDiffYears #[kg/a] --> Complete ice mass change (per year)!!!
        
        
        
        #-------------------------------------------------------------------------------
        #Calculating thickness change dH = [m/a] and mass balances B
        
        #mass = volume * density
        #M [Massenbilanz] = (V [Volumen] * p [Dichte]) / A [surface area] --> kg/m² ~ m. w. e. 
        #--> ([kg]/[m²] for ice) / (density of water = 999.972 (kg/m3) = ~ 1000 (kg/m3)) --> meter water equivalent
        #--> Depth in kg of ice over one square meter [kg]/[m²]. The unit corresponds to about 1 m of water equivalent.
        
        #SNOW WATER EQUIVALENT
        #Snow water equivalent can be presented in units of kg/m2 or meters of depth of liquid water that would result from melting the snow. SWE is the product of depth and density
        #SWE = depth (m) x density (kg/m3) --> (units: kg/m2)
        #SWE = depth (m) x density (kg/m3) / density of water (kg/m3) --> (units: m)
        
        
        #B_i_mean (dependent on ice density) = average net balance = specific net balance = Dividing total mass balance by glacier surface [m]
        meanThicknessChange = np.mean(pNpGlacier_[:,3]) # [m]  
        #B_i_mean_dot (dependent on ice density) = average net balance rate # [m/a] 
        perYearMeanThicknessChange = meanThicknessChange / dateDiffYears # [m/a] 
        
        #B_mean = average net balance = specific net balance = Dividing total mass balance by glacier surface [kg/m²]
        massChanges = (meanThicknessChange * float(ICE_DENSITY)) / (float(WATER_DENSITY)) # [m w.e./a]
        #B_mean_dot = average net balance rate = pecific net balance [kg/m²]
        perYearMassChanges = massChanges  / dateDiffYears #[m w.e./a] --> Mass balance (per year)!!!!
        
        
        
        #-------------------------------------------------------------------------------
        #Checks of calculation by longer, but standard means according to definitions:
        
        #B_i_mean (dependent on ice density) = average net balance = specific net balance = Dividing total mass balance by glacier surface [m]
        perAreaVolumeChange = sumVolumeChange / glacierSurfaceArea #[m³] / [m²] = [m] --> to check, same as meanThicknessChange
        #B_i_mean_dot (dependent on ice density) = average net balance rate #[m/a]
        perAreaYearVolumeChange = sumVolumeChange / (dateDiffYears * glacierSurfaceArea) #[m/a] --> to check, same as perYearMeanThicknessChange
        
        #B_mean = average net balance = specific net balance = Dividing total mass balance by glacier surface [kg/m²]
        perAreaMassBalance = sumMassBalance / glacierSurfaceArea # [kg] / [m²] ~ [m w.e.] --> to check, same as massChanges
        #B_mean_dot = average net balance rate #[kg]/[m²]*[1/a]
        perAreaYearMassBalance = sumMassBalance / (dateDiffYears * glacierSurfaceArea) #[kg]/[m²]*[1/a] ~ [m w.e./a] --> to check, same as perYearMassChanges
        
        
        
        #-------------------------------------------------------------------------------
        #Error calculation
        #Gaußsche Fehlerfortpflanzungsgesetz beachten: http://de.wikipedia.org/wiki/Fehlerfortpflanzung
        
        #Thickness change uncertainty
        dtmErrorGlacier = math.sqrt(pow(float(dtmError_),2) + pow(float(srtmPenetrationError_),2)) # [m]
        perYearDtmErrorGlacier = float(dtmErrorGlacier) / dateDiffYears #[m/a]
        
        #Glacier mass balance uncertainty
        dtmErrorGlacierWE = (dtmErrorGlacier * float(ICE_DENSITY)) / (float(WATER_DENSITY)) # [m w.e.]
        iceDensityErrorMBWE = (float(meanThicknessChange) * float(ICE_DENSITY_ERROR)) / (float(WATER_DENSITY)) # [m w.e.]
        
        dtmDensityMassError = math.sqrt(pow(dtmErrorGlacierWE,2) + pow(iceDensityErrorMBWE,2)) # [m w.e.]
        perYearDtmDensityMassError = dtmDensityMassError / dateDiffYears # [m w.e.]
        
        
        #---
        #Volume change uncertainty
        sumVolumeChangeError = math.sqrt(pow((float(dtmError_) * glacierSurfaceArea),2) + pow((float(srtmPenetrationError_) * glacierSurfaceArea),2)) #[m³]
        perYearSumVolumeChangeError = sumVolumeChangeError / dateDiffYears #[m³/a]
        
        sumVolumeChangeErrorIce = sumVolumeChangeError * float(ICE_DENSITY) #[m³]*[kg/m³] = [kg]
        iceDensityErrorVCIce = float(sumVolumeChange) * float(ICE_DENSITY_ERROR) #[m³]*[kg/m³] = [kg]
        
        sumMassBalanceErrorIce= math.sqrt(pow(sumVolumeChangeErrorIce,2) + pow(iceDensityErrorVCIce,2)) #[m³]*[kg/m³] = [kg]
        perYearSumMassBalanceErrorIce = sumMassBalanceErrorIce / dateDiffYears #[kg/a]
        
        
        
        #Print Summary  
        #-------------------------------------------------------------------------------
        print("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
        print("GEODETIC GLACIER MASS BALANCE CALCULATION for glacier '" + str(glacierName_) + "':")
        print("Penetration of the SRTM C-band signal with an error of '" + str(srtmPenetrationError_) + "'m considered in calculation!")
        print("\nTime interval calculation values:")
        print(" Old DEM is of date '" + str(pDateOld_) + "', new DEM is of date '" + str(pDateNew_) + "', thereof the time intervall is '" + str(pDateDiff) + "' or '" + str(dateDiffYears) + "' integer years.") 
       
        print("\nGlacier area calculation values:")
        print(" Pixelsize is '" + str(PIXELSIZE) + "'m, thereof the area per pixel is '" + str(areaPerPixel) + "'m². With '" +str(pNpGlacier_.shape[0]) \
              + "' dh-pixel existing, the complete surface area of glacier '" + str(glacierName_) + "' is '" + str(glacierSurfaceArea) + "'m² (~'" + str(round(glacierSurfaceArea/pow(10,6),2)) + "'km²).")
        
        print("\nThickness change calculation:")
        print(" B_i_mean = average net balance = specific net balance = mean thickness change for '" + str(dateDiffYears) + "' years: '"+ str(round(meanThicknessChange,4)) + "'m; (Sum of all thickness change values: '" + str(round(sumThicknessChange,4)) + "'m.)")
        print(" B_i_mean_dot = average net balance rate = specific net balance rate = mean thickness change per year: '" + str(round(perYearMeanThicknessChange,4)) + "'m/a (~'" + str(round(perYearMeanThicknessChange/pow(10,-3),2)) + "'mm/a) ") 
        print(" !!! DTM Error = mean thickness change error: '" + str(round(dtmErrorGlacier,4)) + "'m; Mean thickness change error per year: '" + str(round(perYearDtmErrorGlacier,4)) + "'m/a ")
        
        print("\nVolume change calculation:")
        print("--> For sum of ice thickness difference values of '" + str(round(sumThicknessChange,4)) + "'m and and area per pixel of '" + str(areaPerPixel) + "'m²:")
        print(" B_i = glacier net balance = total net balance = glacier volume change for '" + str(dateDiffYears) + "' years: '" + str(round(sumVolumeChange,4)) + "'m³ (~'" + str(round(sumVolumeChange/pow(10,9),4)) + "'km³, ~'" + str(round(sumVolumeChange/pow(10,6),4)) + "'km³*10^-3).")
        print(" B_i_dot = glacier net balance rate = total net balance rate = glacier volume change per year: '" + str(round(perYearVolumeChange,4)) + "'m³/a (~'" + str(round(perYearVolumeChange/pow(10,9),4)) + "'km³/a, ~'" + str(round(perYearVolumeChange/pow(10,6),4)) + "'km³*10^-3/a).")
        print(" !!! Glacier volume change error: '" + str(round(sumVolumeChangeError,4)) + "'m³ (~'" + str(round(sumVolumeChangeError/pow(10,9),4)) + "'km³, ~'" + str(round(sumVolumeChangeError/pow(10,6),4)) + "'km³*10^-3); Glacier volume change error per year: '" + str(round(perYearSumVolumeChangeError,4)) + \
              "'m³/a (~'" + str(round(perYearSumVolumeChangeError/pow(10,9),4)) + "'km³/a, ~'" + str(round(perYearSumVolumeChangeError/pow(10,6),4)) + "'km³*10^-3/a)")
        print(" B = glacier net balance = total net balance = glacier ice mass change for '" + str(dateDiffYears) + "' years: '" + str(round(sumMassBalance,4)) + "'kg (~'" + str(round(sumMassBalance / pow(10,12),4)) + "'Gt, ~'" + str(round(sumMassBalance / pow(10,9),4)) + "'Gt*10^-3)") #1 kilogram is equal to 1.0E-12 gigatonne.
        print(" B_dot = glacier net balance rate = total net balance rate = glacier ice mass change per year: '" + str(round(perYearSumMassBalance,4)) + "'kg/a (~'" + str(round(perYearSumMassBalance / pow(10,12),4)) + "'Gt/a, ~'" + str(round(perYearSumMassBalance / pow(10,9),4)) + "'Gt*10^-3/a)")
        print(" !!! Glacier ice mass change error: '" + str(round(sumMassBalanceErrorIce,4)) + "'kg (~'" + str(round(sumMassBalanceErrorIce / pow(10,12),4)) + "'Gt, ~'" + str(round(sumMassBalanceErrorIce / pow(10,9),4)) + "'Gt*10^-3); Glacier ice mass change error per year: '" + str(round(perYearSumMassBalanceErrorIce,4)) + \
              "'Gt/a (~'" + str(round(perYearSumMassBalanceErrorIce / pow(10,12),4)) + "'Gt/a, ~'" + str(round(perYearSumMassBalanceErrorIce / pow(10,9),4)) + "'Gt*10^-3/a)")
               
        print("\nMass balance calculation (1 kg/m² of ice ~ 1 m w. e. (with water density of 999.972 kg/m³)):")
        print("--> With mean thickness change of '"+ str(round(meanThicknessChange,4)) + "'m and ice density of '" + str(ICE_DENSITY) + "'kg/m³ divided by water density of '" + str(WATER_DENSITY) + "'kg/m³:")
        print(" B_mean = average net balance = specific net balance = glacier mass changes for '" + str(dateDiffYears) + "' years: '"+ str(round(massChanges,4)) + "'m w. e. ('" + str(round(massChanges * pow(10,3),4)) + "' mm w. e.)")
        print(" B_mean_dot = average net balance rate = specific net balance = glacier mass changes per year: "+ str(round(perYearMassChanges,4)) + "'(m w. e.)/a ('" + str(round(perYearMassChanges * pow(10,3),4)) + "' (mm w. e.)/a)")
        print(" !!! Mass balance error (ice density uncertity: '" + str(ICE_DENSITY_ERROR) + "'kg/m³): '" + str(round(dtmDensityMassError,4)) + "'m w. e.  ('" + str(round(dtmDensityMassError * pow(10,3),4)) + "' mm w. e.) ; Mass balance error per year: '" + str(round(perYearDtmDensityMassError,4)) + "'(m w. e.)/a ('" + str(round(perYearDtmDensityMassError * pow(10,3),4)) + "' (mm w. e.)/a)")
        
        print("\nCross-check with other method:")
        print(" Per area of '"  + str(glacierSurfaceArea) + "'m² for '" + str(dateDiffYears) + "' years: '" + str(round(perAreaVolumeChange,4)) + "'m --> Same as mean thickness change?")
        print(" Per area of '" + str(glacierSurfaceArea) + "'m² and per year: '" + str(round(perAreaYearVolumeChange,4)) + "'m/a (~'" + str(round(perAreaYearVolumeChange/pow(10,-3),2)) + "'mm/a) --> Same as mean thickness change per year?")
        print(" Glacier mass changes for '" + str(dateDiffYears) + "' years: '" + str(round(perAreaMassBalance,4)) + "'kg/m² (~'" + str(round(perAreaMassBalance/pow(10,3),4))+ "'m w. e.) --> ~Same as glacier mass changes?")
        print(" Glacier mass changes per year: '" + str(round(perAreaYearMassBalance,4)) + "'(kg/m²)/a (~'" + str(round(perAreaYearMassBalance/pow(10,3),4)) + "'(m w. e.)/a)  --> ~Same as glacier mass changes per year?")
        
        print("\n--> Summary mass balance calculations for glacier '" + str(glacierName_) + "' ('" + str(pDateOld_) + "' to '" + str(pDateNew_) + "'):")
        print(" Mean glacier thickness change for '" + str(dateDiffYears) + "' years: '"+ str(round(meanThicknessChange,2)) + "'+-'" + str(round(dtmErrorGlacier,2)) +"'m (B_i_mean)")
        print(" Glacier ice mass change for '" + str(dateDiffYears) + "' years: '" + str(round(sumMassBalance / pow(10,9),2)) + "'+-'" + str(round(sumMassBalanceErrorIce / pow(10,9),2)) + "'Gt*10^-3 ('" + str(round(sumMassBalance / pow(10,12),2)) + "'+-'" + str(round(sumMassBalanceErrorIce / pow(10,12),2)) + "'Gt) (B = glacier net balance = total net balance)") 
        print(" Glacier mass balance per year: '" + str(round(perYearMassChanges,2)) + "'+-'" + str(round(perYearDtmDensityMassError,2)) +"'(m w. e.)/a (ice density of'" + str(ICE_DENSITY) + "'+-'" +str(ICE_DENSITY_ERROR) + "'kg/m³ and water density of '" + str(WATER_DENSITY) + "'kg/m³) (B_mean_dot)")
        
        
        
        return pNpGlacier_





    def calculateDtmAccuracy(self, pNpNoGlacier_, dtmName_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        print("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
        print("DTM ACCURACY MEASURES for DTM '" + str(dtmName_) + "':")
        
        
        #Calculating area per pixel and complete valid DTM surface area 
        #-------------------------------------------------------------------------------
        areaPerPixel = pow(float(PIXELSIZE),2) #Area per Pixel, depending on pixel size --> (pixel_size * pixel_size) [m²]
        dtmSurfaceArea = pNpNoGlacier_.shape[0] * areaPerPixel  #Complete surface area of glacier --> (number of amount of dh-pixels * area per pixel) [m²]
        
        print(" Pixelsize is '" + str(PIXELSIZE) + "'m, thereof the area per pixel is '" + str(areaPerPixel) + "'m². With '" +str(pNpNoGlacier_.shape[0]) \
              + "' valid dh-pixel existing, the complete area of the DTM '" + str(dtmName_) + "' is '" + str(dtmSurfaceArea) + "'m² (~'" + str(round(dtmSurfaceArea/pow(10,6),2)) + "'km²).")
        
        
        
        #Calculate proposed accuracy measures for DEMs according to hoehleetal2009 
        #-------------------------------------------------------------------------------
        pDhValues = np.copy(pNpNoGlacier_[:,3])
        medianDh = np.median(pDhValues)
        nmadDh = 1.4826 * np.median(abs(pDhValues - medianDh)) #Calculate Normalized Median Absolute Deviation (NMAD) according to Höhle and Höhle (2009)
        
        stdQuantileDh = stats.mstats.mquantiles(abs(pDhValues), prob=[0.683]) #Use absolute value of quantile according to Höhle and Höhle (2009)
        doubleStdQuantileDh = stats.mstats.mquantiles(abs(pDhValues), prob=[0.95]) #Use absolute value of quantile according to Höhle and Höhle (2009)
        
        print(" DTM accuracy measures handling outliers and non-normal distribution: Normalized Median Absolute Deviation (NMAD) -> dh '" + str(round(nmadDh,2)) + "'m; Median (50% quantile) -> dH: '" + str(round(medianDh,2)) + \
              "'m; 68.3% quantile: -> abs(dH): '" + str(round(stdQuantileDh,2)) + "'m; 95% quantile: -> abs(dH): '" + str(round(doubleStdQuantileDh,2)) + "'m.")
        
        
        #Calculate accuracy measures for DEMs presenting normal distribution of errors
        #-------------------------------------------------------------------------------
        rmseDh = math.sqrt(np.mean(pow(pDhValues,2)))
        #rmseDh = math.sqrt((np.sum(pow(pDhValues,2)))/len(pDhValues)) #RMSE alternative
        meanDh = np.mean(pDhValues)
        stdDh = np.std(pDhValues)
        
        print (" DTM accuracy measures in case of normal distribution of errors: Mean: '" + str(round(meanDh,2)) + "'m; RMSE: '" + str(round(rmseDh,2)) + "'m; Standard deviation: '" + str(round(stdDh,2)) +"'m.") 
        
        
        #Calculate further accuracy measures 
        #-------------------------------------------------------------------------------
        #Find mode value 
        pNpNoGlacierCalcInt = pNpNoGlacier_.astype(int)
        [modeAll], [cntsModeAll] = stats.mode(pNpNoGlacierCalcInt[:,3])#dH
        
        print (" DTM accuracy further measures:  Min: '" + str(round(np.min(pNpNoGlacier_[:,3]),2)) + "'m; Max: '" + str(round(np.max(pNpNoGlacier_[:,3]),2)) \
               + "'m; Integer mode: '" + str(modeAll) + "'m with '" + str(cntsModeAll) + "' counts.") 
        
        
        return round(nmadDh,4)
        
        
        
        
    #def fillNoDataAreasWithZeros(self, pNpGlacier_):
    #    """
    #    Program description: #
    #
    #    INPUT_PARAMETERS:
    #    inputValue_      - 
    #
    #    COMMENTS:
    #    """
    #   
    #    ################TEMP ADAPTION OF NUMPY ARRAY
    #    i = 0
    #    for value in pNpGlacier_[:,3]:
    #        i+=1
    #        if (value<-35) or (value>35):
    #            pNpGlacier_[i-1,3] = 0
    #     
    #    print pNpGlacier_.shape, np.min(pNpGlacier_[:,3]), np.max(pNpGlacier_[:,3]), np.mean(pNpGlacier_[:,3])
    #    
    #    return pNpGlacier_
    #    ############################
        
        
    def __calcMeanDhNodataAbl(self, pNpGlacierAbl_, ela_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        #Sort all entries along the zValues
        #-------------------------------------------------------------------------------
        pNpGlacierAblSortZ = pNpGlacierAbl_[pNpGlacierAbl_[:,2].argsort()] 
        
        
        #Calculate mean values in sorted zValues array for each Z_ABL_INTERVAL section
        #-------------------------------------------------------------------------------
        
        #Initialization
        zLocalMax = pNpGlacierAblSortZ[0,2] + float(Z_ABL_INTERVAL) #First local maximal z Value (deepest z-Value + half offset)
        pListGlacierAblSortZMean = [] #List for first saving mean dh values
        
        cntRow = 0 #Counter all rows
        cntDhValues = 0 #Counter of all dh-Values per section
        dhSumValue = 0 #Sum of all dh-Values per section
        
        
        #Calculate all mean dh-Values per section until final Z-Max = ELA
        for zValue in pNpGlacierAblSortZ[:,2]: #For all dh-Values within a section
            
            dhSumValue = dhSumValue + pNpGlacierAblSortZ[cntRow,3] #Sum of all dh-Values
            cntDhValues += 1
            cntRow += 1
            
            if zValue > zLocalMax: #Z-Value reaches next section
            
                dhSumValueMean = dhSumValue / cntDhValues #Calculate mean per section
                meanSectionZ = zLocalMax - float(Z_ABL_INTERVAL) / 2 #Z-Value in middle of section as closest z-Value to z-Value of nodata-Dh entry
                pListGlacierAblSortZMean.append([meanSectionZ, dhSumValueMean]) #Save local Z-max value of section of mean per section to list
                
                zLocalMax = zLocalMax + float(Z_ABL_INTERVAL) #Next local Z-max value
                dhSumValue = 0
                cntDhValues = 0
        
        if dhSumValue != 0: #If there are remaining values 
            
            dhSumValueMean = dhSumValue / cntDhValues #Calculate mean per section
            meanSectionZ = float(ela_) #((zLocalMax - float(Z_ABL_INTERVAL))+float(ELA))/2
            pListGlacierAblSortZMean.append([meanSectionZ, dhSumValueMean]) #Save local Z-max value of section of mean per section to list
            
        #Save list with dh-Mean values per section to Numpy
        pNpGlacierAblSortZMean = np.empty((len(pListGlacierAblSortZMean), 2), dtype = pNpGlacierAblSortZ.dtype) #zLocalMax, dhSumValueMean
        pNpGlacierAblSortZMean = np.asarray(pListGlacierAblSortZMean)
        
        
        print (" Calculated '"+str(pNpGlacierAblSortZMean.shape[0])+"' regional nodata replacement values for ablation area from '" + str(pNpGlacierAblSortZ[0,2]) + "'m to '" + str(pNpGlacierAblSortZ[-1,2]) \
               + "'m in height intervalls of every '"+str(Z_ABL_INTERVAL)+"' meter.")
        
        
        return  pNpGlacierAblSortZMean #zLocalMax, dhSumValueMean
        
        
        
        
        
        
        
        
        
        
    def __openCsvFile(self, csvFileName_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        
        if not csvFileName_.endswith('.csv'): #Add filename suffix '.csv' if this is missing
            self.csvFileName = self.csvFileName + '.csv'
        
        #csv.reader(csvfile[, dialect='excel'][, fmtparam])
        try:
            pDocCsv = csv.reader(open(csvFileName_, 'r') , dialect= 'excel')
        except:
            raise Exception ("Opening of file '" + str(csvFileName_) + "' failed. Check if it exists and if filename suffix is set.")
       
        return pDocCsv
    
    
    

    def __manualDhNodata(self, pNpGlacierCalc, cntRow, cntManReplc, methodNodata_, ela_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        #Mustag Ata
        #-------------------------------------------------------------------------------
        
        if ('manualDhMustagPLEKH9'  in methodNodata_) or ('manualDhMustagALOSKH9' in methodNodata_):
            
            #Nodata area Kekesayi glacier southern part PLE to KH9
            if ((float(pNpGlacierCalc[cntRow-1,2]) <= float(ela_)) and
                 (float(513000) < float(pNpGlacierCalc[cntRow-1,0]) and float(pNpGlacierCalc[cntRow-1,0]) < float(516500)) and 
                 (float(4229800) < float(pNpGlacierCalc[cntRow-1,1]) and float(pNpGlacierCalc[cntRow-1,1]) < float(4232000))):
                cntManReplc+=1
                pNpGlacierCalc[cntRow-1,3] = float(5)
                
            #Nodata area Kekesayi glacier north-western part PLE to KH9
            if ((float(pNpGlacierCalc[cntRow-1,2]) <= float(ela_)) and
                 (float(512200) < float(pNpGlacierCalc[cntRow-1,0]) and float(pNpGlacierCalc[cntRow-1,0]) < float(515400)) and 
                 (float(4235400) < float(pNpGlacierCalc[cntRow-1,1]) and float(pNpGlacierCalc[cntRow-1,1]) < float(4237600))):
                cntManReplc+=1
                pNpGlacierCalc[cntRow-1,3] = float(-5)
                
            
            #Nodata ablation area south-western glacier G075075E38189N
            if ((float(pNpGlacierCalc[cntRow-1,2]) <= float(ela_)) and
                 (float(505800) < float(pNpGlacierCalc[cntRow-1,0]) and float(pNpGlacierCalc[cntRow-1,0]) < float(507500)) and 
                 (float(4226630) < float(pNpGlacierCalc[cntRow-1,1]) and float(pNpGlacierCalc[cntRow-1,1]) < float(4227100))):
                cntManReplc+=1
                pNpGlacierCalc[cntRow-1,3] = float(-8)
            
            #Nodata positive and negative glacier area south-western glacier Kuosikulake
            if ((float(503400) < float(pNpGlacierCalc[cntRow-1,0]) and float(pNpGlacierCalc[cntRow-1,0]) < float(505000)) and 
                 (float(4226500) < float(pNpGlacierCalc[cntRow-1,1]) and float(pNpGlacierCalc[cntRow-1,1]) < float(4228000))):
                cntManReplc+=1
                if (float(pNpGlacierCalc[cntRow-1,2]) >= float(4640)):
                    pNpGlacierCalc[cntRow-1,3] = float(8)
                else:
                    pNpGlacierCalc[cntRow-1,3] = float(-30)
            if ((float(pNpGlacierCalc[cntRow-1,2]) <= float(ela_)) and
                 (float(505100) < float(pNpGlacierCalc[cntRow-1,0]) and float(pNpGlacierCalc[cntRow-1,0]) < float(507000)) and 
                 (float(4228000) < float(pNpGlacierCalc[cntRow-1,1]) and float(pNpGlacierCalc[cntRow-1,1]) < float(4229500))):
                cntManReplc+=1
                pNpGlacierCalc[cntRow-1,3] = float(5)
            
        
        if ('manualDhMustagPLESRTM'  in methodNodata_) or ('manualDhMustagALOSSRTM' in methodNodata_):
            
            #Nodata ablation area south-western glacier G075075E38189N
            if ((float(pNpGlacierCalc[cntRow-1,2]) <= float(ela_)) and
                 (float(505800) < float(pNpGlacierCalc[cntRow-1,0]) and float(pNpGlacierCalc[cntRow-1,0]) < float(507500)) and 
                 (float(4226630) < float(pNpGlacierCalc[cntRow-1,1]) and float(pNpGlacierCalc[cntRow-1,1]) < float(4227100))):
                cntManReplc+=1
                pNpGlacierCalc[cntRow-1,3] = float(-15)
            
            #Nodata positive and negative glacier area south-western glacier Kuosikulake
            if ((float(503400) < float(pNpGlacierCalc[cntRow-1,0]) and float(pNpGlacierCalc[cntRow-1,0]) < float(505000)) and 
                 (float(4226500) < float(pNpGlacierCalc[cntRow-1,1]) and float(pNpGlacierCalc[cntRow-1,1]) < float(4228000))):
                cntManReplc+=1
                if (float(pNpGlacierCalc[cntRow-1,2]) >= float(4640)):
                    pNpGlacierCalc[cntRow-1,3] = float(15)
                else:
                    pNpGlacierCalc[cntRow-1,3] = float(-10)
            if ((float(pNpGlacierCalc[cntRow-1,2]) <= float(ela_)) and
                 (float(505100) < float(pNpGlacierCalc[cntRow-1,0]) and float(pNpGlacierCalc[cntRow-1,0]) < float(507000)) and 
                 (float(4228000) < float(pNpGlacierCalc[cntRow-1,1]) and float(pNpGlacierCalc[cntRow-1,1]) < float(4229500))):
                cntManReplc+=1
                pNpGlacierCalc[cntRow-1,3] = float(-5)
                
            
            #Nodata area Kekesayi glacier northern part PLE to KH9
            #if ((float(515600) < float(pNpGlacierCalc[cntRow-1,0]) and float(pNpGlacierCalc[cntRow-1,0]) < float(517400)) and 
            #    (float(4235000) < float(pNpGlacierCalc[cntRow-1,1]) and float(pNpGlacierCalc[cntRow-1,1]) < float(4237800))):
            #    cntManReplc+=1
            #    if (float(pNpGlacierCalc[cntRow-1,2]) <= float(ela_)):
            #        pNpGlacierCalc[cntRow-1,3] = float(-20)
            #    elif (float(pNpGlacierCalc[cntRow-1,2]) <= (float(ela_)) and float(pNpGlacierCalc[cntRow-1,2]) <= (float(ela_)+75)):
            #        pNpGlacierCalc[cntRow-1,3] = float(-15)
            #    elif (float(pNpGlacierCalc[cntRow-1,2]) <= (float(ela_)+75) and float(pNpGlacierCalc[cntRow-1,2]) <= (float(ela_)+150)):
            #        pNpGlacierCalc[cntRow-1,3] = float(-10)
            #    elif (float(pNpGlacierCalc[cntRow-1,2]) <= (float(ela_)+150) and float(pNpGlacierCalc[cntRow-1,2]) <= (float(ela_)+250)):
            #        pNpGlacierCalc[cntRow-1,3] = float(-5)
            #    else:
            #        cntManReplc-=1
        
        
        #Gurla Mandhata
        #-------------------------------------------------------------------------------
        
        if "manualDhGurla" in methodNodata_:
           
            #Nodata area at NamuNaNiWestWest Glacier TDX DTM in specific accumulation zone
            if ((float(521870) < float(pNpGlacierCalc[cntRow-1,0]) and float(pNpGlacierCalc[cntRow-1,0]) < float(523870)) and 
                 (float(3367680) < float(pNpGlacierCalc[cntRow-1,1]) and float(pNpGlacierCalc[cntRow-1,1]) < float(3368800))):
                cntManReplc+=1
                pNpGlacierCalc[cntRow-1,3] = float(8)
            
            #Nodata area at NumaNaNi, NamuNaNiWest and NamuNaNiEast Glacier TDX DTM in specific accumulation zone
            if ((float(523900) < float(pNpGlacierCalc[cntRow-1,0]) and float(pNpGlacierCalc[cntRow-1,0]) < float(529270)) and 
                 (float(3367200) < float(pNpGlacierCalc[cntRow-1,1]) and float(pNpGlacierCalc[cntRow-1,1]) < float(3368600))):
                cntManReplc+=1
                pNpGlacierCalc[cntRow-1,3] = float(10)
            
            #Nodata area at NanNumaNaNi Glacier TDX DTM in specific accumulation zone
            if ((float(529700) < float(pNpGlacierCalc[cntRow-1,0]) and float(pNpGlacierCalc[cntRow-1,0]) < float(530800)) and 
                 (float(3366600) < float(pNpGlacierCalc[cntRow-1,1]) and float(pNpGlacierCalc[cntRow-1,1]) < float(3367400))):
                cntManReplc+=1
                pNpGlacierCalc[cntRow-1,3] = float(6)
            
            #Nodata area at Gunala Glacier TDX DTM in specific accumulation zone
            if ((float(533000) < float(pNpGlacierCalc[cntRow-1,0]) and float(pNpGlacierCalc[cntRow-1,0]) < float(535100)) and 
                 (float(3369500) < float(pNpGlacierCalc[cntRow-1,1]) and float(pNpGlacierCalc[cntRow-1,1]) < float(3370750))):
                cntManReplc+=1
                pNpGlacierCalc[cntRow-1,3] = float(10)
            
            #Nodata area at GunalaEastEast Glacier TDX DTM in specific accumulation zone
            if ((float(537300) < float(pNpGlacierCalc[cntRow-1,0]) and float(pNpGlacierCalc[cntRow-1,0]) < float(538000)) and 
                 (float(3374000) < float(pNpGlacierCalc[cntRow-1,1]) and float(pNpGlacierCalc[cntRow-1,1]) < float(3374500))):
                cntManReplc+=1
                pNpGlacierCalc[cntRow-1,3] = float(-8)
        
        
        return pNpGlacierCalc, cntManReplc
    
    
    
    
    def __manualZForElaShift(self, pNpGlacier_, methodManual_, flag_):
        """
        Program description: 

        INPUT_PARAMETERS:
        inputValue_      - 

        COMMENTS:
        """
        
        #Manually adapting z-values at specific areas to manipulate ELA-altitude in regard of later nodata / gap filling operations
        #-------------------------------------------------------------------------------
        pNpGlacierAdapt = np.copy(pNpGlacier_)
        cntRow = 0
        cntManReplc = 0
        
        for zValue in pNpGlacier_[:,2]: #Z
            
            if "manualZMustag" in methodManual_:
                
                #Make ELA higher et Kekesayi glacier northern area
                if ((float(515600) < float(pNpGlacierAdapt[cntRow,0]) and float(pNpGlacierAdapt[cntRow,0]) < float(517400)) and 
                (float(4235000) < float(pNpGlacierAdapt[cntRow,1]) and float(pNpGlacierAdapt[cntRow,1]) < float(4237800))):
                    if flag_ == "before":
                        pNpGlacierAdapt[cntRow,2] -= 190
                        cntManReplc+=1
                    elif flag_ == "after":
                        pNpGlacierAdapt[cntRow,2] += 190
                        cntManReplc+=1
                
            cntRow+=1
        
        print("\nAdapted '" + str(cntManReplc) + "' z-values '" + str(flag_) + "' outlier detection of accumulation / ablation area for method '" + str(methodManual_) + "'.") 
        
        return pNpGlacierAdapt
       
       
       
       
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
    pInterfaceMain = MainInterface() #Initialize
    
    
    #study = "test"
    #study = "halji"
    study = "mustag"
    #study = "gurla"
    
    try:
        
        if study == "test":
            
            ela = 5660 
            
            #pNpNoGlacierTest = pInterfaceMain.testDtmAccuracy()
            #pInterfaceMain.testPlots(pNpNoGlacierTest, "input_DEM", "SRTM", "stable terrain", ela) #pNpGlacier_, demNameNew_, demNameOld_, name_, ela_
            
            pNpGlacierTest = pInterfaceMain.testMassBalance("test glacier", ela)
            pInterfaceMain.testPlots(pNpGlacierTest, "input_DEM", "SRTM-3", "test glacier", ela) #pNpGlacier_, demNameNew_, demNameOld_, name_, ela_
            
            
            
        if study == "halji":
            
            ela = 5660 #Halji ELA for glacier G081388E30289N (Shi, 2008) : 5660m --> Gletscher ID Unterschiede GLIMS zu RGI, G081388E30289N in Glims ist G081393E30265N in RGI
            
            #pNpNoGlacierTest = pInterfaceMain.testDtmAccuracy()
            #pInterfaceMain.testPlots(pNpNoGlacierTest, "input_DEM", "SRTM", "stable terrain", ela) #pNpGlacier_, demNameNew_, demNameOld_, name_, ela_
            
            glacierNamesList = ["G081470E30264N", "G081437E30281N", "G081393E30265N"] #Bei Gletscher IDs existieren Unterschiede von GLIMS zu RGI! RGI ist genauer kartiert, daher GLIMS IDs aus RGI v3 entnommen
            
            for glacierName in glacierNamesList: 
                pNpGlacierDhDem = pInterfaceMain.haljiMassBalance(glacierName, ela)
                pInterfaceMain.testPlots(pNpGlacierDhDem, "PLE", "SRTM-3", glacierName, ela) #pNpGlacier_, demNameNew_, demNameOld_, name_, ela_
            
            
            
        if study == "mustag":
            
            #Glacier specific ELAs
            #Pamir ELA Gardelleetal: 4580+-250m
            #"Muztag Ata Glacier" ELA nach yaoetal2012: 2001/02: 5275m; 2002/03: 5450m; 2005/06: 5750m; 2007/08: 5450m; 2009/10: 5425m --> Mean: (5275+5450+5750+5450+5425)/5 = 5470
            #Mean of all ELAs of Shi(2008) in study area: 5316.25 = ~5316m (Adapted ELA: 5296m for only investigated glaciers; 5285m for all glaciers))
            #elasMustagDictOrig =  {"allglacier": "5316", "kekesayi": "4900", "G075233E38272N": "4800", "G075175E38297N": "4920", "G075153E38307N": "4940", "G075101E38308N": "5120", "kematulejia": "5940", \
            #       "G075060E38278N": "5720", "kalaxiong": "5760", "muztag_ata": "5470", "G075071E38240N": "5660", "kuosikulake": "5510", "G075101E38180N": "5510", \
            #       "kuokuosele": "5390", "G075171E38163N": "5310", "G075195E38150N": "5090", "G075219E38156N": "5020"} #Shi(2008)
            elasMustagDict =  {"allglacier": "6000", "kekesayi": "4900", "G075233E38272N": "4770", "G075175E38297N": "4820", "G075153E38307N": "4940", "G075101E38308N": "4970", "kematulejia": "5940", \
                   "G075084E38279N": "5940", "G075060E38278N": "5720", "kalaxiong": "5460", "muztagata": "5470", "G075071E38240N": "5460", "kuosikulake": "5410", "G075075E38189N": "5410", "G075101E38180N": "5510", \
                   "kuokuosele": "5190", "G075171E38163N": "5110", "G075195E38150N": "5090", "G075219E38156N": "5020"} #Shi(2008), modified according to difference images
            
            #Correction value for SRTM C-band radar penetration, do NOT apply on debris covered glaciers and differencing from optical DEMs
            #Mustag Ata following gardelleetal2013 SRTM penetration for Pamir: 1.8+-1.5m 
            #Mustag Ata following kaeaebetal2012, SRTM penetration as mean of (HK+KK+JK)/3 regions: glacier 2.1+-0.6m; firn/snow: 4.3+-0.9m; clean ice: 1.5+-0.7m
            
            #Do not average uncertainties, use the highest uncertainty! Thereof 0.9m instead of 0.6m
            srtmPenetration = [-4.3, -1.5, 0.9] # penetrationDepthAcc, penetrationDepthAbl, srtmPenetrationError --> consider math sign!
            
            debrisBoundary = {"allglacier": "0", "kekesayi": "4600", "G075233E38272N": "4570", "G075175E38297N": "4420", "G075153E38307N": "4540", "G075101E38308N": "4870", "kematulejia": "0", \
                   "G075084E38279N": "0", "G075060E38278N": "0", "kalaxiong": "4660", "muztagata": "0", "G075071E38240N": "0", "kuosikulake": "4810", "G075075E38189N": "0", "G075101E38180N": "0", \
                   "kuokuosele": "4790", "G075171E38163N": "5010", "G075195E38150N": "4990", "G075219E38156N": "0"} #Height of debris boundary, glacier is debris covered below this altitude (0 for clean-ice glaciers)
            
            #DEM difference image combinations
            #dhDemNamesList = ["ple1_m_srtm"]
            #dhDemNamesList = ["ple1_m_srtm", "alos2_m_srtm"]
            #dhDemNamesList = ["ple1_m_srtm", "ple1_m_kh91", "srtm_m_kh91"]
            #dhDemNamesList = ["alos2_m_srtm", "alos2_m_kh91", "srtm_m_kh91"]
            #dhDemNamesList = ["ple1_m_alos2", "ple1_m_srtm", "alos2_m_srtm", "ple1_m_kh91", "alos2_m_kh91", "srtm_m_kh91"]
            dhDemNamesList = ["ple1_m_alos1", "ple1_m_alos2", "ple1_m_srtm", "alos1_m_srtm", "alos2_m_srtm", "ple1_m_kh91", "alos1_m_kh91", "alos2_m_kh91", "srtm_m_kh91"]
            
            
            #DEM accuracy
            for dhDemName in dhDemNamesList: 
                
                pNpNoGlacierDhDem, dtmAccuracyNmad = pInterfaceMain.mustagDtmAccuracy(dhDemName)
                pInterfaceMain.testPlots(pNpNoGlacierDhDem, dhDemName.split('_')[0], dhDemName.split('_')[-1], "stable_terrain", elasMustagDict["allglacier"]) #pNpGlacier_, demNameNew_, demNameOld_, name_, ela_
                
                #DEM mass balance for specific glaciers
                #glacierNamesList = ["muztagata", "G075071E38240N"] #only ["ple1_m_srtm"]
                #glacierNamesList = ["G075084E38279N"] #only ["ple1_m_srtm", "alos2_m_srtm"]
                #glacierNamesList = ["G075101E38308N", "kematulejia", "kalaxiong"] #only ["ple1_m_srtm", "ple1_m_kh91", "srtm_m_kh91"]
                #glacierNamesList = ["kuosikulake"] #only ["alos2_m_srtm", "alos2_m_kh91", "srtm_m_kh91"]
                #glacierNamesList = ["G075233E38272N", "G075175E38297N", "G075075E38189N", "kuokuosele", "G075171E38163N"] #only ["ple1_m_alos2", "ple1_m_srtm", "alos2_m_srtm", "ple1_m_kh91", "alos2_m_kh91", "srtm_m_kh91"]
                #glacierNamesList = ["kekesayi", "allglacier"] #only ["ple1_m_alos1", "ple1_m_alos2", "ple1_m_srtm", "alos1_m_srtm", "alos2_m_srtm", "ple1_m_kh91", "alos1_m_kh91", "alos2_m_kh91", "srtm_m_kh91"]
                                
                glacierNamesList = ["allglacier"]
                
                #Create empty numpy array which will contain the sum of data for all glaciers from the DEM
                #pNpGlacierAll = np.ndarray(shape = (0, 4))
                
                
                for glacierName in glacierNamesList: 
                    
                    pNpGlacierDhDem = pInterfaceMain.mustagMassBalance(glacierName, dhDemName, dtmAccuracyNmad, elasMustagDict[glacierName], srtmPenetration, debrisBoundary)
                    pInterfaceMain.testPlots(pNpGlacierDhDem, dhDemName.split('_')[0], dhDemName.split('_')[-1], glacierName, elasMustagDict[glacierName]) #pNpGlacier_, demNameNew_, demNameOld_, name_, ela_
                    
                    #pNpGlacierAll = np.vstack((pNpGlacierAll, pNpGlacierDhDem)) #Sum up all glacier data from DEM
                
                #Plot again for all glacier data
                #pInterfaceMain.testPlots(pNpGlacierAll, dhDemName.split('_')[0], dhDemName.split('_')[-1], "all DEMs " + glacierName, 5316) #pNpGlacier_, demNameNew_, demNameOld_, name_, ela_
            
            
            
            
        if study == "gurla":
            
            
            #Glacier specific ELAs
            #West Nepal ELA nach Gardelletal: 5590+-138m
            #"Gurla Glacier" ELA nach yaoetal2012: 2004/06: 6300m; 2007/08: 6480; 2008/09: 6200; 2009/10: 6150 --> Mean: (2004 - 2009: 6300 + 6300 + 6480 + 6200 + 6150) / 5 = 6286 
            #Mean of all ELAs of Shi(2008) in study area: 6024.25 = ~6025m
            #elasGurlaDictOrig =  {"NamuNaNiWestWest": "5950", "NamuNaNiWest": "6200", "NamuNaNi": "6050", "NamuNaNiEast": "6230", "GurlaWest": "6160", \
            #    "Gurla": "6075", "NanNamuNaNi": "5915", "Gunala": "6070", "GurlaEast": "", "GurlaSouthEast": "5800" , "GunalaEast": "5960", "GunalaEastEast": ""} #Shi(2008)
            elasGurlaDict = {"allglacier": "6025", "NamuNaNiWestWest": "5950", "NamuNaNiWest": "6200", "NamuNaNi": "6050", "NamuNaNiEast": "6000", "GurlaWest": "6160", \
                "Gurla": "6286", "NanNamuNaNi": "5850", "Gunala": "6070", "GurlaEast": "6075", "GurlaSouthEast": "5900" , "GunalaEast": "5960", "GunalaEastEast": "5825"} #Shi(2008), modified according to difference images
            
            
            #Correction value for SRTM C-band radar penetration, do NOT apply on debris covered glaciers and differencing from optical DEMs
            #Gurla Mandhata following gardelleetal2013 SRTM penetration for Everest (West Nepal): -1.4m, uncertainty +-1.5m but adapted to +-1.0m; Further employed penetration depth: -2.8m  
            #Gurla Mandhata following kaeaebetal2012, SRTM penetration for HP study site: glacier 1.5+0.4m; firn/snow: 2.3+-0.6m; clean ice: 1.7+-0.6m
            srtmPenetration = [-2.3, -1.7, 0.6] # penetrationDepthAcc, penetrationDepthAbl, srtmPenetrationError --> consider math sign!
            debrisBoundary = {"allglacier": "0", "NamuNaNiWestWest": "0", "NamuNaNiWest": "0", "NamuNaNi": "0", "NamuNaNiEast": "0", "GurlaWest": "0", \
                "Gurla": "0", "NanNamuNaNi": "5850", "Gunala": "0", "GurlaEast": "0", "GurlaSouthEast": "0" , "GunalaEast": "0", "GunalaEastEast": "0"} #Height of debris boundary, glacier is debris covered below this altitude (0 for clean-ice glaciers)
            
            
            #DEM accuracy #["tdx", "ple", "spot"]
            pNpNoGlacierDhTdx, dtmAccuracyNmadTdx = pInterfaceMain.gurlaDtmAccuracy("tdx")
            pNpNoGlacierDhPle, dtmAccuracyNmadPle = pInterfaceMain.gurlaDtmAccuracy("ple")
            pNpNoGlacierDhSpot, dtmAccuracyNmadSpot = pInterfaceMain.gurlaDtmAccuracy("spot")
            
            pInterfaceMain.gurlaPlots(pNpNoGlacierDhTdx, pNpNoGlacierDhPle, "TDX", "PLE", "stable terrain", elasGurlaDict["allglacier"]) #pNpGlacier1_, pNpGlacier2_, demName1_, demName2_, name_, ela_
            
            
            #DEM mass balance for specific glaciers
            #glacierNamesList = ["Gurla"]
            glacierNamesList = ["NamuNaNiWestWest", "NamuNaNiWest", "NamuNaNi", "NamuNaNiEast", "GurlaWest", \
                "Gurla", "NanNamuNaNi", "Gunala", "GurlaEast", "GurlaSouthEast" , "GunalaEast", "GunalaEastEast"]
            glacierNamesListDict = {"NamuNaNiWestWest": "G081237E30450N", "NamuNaNiWest": "G081255E30450N", "NamuNaNi": "G081267E30456N", "NamuNaNiEast": "G081296E30455N", "GurlaWest": "G081293E30477N", \
                "Gurla": "G081317E30454N", "NanNamuNaNi": "G081307E30424N", "Gunala": "G081351E30471N", "GurlaEast": "G081365E30442N", "GurlaSouthEast": "G081349E30428N", "GunalaEast": "G081379E30472N", "GunalaEastEast": "G081387E30495N"}
            
            
            #Create empty numpy arrays which will contain the sum of data for all glaciers from DEM #1 and DEM #2
            pNpGlacierAll1 = np.ndarray(shape = (0, 4))
            pNpGlacierAll2 = np.ndarray(shape = (0, 4))
            
            
            #For each glacier calculate mass balance and make plots
            for glacierName in glacierNamesList: 
                
                pNpGlacierDhTdx = pInterfaceMain.gurlaMassBalance(glacierName, "tdx", dtmAccuracyNmadTdx, elasGurlaDict[glacierName], srtmPenetration, debrisBoundary)
                pNpGlacierDhPle = pInterfaceMain.gurlaMassBalance(glacierName, "ple", dtmAccuracyNmadPle, elasGurlaDict[glacierName], srtmPenetration, debrisBoundary)
                pNpGlacierDhSpot = pInterfaceMain.gurlaMassBalance(glacierName, "spot", dtmAccuracyNmadSpot, elasGurlaDict[glacierName], srtmPenetration, debrisBoundary)
                                
                pNpGlacierAll1 = np.vstack((pNpGlacierAll1, pNpGlacierDhTdx)) #Sum up all glacier data from DEM #1
                pNpGlacierAll2 = np.vstack((pNpGlacierAll2, pNpGlacierDhPle)) #Sum up all glacier data from DEM #2
                
                pInterfaceMain.gurlaPlots(pNpGlacierDhTdx, pNpGlacierDhPle, "TDX", "PLE", "glacier "+glacierNamesListDict[glacierName], elasGurlaDict[glacierName]) #pNpGlacier1_, pNpGlacier2_, demName1_, demName2_, name_, ela_
                
            #Plot again for all glacier data
            pInterfaceMain.gurlaPlots(pNpGlacierAll1, pNpGlacierAll2, "TDX", "PLE", "all glaciers", 6025) #pNpGlacier1_, pNpGlacier2_, demName1_, demName2_, name_, ela_
                
                
                
    except Exception: #If Exception occurred in this module or all connected sub-modules
        print('Exception error occurred (see below)! ')
        raise

    finally:
        print("Finished. Total processing time [s]: '" + str(time.time() - startTime) + "'.")
        #plt.show() #Plot all 
        print("_____________________________________________________________________________________________")
      

if __name__ == "__main__":
    main()
