#!/usr/bin/env python
#!/usr/anaconda3/envs/my_env_py3/bin python
#This script is to semi-automate basic data reduction of TANSPEC low-resolution spectroscopic data.
#Highlight ---> No 'IRAF' task has been used. Solely in Python.
#---------------------------------------------supriyo+indiajoe+varghesereji
import os
import os.path
import glob
import argparse
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import sys
import numpy as np
import warnings
import re
import shlex
import readline
import shutil
import subprocess
import time
import pkgutil
from astropy.stats import mad_std
from scipy import signal, ndimage

import datetime
from astropy.time import Time
from astropy.io import ascii
import astropy.table as table 
import astroscrappy
import configparser
from .libs import imarith
import WavelengthCalibrationTool.reidentify as reident
import SpectrumExtractor.spectrum_extractor as specextractor

import warnings
warnings.filterwarnings("ignore")


from .reduce_specXD import fxn
from .reduce_specXD import WriteSpecToFitsFile
from .reduce_specXD import SpecMake
from .reduce_specXD import MakeMasterFlat
from .reduce_specXD import DivideSmoothGradient
from .reduce_specXD import RemoveCosmicRays
from .reduce_specXD import SubtractSmoothGradient
from .reduce_specXD import CombDith_Flatcorr_subrout
from .reduce_specXD import Manual_InstectCalframes_subrout
from .reduce_specXD import DitherDetection
from .reduce_specXD import ImagePlot
from .reduce_specXD import Manual_InspectObj_subrout
from .reduce_specXD import is_number
from .reduce_specXD import ListOfBiasesForGivenImg
from .reduce_specXD import SelectionofFrames_subrout
from .reduce_specXD import CreatLogFilesFromFits_subrout
from .reduce_specXD import LoadDirectories
from .reduce_specXD import KeyboardInterrupt_handler
from .reduce_specXD import InitialTest


def LrSpectralExtraction_subrout(PC):
    """ Extracts spectra from 2d image and also applies wavelength calibration """

    directories = LoadDirectories(PC,CONF=False)
    for night in directories:
        PC.currentnight = night # Updating the night directory for using PC.GetFullPath()
        print('Working on night: '+night)
        try:
            #First we load the Spectrastoextract_Lamp_BandFilter.txt
            SpecslistFILE=open(PC.GetFullPath('Spectrastoextract_Lamp_BandFilter.txt'),'r')
        except IOError:
            print('Could not find Spectrastoextract_Lamp_BandFilter.txt in this directory. Hence skipping %s directory'%(night))
            continue
        #Image to Lamp file dictionary
        Img2Lamp = []  #They are first born as lists...
        Img2NeLamp = []
        Img2Filt = []
        AllImglist = []
        for line in SpecslistFILE:
            line = shlex.split(line.rstrip())
            Img2Lamp.append((line[0],line[1]))
            Img2NeLamp.append((line[0],line[2]))
            Img2Filt.append((line[0],line[3]))
            AllImglist.append(line[0])
        SpecslistFILE.close()
        Img2Lamp = dict(Img2Lamp)
        Img2NeLamp = dict(Img2NeLamp)
        Img2Filt = dict(Img2Filt)
        #Initialising an empty dictionary for storing spectra of each filter
        Filt2finalspecs = dict()
        for filt in set(Img2Filt.values()) : Filt2finalspecs[filt] = []

        OutputObjSpecWlCaliList = [] # Initialising an empty array for storing wl calibrated spectra
        for img in AllImglist:
            print("Working on image "+night+" "+img)
            leftover = glob.glob(os.path.splitext(img)[0]+'_*.fits')  #Clean up of some previous attempts if any..
            leftover += glob.glob(os.path.splitext(img)[0]+'.ms.fits')
            if len(leftover) > 0 :
                if not os.path.isdir('Leftover'): os.makedirs('Leftover') 
                for lft in leftover :
                    shutil.move(lft,'Leftover/')

            # Spectrum extraction for object
            SpectrumFile = PC.GetFullPath(img)
            OutputObjSpec = SpectrumFile.rstrip('.fits')+'.ms.fits'
            #ReFitApertureInXD = 'p1'
            APERTUREWINDOW = PC.APERTUREWINDOW
            BKGWINDOWS = PC.BKGWINDOWS            
                        
            #finding the configuration file for spectrum extraction as well as ContinuumFile and ApertureLabel as references
            pkgpath = os.path.split(pkgutil.get_loader('pyTANSPEC').get_filename())[0]
            ContinuumFile = os.path.join(pkgpath,'data/Avg_continuum2_lr_s4.0_hg.fits')
            ApertureLabel = os.path.join(pkgpath,'data/ApLabel_bg_removed_star.npy')
            ApertureTraceFilename = os.path.join(pkgpath,'data/TANSPEC_lr_trace.pkl')
            ConfigFileSpecExt = os.path.join(pkgpath,'config/spectrum_extractor_TANSPEC_lr.config') #Default config file
            
            #Finding if any 'spectrum_extractor_TANSPEC.config' file exist in working directory
            ConfigFileSpecExt_exists = os.path.exists(os.path.join(PC.RAWDATADIR,PC.SPECCONFIGFILE))

            if ConfigFileSpecExt_exists:
                ConfigFileSpecExt = PC.SPECCONFIGFILE #customized config file if the users want
                print('\n \033[1m' + 'Uses customized config file at {} users for spectrum extraction'.format(os.path.join(PC.RAWDATADIR,PC.SPECCONFIGFILE))  + '\033[0m')
            else:
                ConfigFileSpecExt = ConfigFileSpecExt
                print('\n \033[1m' + 'Uses default config file (https://github.com/astrosupriyo/pyTANSPEC/blob/main/pyTANSPEC/config/spectrum_extractor_TANSPEC_lr.config) for spectrum extraction' + '\033[0m')
            
            #creating a new configuration file 
            config = configparser.ConfigParser()
            config.optionxform = str  # preserve the Case sensitivity of keys
            config.read(ConfigFileSpecExt)
            
            #set these two variables, which are coming from the main TANSPEC confg file, into the spectrum_extractor_TANSPEC.config
            config.set("tracing_settings","ContinuumFile",str(ContinuumFile))
            config.set("tracing_settings","ApertureLabel",str(ApertureLabel))
            config.set("tracing_settings","ApertureTraceFilename",str(ApertureTraceFilename))
            config.set("extraction_settings","ApertureWindow",str(APERTUREWINDOW))
            config.set("extraction_settings","BkgWindows",str(BKGWINDOWS))
            
            #write new config file into the output directory.
            new_config_file_star = SpectrumFile.rstrip('.fits')+'.config'
            with open(new_config_file_star, 'w') as configfile:
                config.write(configfile)            
            

            OutputObjSpec, Avg_XD_shift, PixDomain = specextractor.main([SpectrumFile, new_config_file_star, OutputObjSpec])
            
            #plot the spectra for the first look
            hdulist = fits.open(OutputObjSpec)
            #Scale factor to scale the spectra
            ScaleFac = (APERTUREWINDOW[1] - APERTUREWINDOW[0]) / (BKGWINDOWS[1][1] - BKGWINDOWS[1][0])
            SigWithBkg = hdulist[0].data
            bkg = np.mean([hdulist[1].data, hdulist[2].data], axis=0)
            SigWithoutBkg = SigWithBkg - (bkg*ScaleFac)
            
            plt.plot(SigWithoutBkg.flatten(), alpha=0.5)
            plt.title("Extracted Spectra from all {} orders".format((SigWithBkg.shape)[0]))
            plt.xlabel("Flattened pixels")
            plt.ylabel("Counts")
            plt.show()

            # Writing a configuration file for the Lamp
            ReFitApertureInXD = [tuple(Avg_XD_shift), tuple(PixDomain)]
            config.set("tracing_settings","ReFitApertureInXD",str(ReFitApertureInXD))
                       
            #write new config file into the output directory.
            new_config_file_lamp = SpectrumFile.rstrip('.fits')+'.Lamp.config'
            with open(new_config_file_lamp, 'w') as configfile:
                config.write(configfile)
                
             # Spectrum extraction for Ar lamp             
            LampFile = os.path.join(PC.RAWDATADIR,night,Img2Lamp[img])
            OutputArLampSpec = os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,img[:-5]+'_arc1.fits')            
            OutputArLampSpec, Avg_XD_shift, PixDomain = specextractor.main([LampFile, new_config_file_lamp,           
                                                                                                                                OutputArLampSpec])
             
            # Spectrum extraction for Ne lamp             
            NeLampFile = os.path.join(PC.RAWDATADIR,night,Img2NeLamp[img])
            OutputNeLampSpec = os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,img[:-5]+'_arc2.fits')
            OutputNeLampSpec, Avg_XD_shift, PixDomain = specextractor.main([NeLampFile, new_config_file_lamp, 
                                                                                                                                 OutputNeLampSpec])
               
            #Combine Ar (7 orders from the redder end) and Ne (3 orders from the bluer end) lamps
            hdularc1 = fits.open(OutputArLampSpec)
            hdularc2 = fits.open(OutputNeLampSpec)
            
            #hdularc1data = hdularc1[0].data
            #hdularc2data = hdularc2[0].data
            hdularc1data = fits.getdata(OutputArLampSpec, ext=0)
            hdularc1header = fits.getheader(OutputArLampSpec, ext=0)
            hdularc2data = fits.getdata(OutputNeLampSpec, ext=0)
            hdularc2header = fits.getheader(OutputNeLampSpec, ext=0)
            #lamphdr = hdularc2[0].header
            hdularc1_slit = hdularc1header['Slit']
            hdularc2_slit = hdularc2header['Slit']

            if PC.INSTRUMENT == 'TANSPEC':
                if hdularc1_slit == hdularc2_slit: # matching slits for Ar and Ne lamps
                    slit = hdularc1_slit 
                else:
                    sys.exit("Slit width of lamps does not match: check header of both lamps. ")
                for i in range(3):
                    hdularc1data[i] = hdularc2data[i] # 1st 3-orders (from bluer end) of Ar replaced by Ne
            else:
                   sys.exit("Script is for TANSPEC data reduction")     

#            if PC.INSTRUMENT == 'TANSPEC':
#                for i in range(3):
#                    hdularc1data[i] = hdularc2data[i] # 1st 3-orders (from bluer end) of Ar replaced by Ne 
            
            hdularc1[0].data = hdularc1data
            OutputCombLampSpec = OutputArLampSpec[:-9]+'combarc.fits'
            hdularc1.writeto(OutputCombLampSpec)

            #WL calibration
            pkgpath = os.path.split(pkgutil.get_loader('pyTANSPEC').get_filename())[0]
            RefDispTableFile = os.path.join(pkgpath,'data','LAMPIDENTDIR',slit,'tanspecArNe' + '{}' + '.txt')
            OutDispTableFile =  os.path.splitext(OutputCombLampSpec)[0] + '.OutDispTableFile' + '{}'
            OutputWavlFile =  os.path.splitext(OutputCombLampSpec)[0] + '.OutputWavlFile' + '{}' + '.npy'
  
            ModelForDispersion =  PC.WLFITFUNC
            _ = reident.main([OutputCombLampSpec, RefDispTableFile, OutDispTableFile, "--OutputWavlFile", OutputWavlFile,
                                         "--ModelForDispersion", ModelForDispersion, "--SavePlots", "--StackOrders"])            
            
            AllOutWlSolFile = OutputWavlFile.format('all')
            OutputObjSpechdul = fits.open(OutputObjSpec)
            AllOutWlSol = np.load(AllOutWlSolFile)
            #to make wavelength data as ImageHdu
            AllOutWlSolImgHDU = fits.ImageHDU(AllOutWlSol, name = 'Wavelength')
            OutputObjSpechdul.append(AllOutWlSolImgHDU)
            
            #wl calibrated spectra
            OutputObjSpecWlCali =  OutputObjSpec.rstrip('fits') + 'wlc.fits'
            print(OutputObjSpecWlCali)
            OutputObjSpechdul.writeto(OutputObjSpecWlCali)            
            OutputObjSpecWlCaliList.append(OutputObjSpecWlCali)                        
       
        ## to combine images or not
        N = len(OutputObjSpecWlCaliList)

        if PC.SCOMBINE == 'YES' and N> 1:
            OutputObjSpecWlCaliFinal = OutputObjSpecWlCaliList[0].rstrip('.fits') + '.final' + '.avg.fits'
            OutputObjSpecWlCaliFinalhdul = SpecMake(OutputObjSpecWlCaliList, method = 'mean',ScaleF=ScaleFac)
            OutputObjSpecWlCaliFinalhdul.writeto(OutputObjSpecWlCaliFinal, overwrite=True)
        elif PC.SCOMBINE == 'NO' and N> 1:
            for i in range(N):
                OutputObjSpecWlCaliFinal = OutputObjSpecWlCaliList[i].rstrip('.fits') + '.final.fits'
                OutputObjSpecWlCaliFinalhdul = SpecMake([OutputObjSpecWlCaliList[i]], method = None,ScaleF=ScaleFac)
                OutputObjSpecWlCaliFinalhdul.writeto(OutputObjSpecWlCaliFinal, overwrite=True)
        elif N == 1:
            OutputObjSpecWlCaliFinal = OutputObjSpecWlCaliList[0].rstrip('.fits') + '.final.fits'
            OutputObjSpecWlCaliFinalhdul = SpecMake(OutputObjSpecWlCaliList, method = None, ScaleF=ScaleFac)
            OutputObjSpecWlCaliFinalhdul.writeto(OutputObjSpecWlCaliFinal, overwrite=True)
        else:
            raise NotImplementedError('Unknown combine {0}'.format(PC.SCOMBINE))

    print('All nights over...')



class PipelineConfig(object):
    """ This class is just a collection of variables required to run the pipeline """
    def __init__(self,ConfigFilename = None):
        if ConfigFilename is not None:
            self.LoadFromFile(ConfigFilename)

        self.RAWDATADIR = os.getcwd()
        # variables for helper functions
        self._FrameDatabaseInitialised = False
        self.currentnight = None

    def LoadFromFile(self,ConfigFilename):
        """ Loads the configuration from the input file """
        with open(ConfigFilename,'r') as configfile:
            for con in configfile:
                con = con.rstrip()
                if len(con.split()) >= 2 :
                    
                    if con.split()[0] == "INSTRUMENT=" :
                        self.INSTRUMENT = con.split()[1]
                    elif con.split()[0] == "TODO=" :
                        self.TODO = con.split()[1]
                    elif con.split()[0] == "OUTPUTDIR=" :
                        self.OUTDIR = con.split()[1]
                    #elif con.split()[0] == "VERBOSE=" :
                        #self.VER = con.split()[1]
                    elif con.split()[0] == "TEXTEDITOR=" :
                        self.TEXTEDITOR = shlex.split(con)[1]
                    #elif con.split()[0] == "OVERSCAN=" :
                    #    self.OVERSCAN = con.split()[1]
                    elif con.split()[0] == "Selected_Window_Dither_Find=" :
                        self.WindowSelection = (int(con.split()[1]), int(con.split()[2]))
                    
                    elif con.split()[0] == "COMBINEIMGS=" :
                        self.IMGCOMBINE = con.split()[1][0].upper()
                    #elif con.split()[0] == "IMGCOMBMETHOD=" :
                   #     self.IMGCOMBMETHOD = con.split()[1]
                    #elif con.split()[0] == "FLATSTATSECTION=" :
                    #    self.FLATSTATSECTION = con.split()[1]
                    #elif con.split()[0] == "USEALLFLATS=" :
                    #    self.USEALLFLATS = con.split()[1][0].upper()

                    #elif con.split()[0] == "SEPARATE_SKY=" :
                    #    self.SEPARATESKY = con.split()[1][0].upper()
                    elif con.split()[0] == "REMOVE_COSMIC=" :
                        self.REMOVECOSMIC = con.split()[1][0].upper()

                    elif con.split()[0] == "GRADIENT_REMOVE=" :
                        self.GRADREMOVE = con.split()[1][0].upper()
                    elif con.split()[0] == "GRAD_FILT_SIZE=" :
                        self.MEDSMOOTHSIZE = (int(con.split()[1]), int(con.split()[2]))

                    #elif con.split()[0] == "BPMPHOTO=" :
                        #self.PhotBadPixelMaskName = con.split()[1]
                    #elif con.split()[0] == "BPMSPEC=" :
                        #self.SpecBadPixelMaskName = con.split()[1]
        
                    elif con.split()[0] == "THRESHOLD=" :
                        self.THRESHOLD = float(con.split()[1])
                    elif con.split()[0] == "EPADU=" :
                        self.EPADU = float(con.split()[1])
                    elif con.split()[0] == "READNOISE=" :
                        self.READNOISE = float(con.split()[1])
                    elif con.split()[0] == "DATAMAX=" :
                        self.DATAMAX = con.split()[1]
                    elif con.split()[0] == "XYMATCHMIN=" :
                        self.XYMATCHMIN = int(con.split()[1])
            
                    elif con.split()[0] == "APERTURE=" :
                        self.APERTURE = con.split()[1]
                    elif con.split()[0] == "ANNULUS=" :
                        self.ANNULUS = con.split()[1]
                    elif con.split()[0] == "DANNULUS=" :
                        self.DANNULUS = con.split()[1]

                    elif con.split()[0] == "EXPTIME=" :
                        self.EXPTIMEHDR = con.split()[1]
                    elif con.split()[0] == "FILTER=" :
                        self.FILTERHDR = con.split()[1]
                    elif con.split()[0] == "GRATING=" :
                        self.GRATINGHDR = con.split()[1]
                    elif con.split()[0] == "SLIT=" : 
                        self.SLITHDR = con.split()[1]
                    elif con.split()[0] == "CALMIR=" : 
                        self.LAMPHDR = con.split()[1]
                    elif con.split()[0] == "GAIN=" : 
                        self.GAINHDR = con.split()[1]          
                    elif con.split()[0] == "ARGONL=" : 
                        self.ARGONLHDR = con.split()[1]  
                    elif con.split()[0] == "NEONL=" : 
                        self.NEONLHDR = con.split()[1]
                    elif con.split()[0] == "CONT1L=" : 
                        self.CONT1LHDR = con.split()[1]
                    elif con.split()[0] == "CONT2L=" : 
                        self.CONT2LHDR = con.split()[1]
                            

                    elif con.split()[0] == "UT=" :
                        self.UTHDR = con.split()[1]
                    elif con.split()[0] == "ODATE=" :
                        self.DATEHDR = con.split()[1]

                    elif con.split()[0] == "OBJECT=" : 
                        self.OBJECTHDR = con.split()[1]
                    elif con.split()[0] == "COMMENT=" :
                        self.COMMENTHDR = con.split()[1]
                    elif con.split()[0] == "XSTART=" : 
                        self.XSTARTHDR = con.split()[1]
                    elif con.split()[0] == "XEND=" : 
                        self.XENDHDR = con.split()[1]
                    elif con.split()[0] == "YSTART=" : 
                        self.YSTARTHDR = con.split()[1]
                    elif con.split()[0] == "YEND=" : 
                        self.YENDHDR = con.split()[1]
                    elif con.split()[0] == "XBIN=" : 
                        self.XBINHDR = con.split()[1]
                    elif con.split()[0] == "YBIN=" : 
                        self.YBINHDR = con.split()[1]

                    elif con.split()[0] == "RA_HDR=" : 
                        self.RAHDR = con.split()[1]
                    elif con.split()[0] == "DEC_HDR=" :
                        self.DECHDR = con.split()[1]


                    elif con.split()[0] == "OUTPUT=" :
                        self.OUTPUTFILE = con.split()[1]
                    #elif con.split()[0] == "BACKUP=" :
                        #self.BACKUPDIR = con.split()[1]

                    #elif con.split()[0] == "CONVOLVEIMG=" :
                     #   self.CONVOLVEIMG = con.split()[1]
                    #elif con.split()[0] == "DOPSF=" :
                    #    self.DOPSF = con.split()[1]

                    #elif con.split()[0] == "LAMPDIRECTORY=" :
                    #    self.LAMPREPODIR = con.split()[1]
                    elif con.split()[0] == "SPECAPERTURE=" :
                        self.SPECAPERTURE = con.split()[1]
                    elif con.split()[0] == "SPECAPERTURE_LLIMIT=" :
                        self.SPECAPERTURE_LLIMIT = con.split()[1]
                    elif con.split()[0] == "SPECAPERTURE_ULIMIT=" :
                        self.SPECAPERTURE_ULIMIT = con.split()[1]


                    elif con.split()[0] == "BACKGROUND=" :
                        self.BACKGROUND = con.split()[1]
                    elif con.split()[0] == "TRACEFUNC=" :
                        self.TRACEFUNC = con.split()[1]
                    elif con.split()[0] == "TRACEORDER=" :
                        self.TRACEORDER = con.split()[1]
                    elif con.split()[0] == "NORMFUNC=" :  # Not used yet
                        self.NORMFUNC = con.split()[1]
                    elif con.split()[0] == "NORMORDER=" :  # Not used yet
                        self.NORMORDER = con.split()[1]
                    elif con.split()[0] == "SCOMBINE=" :
                        self.SCOMBINE = con.split()[1]
                    elif con.split()[0] == "DISPAXIS=" :
                        self.DISPAXIS = con.split()[1]
                    elif con.split()[0] == "REMOVE_CONTINUUM_GRAD=" :
                        self.CONTINUUMGRADREMOVE = con.split()[1][0].upper()
                    elif con.split()[0] == "CONT_GRAD_FILT_SIZE=" :
                        self.DVDMEDSMOOTHSIZE = (int(con.split()[1]), int(con.split()[2]))
                    elif con.split()[0] == "ApertureWindow=" :    
                        self.APERTUREWINDOW = (int(con.split()[1]), int(con.split()[2]))
                    elif con.split()[0] == "BkgWindows=" :
                        self.BKGWINDOWS = [(int(con.split()[1]),int(con.split()[2])), (int(con.split()[3]), int(con.split()[4]))]
                    elif con.split()[0] == "SPECCONFIGFILE=" :
                        self.SPECCONFIGFILE = con.split()[1]
                    elif con.split()[0] == "WLFitFunc=" :
                        self.WLFITFUNC = con.split()[1]                        

                    elif con.split()[0] == "NIGHTLOGFILE=" :
                        self.NIGHTLOGFILE = con.split()[1]


    # Some helper functions for obtaining the full path to files.
    # They are not necessory for many steps.
    def GetFullPath(self,filename):
        """ Returns full path to filename.
        If the filename is a raw image, then it will give adress to the raw data directory.
        If the filename is new, then it will give address to the output directory of PC.currentnight""" 
        if not self._FrameDatabaseInitialised:
            #Initialise the database. # Don't call this if the Logfiles are not already created (in 0th step)
            self.FrameDatabase = {}
            for night in LoadDirectories(self,CONF=False):
                with open(os.path.join(self.RAWDATADIR,night,self.NIGHTLOGFILE),'r') as imglogFILE:
                    # Skip blank lines and Commented out lines with #
                    images = [imageLINE.split()[0] for imageLINE in imglogFILE if \
                              ((imageLINE.strip() != '') and (imageLINE[0] !='#'))]
                self.FrameDatabase[night] = set(images)
            self._FrameDatabaseInitialised = True

        if self.currentnight is None:
            # Raise error
            print('ERROR: Current directory not defined PC.currentnight')
            raise ValueError('ERROR: PC.currentnight not defined')

        # All fine now,  give full path intelegently.
        if filename in self.FrameDatabase[self.currentnight]:
            return os.path.join(self.RAWDATADIR,self.currentnight,filename)
        else:
            return os.path.join(self.RAWDATADIR,self.OUTDIR,self.currentnight,filename)
            

################################# Instrument Definitions  #################################################################
# Details as well as the name of each function calls of each instrument
# each step in pipeline is a function call with one input argument PC

# Ascii art courtesy : http://patorjk.com/software/taag/

class InstrumentObject(object):
    """ All the definitions and procedures of each instrument is defined inside this class."""
    def __init__(self,PC):
        self.PC = PC
        self.Name = self.PC.INSTRUMENT
        self.Banner = InstrumentDictionaries[self.Name]['Banner']
        self.About = InstrumentDictionaries[self.Name]['About']
        self.Steps = InstrumentDictionaries[self.Name]['Steps']
        self.StepsToRun = self.GetListofStepsToRun()

    def GetListofStepsToRun(self):
        """ Returns the list of steps to run based onPC, the PipelineConfiguration object """
        StepsToRun = []
        if self.Name in ['TANSPEC']:
#            StepsToRun += [0, 1, 2, 3, 4, 5, 11, 12]
            StepsToRun += [0, 1, 2, 3, 4, 5, 6]
#            if self.PC.IMGCOMBINE == 'Y':
#                StepsToRun += [6]
#            if self.PC.TODO == 'P':
#                StepsToRun += [7, 8, 9, 10]
#            elif self.PC.TODO == 'S':
#                StepsToRun += [11, 12]

        else:
            print('Unknown Instrument')

        return StepsToRun

    def StandardiseHeader(self,prihdr,filename=None):
        """ Return the header object after standardising the values in it, specific to instrument."""
        
        if self.Name == 'TANSPEC':
            prihdr[self.PC.COMMENTHDR] = prihdr['COMMENT'][-1].split('/')[0][1:].strip()
            del prihdr['COMMENT']
            try:
                del prihdr['']
            except KeyError:
                pass
            if  prihdr[self.PC.GRATINGHDR] == '-1':
                prihdr[self.PC.GRATINGHDR]= 'grating1'
        
        elif self.Name == 'HFOSC':
            ut_sec = float(prihdr[self.PC.UTHDR])
            prihdr[self.PC.UTHDR] = str(datetime.timedelta(seconds=ut_sec))  # Converting to HH:MM:SS
            # XSTART, XEND, YSTART, YEND, XBIN, YBIN etc not in HFOSC header. so we shall add dummy values
            prihdr[self.PC.XSTARTHDR] = 1
            prihdr[self.PC.XENDHDR] = prihdr['NAXIS1']
            prihdr[self.PC.YSTARTHDR] = 1
            prihdr[self.PC.YENDHDR] = prihdr['NAXIS2']
            # If the Xsize is 250 then it must be most likely binned by 2 in Xaxis for HFOSC!!
            prihdr[self.PC.XBINHDR] = 2 if prihdr['NAXIS1'] == 250 else 1
            prihdr[self.PC.YBINHDR] = 1

        elif self.Name == 'IFOSC':
            prihdr[self.PC.EXPTIMEHDR] = float(prihdr[self.PC.EXPTIMEHDR])/1000.0  # Convert ms to sec for IFOSC

        elif self.Name == 'ARIES1.3m512':
            date_ut = prihdr[self.PC.UTHDR]  # Is in the form YYYY-MM-DDTHH:MM:SS
            date = date_ut.split('T')[0]
            ut = date_ut.split('T')[1]
            prihdr[self.PC.UTHDR] = ut  # in HH:MM:SS
            prihdr[self.PC.DATEHDR] = date
            
            # Filter and object name has to be read form filename
            fsearch = re.search('(.*)([ubvriUBVRI]).*',os.path.splitext(filename)[0])
            if fsearch:
                prihdr[self.PC.OBJECTHDR] = fsearch.group(1)
                prihdr[self.PC.FILTERHDR] = fsearch.group(2).upper()
            # Add the X start and Ystart header objects
            subimageformat = prihdr['SUBRECT']
            prihdr[self.PC.XSTARTHDR] = int(subimageformat.split(',')[0])
            prihdr[self.PC.XENDHDR] = int(subimageformat.split(',')[1])
            prihdr[self.PC.YSTARTHDR] = int(subimageformat.split(',')[2])
            prihdr[self.PC.YENDHDR] = int(subimageformat.split(',')[3])
        return prihdr

    def IdentifyFrame(self,ObjectLine):
        """ Returns what object the input frame is based on the ObjectLine 
        'Filename, PC.OBJECTHDR, PC.EXPTIMEHDR, PC.FILTERHDR, PC.GRATINGHDR, PC.LAMPHDR, PC.GAINHDR, PC.SLITHDR, PC.ARGONLHDR, PC.NEONLHDR, PC.CONT1LHDR, PC.CONT2LHDR, PC.DATEHDR, PC.UTHDR, PC.RAHDR, PC.DECHDR, PC.COMMENTHDR, PC.XSTARTHDR, PC.XENDHDR, PC.YSTARTHDR, PC.YENDHDR, PC.XBINHDR, PC.YBINHDR, FileNumber' """
        Frame = 'UNKNOWN'
        if self.Name == 'HFOSC':
            if float(shlex.split(ObjectLine)[2]) == 0:
                Frame = 'BIAS'
            elif 'GRISM' in shlex.split(ObjectLine)[4].upper():
                if 'FE-' in shlex.split(ObjectLine)[5].upper():
                    Frame = 'LAMP_SPEC'
                elif 'HALOGEN' in shlex.split(ObjectLine)[5].upper():
                    Frame = 'FLAT_SPEC'
                else:
                    Frame = 'OBJECT_SPEC'
            else:
                Frame = 'OBJECT_IMG'
            

        return Frame



InstrumentDictionaries = {'TANSPEC':{
    'Banner' :"""

                     )  (    (                
  *   )    (      ( /(  )\ ) )\ )        (    
` )  /(    )\     )\())(()/((()/( (      )\   
 ( )(_))((((_)(  ((_)\  /(_))/(_)))\   (((_)  
(_(_())  )\ _ )\  _((_)(_)) (_)) ((_)  )\___  
|_   _|  (_)_\(_)| \| |/ __|| _ \| __|((/ __| 
  | |     / _ \  | .` |\__ \|  _/| _|  | (__  
  |_|    /_/ \_\ |_|\_||___/|_|  |___|  \___| 
                                              Data Reduction Pipeline for LR mode...
""",
    'About' : "TANSPEC on 3.6-m DOT, ARIES ...",
    'Steps' : {
        0 : {'Menu': 'Generate Log files of fits files in each directory.',
             'RunMessage': "RUNNING TASK:0  Generating log files of fits files in each directory..",
             'function': CreateLogFilesFromFits_subrout },
        1 : {'Menu': 'Selection of object frames, Flats/Sky/Lamps etc. to reduce',
             'RunMessage': "RUNNING TASK:1 Selecting object frames, Flats/Sky/Lamps etc..",
             'function': SelectionofFrames_subrout },
        2 : {'Menu': 'Visually inspect and/or reject object images one by one.',
             'RunMessage': "RUNNING TASK:2  Visual inspection and/or rejection of object frames..",
             'function': Manual_InspectObj_subrout},
        3 : {#'Menu': 'Visually inspect and/or reject Bias/Flats/Sky/Lamps one by one.',
             'Menu': 'Visually inspect and/or reject Flats/Lamps one by one.',
             'RunMessage': "RUNNING TASK:3  Visual inspection and/or rejection of Flats/Lamps frames..",
             'function': Manual_InspectCalframes_subrout},
#        4 : {'Menu': 'Subtract Overscan,Bias/sky',
#             'RunMessage': "RUNNING TASK:4  Subtracting biases or sky..",
#             'function': Bias_Subtraction_subrout},
       4  : {#'Menu': 'Apply Flat Correction and/or Bad pixel interpolation, and/or CR removal',
             'Menu': 'Apply Flat Correction and/or CR removal',
             'RunMessage': "RUNNING TASK:4  Flat correction etc..",
             'function': CombDith_FlatCorr_subrout},
#        6 : {'Menu': 'Align and combine images.',
#             'RunMessage': "RUNNING TASK:6  Aligning and combining images..",
#             'function': AlignNcombine_subrout },
#        7 : {'Menu': 'Make the list of images, Images4Photo.in to do Photometry.',
#             'RunMessage': "RUNNING TASK:7  Makeing list of images, Images4Photo.in to do Photometry..",
#             'function': Createlist_subrout },
#        8 : {'Menu': 'Select Stars and Sky region of the field on first image',
#             'RunMessage': "RUNNING TASK:8  Selecting Stars and Sky region of the field from first image..",
#             'function': Star_sky_subrout },
#        9 : {'Menu': 'Create Sextracter config file & coordinate output of first image.',
#             'RunMessage': "RUNNING TASK:9  Create Sextracter config file & coordinate output of first image..",
#             'function': Sextractor_subrout },
#        10 : {'Menu': 'Do Photometry',
#             'RunMessage': "RUNNING TASK:10 Doing Photometry..",
#             'function': Photometry },
        5: {'Menu': 'Input Spectrum pair subtraction and filenames.',
             'RunMessage': "RUNNING TASK:5  Inputing Spectrum pair subtraction and filenames...",
             'function': SpectralPairSubtraction_subrout },
        6: {'Menu': 'Extract wavelength calibrated 1D spectra from image.',
             'RunMessage': "RUNNING TASK:6  Extracting wavelength calibrated 1D spectra..",
             'function': LrSpectralExtraction_subrout }
    }
}
} # END of Instrument Definitions..#######
 
def parse_args(raw_args=None):
    """ Parses the command line input arguments """
    parser = argparse.ArgumentParser(description="Spectrum Extraction Pipeline")
    parser.add_argument('ConfigFile', type=str,
                         help="Configuration file which contains settings for TANSPEC XD-mode extraction. Hint: You can download sample config file from https://github.com/astrosupriyo/pyTANSPEC/blob/main/pyTANSPEC/config/TANSPEC_lrscript.config")
                         
    args = parser.parse_args(raw_args)  
    return args     
          
#-----Main Program Calls Begins here........................
# If you are using the functions seperately from module, see the Class PipelineConfig 
# to create the PipelineConfiguration instance which has to be passed on to many function calls.

def main(raw_args=None):
    """ The Pipeline main funtion to run while run from terminal as a standalone pipeline . """
    """ Extracts 2D spectrum image into 1D spectrum """
    args = parse_args(raw_args)
    
#    if len(sys.argv) < 2 :
#        print('-'*10)
#        print('Usage : {0} ScriptSettings.conf'.format(sys.argv[0]))
#        print('where,')
#        print('     ScriptSettings.conf is the configuration file for this run of reduction pipeline')
#        print(' ')
#        print("Note: This script should be run from the directory containing all the night's data directories.")
#        print('-'*10)
#        sys.exit(1)

    #pkgpath = os.path.split(pkgutil.get_loader('pyTANSPEC').get_filename())[0]
    #config_file_Main = os.path.join(pkgpath,'config/TANSPECscript.conf')
    
    try : 
#        PC = PipelineConfig(ConfigFilename = sys.argv[1])
        PC = PipelineConfig(ConfigFilename = args.ConfigFile)
    except IOError :
        print ("Error: Cannot open the file "+args.ConfigFile+". Setup the config file based on ScriptSettings.conf file correctly, before running the script.")
        sys.exit(1)

    # Load Instrument
    Instrument = InstrumentObject(PC)
    # Print Welcome Banner.. 
    print("Time: {0}".format(time.strftime("%c")))
    print(Instrument.Banner)
    print(Instrument.About)

    print('Running Initial Tests...')
    try :
        InitialTest(PC)  # Some important Initial tests for Sanity before running the code.
    except Exception as e :
        print(e)
        sys.exit(1)

    print("\n *** Very Very Important: Backup your RAW data first. Don't proceed without backup. ***\n")
    input("  Hit Enter to continue...")
    print('-'*10)
    print('IMP: All outputs of this run will be written to the directory '+os.path.join(PC.RAWDATADIR,PC.OUTDIR)+'\n')
    try:
        os.makedirs(os.path.join(PC.RAWDATADIR,PC.OUTDIR))
    except OSError:
        if os.path.isdir(os.path.join(PC.RAWDATADIR,PC.OUTDIR)) :
            try:
                with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,'StepsFinished'),'r') as stepsoverFLS:
                    StepsOver = stepsoverFLS.read()
            except IOError:
                    StepsOver = 'Nothing...'
            print('Steps you have already finished : '+StepsOver)        
        
            print("\n WARNING : Output directory "+os.path.join(PC.RAWDATADIR,PC.OUTDIR)+" already exists.\n\n         If you run the steps that have already been finished, everything inside it will be overwritten. Be warned...")
            while True:
                UserInput = input('\n Wanna proceed (hit Y/N): ').strip(' ')
                if UserInput.upper() == 'Y':
                    break 
                elif UserInput.upper() == 'N':
                    sys.exit()
                else:
                    continue
        else:
            raise
    else:
        print('Created directory :'+os.path.join(PC.RAWDATADIR,PC.OUTDIR))
        
        
    if PC.TODO == 'P' : todoinwords = 'Photometry' # For future
    elif PC.TODO == 'S' : todoinwords = 'Spectroscopy'


    print(" ---------------- Welcome to \033[91m {0} {1} \033[0m Pipeline --------------- \n".format(PC.INSTRUMENT,todoinwords))
    print("Enter the Serial numbers (space separated if more than one task in succession)")
    for i in Instrument.StepsToRun:
        print("\n {0}  {1}".format(i,Instrument.Steps[i]['Menu']))
    print("--------------------------------------------------------------- \n")

    try:
        todo = input('Enter the list : ')
        todo = todo.split()

        for task in todo :
            print("Time now: {0}".format(time.strftime("%c")))
            CalledTheTask=True
            try:
                task = int(task)
                if task not in Instrument.StepsToRun : 
                    raise IndexError('Task asked to run is not consistent with .conf file') # Raise Index error if a task not supposed to be run was given.
                print(Instrument.Steps[task]['RunMessage'])
            except (ValueError, IndexError) :
                print('Cannot understand the input task: {0}'.format(task))
                print('Skipping task :{0}'.format(task))
                CalledTheTask=False
            else: # Run the task..
                Instrument.Steps[task]['function'](PC)

            if CalledTheTask :
                with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,'StepsFinished'),'a') as stepsoverFLS:
                    stepsoverFLS.write(str(task)+' ')
    except KeyboardInterrupt :
        KeyboardInterrupt_handler()

    print("All tasks over... Enjoy!!!")
    print('Reduced data are at '+os.path.join(PC.RAWDATADIR,PC.OUTDIR))       
if __name__ == "__main__":
    main()           


# End of code
