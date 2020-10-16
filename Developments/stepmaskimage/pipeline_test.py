import os
import sys
import urllib
import logging
from darepype.drp import DataParent, PipeLine, DataFits
from stepwebastrometry import StepWebAstrometry

filenames = ['M5_r-band_60s_bin2_200711_053415_itzamna_seo_0_RAW.fits']
            #['M5_r-band_60s_bin2_200711_053415_itzamna_seo_0_RAW_TABLE.fits']
            #['M5_r-band_60s_bin2_200711_053415_itzamna_seo_0_RAW.fits']
            #['M5_g-band_60s_bin2_200711_053258_itzamna_seo_0_RAW.fits']#,
            #   'M5_i-band_60s_bin2_200711_053538_itzamna_seo_0_RAW.fits',
            #   'M5_r-band_60s_bin2_200711_053415_itzamna_seo_0_RAW.fits']
logfilename = '/Users/josh/Desktop/pipeline_test/pipelog.txt'
# Location where the code is
codefolder = '/Users/josh/pipeline'
# Location where you want the data to be
datafolder = '/Users/josh/Desktop/pipeline_test/data'
# Location of config files
baseconfig = os.path.join(codefolder,'pipeline', 'Developments', 'stepwebastrometry', 'pipeconf_stonedge_auto.txt')
logging.basicConfig(filename = logfilename, level = logging.DEBUG,
                    format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s")

os.chdir(codefolder)
if not os.path.exists(datafolder):
    os.mkdir(datafolder)
infilenames = [os.path.join(datafolder,f) for f in filenames]

### Load a file into a DataFits object
dparent = DataParent(config = baseconfig)
dfits = dparent.load(infilenames[0])
print(dfits.filename)

### Look at the FITS header of the loaded object
# print(repr(dfits.header))

### OPTIONAL BUT RECOMMENDED: Check if all necessary files exist
error_flag = False
# Check if configuration file exists
if not os.path.exists(baseconfig):
    print('ERROR: The config file you specified, %s,\n  does NOT exist on your computer, fix "config" above' % baseconfig)
    error_flag = True
# Check if input files exist
for name in infilenames:
    if not os.path.exists(name):
        print('ERROR: The input file you specifed, %s,\n  does NOT exist on your computer, fix "inputnames" above' % name)
        error_flag = True
if not error_flag:
    print("All Good")

os.chdir('/Users/josh/pipeline/pipeline/Developments/stepwebastrometry')
step = StepWebAstrometry()
indata = []
for f in infilenames:
    fits = DataFits(config=baseconfig)
    fits.load(f)
    indata.append(fits)

outdata = step(indata[0])
print('Done')
