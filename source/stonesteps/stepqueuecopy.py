#!/usr/bin/env python
''' PIPE STEP QUEUE COPY

    Pipe step to copy queue data. Input to the step are all selected
    headers in datain. Simplest is to pipe this step after StepLoadInput.
    
    The step:
    * identifies the necessary parameters
      (high/lowgain,  username, date/time, object)
    * Sorts the files by observation (user / object)
    * For each observation;
      * Makes the folder
      * For each file: makes the filename and copies the file
      * Makes a piperun file
      
    File format
    * Low gain file must have GAIN > 2.0 e/ADU and '.fitsL.fits' at
      the end of filename. Also data in second header.
    * High gain file must have GAIN < 2.0 e/ADU and '.fits' at end of filename
    * The observer is the last '_' separated word in the filename
    * Input filename is of the form
      m16_g-band_256s_bin1_2022-06-08_seo_apagul.fitsL.fits
    * Output filename is of the form
      m16_g-band_256s_bin1L_220608_104633_apagul_seo_0_RAW.fits
    * Observation folder stays unchanged and should have observer at the end
      m16_256s_bin2_2022-06-08_seo_apagul
      * Science files are in 
        observation/raw/science
      
    @author: Marc Berthoud
'''

import logging # logging object library
from darepype.drp.stepmiparent import StepMIParent
from datetime import datetime
import os
import shutil
from astropy.io import fits

class StepQueueCopy(StepMIParent):
    """ DarePype Step Queue Copy Object
        The object is callable. It requires a valid configuration input
        (file or object) when it runs.
    """
    stepver = '0.1' # pipe step version
    
    def setup(self):
        """ ### Names and Parameters need to be Set Here ###

            Sets the internal names for the function and for saved files.

            Defines the input parameters for the current pipe step.
            Setup() is called at the end of __init__
            The parameters are stored in a list containing the following
            information:
            - name: The name for the parameter. This name is used when
                    calling the pipe step from command line or python shell.
                    It is also used to identify the parameter in the pipeline
                    configuration file.
            - default: A default value for the parameter. If nothing, set
                       '' for strings, 0 for integers and 0.0 for floats
            - help: A short description of the parameter.
        """
        ### Set Names
        # Name of the pipeline reduction step (same case as config file)
        self.name='queuecopy'
        # Identifier for saved file names.
        self.procname = 'QCP'
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['outpath', '/data/queue',
            'Output path will be OUTPATH/observer/observation/file.fits'])
        self.paramlist.append(['piperunpath', '/data/queue/piperuns',
            'Folder for the piperun files'])

    def run(self):
        """ Runs the data reduction algorithm. The self.datain is run
            through the code, the result is in self.dataout.
        """
        # For each file add gain, observer to header, then make filename
        observers = [] # list of observers
        observations = {}
        for dat in self.datain:
            # Get filepath / name
            filepath, filename = os.path.split(dat.filename)
            obspath = os.path.split(filepath)[0]
            obspath = os.path.split(obspath)[0]
            obspath = os.path.split(obspath)[1]

            # Check gain / remove end of filename
            # fitsL.fits i.e. bin1L files have gain 18.31
            # fits i.e. bin1H files have gain 0.86
            try:
                gain = dat.getheadval('GAIN')
            except KeyError:
                # all data is in the second header, load it
                hdus = fits.open(dat.filename)
                dat.header = hdus[1].header
            if not (dat.getheadval('GAIN') > 2.0 and '.fitsL.fits' in filename):
                self.log.warn("Image %s has mismatched gain/filename" % filename)
            # Remove end of filename
            fname = filename
            if fname[-5:] == '.fits':
                fname = fname[:-5]
            if fname[-6:] == '.fitsL':
                fname = fname[:-6]
            # Get observer into OBSERVER keyword
            observer = fname.split('_')[-1]
            dat.setheadval('OBSERVER',observer.capitalize())
            # Make output filename
            outfname = fname[:fname.index('_bin')+5] # get all to binN
            try:
                if dat.getheadval('GAIN') > 2.0: outfname += 'L_' # add high/low gain
                else: outfname += 'H_'
            except:
                outfname += 'L_'
            fdate = datetime.strptime(dat.getheadval('DATE-OBS'),'%Y-%m-%dT%H:%M:%S')
            outfname += fdate.strftime('%y%m%d_%H%M%S_') # add date_time
            outfname += observer + '_seo_0_RAW.fits' # and end of name
            # Get output path
            outpath = os.path.join(self.getarg('outpath'), 
                                   observer.capitalize(), obspath)
            outfname = os.path.join(outpath,outfname)
            dat.setheadval('OUTFNAME', outfname)
            self.log.debug(f'{dat.filename} -> {outfname}')
            # Add observer and observation
            if observer not in observers:
                observers.append(observer)
                observations[observer] = []
            # Add observation
            if obspath not in observations[observer]:
                observations[observer].append(obspath)
        # Print observers and observations
        self.log.info(f'Observer {observer}: {" ".join(observations[observer])}')
        # Make folders and piperun objects
        # Loop through observers
        for observer in observers:
            # Make folder
            outpath = os.path.join(self.getarg('outpath'), 
                                   observer.capitalize())
            if not os.path.exists(outpath):
                os.mkdir(outpath)
            # Loop through observations
            for obspath in observations[observer]:
                # Make folder
                outpath = os.path.join(self.getarg('outpath'),
                                       observer.capitalize(), obspath)
                if not os.path.exists(outpath):
                    os.mkdir(outpath)
                # Make piperun filepathname
                piperun = obspath + '.txt'
                piperun = os.path.join(self.getarg('piperunpath'), piperun)
                # Make file content
                text = (
                    f'# === Piperun file for {observer} ===\n'
                    f'\n'
                    f'# !!! Auto-generated Piperun file - may be overwritten !!!\n'
                    f'pythonpath = /data/scripts/pipeline/source\n'
                    f'pipeconf = /data/scripts/pipeline/config/pipeconf_SEO.txt\n'
                    f'/data/scripts/pipeline/config/dconf_stars.txt\n'
                    f'pipemode = seo_server_hdr\n'
                    f'loglevel = DEBUG\n'
                    f'logfile = /data/scripts/pipeline/PipeLineLog.txt\n'
                    f'{os.path.join(outpath, obspath+"_log.txt")}\n' # Additional logfile
                    f'inputfiles = {os.path.join(outpath, "*RAW.fits")}\n'
                    f'outputfolder = {outpath}\n'
                    )
                # Save file
                self.log.info('Writing piperun at %s' % piperun)
                with open(piperun,'wt') as outf:
                    outf.write(text)
        # Copy the files
        for dat in self.datain:
            # Set up command
            shutil.copy(dat.filename, dat.getheadval('OUTFNAME'))               
        # Populate dataout (just so there's something in it)
        self.dataout = self.datain[0]
    
if __name__ == '__main__':
    """ Main function to run the pipe step from command line on a file.
        Command:
          python stepparent.py input.fits -arg1 -arg2 . . .
        Standard arguments:
          --config=ConfigFilePathName.txt : name of the configuration file
          -t, --test : runs the functionality test i.e. pipestep.test()
          --loglevel=LEVEL : configures the logging output for a particular level
          -h, --help : Returns a list of parameters
    """
    StepQueueCopy().execute()

""" === History ===
* 2022-6-27: Marc Berthoud - First Version
  * Development thoughts:
    * Stepdatagroup can not be used since header keys may be missing
    * Only the headers are input to avoid having to load all files into
      memory. This is the reason why the files are copied by StepQueueCopy
      instead of using the regular darepype save.
    
"""
