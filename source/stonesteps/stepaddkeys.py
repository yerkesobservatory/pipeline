#!/usr/bin/env python
""" PIPE STEP ADD KEYS - Version 1.0.0

    This pipe step adds FITS keywords to the file based on information
    in the file name.
    
    The step also copies selected keywords which are not available in the 
    primary HDU but are present in other HDUs.
    
    @author: Joe Polk
"""

import re
import logging # logging object library
import os
from darepype.drp import StepParent

class StepAddKeys(StepParent):
    """ HAWC Pipeline Step Parent Object
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
        # Name of the pipeline reduction step
        self.name='addkeys'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'keys'
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['filternames', ['unknown'], 'List of valid strings for filter names'])
        self.paramlist.append(['keystocopy',[],'List of Keywords to copy to primary HDU if they are missing there (default is empty [] list)'])

    def run(self):
        """ Runs the data reduction algorithm. The self.datain is run
            through the code, the result is in self.dataout.
        """
        ### Copy input file to output file
        self.dataout = self.datain.copy()
        ### Add keywords to primary HDU
        keystocopy = self.getarg('keystocopy')
        for key in keystocopy:
            # See if the key is in primary header
            try:
                val = self.dataout.getheadval(key, errmsg = False)
            except KeyError:
                # Key is not in primary header - look in all headers
                try: 
                    val = self.dataout.getheadval(key, dataname = 'allheaders', errmsg = False)
                    self.dataout.setheadval(key, val)
                    self.log.debug('Setting primary header[%s] to %s' % (key,repr(val)))
                except KeyError:
                    self.log.debug('Header key %s not found in data' % key)
        ### Get file name only (no path)
        fileonly = os.path.split(self.dataout.filename)[1]
        ### Add OBSERVER
        # Check if observer keyword exists and is valid
        got_observer = False
        try:
            observer = self.dataout.getheadval('OBSERVER', errmsg = False)
            # Make sure it's not invalid entry
            if not observer.lower() in ['', 'unk', 'unknown', 'remote', 'sirius'] :
                got_observer = True
        except KeyError:
            # if there's a key error -> there's no OBSERVER
            pass # as got_observer is already false
        if not got_observer:
            ## File patterns:
            # filepatt[0] <==> OBJECT_BAND_EXPOSURE_BINNING_YYMMDD_HHMMSS_OBSERVER_seo_FILENUM_RAW.fits (format starting spring 2018)
            # filepatt[1] <==> OBJECT_BAND_EXPOSURE_BINNING_OBSERVER_DATE(YYYYMmmDD)_TIME(11h11m11s)_num0000_HjfyiYt5_seo.fits
            # filepatt[2] <==> OBJECT_BAND_EXPOSURE_BINNING_YYYYmmmDD_OBSERVER_OBSNUM_seo.fits
            # filepatt[3] <==> YYYY-MM-DD_OBSERVER_OBJECT_Ez5/ (astroclass)

            filepatt = ['(_\d{6}){2}', # two sets of 6 numbers (date and time)
                        'bin\d_[a-zA-Z]', # Binning followed by alphabetic character (observer)
                        'bin\d_\d{4}[a-z]{3}\d{2}', # Binning followed by Date in YYYYmmmDD format
                        '^20\d{2}\-[0-1]\d\-[0-3]\d'] # Starts with date in YYYY-MM-DD format

            obsbefore = ['\d{6}_','bin\d_','20\d{2}[a-z]{3}\d{2}_','^.{11}']
            obsafter = ['_seo_','_20\d{2}','_num','_[a-zA-Z]|_\d']

            fileindex = 0
            for i in range(len(filepatt)):
                if re.compile(filepatt[i]).search(self.dataout.filename):
                    fileindex = i
                    break
            # splits file name before observer name
            obsTemp = re.split(obsbefore[fileindex],self.dataout.filename)[-1]
            #splits file name after observer name
            observer = re.split(obsafter[fileindex],obsTemp)[0]
            self.log.debug('File name is type %d, which fits regexp pattern \'%s\', found observer = %s'
                           % (fileindex, filepatt[fileindex], observer ) )
        else:
            self.log.debug('Observer in header = ' + observer)
        ### Add OBJECT name
        got_object = False # assume it's not there
        try:
            objname = self.dataout.getheadval('OBJECT', errmsg = False)
            # Make sure it's not invalid entry
            if not objname.lower() in ['', 'unk', 'unknown'] :
                got_object = True
                self.log.debug('Object in header = ' + objname)
        except KeyError:
            pass # b/c got_object is already false
        if not got_object:
            # Getting the object from the file name (all before first "_")
            objname = fileonly.split('_')[0]
            self.log.debug('Object from Filename = ' + fileonly.split('_')[0])
        if objname.lower() in  ['', 'unk', 'unknown']:
            # Finds and formats RA/DEC values from the fits header if present,
            # if not sets objname to unknown
            objname = ''
            self.log.info('Object unknown, changeing object name to RA/DEC')
            try:
                ra = self.dataout.getheadval('RA')
                dec = self.dataout.getheadval('DEC')
                # Make sure it's not invalid entry
                if ra.lower() in ['', 'unk', 'unknown'] or dec.lower in ['', 'unk', 'unknown']:
                    got_radec = False
                else:
                    got_radec = True
            except KeyError:
                # if there's a key error -> there's no RA or DEC
                got_radec = False
            if got_radec==True:
                regexRA=r'(\d{1,2}):(\d{2}):(\d{2}).?(\d{1,2})?'
                RAre=re.findall(regexRA, self.dataout.getheadval('RA'))
                if RAre[0][3]=='':
                    obsRA=RAre[0][0]+"h"+RAre[0][1]+"m"+RAre[0][2]+".0"+"s"
                else:
                    obsRA=RAre[0][0]+"h"+RAre[0][1]+"m"+RAre[0][2]+"."+RAre[0][3]+"s"
                regexDEC=r'(\d{1,2}):(\d{2}):(\d{2}).?(\d{1,2})?'
                DECre=re.findall(regexDEC, self.dataout.getheadval('DEC'))
                DECsign=self.dataout.getheadval('DEC')[0]
                if DECre[0][3]=='':
                    obsDEC=DECsign+DECre[0][0]+"d"+DECre[0][1]+"m"+DECre[0][2]+".0"+"s"
                else:
                    obsDEC=DECsign+DECre[0][0]+"d"+DECre[0][1]+"m"+DECre[0][2]+"."+DECre[0][3]+"s"
                objname = obsRA+obsDEC
                self.log.debug('Object RA/DEC = ' + objname)
                got_object=True
            else:
                objname='unknown'
        ### Add FILTER
        got_filter = False # assume it's not there
        try:
            filtername = self.dataout.getheadval('FILTER', errmsg = False)
            if not filtername.lower() in ['', 'unk', 'unknown'] :
                got_filter = True
        except KeyError:
            pass # b/c got_filter is already false
        if not got_filter:
            # Getting the filter from the file name (second '_' separated part)
            filtername = fileonly.split('_')[1]
#             filtername = 'unknown' # in case no filter name is found
#             for f in self.getarg('filternames'):
#                 if f in fileonly:
#                     filtername = f
#                     break # exit the for loop
            self.log.debug('Filter = ' + filtername)
        ### Make changes to file
        # Put keyword into the output file
        # (need: OBSERVER and OBJECT keywords with values from the filename)
        self.dataout.setheadval('OBSERVER', observer )
        self.dataout.setheadval('OBJECT', objname )
        self.dataout.setheadval('FILTER', filtername )
        self.log.debug('Keys Updated')
    
    def test(self):
        """ Test Pipe Step Parent Object:
            Runs a set of basic tests on the object
        """
        # log message
        self.log.info('Testing pipe step %s' %self.name)

        # log message
        self.log.info('Testing pipe step %s - Done' %self.name)
    
if __name__ == '__main__':
    """ Main function to run the pipe step from command line on a file.
        Command:
          python stepparent.py input.fits -arg1 -arg2 . . .
        Standard arguments:
          --config=ConfigFilePathName.txt : name of the configuration file
          -t, --test : runs the functionality test i.e. pipestep.test()
          --loglevel=LEVEL : configures the logging output for a particular level
          -h, --help : Returns a list of 
    """
    StepAddKeys().execute()

""" === History ===
2016-12-20: Joe Polk, Marc Berthoud: First Version
"""
