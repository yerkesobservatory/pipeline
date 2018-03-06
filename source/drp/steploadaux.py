#!/usr/bin/env python
""" PIPE STEP LOAD AUXILIARY FILE - Version 1.0.0

    This module provides a pipe step with a simple mechanism to search
    auxiliary files and match them to input data. An initial list of
    potential aux files is determined from a path/filename string
    (a glob), the final file(s) is determined by matching header
    keywords between the auxiliary files and the step input data.
    
    @author: berthoud
    
    To check:
    - check from MI step and SI step
    - make sure it works with any pipedata type
    To Change:
    OK Test function with environment: 1 input file, 3 bias and 3 dark files (fudge hdr keys)
      OK make a test object to play with, feed it input file and look what it finds to fit it
    OK Loadauxsetup: run it multiple time and loadauxfile have auxfilepar option
       - Also update fitkeys to %s(fit)keys %sbkup %sfile
         OPEN QUESTION: change all or just fitkeys? ANSWER: %sfile, bkup%s, %sfitkeys (with fitskeys as fallback)
    OK Loadauxfile: set up with auxfilepar option
    OK loadauxfile() and loadauxname() for returning filenames only
      - copy all code and update
      - update loadauxfile
      - test both
    - Loadauxfile: set up with multiple option
      - Changes to loadauxname: namelist/single option
        - idea: load all hdus then loop through keys and shorten list with each run
                from back to front and don't delete last file.
        - first load all headers
        - loop through all keywords (oldlist and newlist)
          end if all keys checked or break if len(newlist) == 0
      - loadauxfile: for filelist/single option
    - loadauxname() with option to give range for auxdaterange (fraction of day)
      
    Future Development:
    - Function to load aux files w/o parameters (if needed)
"""

import os # os library
import glob # glob library
import string # for string.join
import time # time library
from drp.dataparent import DataParent # Pipeline Data object
from drp.stepparent import StepParent
from drp.stepmiparent import StepMIParent # To check if we have datain or [datain, datain, ...]
#from docutils.parsers import null

class StepLoadAux(StepParent):
    """ HAWC Pipeline Step Parent Object
        The object is callable. It requires a valid configuration input
        (file or object) when it runs.
    """

    def loadauxsetup(self, auxpar = 'aux'):
        """ This object does not have a conventional setup function.
            Since it is intended to be inherited by other steps. This
            code should be used to set up the child step with the
            parameters auxfile and auxfitkeys.
            
            auxpar: Parameter name for the parameter containing
                the filepathname glob to search for auxiliary
                files.
            
            This function should be called in the setup function of
            the child step after self.paramlist has been initiated.
        """
        # Set name of the auxfile parameter
        self.auxpar = auxpar
        # Append parameters
        self.paramlist.append([auxpar + 'file', '%sfolder/*.fits' % auxpar,
            'Filename for auxiliary file(s). Can contain * and ? ' +
            'wildcards to match multiple files to be selected using fitkeys ' +
            '(default = %sfolder/*.fits)' % auxpar])
        self.paramlist.append(['bkup'+auxpar, 'bkup%sfolder/*.fits' % auxpar,
            'Back up filename for auxiliary file(s). Can contain * and ? ' +
            'wildcards to match multiple files to be selected using fitkeys ' +
            '(default = bkup%sfolder/*.fits)' % auxpar])
        self.paramlist.append([auxpar + 'fitkeys', [],
            'List of header keys that need to match auxiliary data file ' +
            '(default = []) - only used if multiple files ' +
            'match %sfile' % auxpar])
        if 'daterange' not in [ par[0] for par in self.paramlist] :
            self.paramlist.append(['daterange',1.0,
                'If DATE-OBS is in fitkeys, files are matched within this many days.'])

    def loadauxname(self, auxpar = '', data = None, multi = False):
        """ Searches for files matching auxfile. If only one match is
            found, that file is returned. Else the header
            keywords listed in auxfitkeys are matched between the
            data and the auxfiles which were found. The first auxfile
            for which these keywords values best match the ones
            from data is selected. The filename of the best match
            is returned.
            
            auxpar: A name for the aux file parameter to use. This
                    allows loadauxfiles to be used multiple times
                    in a given pipe step (for example for darks and
                    flats). Default value is self.auxpar which is set
                    by loadauxsetup().
            data: A pipedata object to match the auxiliary file to.
                  If no data is specified self.datain is used (for
                  Multi Input steps self.datain[0]).
        """
        ### Setup
        # Set auxpar
        if len(auxpar) == 0:
            auxpar = self.auxpar
        # Get parameters
        auxfile = os.path.expandvars(self.getarg(auxpar + 'file'))
        fitkeys  = self.getarg(auxpar + 'fitkeys')
        if len(fitkeys) == 1 and len(fitkeys[0]) == 0:
            fitkeys = []
        ### Look for files - return in special cases
        # Glob the list of files
        auxlist = glob.glob(auxfile)
        # Throw exception if no file found
        if len(auxlist) < 1:
            self.log.warn('No files found under %s - looking in backup' % auxfile)
            auxfile = os.path.expandvars(self.getarg('bkup'+auxpar))
            auxlist = glob.glob(auxfile)
            if len(auxlist) < 1:
                msg = 'No %s files found under %s' % (auxpar, auxfile)
                self.log.error(msg)
                raise ValueError(msg)
        # Get datain object (depends on step being SingleInput or MultiInput)
        if data == None:
            if issubclass(self.__class__, StepMIParent):
                data = self.datain[0]
            else:
                data = self.datain 
        # Return unique file, or all files if fitkeys is empty
        if len(auxlist) == 1 or len(fitkeys) == 0:
            if len(auxlist) == 1:
                self.log.info('LoadAuxName: Found unique file = %s' % auxlist[0])
            else:
                self.log.info('LoadAuxName: No fitkeys: Return first %sfile match = %s' %
                              (self.auxpar, auxlist[0]) )
            data.setheadval('HISTORY','%s: Best %sfile = %s' % 
                            (self.name, self.auxpar, os.path.split(auxlist[0])[1],))
            if multi:
                return auxlist
            else:
                return auxlist[0]
        ### Select files with Fitkeys
        # check format (make first element uppercase)
        try:
            _ = fitkeys[0].upper()
        except AttributeError:
            # AttributeError if it's not a string
            self.log.error('LoadAuxFile: fitkeys config parameter is ' +
                           'incorrect format - need list of strings')
            raise TypeError('fitkeys config parameter is incorrect format' +
                            ' - need list of strings')
        # Load all headers from auxlist into a auxheadlist (pipedata objects)
        auxheadlist = []
        for auxnam in auxlist:
            auxheadlist.append(DataParent(config = self.config).loadhead(auxnam))
        # Look through keywords, only keep auxfiles which fit keys
        for key in fitkeys:
            newheadlist = []
            # Look through auxfiles, transfer good ones
            if key in 'DATE-OBS': # SPECIAL CASE DATE-OBS: 
                # get time for data
                datime = time.mktime(time.strptime(data.getheadval('DATE-OBS'), 
                                                   '%Y-%m-%dT%H:%M:%S'))
                # get time offset (from data) for each auxfile
                auxtimes = []
                for auxhead in auxheadlist:
                    auxtime = time.mktime(time.strptime(auxhead.getheadval('DATE-OBS'), 
                                                        '%Y-%m-%dT%H:%M:%S'))
                    auxtimes.append(abs(auxtime-datime))
                # only keep auxfiles which are within daterange of closest auxfile
                mindiff = min(auxtimes)
                timerange = self.getarg('daterange') * 86400
                for auxi in range(len(auxheadlist)):
                    if auxtimes[auxi] - mindiff < timerange:
                        newheadlist.append(auxheadlist[auxi])
            else: # Normal Keyword compare
                for auxhead in auxheadlist:
                    # Check if the auxfile fits (compare with data)
                    if auxhead.getheadval(key) == data.getheadval(key) :
                        # it fits -> add to newheadlist
                        newheadlist.append(auxhead)
            # break key loop if no files left
            if len(newheadlist) == 0:
                break
            else:
                auxheadlist = newheadlist
        
        ### Select file to return
        if multi:
            # Return all filenames
            auxname = [aux.filename for aux in auxheadlist]
            # Return message
            if len(auxname) > 3:
                listnames = "%d files: %s to %s" % (len(auxname),auxname[0],auxname[-1])
            else:
                listnames = string.join(auxname)
            if len(newheadlist) > 0:
                self.log.info('LoadAuxName: Matching %s found are <%s>' % 
                              (auxpar, listnames) )
            else:
                self.log.warn('LoadAuxName: NO MATCH finding aux files')
                self.log.warn('Returning files <%s>' % listnames )
        else:
            # Return first filename
            auxname = auxheadlist[0].filename
            # Select best file
            if len(newheadlist) > 0:
                self.log.info('LoadAuxName: Matching %s found is <%s>' % 
                              (auxpar, auxname) )
            else:
                self.log.warn('LoadAuxName: NO MATCH finding aux file')
                self.log.warn('Returning first file <%s>' % auxname )
            listnames = auxname # just so we can use it below
        data.setheadval('HISTORY','%s: Best %s = %s' % 
                        (self.name, auxpar, listnames))
        # Return selected file
        return auxname

    def loadauxfile(self, auxpar = '', data = None, multi = False):
        """ Uses loadauxname to search for files matching auxfile.
            See loadauxname for parameter description.
            
            A pipedata object with the best match is returned.
        """
        # Get auxname
        auxname = self.loadauxname(auxpar,data, multi)
        # Load auxdata
        if multi:
            auxdata = [ DataParent(config=self.config).load(auxnam) for auxnam in auxname ]
        else:
            auxdata = DataParent(config=self.config).load(auxname)
        # Return selected file
        return auxdata
    
    def test(self):
        """ Test Pipe Step Parent Object:
            Runs a set of basic tests on the object
        """
        # log message
        self.log.info('Testing pipe step %s' %self.name)
        # Set up the step
        self.name = 'loadaux' # this is not used normally as loadaux is normally used as parent
        self.loadauxsetup('test1')
        self.loadauxsetup('test2')
        for par in self.paramlist:
            print(par)
        # Load input data
        self.datain = DataParent(config=self.config).load('IN_a0_1.fits')
        # Get test1 auxfile
        auxf = self.loadauxname('test1',multi=True)
        print('********** ' + repr(auxf))
        # Get test2 auxfile
        auxf = self.loadauxname('test2')
        print('********** ' + repr(auxf))
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
    StepLoadAux().execute()

""" === History ===
"""
