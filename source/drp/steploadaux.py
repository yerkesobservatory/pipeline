#!/usr/bin/env python
""" PIPE STEP LOAD AUXILIARY FILE - Version 1.0.0

    This module provides a pipe step with a simple mechanism to search
    auxiliary files and match them to the loaded data.
    
    @author: berthoud
    
    To check:
    - check from MI step
"""

import os # os library
import glob # glob library
from drp.dataparent import DataParent # Pipeline Data object
from drp.stepparent import StepParent
from drp.stepmiparent import StepMIParent # To check if we have datain or [datain, datain, ...]
#from docutils.parsers import null

class StepLoadAux(StepParent):
    """ HAWC Pipeline Step Parent Object
        The object is callable. It requires a valid configuration input
        (file or object) when it runs.
    """

    def loadauxsetup(self, auxfilepar = 'auxfile'):
        """ This object does not have a conventional setup function.
            Since it is intended to be inherited by other steps. The
            code below is optional if the child step wants to be set
            up with the parameters auxfile and fitkeys. Alternatively
            (or as an option) the parameter auxfile and fitkeys can
            be passed to loadauxfile().
            
            auxfilepar: Parameter name for the parameter
                containing the filepathname to search for auxiliary
                files.
            
            This function should be called in the setup function of
            the child step after self.paramlist has been initiated.
        """
        # Set name of the auxfile parameter
        self.auxfilepar = auxfilepar
        # Append parameters
        self.paramlist.append([auxfilepar, '%s/*.fits' % auxfilepar,
            'Filename for auxiliary file(s). Can contain * and ? ' +
            'wildcards to match multiple files to be selected using fitkeys ' +
            '(default = %s/*.fits)' % auxfilepar])
        self.paramlist.append(['bkup'+auxfilepar, 'bkup%s/*.fits' % auxfilepar,
            'Back up filename for auxiliary file(s). Can contain * and ? ' +
            'wildcards to match multiple files to be selected using fitkeys ' +
            '(default = bkup%s/*.fits)' % auxfilepar])
        self.paramlist.append(['fitkeys', [],
            'List of header keys that need to match auxiliary data file ' +
            '(default = []) - only used if multiple files ' +
            'match %s' % auxfilepar])

    def loadauxfile(self, auxfile = '', fitkeys = [], data = None):
        """ Searches for files matching auxfile. If only one match is
            found, that file is returned. Else the header
            keywords listed in fitkeys are matched between the
            data and the auxfiles which were found. The first auxfile
            for which these keywords values best match the ones
            from data is selected. The loaded pipe data (child of
            dataparent) is returned.
            
            auxfile: A filepathname (with glob wildcards) for the
                     auxfiles to consider. If not specified, 
                     self.getarg('auxfile') is used.
            fitkeys: A list of header keywords to be matched.
                     If no list is specified, self.getarg('fitkeys')
                     is used.
            data: A pipedata object to match the auxiliary file to.
                  If no data is specified self.datain is used (for
                  Multi Input steps self.datain[0]).
        """
        # Get parameters
        auxfile = os.path.expandvars(self.getarg(self.auxfilepar))
        fitkeys  = self.getarg('fitkeys')
        if len(fitkeys) == 1 and len(fitkeys[0]) == 0:
            fitkeys = []
        # Glob the list of files
        auxlist = glob.glob(auxfile)
        # Throw exception if no file found
        if len(auxlist) < 1:
            auxfile = os.path.expandvars(self.getarg('bkup'+self.auxfilepar))
            auxlist = glob.glob(auxfile)
            if len(auxlist) < 1:
                msg = 'No %s files found under %s' % (self.auxfilepar, auxfile)
                self.log.error(msg)
                raise ValueError(msg)
        # Get datain object (depends on step being SingleInput or MultiInput)
        if data == None:
            if issubclass(self.__class__, StepMIParent):
                data = self.datain[0]
            else:
                data = self.datain 
        # Return unique file
        if len(auxlist) == 1:
            self.log.info(' LoadAuxFile: Found unique file = %s' % auxlist[0])
            auxdata = DataParent(config=self.config).load(auxlist[0])
            data.setheadval('HISTORY','%s: Best %s = %s' % 
                            (self.name, self.auxfilepar, os.path.split(auxdata.filename)[1],))
            return auxdata
        # check format (make first element uppercase)
        try:
            _ = fitkeys[0].upper()
        except AttributeError:
            # AttributeError if it's not a string
            self.log.error('LoadAuxFile: fitkeys config parameter is ' +
                           'incorrect format')
            raise TypeError('fitkeys config parameter is incorrect format')
        # get keywords from data
        datakeys=[]
        for fitkey in fitkeys:
            datakeys.append(data.getheadval(fitkey))
        # match aux files, return best aux file
        bestind = -1 # index of file with best match in filelist
        bestfitn = 0 # number of keywords that match in best match
        auxind = 0 # index for going through the list
        while auxind < len(auxlist) and bestfitn < len(fitkeys):
            # load keys of aux file
            auxhead = DataParent(config=self.config).loadhead(auxlist[auxind])
            auxkeys=[]
            for fitkey in fitkeys:
                try:
                    auxkeys.append(auxhead.getheadval(fitkey))
                except KeyError:
                    self.log.warning('LoadAuxFile: missing key [%s] in auxfile <%s>'
                                     % (fitkey, auxlist[auxind] ) )
                    auxkeys.append('')
            # count number of fitting keywords
            keyfitn=0
            while ( keyfitn < len(fitkeys) and
                    datakeys[keyfitn] == auxkeys[keyfitn] ):
                keyfitn = keyfitn + 1
            # compare with previous best find
            if keyfitn > bestfitn:
                bestind = auxind
                bestfitn = keyfitn
            auxind=auxind+1
        # Select best file
        if bestfitn == len(fitkeys):
            self.log.info('LoadAuxFile: Matching %s found is <%s>' % 
                          (self.auxfilepar, auxlist[bestind]) )
            auxdata = DataParent(config=self.config).load(auxlist[bestind])
        elif bestind > -1:
            self.log.warn('LoadAuxFile: Could not find perfect aux file match')
            self.log.warn('Best flat file found is <%s>' % auxlist[bestind] )
            auxdata = DataParent(config=self.config).load(auxlist[bestind])
        else:
            self.log.warn('LoadAuxFile: NO MATCH finding aux file')
            self.log.warn('Returning first file, <%s>' % auxlist[bestind] )
            auxdata = DataParent(config=self.config).load(auxlist[bestind])
        # Return selected file
        data.setheadval('HISTORY','%s: Best %s = %s' % 
                        (self.name, self.auxfilepar, os.path.split(auxdata.filename)[1],))
        return auxdata
    
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
    StepLoadAux().execute()

""" === History ===
"""
