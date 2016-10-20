#!/usr/bin/env python
""" PIPE STEP DATA GROUP - Version 1.0.0

    This module defines the pipeline step data group object. This pipe
    step divides the list of input files into groups according to given
    header data keywords. Each group is then reduced in a MISO or MIMO
    reduction step (the redstep) and all the output files are returned by
    StepDataGroup.
    
    The step to reduce the files with is specified in the configuration
    file and is found using the same process the pipeline uses.
    
    @author: berthoud
"""

import logging # logging object library
from drp.pipedata import PipeData
from drp.stepmoparent import StepMOParent

class StepDataGroup(StepMOParent):
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
        self.name='datagroup'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'unk' # is taken from the redstep below
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        # Reduction step, unset if == None
        self.redstep = None
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['redstepname', 'stepmiparent',
            'Name of the Pipestep to use'])
        self.paramlist.append(['groupkeys','OBJECT|DATASRC',
            'List of header keywords to decide data group membership (| separated)'])        
        self.paramlist.append(['groupkfmt','',
            'List of group key formats for string comparison (unused if "", | separated)'])

    def run(self):
        """ Runs the data reduction algorithm. The self.datain is run
            through the code, the result is in self.dataout.
        """
        ### Get redstep if it's not loaded
        if self.redstep == None:
            # Get list of step packages
            steppacks = self.config['general']['steppacks']
            if isinstance(steppacks,str): # ensure we have a list if only
                steppacks = [steppacks]    # 1 item, steppacks is str
            # get stepname
            stepname = self.getarg('redstepname')
            for pack in steppacks:
                if pack == '.': # look in local directory
                    mod = stepname.lower()
                else:
                    mod = '%s.%s' % (pack,stepname.lower())
                self.log.debug('looking for %s' % mod)
                # import the module
                try:
                    stepmodule = __import__(mod, globals(), locals(),
                                            [stepname])
                    self.log.debug('Pipe step %s found in %s' %
                                   (stepname,pack))
                    break
                except ImportError,msg:
                    tmp = 'No module named %s' % stepname.lower()
                    if str(msg).startswith(tmp): # module not present in directory
                        self.log.debug('Pipe step %s not found in %s' %
                                       (stepname, pack))
                    else: # module tries to import a missing package
                        raise
                except:   # print out import errors not due to missing
                    raise # modules (e.g., typos in code)
            # If not found -> Raise error
            if stepmodule == None:
                msg = 'Could not find step=%s' % stepname
                msg += ' in pipemode=%s' % self.pipemode
                self.log.error('Setup: ' + msg)
                raise ImportError(msg)
            # Make a step instance and add to step list
            try:
                self.redstep=stepmodule.__dict__[stepname]()
            except KeyError, error:
                msg = 'Pipe step=%s' % stepname
                msg+= ' not found in module=%s' % mod
                self.log.error('Setup: %s' % msg)
                raise error                
        ### Group the input files
        # Setup datagroups, get keys and key formats
        datagroups=[]
        groupkeys = self.getarg('groupkeys').split('|')
        groupkfmt = self.getarg('groupkfmt')
        if len(groupkfmt) == 0:
            groupkfmt = None
        else:
            groupkfmt = groupkfmt.split('|')
        # Loop over files
        for data in self.datain:
            groupind = 0
            # Loop over groups until group match found or end reached
            while groupind < len(datagroups) :
                # Check if data fits group
                found = True
                gdata = datagroups[groupind][0]
                for keyi in range(len(groupkeys)):
                    # Get key from group and new data - format if needed
                    key = groupkeys[keyi]
                    dkey = data.getheadval(key)
                    gkey = gdata.getheadval(key) 
                    if groupkfmt != None:
                        dkey = groupkfmt[keyi] % dkey
                        gkey = groupkfmt[keyi] % gkey
                    # Compare
                    if dkey != gkey :
                        found = False
                # Found -> add to group
                if found:
                    datagroups[groupind].append(data)
                    break
                # Not found -> increase group index
                groupind += 1
            # If not in any group -> make new group
            if groupind == len(datagroups):
                datagroups.append([data,])
        # info messages
        self.log.debug(" Found %d data groups" % len(datagroups))
        for groupind in range(len(datagroups)):
            group = datagroups[groupind]
            msg = "  - Group %d len=%d" % (groupind, len(group) )
            for key in groupkeys:
                msg += " %s = %s" % (key, group[0].getheadval(key))
            self.log.debug(msg)
        ### Reduce input files - collect output files
        self.dataout = []
        # Loop over groups -> save output in self.dataout
        for group in datagroups:
            dataout = self.redstep(group)
            # add output to dataout
            if isinstance(dataout,PipeData):
                self.dataout.append(dataout)
            else:
                for data in dataout:
                    self.dataout.append(dataout)
        # Set procname to redstep.procname
        self.procname = self.redstep.procname
    
    def reset(self):
        """ Resets the step to the same condition as it was when it was
            created. Internal variables are reset, any stored data is
            erased.
        """
        self.log.debug('Reset: done')
        
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
    StepDataGroup().execute()

""" === History ===
"""
