#!/usr/bin/env python
""" PIPE STEP DATA GROUP - Version 1.0.0

    This module defines the pipeline step data group object. This pipe
    step divides the list of input files into groups according to given
    header data keywords. Each group is then reduced in a MISO or MIMO
    reduction step (the redstep) and all the output files are returned by
    StepDataGroup.
    
    To avoid re-reducing the same groups of files, the outputs are stored
    and if group elements match earlier group results are returned.
    
    @author: berthoud
"""

import logging # logging object library
from drp.dataparent import DataParent
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
        # Groups variables: store information about reduced data to avoid
        #   re-reducing unchanged data
        self.groupoutputs = [] # List of output data objec for each group
        self.groupidkeys = [] # For each group a list of the dataidkey for
                              # the data objects that were used for this group
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['redstepname', 'StepMIParent',
            'Name of the Pipestep to use'])
        self.paramlist.append(['groupkeys','OBJECT|DATASRC',
            'List of header keywords to decide data group membership (| separated)'])        
        self.paramlist.append(['groupkfmt','',
            'List of group key formats for string comparison (unused if "", | separated)'])
        self.paramlist.append(['fileidkey','',
            'Header keyword to re-identify files to avoid re-reducing the same groups' +
            ' (default is '' indicating all data has to be re-reduced)'])

    def run(self):
        """ Runs the data reduction algorithm. The self.datain is run
            through the code, the result is in self.dataout.
        """
        ### Get redstep if it's not loaded
        if self.redstep == None:
            # Get the step
            datap = DataParent(config=self.config)
            self.redstep = datap.getobject(self.getarg('redstepname'))
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
            msg = "  Group %d len=%d" % (groupind, len(group) )
            for key in groupkeys:
                msg += " %s = %s" % (key, group[0].getheadval(key))
            self.log.debug(msg)
        ### Reduce input files - collect output files
        self.dataout = []
        # Make new variables for groupidkeys and groupoutputs
        groupidkeys = []
        groupoutputs = []
        # Loop over groups -> save output in self.dataout
        for groupi in range(len(datagroups)):
            group = datagroups[groupi]
            # Get fileidkeys to see if unchanged groups should be re-reduced
            fileidkey = self.getarg('fileidkey')
            if len(fileidkey):
                # Get fileidkeys for the current new group
                newkeys = [dat.getheadval(fileidkey) for dat in group]
                copykeys = ['x']
                # Search for fit in existing groups: fit is index
                fit = -1
                for fit in range(len(self.groupidkeys)):
                    # Make copy of new keys
                    copykeys = list(newkeys)
                    # For each key in group[fit]
                    for val in self.groupidkeys[fit]:
                        if val in copykeys:
                            # Remove key from copykeys if found
                            del copykeys[copykeys.index(val)]
                        else:
                            # Else: group[fit] does not match, go to next group
                            copykeys = ['x'] # if all have been removed
                            break
                    # Check if any values left in copykeys
                    if len(copykeys) == 0:
                        # No values left in copykeys, group[fit] is valid match
                        break
                # Any values left in copykeys -> no match found
                if len(copykeys):
                    fit = -1
                    self.log.debug('New datagroup # %d has no previous match'
                                   % groupi)
                else:
                    self.log.debug('New datagroup # %d matches previous group # %d' 
                                   % (groupi, fit))
            else:
                fit = -1
            # Reduce the data
            if fit < 0:
                dataout = self.redstep(group)
                # Add groupoutputs and groupidkeys
                if len(fileidkey):
                    groupoutputs.append(dataout)
                    idkeys = [dat.getheadval(fileidkey) for dat in group]
                    groupidkeys.append(idkeys)
            else:
                groupoutputs.append(self.groupoutputs[fit])
                groupidkeys.append(self.groupidkeys[fit])
                dataout = self.groupoutputs[fit]
            # add output to dataout
            if issubclass(dataout.__class__,DataParent):
                self.dataout.append(dataout)
            else:
                for data in dataout:
                    self.dataout.append(dataout)
        # Copy groupidkeys and groupoutputs
        self.groupoutputs = groupoutputs
        self.groupidkeys = groupidkeys
        # Set procname to redstep.procname
        self.procname = self.redstep.procname
    
    def reset(self):
        """ Resets the step to the same condition as it was when it was
            created. Internal variables are reset, any stored data is
            erased.
        """
        self.redstep = None
        self.groupoutputs = []
        self.groupidkeys = []
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
