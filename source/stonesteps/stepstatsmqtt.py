#!/usr/bin/env python
""" PIPE STEP STATS MQTT- Version 1.0.0

    This pipe step collects reduction information and sends it to an MQTT
    server.
    
    @author: Marc Berthoud

    2DO:
    - Make MIMO step tso get names then use filenamebegin*filenameend to get info on files by glob
    - ./ Set up step
    - Program get statistics:
      - ./ Keep list of all good filepathnames (to exclude double counting)
      - ./ get list (from config) which steps to log and collect that data
      - ./ check if self.pipeline.outfiles and infiles are populated. AttributeError if using old darepype
        - For each filename in list get file type and add to statistics
      - ./ else get original input files
        - for each file look for each file type to add to statistics 
    - ./ Set up local single step and pipeline with 2 steps and test printing of file list
    - ./ Send statistics messages for each step to MQTT  
    - ./ Add list of objects, bands and exposure times from files in self.datain
    - ./ Test on local - sys.path.insert(0,'/Users/berthoud/instruments/software/DarePype') 
    - Test on stars in a folder as steponly, as part of pipe, with recent SEO data
      - pull on stars - set pythonpath then run the step, set up server pipeline for itzamna
        and queue reduction (test with running last observation day and any queue)
      - ALSO: make sure write system details about how reduction works (modes, config and such)

"""

import os
import logging
import paho.mqtt.client as mqtt
from darepype.drp import DataParent
from darepype.drp import StepMOParent

class StepStatsMqtt(StepMOParent):
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
        self.name='statsmqtt'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'STS'
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['steplist', 'RAW WCS FCAL',
                               'List of step identifiers to report on'])
        self.paramlist.append(['keylist', 'OBJECT',
                               'List of header keywords to report on'])
        self.paramlist.append(['broker','mqtt.server.com',
                               'url for MQTT broker'])
        self.paramlist.append(['port',8883, 'port for MQTT'])
        self.paramlist.append(['userpass','user|password',
                               "username and password, separated by '|'. " +
                               "alternatively full path of a file containing this information"])
        self.paramlist.append(['topic','server/stepstatsmqtt',
                               'mqtt topic to publish under'])
        # confirm end of setup
        self.log.debug('Setup: done')

    def run(self):
        """ Runs the data reduction algorithm. The self.datain is run
            through the code, the result is in self.dataout.
        """
        ### Get number of each step
        # Make dictionary numbdict[stepname]=filenamelist
        numbdict = {}
        steplist = [ s.strip() for s in self.getarg('steplist').split(' ') if len(s.strip())]
        for stepname in steplist:
            numbdict[stepname] = []
        # Get filelist from pipeline
        dp = DataParent(config = self.config)
        if self.pipeline:
            for fname in self.pipeline.outfiles+self.pipeline.infiles:
                dp.filename = fname
                fstep = fname[len(dp.filenamebegin):-len(dp.filenameend)]
                if fstep in steplist:
                    if not fname in numbdict[fstep]:
                        numbdict[fstep].append(fname)
        # Get header keyword information
        keylist = [ s.strip() for s in self.getarg('keylist').split(' ') if len(s.strip())]
        keydict = {}
        for key in keylist:
            keydict[key] = []
        # For each input file, find all files and look through header keys
        for inf in self.datain:
            fstep = inf.filename[len(inf.filenamebegin):-len(inf.filenameend)]
            # Add file to dictionary
            if fstep in steplist:
                if not inf.filename in numbdict[fstep]:
                    numbdict[fstep].append(inf.filename)
            # Look for files with that filename from other pipesteps
            for step in steplist:
                sfname = inf.filenamebegin+step+inf.filenameend
                if os.path.exists(sfname):
                    if not sfname in numbdict[step]:
                        numbdict[step].append(sfname)
            # Get header keywords
            for key in keylist:
                keyval = inf.getheadval(key)
                if not keyval in keydict[key]:
                    keydict[key].append(keyval)
        #for step in steplist: print(f'{step} {numbdict[step]}')
        ### Send statistics messages to mqtt
        # Get username password
        userpass = self.getarg('userpass')
        if os.path.exists(userpass):
            userpass = open(userpass,'rt').read().strip()
        userpass = userpass.split('|')
        # Set up client
        client =  mqtt.Client(client_id='mqtt_py_publisher')
        client.username_pw_set(userpass[0],userpass[1])
        client.tls_set()
        client.connect(self.getarg('broker'),self.getarg('port'))
        client.loop_start()
        stepnumbs = ''
        for step in steplist:
            stepnumbs += f', {len(numbdict[step])} {step} files'
        client.publish(self.getarg('topic'), stepnumbs[2:] )
        keymsg = ''
        for key in keylist:
            keymsg += f', {key}:'
            for keyval in keydict[key]:
                keymsg += f' {keyval}'
        client.publish(self.getarg('topic'), keymsg[2:] )
        client.disconnect()
        client.loop_stop()
        ### Copy datain to dataout
        self.dataout = self.datain
        self.log.debug('Run: Done')
        
        
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
    StepStatsMqtt().execute()

""" === History ===
2024-09-12 MGB: First version
"""