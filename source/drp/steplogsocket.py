#!/usr/bin/env python
""" PIPE STEP LOG SOCKET - Version 1.0.0

    This module runs a background socket connection that accepts log messages to
    be added to the pipeline log. The socket should be parent of all pipeline
    steps which use that socket. To start the socket connection logsocketstart
    has be be called in the setup() function of these steps, repeated calls to
    logsocketstart are expected and ignored. The socket connection runs in a
    separate thread. Configuration settings for logsocket are expected in the
    [logsocket] section of the configuration file.

    This object is a child of StepParent so that self.log is available and the
    -t option can be used for testing.

    @author: berthoud
"""

import os # os library
import socket # socket library
import thread # thread library
import time # time libary
import logging # logging object library
from drp.stepparent import StepParent

class StepLogSocket(StepParent):
    """ Pipeline LogSocket Object
        Should be parent object of pipesteps that require
        the socket.
    """
    socketactive = [False,]
    activeloghost = [None,]
    activelogport = [None,]

    def logsocketstart(self):
        """ Opens the log socket if it is not yet active
        """
        # Check if socket is active
        if not self.socketactive[0]:
            # Set up socket
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            try:
                conf = self.config['logsocket']
            except:
                self.log.warn('Missing [logsocket] configuration ' +
                              '- using default setup: localhost:50747')
                conf={'host':'localhost','port':50747}
            # Try specified port + up to 20
            cport = int(conf['port'])
            err = False
            for i in range(20):
                try:
                    self.log.debug('Attempting port ' + str(cport))
                    sock.bind((conf['host'],cport))
                    # Start listening
                    sock.listen(5)
                    # Start thread
                    thread.start_new_thread(self.logsocketrun,(sock,))
                    err = False
                    break
                except Exception:
                    err = True
                    cport += 1
                    continue
            # If still failed, then raise error
            if err: raise

            # Save host and port
            self.activeloghost[0] = socket.gethostbyname(conf['host'])
            self.activelogport[0] = cport

            self.log.info('Started Logging Socket on port %d' % cport)
            self.socketactive[0]=True
        else:
            self.log.debug('Multiple calls to LogSocketStart - socket already running')

        self.loghost = self.activeloghost[0]
        self.logport = self.activelogport[0]

    def logsocketrun(self, sock):
        """ Log socket subprocess loop
        """
        try:
            conf = self.config['logsocket']
        except:
            self.log.warn('Missing [logsocket] configuration ' +
                          '- using default backup logger: pipe.logserver')
            conf={'deflogger':'pipe.logserver'}
        while 1:
            # Get new connection
            conn, addr = sock.accept()
            #print('Conected with %s at address %s' % (addr[0],str(addr[1])))
            # Get the message
            reply = conn.recv(1024)
            #print('  Got message: %s' % reply)
            # Close the connection
            conn.close()
            #print('  Connection Closed')
            ### Make log message
            #print('  Making log message:')
            # Split incoming message
            split = reply.split('\t')
            if len(split) < 2:
                lvl = 'INFO'
                lgr = conf['deflogger']
                msg = split[0]
            elif len(split) < 3:
                lvl = split[0]
                lgr = conf['deflogger']
                msg = split[1]
            else:
                lvl = split[0]
                lgr = split[1]
                msg = '\t'.join(split[2:])
            # Get logging level
            log=logging.getLogger(lgr)
            if lvl.upper() in ['DEBUG','DEBUGGING']: log.debug(msg)
            elif lvl.upper() in ['INFO','INFORMATION']: log.info(msg)
            elif lvl.upper() in ['WARN','WARNING']: log.warn(msg)
            elif lvl.upper() in ['ERR','ERROR']: log.error(msg)
            elif lvl.upper() in ['CRIT','CRITICAL']: log.critical(msg)
            else: log.info(msg)
            #if lvl.upper() in ['DEBUG','DEBUGGING']: lvl=logging.DEBUG
            #elif lvl.upper() in ['INFO','INFORMATION']: lvl=logging.INFO
            #elif lvl.upper() in ['WARN','WARNING']: lvl=logging.WARN
            #elif lvl.upper() in ['ERR','ERROR']: lvl=logging.ERROR
            #elif lvl.upper() in ['CRIT','CRITICAL']: lvl=logging.CRITICAL
            #else: lvl=logging.INFO
            #rec=logging.LogRecord(lgr,lvl,'','',msg,{},{})
            #log.handle(rec)

    def test(self):
        """ Test Pipe Step Parent Object:
            Runs a set of basic tests on the object
        """
        # log message
        self.log.info('Testing pipe step %s' %self.name)
        # set up logger
        self.logsocketstart()
        # wait
        time.sleep(1)
        # send log message via python socket
        self.log.debug('Python Socket Test:')
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect(('localhost', 50747))
        message = 'error\tpipe.logtest.client\tTest Message'
        s.sendall(message)
        s.close()
        # send log message via c socket
        self.log.debug('C-Program Test')
        os.system("/Users/berthoud/tech/tcpip/simpleclient $'warn\tpipe.logtest.C\tWarn Message'")
        os.system("/Users/berthoud/tech/tcpip/simpleclient $'info\tpipe.logtest.C\tInfo Message'")
        os.system("/Users/berthoud/tech/tcpip/simpleclient $'debug\tpipe.logtest.C\tDebug Message'")
        # wait
        time.sleep(1)
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
    StepLogSocket().execute()

""" === History ===
    2015-4-5 Marc Berthoud: Ver 0.1.0: Written and tested
"""
