""" PIPE DATA TXT - Verison 1.0.0

    This module defines the pipeline TXT file object.
    The object is used by the pipeline to store raw, intermediate and
    reduced data.

"""

# Imports
import os
import logging
import re
import time
from dataparent import DataParent # Pipe Data Parent

class DataText(DataParent):
    """ Pipeline Data Text Object
        The data is stored as a list of strings
    """
    
    # File Name Fit: Regexp expression that fits valid fileanmes
    filenamefit = r'\.txt\Z' # fill files .txt and end of name
    
    def __init__(self, filename='', config=None):
        """ Constructor: Initialize data object and variables
        """
        # Run Parent Constructor
        super(DataText,self).__init__(config=config) # loads config, set up header
        # Set up internal Variables
        self.log=logging.getLogger('pipe.data.text') # Logger
        # Match for header: a '#' at start of line followed by space or end of line
        self.headermatch = r'\A#(\Z|\s)'
        # Match for header keyword: Alphanumeric values (with spaces) followed by
        #     ':' or '=' and any number of spaces
        self.keymatch = r'\A#[A-Za-z0-9\s]{0,20}(:|=)\s*'
        # if file exists, load it
        if os.path.exists(filename):
            self.load(filename)
        elif filename != '': # user specified a non-existant filename
            msg = 'No such file %s' %filename
            self.log.error(msg)
            raise IOError(msg)
        self.log.debug('Init: done')
        
    def sethead(self, line):
        """ Adds a text line that has been identified as a header line to
            the header.
            Header keyword lines are identified with self.keymatch.
            - line: the header line
        """
        # Identify keywords
        match = re.search(self.keymatch, line)
        if match:
            # Get key
            key = match.group()
            # Get value (rest of line after match)
            value = line[len(key):].strip()
            # remove header match and clearnup key
            match = re.search(self.headermatch, key)
            key = key[len(match.group()):].strip()[:-1].strip()
            # Set header keyword
            self.setheadval(key,value)
        else:
            # Not a key, it's a comment
            comment = line
            # remove header match and clearnup comment
            match = re.search(self.headermatch, comment)
            comment = comment[len(match.group()):].strip()
            # Set header keyword
            self.setheadval('COMMENT',comment)
        
    
    def loadhead(self,filename=''):
        """ Loads the header for text file given.
            Header lines are identified with self.headermatch.
            - filename: a string with the file name to load
        """
        # set self.filename and filename
        if len(filename) > 0:
            self.filename=filename
        else:
            filename=self.filename
        # Open file
        inf = file(filename,'rt')
        # Load header
        for line in inf:
            # identify header lines
            match = re.search(self.headermatch, line)
            if match: self.sethead(line)
        # Close file
        inf.close()
        self.log.debug('LoadHead: loaded text file')
        
    def load(self,filename=''):
        """ Loads the data and the header for a given text file.
            - filename: a string with the file name to load
        """
        # set self.filename and filename
        if len(filename) > 0:
            self.filename=filename
        else:
            filename=self.filename
        # Open file
        inf = file(filename,'rt')
        # Load header and data
        self.data = []
        for line in inf:
            # identify header lines
            match = re.search(self.headermatch, line)
            if match: self.sethead(line)
            else: self.data.append(line.strip())
        # Close file
        inf.close()
        self.log.debug('Load: loaded text file %s' % filename)
        
    def save(self,filename=''):
        """ Save the data to the specified file. Existing files are
            overwritten.
            - filename: a string with the file name to save
        """
        # get file name
        if filename == None:
            filename = self.filename
        # update pipeline keywords
        self.setheadval('Pipeline Version','Pipe v'+DataParent.pipever.replace('.','_'))
        self.setheadval('This filename',os.path.split(filename)[-1])
        self.setheadval('File Date',time.strftime('%Y-%m-%dT%H:%M:%S'))
        # Open file
        outf = file(filename, 'wt')
        # Save header
        for key in self.header:
            if 'COMMENT' in key:
                for comm in self.header['COMMENT']:
                    outf.write('# %s\n' % comm)
            elif 'HISTORY' in key:
                for hist in self.header['HISTORY']:
                    outf.write('# HISTORY: %s\n' % hist)
            else:
                outf.write('# %s: %s\n' % (key, self.header[key]))
        # Save data
        for line in self.data:
            outf.write('%s\n' % line)
        # Close file
        outf.close()
        self.log.debug('Save: saved text file %s' % filename)
        
""" === History ===
    2016-9-26 Marc Berthoud: Wrote first version
"""
