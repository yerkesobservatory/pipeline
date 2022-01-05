#!/usr/bin/env python
""" 
    Pipestep Flat HDR (High Dynamic Range)

    This module defines the pipeline step that corrects a pair raw image low
    and high gain images files for detector dark and flat effects.
    
    The step requires as input two files, a low gain file and a high gain file.
    It produces one output file.
    
    It uses StepLoadAux functions to call the following files:
        - masterpfit: ADD HERE WHAT THIS FILE IS
        - masterdark: ADD HERE WHAT THIS FILE IS
        - masterflat: ADD HERE WHAT THIS FILE IS
    
    Authors: Carmen Choza, Al Harper, Marc Berthoud
"""

from darepype.drp import DataFits # pipeline data object class
from darepype.drp.stepmiparent import StepMIParent # pipestep Multi-Input parent
from darepype.tools.steploadaux import StepLoadAux # pipestep steploadaux object class

class StepFlatHdr(StepLoadAux, StepMIParent):
    """ Pipeline Step Object to calibrate Flatfield High Dynamic Range files
    """
    
    stepver = '0.1' # pipe step version
    
    