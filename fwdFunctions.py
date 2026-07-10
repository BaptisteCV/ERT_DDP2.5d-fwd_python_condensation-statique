#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    FORWARD FUNCTIONS

    Contain the forward procedures for the Equivalent Porous Media (EPM) model (homogeneous) and the Discrete Dual Porosity (DDP) model (matrix+
    Discrete Fracture Network (DFN)).
    The forward models are adapted from Caballero-Sanz et al. (2017) to data acquired in time-lapse cross-borehole ERT during aquifer remediation.

    Last modification: 25/02/12 - Lea

"""

import os
from os.path import join as pjoin
import sys
import numpy as np
import subprocess
from post_processing.mainEPM import mainEPM
from post_processing.mainDFN import mainDFN
try:
    sys.path.append(pjoin(os.getcwd(), '..', '..','src','electrical_resistivity', 'Release'))
    from  PERFORM import perform
except:
    print('PERFORM.so not found. C++ code may be run in executable mode.')


def fwdEPM(tabVar):

    # Change the name of the test in SimuType.txt file
    with open(pjoin('input','SimuType.txt'), 'w') as f:
        f.write(str(tabVar[0]))
    
    # C++ (simulate current injection in the Fourier domain)
    if tabVar[6]==2:
        res_cxx=0
        subprocess.run([pjoin("electrical_resistivity","Release", "ElectricalResistivity")])
    else:
        res_cxx = perform()

    # Calculate the electrical resistivity from the C++ outputs
    rho_a_EPM = mainEPM(tabVar,res_cxx)

    return rho_a_EPM


def fwdDFN(tabVar):

    # Change the name of the test in SimuType.txt file
    with open(pjoin('input','SimuType.txt'), 'w') as f:
        f.write(str(tabVar[1]))

    # C++ (simulate current injection in the Fourier domain)
    if tabVar[6]==2:
        res_cxx=0
        subprocess.run([pjoin("electrical_resistivity","Release", "ElectricalResistivity")])
    else:
        res_cxx = perform()

    # Calculate the electrical resistivity from the C++ outputs
    rho_a_DFN = mainDFN(tabVar,res_cxx)

    return rho_a_DFN

