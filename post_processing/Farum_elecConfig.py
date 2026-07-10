#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np


'''
    FARUM ELECTRODE CONFIGURATIONS

    Separate the x and y coordinates of the electrode ordered by quadrupole. This script is adapted to a specific 
    data-set acquired with time-lapse cross-borehole ERT during aquifer remediation.
    Last modification: 24/07/11 - Lea

    :input:
        data_xx_yy: [str] coordinates of the electrodes per quadrupole
    :output: 
        [elec_configx, elec_configy]: [array] coordinates of the electrodes per quadrupole separated between x and y
        coordinates
'''


def Farum_elecConfig(data_xx_yy):

    elecPos = data_xx_yy
    nrows = np.shape(elecPos)[0]
    elec_configx = np.zeros([nrows, 4])
    elec_configy = np.zeros([nrows, 4])

    for j in range(nrows):
        elec_configx[j, 0: 4] = elecPos[j, 0: 4] + 3  # array of x coordinates: xA, xB, xM, xN
        elec_configy[j, 0: 4] = -elecPos[j, 4: 8]  # array of y coordinates: yA, yB, yM, yN

    return [elec_configx, elec_configy]
