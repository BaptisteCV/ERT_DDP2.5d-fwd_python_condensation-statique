#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import os
from os.path import join as pjoin
import fwdFunctions as fwd

'''
    MAIN

    Main script to run the forward of electrical resistivity data in a context of heterogeneous or fractured media.

    Last modification: 26/02/24 - Lea
'''


if __name__ == "__main__":

    # ---------- Get the name of the current folder (name depends on the chain) ----------
    path_folder = os.getcwd()
    name_new_folder = os.path.basename(os.getcwd())
    print("***** Demat da Python main.py adalek " + name_new_folder + " ! *****")

    # ---------- Set parameters ----------
    # To be changed by the user
    data_xx_yy = np.loadtxt(pjoin('fielddata','data_xxyy_R4_noratio.txt')) # tabVar[2]
    nb_electrodes = 128                                                              # tabVar[3]
    indexData = pd.read_csv(pjoin('fielddata','data.dat'), delimiter='\t', header=None,skiprows=1)
    # To not be changed
    expt_name_epm = 'testEPM_CS'                                                        # tabVar[0]
    expt_name_dfn = 'testDFN'                                                        # tabVar[1]
    path_results = 'results'  # tabVar[4]
    indexData = indexData.to_numpy()                                                 # tabVar[5]
    optCompile = 1                                                                   # tabVar[6]: 1 if library, 2 if executable (compilation of the C++ code)
    tabVar = [expt_name_epm, expt_name_dfn, data_xx_yy, nb_electrodes, path_results, indexData, optCompile]

    os.chdir('src')

    # ---------- Forward EPM (calculated only once) ----------

    if not os.path.exists(pjoin(tabVar[4], "rho_a_EPM.txt")):
        bool_error = True
        cpt = 0  # count the number of tests after a bug appeared
        while bool_error == True and cpt < 2:
            cpt += 1
            try:
                rho_a_EPM = fwd.fwdEPM(tabVar)
                bool_error = False
            except SystemExit:
                print(
                    f'An error occurred in {name_new_folder} while solving the forward EPM problem.\n C++ code is ran again in debug mode for the {cpt}th time.')
                if cpt == 1:  # Activate the debug mode of the C++ code
                    lines = open(pjoin(path_folder, 'input', 'testEPM', 'Simu.txt'), 'r').readlines()
                    new_last_line = (str(1))
                    lines[-1] = new_last_line
                    open(pjoin(path_folder, 'input', 'testEPM', 'Simu.txt'), 'w').writelines(lines)

        if bool_error == True:
            print('\nThe error in the forward model could not be solved.')
            print('\n         ***** THE SIMULATION IS INTERRUPTED !!! ***** ')
            exit()
        else:
            np.savetxt(pjoin(tabVar[4], "rho_a_EPM.txt"), rho_a_EPM, delimiter='\t')

    else:
        print("EPM has already been calculated.")

    # ---------- Forward DFN ----------
    bool_error = True
    cpt = 0  # count the number of tests after a bug appeared
    while bool_error == True and cpt < 2:
        cpt += 1
        try:
            rho_a_DFN = fwd.fwdDFN(tabVar)
            bool_error = False
        except SystemExit:
            print(
                f'An error occurred in {name_new_folder} while solving the forward DDP problem.\n C++ code is ran again in debug mode for the {cpt}th time.')
            if cpt == 1:  # Activate the debug mode of the C++ code
                lines = open(pjoin(path_folder, 'input', 'testDFN', 'Simu.txt'), 'r').readlines()
                new_last_line = (str(1))
                lines[-1] = new_last_line
                open(pjoin(path_folder, 'input', 'testDFN', 'Simu.txt'), 'w').writelines(lines)

    if bool_error == True:
        print('\nThe error in the forward model could not be solved.')
        print('\n         ***** THE SIMULATION IS INTERRUPTED !!! ***** ')
        exit()
    else:
        np.savetxt(pjoin(tabVar[4], "rho_a_DFN.txt"), rho_a_DFN, delimiter='\t')

