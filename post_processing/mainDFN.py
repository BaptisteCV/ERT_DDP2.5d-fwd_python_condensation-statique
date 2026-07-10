#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    MAIN DFN

    This function is part of ERT_DDP_2.5D software (Cabllero-Sanz et al., 2017).
    It is the main file for inverting 2.5D data provided in the Fourier domain from the C++ code to the real
    domain. The function was first developed in a context of surface electrode configuration.
    The present version is an adaptation to a cross-borehole configuration.

    Last modification: 25/02/12 - Lea

    :input:
        res_cxx: [array] results from the C++ code
        tabVar: [list] containing all the common parameters and constant to all functions
    :output:
        rho_a_DDP: [array] apparent electrical resistivity values simulated by ERT_DDP_2.5D
"""

import math
import sys
import numpy as np
import pandas as pd
from os.path import join as pjoin
from post_processing.InverseFourierPotential import InverseFourierPotential
from post_processing.Farum_elecConfig import Farum_elecConfig


def mainDFN(tabVar,res_cxx):

    # Electrode configurations
    expt_name = pjoin("input",tabVar[1])
    nb_electrodes = tabVar[3]
    data_xx_yy = tabVar[2]

    # Domain properties
    path_MatrixProperties = expt_name
    file_name_MatrixProperties = 'DomainProperties.txt'
    matrix_properties = np.loadtxt(pjoin(path_MatrixProperties, file_name_MatrixProperties))
    Lx = matrix_properties[0]
    Ly = matrix_properties[1]

    # Results in the Fourier domain (from C++ code)
    if res_cxx==0:
        path_results_DFN = pjoin("output",tabVar[1])
        file_name_DFN= 'Matrix_Potential'

    # Simu
    file_name_simu = 'Simu.txt'
    Simu = pd.read_csv(pjoin(path_MatrixProperties, file_name_simu), delimiter='\t')
    Nx = int(Simu['FOURIER_POTENTIAL_COMPUTATION'][1])
    Ny = int(Simu['FOURIER_POTENTIAL_COMPUTATION'][2])
    delta_x = np.round(Lx / Nx, 4)
    delta_y = np.round(Ly / Ny, 4)

    # Define the electrode configurations
    [configx, configy] = Farum_elecConfig(data_xx_yy)
    matrixIndex=tabVar[5][:,1:4]

    # Define the difference of potential V_MN for each experiment
    nb_exp = np.shape(configx)[0]
    rho_a_DFN = np.zeros((nb_exp, 1))

    # Load all the results
    AllResultsDFN = np.zeros([nb_electrodes, 4 * Ny, Nx])

    for i in range(0, nb_electrodes):
        if res_cxx == 0: # If the executable version of the C++ code is used
            AllResultsDFN[i, :, :] = np.loadtxt(pjoin(path_results_DFN, f'{file_name_DFN}{i+1}.txt'))
        else: # If the library version of the C++ code is used
            # 28/07/25
            try:
                #AllResultsDFN[i, :, :] = res_cxx[i+1].toarray(order='C')
                AllResultsDFN[i, :, :] = res_cxx[i+1].toarray()
            except KeyError:
                print(f'KeyError problem with matrix {i}. Exit function.')
                sys.exit(1)

    for i in range(nb_exp):

        index_elec_A = int(matrixIndex[i, 0])-1
        index_elec_B = int(matrixIndex[i, 1])-1

        position_elec_xA = configx[i, 0]; position_elec_yA = configy[i, 0]
        position_elec_xB = configx[i, 1]; position_elec_yB = configy[i, 1]
        position_elec_xM = configx[i, 2]; position_elec_yM = configy[i, 2]
        position_elec_xN = configx[i, 3]; position_elec_yN = configy[i, 3]

        # Potential distributions in Fourier domain for an injection at electrodes A and B
        Fourier_potential_A_DFN = AllResultsDFN[index_elec_A, :, :]
        Fourier_potential_B_DFN = AllResultsDFN[index_elec_B, :, :]

        # Difference of potential VMN for injection in A and extraction in B
        Fourier_pot_AB_DFN = Fourier_potential_A_DFN - Fourier_potential_B_DFN

        # Inversion of the potential distribution in the real domain
        [x, y, Potential_matrix_AB_DFN] = InverseFourierPotential(Fourier_pot_AB_DFN, Nx, Ny, Lx, Ly, position_elec_xA,
                                                                  position_elec_xB, position_elec_yA, position_elec_yB)
        #[x, y, Potential_matrix_BA_DFN] = InverseFourierPotential(Fourier_pot_AB_DFN, Nx, Ny, Lx, Ly, position_elec_xB,
       #                                                           position_elec_xA, position_elec_yB, position_elec_yA)

        #Potential_matrix_AB_DFN = (Potential_matrix_AB_DFN + Potential_matrix_BA_DFN) / 2

        V_MN_DFN = (Potential_matrix_AB_DFN[math.ceil(position_elec_yM / delta_y) - 1,
                                            math.ceil(position_elec_xM / delta_x) - 1] - Potential_matrix_AB_DFN[
                        math.ceil(position_elec_yN / delta_y) - 1, math.ceil(position_elec_xN / delta_x) - 1])
        #V_MN_DFN_inv = (Potential_matrix_BA_DFN[math.ceil(position_elec_yM / delta_y) - 1,
        #                                        math.ceil(position_elec_xM / delta_x) - 1] - Potential_matrix_BA_DFN[
        #                    math.ceil(position_elec_yN / delta_y) - 1, math.ceil(position_elec_xN / delta_x) - 1])

        #V_MN_DFN = (V_MN_DFN + V_MN_DFN_inv) / 2

        # Geometric factor
        AM = np.sqrt(np.square(np.subtract(position_elec_xM, position_elec_xA)) + np.square(
            np.subtract(position_elec_yM, position_elec_yA)))
        AMinv = 1 / AM
        AN = np.sqrt(np.square(np.subtract(position_elec_xN, position_elec_xA)) + np.square(
            np.subtract(position_elec_yN, position_elec_yA)))
        ANinv = 1 / AN
        BM = np.sqrt(np.square(np.subtract(position_elec_xM, position_elec_xB)) + np.square(
            np.subtract(position_elec_yM, position_elec_yB)))
        BMinv = 1 / BM
        BN = np.sqrt(np.square(np.subtract(position_elec_xN, position_elec_xB)) + np.square(
            np.subtract(position_elec_yN, position_elec_yB)))
        BNinv = 1 / BN

        k_fact = 4 * math.pi / (AMinv - ANinv - BMinv + BNinv)

        # Apparent resistivity
        rho_a_DFN[i] = V_MN_DFN * k_fact


    return rho_a_DFN
