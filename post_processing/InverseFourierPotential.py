# -*- coding: utf-8 -*-

import numpy as np
from post_processing.Inverse_Pot import Inverse_Pot


'''
    INVERSE FOURIER POTENTIAL

    This function is part of ERT_DDP_2.5D software (Cabllero-Sanz et al., 2017).
    It calculates the inverse potential from a matrix using N wavenumbers to do it. It has been created by Victor 
    Caballero in August 2014.
    Last modification: January 2015 - Victor Caballero
    
    NB: 
    c_1 : Position for a lot of experiments
    c_2 : Position of the second source for each experiment

    :input:
        Fourier_pot_AB: [array] potential values between injection electrodes in the Fourier domain
        Nx: [int] number of mesh cells in length
        Ny: [int] number of mesh cells in width
        Lx: [float] length of the domain
        Ly: [float] depth of the domain
        c1x: [array] or [float] x position of electrode A for a given quadrupole
        c2x: [array] or [float] x position of electrode B for a given quadrupole
        c1y: [array] or [float] y depth of electrode A for a given quadrupole
        c2y: [array] or [float] y depth of electrode B for a given quadrupole
    :output: 
        x: [array] coordinates of the electrodes in x
        y: [array] coordinates of the electrodes in y
        Pot_Matrix: [array] potential values between the injection electrodes in the real domain
'''


def InverseFourierPotential(Fourier_pot_AB, Nx, Ny, Lx, Ly, c1x, c2x, c1y, c2y):

    dx = Lx / Nx
    dy = Ly / Ny

    # Matrix block axis
    x = np.linspace(0, Lx - dx, Nx) + dx / 2
    y = np.linspace(0, Ly - dy, Ny) + dy / 2

    # Xu et al. (2000) coefficients
    wk = [0.0217102, 0.2161121, 1.0608400, 5.0765870]
    g = [0.0463660, 0.2365931, 1.0382080, 5.3648010]

    N = np.shape(wk)[0]
    # Bubble method to order wk and g vectors, because in the C++ Code we have them in increasing order
    for i in range(N):
        for j in range(i + 1, N):
            if wk[j] < wk[i]:
                temp_wk = wk[i]
                wk[i] = wk[j]
                wk[j] = temp_wk

                temp_g = g[i]
                g[i] = g[j]
                g[j] = temp_g

    # Optimization method in order to take into account the position of the
    # first and second sources as the same as the x vector, avoiding then changes in
    # resistivity values

    # for r in range(0, np.shape([c1x])[0]): # Don't forget to uncomment this line and indent correctly the following
    # lines from index_c1x to c2y,  if c1x is an array

    index_c1x = np.argmin(abs(c1x - x))
    index_c2x = np.argmin(abs(c2x - x))

    index_c1y = np.argmin(abs(c1y - y))
    index_c2y = np.argmin(abs(c2y - y))

    c1x = x[index_c1x]
    c2x = x[index_c2x]
    c1y = y[index_c1y]
    c2y = y[index_c2y]

    s = [0]

    i_Matrix = np.arange(0, Ny)
    q = np.arange(0, np.shape(s)[0])

    F_Pot_Matrix = np.zeros([Ny, Nx, N, np.shape(s)[0]])

    for p in range(0, np.shape(s)[0]):
        count = 0  # In order to change the sign of the simulation for the second source
        # ----------- Fourier Potential ---------- #
        F_Pot_Matrix1 = np.zeros([Ny, Nx, N])

        for k in range(0, N):
            F_Pot_Matrix1[i_Matrix, 0:Nx, k] = F_Pot_Matrix1[i_Matrix, 0:Nx, k] + Fourier_pot_AB[i_Matrix + k * Ny, 0:Nx]

        # Replacing NaNs and negative and small values for 0s
        F_Pot_Matrix1[np.isnan(F_Pot_Matrix1)] = 0

        if count == 0:  # That means this is the left source
            F_Pot_Matrix[:, :, :, p] = F_Pot_Matrix[:, :, :, p] + F_Pot_Matrix1
        elif count == 1:  # That means this is the right source
            F_Pot_Matrix[:, :, :, p] = F_Pot_Matrix[:, :, :, p] - F_Pot_Matrix1

        count = count + 1

    # ---------- Inverse Potential ---------- #
    Pot_Matrix = Inverse_Pot(F_Pot_Matrix, Ny, Nx, N, c1x, g)  # Maximum at the top

    return [x, y, Pot_Matrix]
