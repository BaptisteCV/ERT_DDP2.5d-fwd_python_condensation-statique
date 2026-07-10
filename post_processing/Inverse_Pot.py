#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import numpy as np
import scipy.fftpack as fftpack


'''
    INVERSE POTENTIAL

    This function is part of ERT_DDP_2.5D software (Cabllero-Sanz et al., 2017).
    It calculates the inverse potential from a matrix using N wavenumbers to do it. It has been created by Victor 
    Caballero in August 2014.
    Last modification: January 2015 - Victor Caballero

    :input:
        F_Pot: [array] potential in the Fourier domain
        Nx: [int] number of mesh cells in length
        Ny: [int] number of mesh cells in width
        c1x: [array] or [float] x position of several experiments
        g: [array] Xu et al. (2000) coefficients
    :output: 
        Potential: [array] inverse potential in the Fourier domain
'''


def Inverse_Pot(F_Pot, Ny, Nx, N, c1x, g):

    if N > 20:  # more than the optimized number

        # Typical inverse potential using at least 1000 wavenumbers
        T = 100
        F = np.ndim(F_Pot)

        if F != 2:
            # Transformation vector between both transforms (see base_discretization_idct)
            trans = np.ones([Ny, Nx, N])
            trans[:, :, 0] = trans[:, :, 0] * math.sqrt(N)
            trans[:, :, 1: N] = trans[:, :, 1: N] * math.sqrt(N / 2)
            F_Pot = F_Pot * trans

            # Multiplying the transformation vector by the transform potential, we can avoid the addition of the idct
            # constant in the next step. We have to do the transform in each vector of the third dimension, that is,
            # taking the vectors in depth.
            Potential = np.zeros([Ny, Nx, N])
            for i in range(Ny):
                for j in range(Nx):
                    c = fftpack.idct(F_Pot[i, j, :])
                    for h in range(N):
                        Potential[i, j, h] = Potential[i, j, h] + c[1, h]

            # In order to normalize the idct with the analytical solution, we have to multiply by 2/T.
            Potential = Potential * 2 / T

        else:
            trans = np.ones([N, Nx])
            trans[0, :] = trans[0, :] * math.sqrt(N)
            trans[1: N, :] = trans[1: N, :] * math.sqrt(N / 2)

            F_Pot = F_Pot * trans
            Potential = np.zeros([N, Nx])

            for i in range(Nx):
                c = fftpack.idct(F_Pot[:, i])
                Potential[:, i] = Potential[:, i] + c

            Potential = Potential * 2 / T
    else:
        # ---------- Optimized inverse potential ----------
        # We are going to develop the inversion using either the Xu,Duan and Zhang or
        # Pidlisecky optimized wavenumbers and Fourier weights. This inverse potential
        # is suitable just when we have an injection point, if not, we need to use the
        # standard one

        # Potential is a 4-D Matrix for inverting matrix potential with lots of
        # experiments, and it is a 3-D Matrix for inverting matrix potential with just one
        # experiment

        F = np.ndim(F_Pot)
        if F == 2:
            Potential = np.zeros([1, np.shape(F_Pot)[1]])
            for k in range(N):
                Potential[0, :] = Potential[0, :] + F_Pot[k, :]

        else:
            #p = np.arange(0, np.shape(c1x)[0]) # if c1x is an array
            p = np.arange(0, 1)
            #Potential = np.zeros([Ny, Nx, 1, np.shape(c1x)[0]])
            Potential = np.zeros([Ny, Nx, 1, 1])
            # F_Pot = F_Pot.reshape(Nx, Ny, N, np.shape(p)[0])
            F_Pot = F_Pot.reshape(Nx, Ny, N, 1)
            for k in range(0, N):
                Potential[:, :, 0, p] = Potential[:, :, 0, p] + F_Pot[:, :, k, p] * g[k]
        # This inversion just takes into account the first plane of the 2.5D
        # solution in the spatial domain, not the whole 3D domain (cube).

        Potential = Potential.reshape(Nx, Ny)

    return Potential
