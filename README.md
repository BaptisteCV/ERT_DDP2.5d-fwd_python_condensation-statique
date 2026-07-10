# ERT_DDP_2.5D
Version Python - February 2026

## Description of the code and changes

This branch contains the forward procedure ERT_DDP_2.5D developed by Caballero-Sanz et
al. (2017) adapted to the context of cross-borehole ERT measurements (Lelimouzin et aL., 2024).

The code is composed of two parts: (i) scripts in C++ that calculate the potential at each electrode
for each experiment in the Fourier domain, (ii) scripts in MATLAB or python that calculate the 
apparent electrical resistivity for each experiment from C++ results in the real world domain.

The post-precessing scripts have been translated from MATLAB to Python by Fabien Rétif 
(Snale company, Aubenas, France) and Léa Lelimouzin. 

The C++ scripts have been parallelized in openMP by Fabien Rétif.

The C++ scripts have been modified in order to compile the scripts as either executable or 
library. The library version works with the PyBind11 package that requires Eigen structures. Part of
the initial ublas structures have been changed for Eigen structures. 
Eigen maps have been used to store the C++ results. A version with Eigen vectors is commented in the scripts.
Changes have been implemented by Delphine Roubinet, Fabien Rétif and Léa Lelimouzin.

Two options to generate synthetic DFNs (no electrical simulations) have been added 
(i) from a fractal definition of the network (like in Gisladottir et al., 2016) and (ii) from a power-law based definition. 

## Project architecture

```text
.
├── ERT_DDP_2.5D
│   ├── documentation
│   ├── example
│   ├── fielddata
│   ├── post_processing
│   │   ├── Farum_elecConfig.py
│   │   ├── Inverse_Pot.py
│   │   ├── InverseFourierPotential.py
│   │   ├── mainDFN.py
│   │   └── mainEPM.py
│   ├── src
│   │   ├── electrical_resistivity
│   │   │   ├── DDPModel
│   │   │   ├── DFNModel
│   │   │   ├── DomainDefintion
│   │   │   ├── EllipseComputation
│   │   │   ├── EPMModel
│   │   │   ├── FourierModel
│   │   │   ├── InputOutput
│   │   │   ├── Release
│   │   │   ├── RngStream
│   │   │   ├── Utilitaries
│   │   │   ├── MAIN.cpp
│   │   │   └── PERFORM.cpp
│   │   ├── input
│   │   │   ├── testDFN
│   │   │   ├── testEPM
│   │   │   │   ├── BoudaryCondition.txt
│   │   │   │   ├── DFN.txt
│   │   │   │   ├── DomainProperties.txt
│   │   │   │   ├── ElectrodesPositions.txt
│   │   │   │   ├── FileName.txt
│   │   │   │   ├── Simu.txt
│   │   │   │   └── WavesValues.txt
│   │   ├── Library
│   │   ├── output
│   │   └── results
│   ├── batch_fwd
│   ├── fwdFunctions.py
│   ├── main.py
│   ├── README.md
│   └── requirements.txt
```

## Input files

### For Python scripts

The file `requirements.txt` contains the list of the packages needed to run Python scripts. They can be installed
in an Anaconda environment or Python virtual environment.

The folder `fielddata` should contain the coordinate text files of the electrode 
configurations relative to the application site. 
By default, the files relative to the study case describes in Lelimouzin et al. (2024) are used.
In `data_xxyy_R4_noratio.txt`, for electrodes marked A, B, M, and N, each column of
the file corresponds, respectively, to the coordinates xA, xB, xM, xN, yA, yB, yM, and yN (x for the 
surface and y for depth). There are as many lines as there are electrical experiments (electrode quadrupole).
In the file `data.dat`, the first row corresponds to the number of experiments (measurements).
Then each row is an experiment. 
The first column provides the index of the experiment.
Columns 2, 3, 4, and 5, correspond to the index of the electrode involve in the measurement.
Column 6 and 7 are not useful in the forward procedure.
The index of the electrode correspond to the row in `electrodes.dat`, which is unused in
the forward procedure but necessary for the inversion procedure. The latter contains the coordinates
of each electrode.


In the script `main.py`, the section ”Set parameters” should be adapted by the user 
according to their study case (number of electrodes, name of the data files, etc.).
The compilation option should be changed: 1 for the library and 2 for the executable.
The script `Farum_elecConfig.py` has been developed for a specific application site and should be adapted to the
user’s new context.

### For C++ scripts

Go to the folder `src/input/`. The folders `testDFN`and `testEPM` contain the same text files with different
parameters. Change ”DomainProperties.txt”, ”ElectrodesPositions.txt”, ”Simu.txt”, and ”BoundaryConditions.txt” files
accordingly to your simulation.
Refer to the documentation of ERT_DDP_2.5D available in the `documentation` folder.

A debug option has been added compared to the original version. It displays more information relative
to matrix or vector dimensions. To activate the option, change 0 to 1
in `src/input/testDFN/` and/or `src/input/testEPM/` at the end of "Simu.txt" file.

## Compile

### Executable (ElectricalResistivity)

To compile the C++ part of ERT_DDP_2.5D as an executable uncomment the writing 
of the matrices in text files in `electrical_resistivity/FourierModel/FourierModelFunctions.cpp`

Do not forget to change the options in `main.py` (2 for the executable).

```bash
cd src/electrical_resistivity/Release
make clean
make exec -j 20
```
The results from the C++ scripts are written in text files stored in the `src/output` folder. They are then
loaded in the Python scripts to calculate the apparent electrical resistivity from the potential difference in the Fourier
domain.

To directly execute the codes from a terminal, go to `src/electrical_resistivty/Release`, then run the code with
the command line: `./ElectricalResistivity`.

### Library (PERFORM.so)

To compile the C++ part of ERT_DDP_2.5D as a library comment the writing of 
the matrices in text files in `electrical_resistivity/FourierModel/FourierModelFunctions.cpp`

Do not forget to change the options in `main.py` (1 for the library).

```bash
cd src/electrical_resistivity/Release
make clean
make lib -j 20
```
Results are directly transfer from C++ to Python thanks to the PyBind11 package and stored in dictionaries.

## Output

If the executable version of the C++ code is selected, the text files containing the predicted potential 
in the Fourier domain are saved in the folder `src/output`. 

For both version of compilation, the folder `src/results` contains the vector of apparent electrical resistivity
calculated at the end of the Python post-processing, i.e. at the end of the forward procedure.

## Test case

A test example is provided in the folder `example`. The latter contains all input files for both C++ and 
Python scripts, as well as the output files (C++ text files from the executable version) and 
the results files (simulated data).
The test corresponds to a cross-borehole ERT experiments with 128 experiments and 2008 measurements.
It is a synthetic case with a DFN containing 18 fractures.

We recommend to use the executable to make tests. 

WARNING: By default, the input files relative to this example are already used in the code.

## New options to generate a DFN (no electrical simulations)

Compared to the version developed by Caballero-Sanz et al. (2017), two options to generate synthetic DFNs have been added
in the C++ code: "FRACTAL" and "POWER_LAW_ANGLE".
An example is provided in the folder `example/power_law_angle/` for the "POWER_LAW_ANGLE" option.

Changes in the input files are required to use these options.
At the end of `Simu.txt`, add a file name for the output DFN.
`Simu.txt` should be set to `DFN_GENERATION`.
The input file `DFN.txt` should be changed as follows:

If DFN_option = "FRACTAL":
- r_min
- C_density
- D_dim
- angle for the first family of fractures
- angle of the second family of fractures 
- fracture aperture 
- conductivity in the fractures
- seed

If DFN_option = "POWER_LAW_ANGLE":
- power-law exponent (a)
- percolation coefficient (p)
- angle for the first family of fractures
- angle of the second family of fractures
- fracture aperture 
- conductivity in the fractures
- seed
- minimal length value (l_min)

For the other options, please refer to the documentation of ERT_DDP_2.5D.

To execute the C++ code in order to generate a DFN, the package OpenGL is needed.

- Uncomment lines 8 to 32 of `src/electrical_resistivity/Utilitaries/Visualization/VisuStructures.cpp`
- Uncomment all lines in `src/electrical_resistivity/Utilitaries/VisualizationVisuStructures.h`
- Uncomment line 8 in `src/electrical_resistivity/Release/objects.mk` and comment line 9
- Uncomment line 185 in `src/electrical_resistivity/InputOutput/Parameters.cpp` and comment line 186
- Uncomment lines 56, 57 and 70 in `src/electrical_resistivity/MAIN.cpp`
- Uncomment lines 17 to 23 of `src/electrical_resistivity/DFNModel/DFNGeneration.cpp`
- Uncomment line 15 of `src/electrical_resistivity/DFNModel/DFNGeneration.h`
- Compile the code in executable mode from the `src/electrical_resistivity/Release` folder

## FAQ

### PERFORM not found

This error means that the C++ library `PERFORM.so` is not find by Python scripts. First, check if the
library `PERFORM.so` has been generated in the folder `src/electrical_resistivty/Release`. If not, follow the steps of the next item
”Error in Parameters.cpp”. Then try to compile the C++ code again.

If the error persists, check the Python version (python3.13 was used on the local computer, and 3.10
on the computing cluster). Check the file PERFORM.d, the version should be the same as the one used in
your Anaconda environment or the version installed on your computer. To check the Python version,
use the command: python3 -V (ubuntu). Depending on the version, some errors may occur with the
Python thread state.

### Error in Parameters.cpp

If compilation fails and an error occurs on the parameters, comment `#include <filesystem>` (write
‘//‘ in front of the line) at the beginning of `src/electrical_resistivty/InputOutput/Parameters.cpp`. Compile the C++ code
again. Another source of error could be the architecture of the repository. Make sure to have kept the
original structure of the folders to avoid file paths issues.

### Fatal error: gl.h: No such file or directory

The OpenGL package may be missing in computing clusters. If it is not possible to install OpenGL,
a solution is to comment all references to the package in the C++ scripts, i.e. all the scripts corresponding to the visualization of the results. By default, TransMontpInv is set with OpenGL references
commented.
To comment, add ”//” in front of the lines in `src/electrical_resistivity/`:

- Comment all ”.cpp” and ”.h” files in the directory Utilitaries/Visualization.
- In DFNModel, comment lines 15 to 22 in DFNGeneration.cpp, and in DFNGeneration.h comment
DFNVisu DFNGeneration(Parameters,Domain);.
- In PERFORM.cpp comment lines 75, 78 to 91, 153, 154, 165, 166, and 171.
- In Release/objects.mk remove -lGL and -lglut.

WARNING: By default all references to OpenGL have already been commented. 

### KeyError and double free() error (WARNING: unsolved issues)

These errors do not affect the results of the procedure. 
They appear in the inversion procedures that called the forward scheme.

A KeyError may occur when loading the C++ results in Python functions using the libraries 
PERFORM.so at a random iteration. There is a problem of dimension between the C++ objects and the
expected format in Python scripts. Simulations may also crash at a random iteration with a message
such as double free() error. Simulations may stay stuck in infinite loops when the latter issues
occur. Regularly check the state of long simulations.

A debug option, has been implemented in the C++ scripts to better understand these
errors. 

These issues have not yet been resolved. We are currently working on them. 
So far, if you notice the bug, stop the simulation and restart it.