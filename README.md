# Discrete Variable Representation (DVR) for Rovibrational structure of a closed-shell diatomic molecule

This repository contains code that computes the rovibrational structure of a diatomic molecule in a given closed-shell electronic potential curve. The eigensystem is obtain by using Discrete Variable Representation (DVR).

## Overview

The DVR Hamiltonian for a given potential energy curve + rotational structure is implemented using a Fourier basis and an uniform grid approach, it follows reference [J. Chem. Phys. 96(3) 1982 (1992)](https://aip.scitation.org/doi/10.1063/1.462100). For radial coordinates, the represetation uses eigenfunctions of a particle-in-a-box for a semi-infinite interval (0,infinity) and the mass is included explicitly in the kinetic energy operator.

### Input
  - **Electronic potential energy curve file**: it should have two columns, 1st: radial coordinate and 2nd: potential energy
  - **Input file**: DVR step, reduced mass of the molecule in amu, maximum values of the rotational and vibrational quantum number to be considered.

### Output
  - Eigensystem written in an external files.
  - Scalar product of the wavefunctions to check orthonormality.
  - Comparition between the Dumham expation series (only for LiCs) and the eigenenergies obtained.

## Main program

The main program *rovib.f90* is set to computed the rovibrational structure of a LiCs in the <img src="https://latex.codecogs.com/gif.latex?X^1\Sigma^+"> electronic state, with <img src="https://latex.codecogs.com/gif.latex?J_{max}=10"> and <img src="https://latex.codecogs.com/gif.latex?\nu_{max}=10">. The potential energy curve for this state is shown in the figure below:

[Figure]

The main program uses other routines to interpolate the curve values (*spline.90*) integrate the eigenfunctions to check orthonormality (*double_integral.f90*) and change measurement units (*unit_convertion.f90*).

### Results folder

This folder has the output files from the LiCs example.

### Use for other examples

In order to run other examples, you must supply the potential energy curve file (*molecule_PEC.in*) and set the input values for your system.
  
### How to run it

It can be run normally as a regular fortran program or use Python and run ``f2py -c rovib.f90 -m rovib`` in the terminal
  
## References

This program was written by Felipe Herrera and modified by Vanessa Olaya as a part of her Master thesis project.
