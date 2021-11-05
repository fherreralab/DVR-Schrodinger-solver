# About
The code contained in this repository solves the following wavefunction which describes the movement of a radial particle. It's a Schrödinger's equation that finds the molecule states when the nuclei and electrons are in specific positions.

<p align="center"><img src="https://latex.codecogs.com/svg.latex?{\left[-\frac{\hbar^2}{2\mu}\left(\frac{d^2}{dr^2%20}\right)+\frac{J\left(J+1\right)\hbar^2}{2\mu%20r^2}+V(r)\right]S(r)=E_{int}S(r)}"></p>

Where:
- <img src="https://latex.codecogs.com/gif.latex?\hbar">'s value is 1.
- <img src="https://latex.codecogs.com/gif.latex?\mu"> is the reduced mass of the system.
- <img src="https://latex.codecogs.com/gif.latex?r"> is the current grid point.
- <img src="https://latex.codecogs.com/gif.latex?J"> is the molecule's rotational number.
- <img src="https://latex.codecogs.com/gif.latex?V(r)"> is the potential energy curve.
- <img src="https://latex.codecogs.com/gif.latex?S(r)"> is the .
- <img src="https://latex.codecogs.com/gif.latex?E_{int}"> are the rotational energies that appear in the code as ``rovib_energies(rot,j) = eigenvalues(j)``.

## Analyzing the equation
In order to understand this equation, we first analyze the following:
### Movement
<p align="center"><img src="https://latex.codecogs.com/gif.latex?-\frac{\hbar^2}{2\mu}\left(\frac{d^2}{dr^2}\right)"></p>

This the kinetic energy solved numerically with the DVR subroutine, which takes into account a one-dimensional quantum system where <img src="https://latex.codecogs.com/gif.latex?r"> is restricted to a given interval as stated in [ paper ].

### Rotation
<p align="center"><img src="https://latex.codecogs.com/gif.latex?\frac{J\left(J+1\right)\hbar^2}{2\mu%20r^2}"></p>

This takes care of the rotational part of our equation and is equivalent to the following Fortran code inside *rovib.f90*: ``NN*(NN+1.d0)/(2.d0*mass*grid(i)*grid(i)``.

### Vibration

<p align="center"><img src="https://latex.codecogs.com/gif.latex?V(r)"></p>

This is, as stated earlier, the potential energy curve which is provided in the *LiCs_PEC.in* file.

## How it works

### Potential Energy Curve
There should be a file called *[molecule]_PEC.in* in the same folder as the rovib.f90 file, where *[molecule]* must be replaced with the molecule name. This file must have two columns: <img src="https://latex.codecogs.com/gif.latex?r">, which should be equally spaced, and the corresponding <img src="https://latex.codecogs.com/gif.latex?V(r)"> for each one of those values.

The Potential Energy Curve is a required input for our code. As an example, we have included the  *LiCs_PEC.in* file which contains these two columns for a LiCs molecule.

<p align="center">[ figure ]</p>

### DVR matrix subroutine

By doing DVR, we numerically obtain the energy of our molecule's Hamiltonian <img src="https://latex.codecogs.com/gif.latex?H_{mol}">. Following [ paper ]'s Appendix A, which goes through calculating a simple generic DVR step by step, we achieve the following:

<p align="center"><img src="https://latex.codecogs.com/gif.latex?-\frac{\hbar^2}{2\mu}\left(\frac{d^2}{dr^2}\right)\longrightarrow\frac{\hbar^2}{2\mu\Delta%20r^2}(-1)^{i-j}\left\{\begin{array}{lr}\frac{\pi^2}{3}-\frac{1}{2i^2},&\text{for%20}i=j\\\frac{2}{(i-j)^2}-\frac{2}{(i+j)^2},&\text{for%20}i\neq%20j\end{array}\right\}"></p>
<p align="center">with <img src="https://latex.codecogs.com/gif.latex?r_i=i\Delta%20r"> given <img src="https://latex.codecogs.com/gif.latex?i=1,%202,%20\dots">.</p>

In our subroutine, we obtain a matrix using this formula, where <img src="https://latex.codecogs.com/gif.latex?i"> and <img src="https://latex.codecogs.com/gif.latex?j"> are the matrix's rows and columns respectively. When <img src="https://latex.codecogs.com/gif.latex?i=j">, we're positioned over our matrix's diagional, where the rotational and electronic energy are added as follows:

```@fortran
hamiltonian(i,i) = (0.5d0/m)*(step**(-2)) * (pi*pi/3.d0 - 0.5d0/dble(i*i)) + potential(i)
```

For the non-diagonal elements where <img src="https://latex.codecogs.com/gif.latex?i\neq%20j">, only kinetic terms are stored:

```@fortran
hamiltonian(i,j) = (0.5d0/m)*(step**(-2)) * (2.d0/(dble(i-j)**2)-2.d0/(dble(i+j)**2)) * (-1.d0)**(i-j)
```

Variables used in this subroutine:
- ``NN`` is <img src="https://latex.codecogs.com/gif.latex?J">, the rotational number of the molecule, which must be a positive number from 0 onwards.
- ``NGP`` are the grid points that define our Hamiltonian's size.

<!--
NGP: Number of Grid Points
Grid: A qué r representa cada punto.
Potential corresponde al valor en el archivo V(r)
J: Número rotacional de la molécula, es ``NN*(NN+1.d0)/(2.d0*mass*grid(i)*grid(i))``.
\hbar es 1 en nuestro código de Fortran.
Cada proyecto debe tener tres escenarios de solución
-->
