module rovib

use unit_conversion_library
implicit none
double precision, allocatable, public :: rovib_wavefunctions(:,:,:), rovib_energies(:,:)
private
save

public :: use_rovib, dvr_radial, dunham, orthonormal_check

CONTAINS

subroutine use_rovib(rmin, rmax, step, mass, NMax, vmax)
!-----------------------------------------------------------------------------------------------------------------	
! 	General program to calculate the rovibrational structure of a diatomic molecule
!	in a given closed-shell electronic potential curve. E(v,J)
!-----------------------------------------------------------------------------------------------------------------
!	programer		date			comments
!-----------------------------------------------------------------------------------------------------------------
!	F. Herrera		April 2, 2013	First working version gives energies and wavefunctions.
!						Computations are done in atomic units. 
!	F. Herrera		April 3, 2013	Accurate for LiCs X state potential. Transition frequencies for
!						1st rotational and 1st vibrational quanta are equal to Dunham
!						expansion in PRA 75,042513 (2007), up to six decimal places.
!						Accuracy is expected to degrade for highest vibrational states 
!						(v > 20). Systematic improvement possible by including short and
!						long range potentials. Molecular species is selected by input.
!						Radial wavefunctions are orthonormal. 
!	F. Herrera		April 4, 2013	Evaluates the dipole integral <v'(J')|d|v(J)> between 
!						selected rovibrational states. 
!-----------------------------------------------------------------------------------------------------------------

! Main variables
integer :: NGP, i, rot, instatus, NV, NR, PEC, j
integer, parameter :: kmax = 5000
double precision :: NN, RMAT(kmax,2)
double precision, intent(in) :: rmin, rmax, step, Nmax
double precision, intent(inout) :: mass
integer, intent(in) :: vmax
double precision, allocatable :: grid(:), potential(:), hamiltonian(:,:), eigenvalues(:), eigenvectors(:,:)
character :: mol*4
	
! Variables for LAPACK diagonalization
double precision, allocatable :: WORK(:), A(:,:)
integer, allocatable :: IWORK(:)
integer :: LWORK, LIWORK, INFO

! Variables for cubic spline interpolation
integer :: mesh
double precision :: x, y
double precision, allocatable :: knots(:,:), moments(:)
	
real*8 threejsymbol

!-----------------------------------------------------------------------------------------------------------------

! Setting initial parameters for LiCs

print*, 'Setting parameters for LiCs'

mol = 'LiCs'

if ((rmax.gt.515.d0).or.(rmin.lt.4.d0)) stop 'Error: grid boundaries outside potential interval!'

! Convertion to atomic units

mass = amu_to_au(mass)

! Define parameters

NGP = nint(dabs(rmax-rmin)/step) + 1			! number of grid points
NR = nint(Nmax) + 1					! number of rotational states (for truncation)
NV = vmax + 1						! number of vibrational states (for truncation)

! Read potential curve

print*, 'Reading potential energy curve'

open(20,file=mol//'_PEC.in')

mesh=0
do
	mesh = mesh + 1
	read(20,*,iostat=instatus) (RMAT(mesh,j), j = 1, 2)
	if (instatus.lt.0) exit
end do

PEC = mesh - 1

if (PEC.eq.0) stop 'Error: Potencial energy file is empty'							

allocate(moments(PEC))
allocate(knots(PEC,2))

do i=1, PEC
	knots(i,1)= RMAT(i,1)				! Already in au
	knots(i,2)= wavenumber_to_hartree(RMAT(i,2))	! Convert to atomic units
end do

call cubicsplines(knots,PEC,0.d0,0.d0,moments)		! Compute moments

! Construction of the DVR Hamiltonian

print*, 'Construction of the DVR Hamiltonian'

allocate(grid(NGP))
allocate(potential(NGP))
allocate(hamiltonian(NGP,NGP))			
allocate(rovib_wavefunctions(NR,NGP,NV))
allocate(rovib_energies(NR,NV))
LWORK = 1 + 6 * NGP + 2 * NGP**2
LIWORK = 3 + 5 * NGP	
allocate(WORK(LWORK))
allocate(IWORK(LIWORK))		
allocate(eigenvalues(NGP))
allocate(eigenvectors(NGP,NGP))
allocate(A(NGP,NGP))
	
! Loop over rotational quantum number

print*, 'Begin DVR'

do rot = 1, NR 									! Begin loop over rotational states

	NN = (rot-1)

	do i = 1, NGP
		x = rmin + dble(i-1)*step
		grid(i) = x							! Set up uniform grid
		call splineinterpol(knots,moments,PEC,x,y) 	
		potential(i) = y + NN*(NN+1.d0)/(2.d0*mass*grid(i)*grid(i))	! Interpolated electronic + rotational
	end do
	
	! Evaluate DVR kinetic + potential energies (in atomic units)

	call DVR_radial(mass,step,NGP,grid,potential,hamiltonian)


	! Diagonalize the Hamiltonian using LAPACK

	A = hamiltonian	
	call DSYEVD ('V','U', NGP, A, NGP,eigenvalues, WORK, LWORK, IWORK,LIWORK,INFO)	
	eigenvectors = A

	! Normalize the eigenvectors by the grid step

	eigenvectors = eigenvectors / dsqrt(step)

	! Write eigenvectors to 3D array
	
	do i = 1, NGP
		do j = 1, NV			
		rovib_wavefunctions(rot,i,j) = eigenvectors(i,j)
		end do
	end do		
		
	do j = 1, NV			
		rovib_energies(rot,j) = eigenvalues(j)
	end do		
	

end do ! End loop over rotational states	

print*, 'End DVR'

deallocate(knots)
deallocate(moments)

end subroutine use_rovib

!==================================================================================================================================================
!==================================================================================================================================================

subroutine DVR_radial(m,step,NGP,grid,potential,hamiltonian)
! 	Program implements a DVR Hamiltonian using a Fourier basis/Uniform grid approach
! 	as described in J. Chem. Phys. 96(3) 1982 (1992). Different boundary conditions are supported.
!	The general representation in terms of particle-in-a-box eigenfunctions is specialized 
!	to the semi-infinite interval (0,infinity), appropriate for a radial coordinate. 
!	The kinetic energy operator includes the mass explicitly. Otherwise the units are hbar = 1
!	
!	INPUT: 
!		 integer "NGP": the number of grid points 
!		 double precision "step": the grid step (step)
!		 double precision "grid(NGP)": array containing the spatial grid points
!		 double precision "potential(NGP)": array containing the potential V(x_i) evaluated at
!								the grid points x_i=grid(i)
!
!	OUTPUT: double precision "hamiltonian(NGP,NGP): the Hamiltonian array!
!-----------------------------------------------------------------------------------------------
!	! Programmer	Date			Comments
!-----------------------------------------------------------------------------------------------
!	F. Herrera		Marh 25, 2013	Solves 1D Harmonic Oscillator. units h=1.
!	F. Herrera		March 26, 2013	Implements intervals #1 and #2.
!	F. Herrera		April 1, 2013	Tested against Morse oscillator. 
!						Generalized for an arbitrary 1D radial potential.
!-----------------------------------------------------------------------------------------------	

implicit none	
integer interval, i, j, indexA,indexB
double precision  pi, xmin, x, omega, De, aa, re
integer :: NGP
double precision, intent(in) :: m, step, grid(NGP), potential(NGP)
double precision, intent(out) :: hamiltonian(NGP,NGP)

pi = dacos(-1.d0)

do i = 1, NGP					
	do j = 1, NGP
		if (i.ne.j) then 
			hamiltonian(i,j) = (0.5d0/m)*(step**(-2)) * (2.d0/(dble(i-j)**2)-2.d0/(dble(i+j)**2)) * (-1.d0)**(i-j)
		end if					
	end do						
	hamiltonian(i,i) = (0.5d0/m)*(step**(-2)) * (pi*pi/3.d0 - 0.5d0/dble(i*i)) + potential(i)
end do	   

! The following loop is left to check if the output Python matrix contains the same values in its diagional
! do i = 1, NGP
!	print*, NGP, hamiltonian(i,i)
! end do
		
end subroutine DVR_radial

!==================================================================================================================================================

subroutine orthonormal_check(indexA,indexB,NGP,grid,eigenvectors)

implicit none
integer i, NGP, indexA,indexB
double precision grid(NGP), eigenvectors(NGP,NGP)

!variables for cubic integration
integer ntab, ia, ib
real norm, error
real, allocatable:: ftab(:), xtab(:)

ntab = NGP
ia = 1
ib = NGP

allocate(ftab(ntab))
allocate(xtab(ntab))

do i = 1, ntab
	xtab(i) = grid(i)
	ftab(i) = eigenvectors(i,indexA)*eigenvectors(i,indexB)
end do
	
call cubint ( ftab, xtab, ntab, ia, ib, norm, error )

write(6,'(2(A,X,F9.6,X))') 'scalar product = ', norm, 'error = ', error
write(6,*)

end subroutine orthonormal_check

!==================================================================================================================================================

subroutine dunham(NVS,energies)
! specific for (7)Li(133)Cs from PRA 75, 042513 (2007)
! Assuming no rotational excitations J = 0, the energy of the vibrational state (in cm-1) is given by 
! E(v,0) = w(v+1/2) -wx(v+1/2)^2 + wy(v+1/2)^3 + ...

implicit none		
integer NVS, v, i
double precision energies(NVS), a, b, c, d,e, f,g, h, x

a = 0.18469891d3 	! Y1,0 = w
b = -0.10002506d1	! Y2,0 = -wx
c = -0.13394d-2 	! Y3,0 = wy
d = -0.6965d-4		! Y4,0
e = -0.65721d-9		! Y7,0
f = 0.163d-13		! Y10,0
g = -0.499d-15		! Y11,0
h = 0.443d-17		! Y12,0
	
do i = 1, NVS
	v = dble(i-1)
	x = v + 0.5d0
	energies(i) = a*x + b*x**2 + c*x**3 + d*x**4 +e*x**7 + f*x**10 + g*x**11 + h*x**12
end do

end subroutine
		
!==================================================================================================================================================
end module rovib