	
	program rovib
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

	use unit_conversion_library

	implicit none

	! Main program variables
	integer :: NGP, i, j, rot, inputstatus, vmax, NV, v, vp, N, Np, NR, Nomega, PEC, PDM, instatus
	integer, parameter :: kmax = 5000
	double precision :: NN, Nmin, Nmax, mass, dissociation_energy, rmin, rmax, step, RMAT(kmax,2)
	double precision, allocatable :: grid(:), potential(:), hamiltonian(:,:), eigenvalues(:), eigenvectors(:,:)
	double precision, allocatable :: dipole_vibration(:,:), rovib_wavefunctions(:,:,:), rovib_energies(:,:)
	double precision :: dipole_integral, DIP(kmax,2), omega,alpha, omega_min, omega_max, omega_step
	double precision :: alpha_imag, alpha_parallel, alpha_perpend, M, Mp, q
	double precision :: w_si, d_si, w_pi, d_pi, w_1, d_1, w_2, d_2, w_11, d_11, w_22, d_22, w_111, d_111, w_222, d_222
	character :: mol*3, rotational*2, magnetic*4, moment*4, moleculestate*16
		
	! Variables for LAPACK diagonalization
	double precision, allocatable :: WORK(:), A(:,:)
	integer, allocatable :: IWORK(:)
	integer :: LWORK, LIWORK, INFO
	
	! Variables for cubic spline interpolation
	integer :: mesh
	double precision :: x, y
	double precision, allocatable :: knots(:,:), moments(:)
		
	real*8 threejsymbol
	
	! Read input

print*, 'Reading'

	open(10,file='rovib_pp.in')
	
	read(10,*)!rmin, rmax, step
	read(10,*) rmin, rmax, step
	read(10,*)!reduced mass (amu) 	Dissociation Energy (cm-1)
	read(10,*)! LiCs
	read(10,*)! RbCs
	read(10,*) mass, dissociation_energy ! KRb
	read(10,*)! LiRb
	read(10,*)! electronic scalar polarizability (au)	electronic tensor polarizability (au)
	read(10,*)! LiCs
	read(10,*)! RbCs
	read(10,*) alpha_parallel, alpha_perpend ! KRb
	read(10,*)! LiRb
	read(10,*)!frequency min (GHz)	Frequency max (GHz)  frequency step (GHz)
	read(10,*) omega_min, omega_max, omega_step
	read(10,*)! Jmax, vmax, q
	read(10,*) Nmax, vmax, q

	! Molecule verification

	if (dabs(dissociation_energy-5875.455).lt.1.d-3) then
	mol = 'LiCs'
	else if (dabs(dissociation_energy-3836.1).lt.1.d-3) then
	mol = 'RbCs'
	else if (dabs(dissociation_energy-4217.91).lt.1.d-3) then
	mol = 'KRb'
	else if (dabs(dissociation_energy-5927.9).lt.1.d-3) then
	mol = 'LiRb'
	else
	STOP 'Error: molecule specification not recognized'
	end if

print*, mol

	if ((rmax.gt.300.d0).or.(rmin.lt.5.d0)) stop 'Error: grid boundaries outside potential interval!'

	! Convertion to atomic units
	
	mass = amu_to_au(mass)
	dissociation_energy = wavenumber_to_hartree(dissociation_energy)
	omega_min = hertz_to_hartree(omega_min*1.d9)	
	omega_max = hertz_to_hartree(omega_max*1.d9)
	omega_step = hertz_to_hartree(omega_step*1.d9)
	
	! Define parameters
	
	NGP = nint(dabs(rmax-rmin)/step) + 1			! number of grid points
	NR = nint(Nmax) + 1					! number of rotational states (for truncation)
	NV = vmax + 1						! number of vibrational states (for truncation)
	Nomega = nint(dabs(omega_max-omega_min)/omega_step) + 1 ! number of omega steps

	! Read potential curve

	open(20,file=mol//'_PEC.in')

	mesh=0
	do
		mesh = mesh + 1
		read(20,*,iostat=inputstatus) (RMAT(mesh,j), j = 1, 2)
		if (inputstatus.lt.0) exit
	end do

	PEC = mesh - 1

	if (PEC.eq.0) stop 'Error: Potencial energy file is empty'							
	
	allocate(moments(PEC))
	allocate(knots(PEC,2))

	do i=1, PEC
		knots(i,1)= RMAT(i,1)			     ! Already in au
		knots(i,2)= wavenumber_to_hartree(RMAT(i,2)) ! Convert to atomic units
	end do

print*, 'Cubicspline potential'

	call cubicsplines(knots,PEC,0.d0,0.d0,moments) ! Compute moments

	! Construction of the DVR Hamiltonian
	
	allocate(grid(NGP))
	allocate(potential(NGP))
	allocate(hamiltonian(NGP,NGP))			
	allocate(dipole_vibration(NV,NV))	
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

	do rot = 1, NR ! Begin loop over rotational states

		Nn = (rot-1)

print*, Nn

		do i = 1, NGP
			x = rmin + dble(i-1)*step
			grid(i) = x							! Set up uniform grid
			call splineinterpol(knots,moments,PEC,x,y) 	
			potential(i) = y + Nn*(Nn+1.d0)/(2.d0*mass*grid(i)*grid(i))	! Interpolated electronic + rotational
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
		
		do i = 1,NGP
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

	! Interpolation of dipole moment function

	open(40,file=mol//'_ETDM.in')

	mesh=0
	do 
		mesh = mesh + 1
		read(40,*,iostat=instatus) (DIP(mesh,j), j = 1, 2) ! All is in au
		if (instatus.lt.0) exit
	end do

	PDM = mesh - 1

	if (PDM.eq.0) stop 'Error: Permanent dipole moment file is empty'							

	allocate(moments(PDM))
	allocate(knots(PDM,2))

	do i=1, PDM
		knots(i,1)= DIP(i,1)			
		knots(i,2)= DIP(i,2)
	end do

print*, 'Cubicspline dipole'

	call cubicsplines(knots,PDM,0.d0,0.d0,moments) 	

print*, 'Polarizability'
			
	! Calculate the rotational polarizability

	! Read parameter for effective polarizabiliaties for frequencies greater than 4500 GHz

	open(13,file='mol_eff_pol.in')
	
	read(13,*)!w_si(cm⁻¹)	 d_si(au)	w_pi(cm⁻¹)	d_pi(au)
	read(13,*)! LiCs
	read(13,*)! RbCs
	read(13,*) w_si, d_si, w_pi, d_pi ! KRb
	read(13,*)! LiRb

	w_1 = 157105	!K
	d_1 = 0.304	!K
	w_2 = 199440	!K
	d_2 = 1.383	!K

	w_11 = 165912	!Rb
	d_11 = 2.728	!Rb
	w_22 = 362918	!Rb
	d_22 = 1.401	!Rb

	w_111 = 142044	!Cs
	d_111 = 4.488	!Cs
	w_222 = 553515	!Cs
	d_222 = 1.860	!Cs
	
	w_si = wavenumber_to_hartree(w_si)
	w_pi = wavenumber_to_hartree(w_pi)
	w_1 = wavenumber_to_hartree(w_1)
	w_2 = wavenumber_to_hartree(w_2)
	w_11 = wavenumber_to_hartree(w_11)
	w_22 = wavenumber_to_hartree(w_22)
	w_111 = wavenumber_to_hartree(w_111)
	w_222 = wavenumber_to_hartree(w_222)

	v = 0
	N = 1
	M = 1
	omega = 0	
	
	write (rotational,'(I2.2)') int(N)
	write (magnetic,'(F4.1)') M
	write (moment,'(F4.1)') q
	
	moleculestate = mol//'_'//rotational//'_'//magnetic//'_'//moment
	
	open(41,file='./molpol1/real_'//moleculestate)
	open(42,file='./molpol1/imag_'//moleculestate)

	do i = 1, Nomega

		omega = omega_min + dble(i-1)*omega_step

print*, omega

		if (omega.le.6.83923d-4) then
			call microwave_polarizability(real(q),real(v),real(N),real(M),NV,NR,NGP,PDM,moments,knots,grid,rovib_energies,&
								rovib_wavefunctions,alpha_parallel,alpha_perpend,omega,alpha,alpha_imag)



			write(41,'(9D25.15)') omega, alpha
			write(42,'(9D25.12)') omega, alpha_imag

		end if

		if (omega.gt.6.83923d-4) then
			alpha = ((2.d0 * w_si * d_si**2) / (w_si**2 - omega**2)) + ((2.d0 * w_pi * d_pi**2) / (w_pi**2 - omega**2)) +&
				((2.d0 * w_1 * d_1) / (w_1**2 - omega**2)) + ((2.d0 * w_2 * d_2) / (w_2**2 - omega**2)) +&
				((2.d0 * w_11 * d_11) / (w_11**2 - omega**2)) + ((2.d0 * w_22 * d_22) / (w_22**2 - omega**2))
			write(41,'(9D25.15)') omega, alpha
		end if

	end do

end program rovib


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
	integer interval, i, j, NGP, indexA,indexB
	double precision  step, pi, xmin, x, omega, m, De, aa, re
	double precision hamiltonian(NGP,NGP),grid(NGP), potential(NGP)	
	
	pi = dacos(-1.d0)
	
	do i = 1, NGP					
		do j = 1, NGP
		  	if (i.ne.j) then 
		  		hamiltonian(i,j) = (0.5d0/m)*(step**(-2)) * (2.d0/(dble(i-j)**2)-2.d0/(dble(i+j)**2)) * (-1.d0)**(i-j)
		  	end if					
		end do						
		hamiltonian(i,i) = (0.5d0/m)*(step**(-2)) * (pi*pi/3.d0 - 0.5d0/dble(i*i)) + potential(i)
	end do	   
	 		
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
	d = -0.6965d-4	! Y4,0
	e = -0.65721d-9	! Y7,0
	f = 0.163d-13	! Y10,0
	g = -0.499d-15	! Y11,0
	h = 0.443d-17	! Y12,0
		
	do i = 1, NVS
		v = dble(i-1)
		x = v + 0.5d0
		energies(i) = a*x + b*x**2 + c*x**3 + d*x**4 +e*x**7 + f*x**10 + g*x**11 + h*x**12
	end do
	
	end subroutine
			
!==================================================================================================================================================

	double precision function dipole(R,mesh,moments,knots)
	! Program interpolates the electric dipole function d(R) for heteronuclear alkali-metal dimer
	! using data from J. Chem. Phys. 122, 204302 (2005)
	! RECEIVES INPUT from subroutine SET_DIPOLESPLINES(mesh,moments,knots)
	! Restricted for LiCs molecules in singlet sigma state, but generalization to other species is trivial.
	! A linear extrapolation is performed for R < 5 a_0 to the short range region and also for R > 20 a_0. 
	
	use unit_conversion_library
	implicit none
	double precision a,b, c,d
	integer mesh
	double precision x, y, R, knots(mesh,2), moments(mesh)	
	
	! parameters for LiCs short range data extrapolation R < 5
	a = (4.7952797d0-4.7773341d0)/(5.1006711d0 - 5.d0)
	b = 4.7773341d0 - a*5.d0
	
	! parameters for LiCs long range extrapolation  R > 20
	c = (-0.036754709d0 + 0.036062174d0)/(20.d0 - 19.899329d0)	
	d = -0.036062174 - c*19.899329
	
	x = R		
		
	if (x.lt.5d0) then	
		dipole = a * x + b
		
	else	if ((x.ge.5d0).and.(x.le.30.d0)) then	
		call splineinterpol(knots,moments,mesh,x,y) 		
		dipole = y
		
	else if (x.gt.30.d0) then
		dipole = c * x + d
	
	end if	
	end function dipole
	
!==================================================================================================================================================

	double precision function dipole_integral(N,Np,v,vp,NR,NV,NGP,grid,mesh,moments,knots,rovib_wavefunctions)		
	! Function evaluates the dipole integral <v(J)|d|v'(J')>
	! assumes Nmin = 0 and vmin = 0.
	
	implicit none
	integer  i, N, Np, v, vp, NR,NV,NGP, rotA,rotB,vibA, vibB, mesh
	double precision grid(NGP),rovib_wavefunctions(NR,NGP,NV), integral, dipole, moments(mesh), knots(mesh,2),dofR
	
	!variables for cubic integration
	integer ntab, ia, ib
	real*8 norm, error
	real*8, allocatable:: ftab(:), xtab(:)
	
	ntab = NGP
	ia = 1
	ib = NGP
	
	allocate(ftab(ntab))
	allocate(xtab(ntab))
	
	rotA = N + 1
	rotB = Np + 1 
	vibA = v + 1
	vibB = vp + 1

	do i = 1, ntab
		xtab(i) = grid(i)
		dofR = dipole(grid(i),mesh,moments,knots)
		ftab(i) = rovib_wavefunctions(rotA,i,vibA)*rovib_wavefunctions(rotB,i,vibB)*dofR
	end do

	call cubint ( ftab, xtab, ntab, ia, ib, norm, error )	
		
	dipole_integral = dble(norm)

	end function dipole_integral

!==================================================================================================================================================

	subroutine microwave_polarizability(q,v,N,M,Nvib,Nrot,NGP,mesh,moments,knots,grid,energies,&
								wavefunctions,alpha_parallel,alpha_perpend,omega,alpha,alpha_imag)
	! Program evaluates the dynamic polarizability at REAL microwave frequencies using 
	! a sum-over-states method. The sumation is restricted for states within the ground electronic state 
	! The polarizability is evaluated for the state (vJM=0) of a diatomic molecule. Generalization is straighforward.
	! Read the static electronic polarizabilities (scalar and tensor) from file and add the contributions
	! to the total polarizability
	
	use unit_conversion_library
	implicit none	
	integer Nvib,Nrot,NGP, mesh, number_terms, indexV,indexR,i,j,k,Mmin,Mmax,Mnum
	double precision grid(NGP), wavefunctions(Nrot,NGP,NVib), energies(Nrot,Nvib),prefactor,prefactor1,omega
	double precision knots(mesh,2), moments(mesh), alpha, integral, energy_factor, numerator,denominator, alpha_rot, alpha_rot_imag
	double precision dipole_integral, alpha_imag, alpha_parallel, alpha_perpend, denominator_imag, alpha_ele
	real v,N,M,vp,Np,Mp,q
	real*8 threejsymbol
		
	alpha = 0.d0
	alpha_imag = 0.d0
	alpha_rot = 0.d0
	alpha_rot_imag = 0.d0
	alpha_ele = 0.d0
	indexV = int(v)+1
	indexR = int(N)+1
		
	do i = 1, Nvib

		vp = real(i-1)
		
		do j = 1, Nrot

			Np = real(j-1)

			Mmin = 	-int(Np)
			Mmax = int(Np)
			Mnum = (Mmax-Mmin) + 1

			integral = dipole_integral(int(N),int(Np),int(v),int(vp),NRot,NVib,NGP,grid,mesh,moments,knots,wavefunctions)

			do k = 1, Mnum
			
				Mp = real(Mmin + (k-1))

				prefactor = (threejsymbol(N,1.0,Np,-M,-q,Mp)**2)*(threejsymbol(N,1.0,Np,0.0,0.0,0.0)**2) &
						*dble(2.0*N+1.0)*dble(2.0*Np+1.0)

				if (dabs(prefactor-0.d0).gt.1.d-12) then

					numerator = 2.d0*(energies(j,i) - energies(indexR,indexV))
					denominator = (energies(j,i) - energies(indexR,indexV))**2 - omega**2 ! real frequencies
					denominator_imag = (energies(j,i) - energies(indexR,indexV))**2 + omega**2 ! imaginary frequencies

					energy_factor = numerator/denominator
					alpha_rot = alpha_rot + prefactor*energy_factor*integral**2
				
					energy_factor = numerator/denominator_imag
					alpha_rot_imag = alpha_rot_imag + prefactor*energy_factor*integral**2

				end if
			end do
		end do
	end do

	do j = 1, Nrot

		Np = real(j-1)

		Mmin = 	-int(Np)
		Mmax = int(Np)
		Mnum = (Mmax-Mmin) + 1

		do k = 1, Mnum
			
			Mp = real(Mmin + (k-1))

			prefactor = (threejsymbol(N,1.0,Np,-M,-q,Mp)**2)*(threejsymbol(N,1.0,Np,0.0,0.0,0.0)**2) &
					*dble(2.0*N+1.0)*dble(2.0*Np+1.0)

			if ((N.eq.1).and.(Np.eq.0)) then
				prefactor1 = 0d0
				else
				prefactor1 = (threejsymbol(N,1.0,Np,-M,-q,Mp)**2)*(threejsymbol(N,1.0,Np,0.0,-1.0,1.0)**2) &
					*dble(2.0*N+1.0)*dble(2.0*Np+1.0)
			end if

			if ((dabs(prefactor-0.d0).gt.1.d-12).or.(dabs(prefactor1-0.d0).gt.1.d-12)) then

				alpha_ele = alpha_ele + prefactor*alpha_parallel + 2.d0*prefactor1*alpha_perpend

			end if	
		end do
	end do

	alpha = alpha_rot + alpha_ele
	alpha_imag = alpha_rot_imag + alpha_ele

	end subroutine microwave_polarizability

