	MODULE unit_conversion_library
		
	CONTAINS

!-----------------------------------------------------

!			ENERGY

!-----------------------------------------------------
	
	function hertz_to_hartree(x)
	
	implicit none	
	double precision :: hertz_to_hartree
	double precision :: x
	
	hertz_to_hartree = x * 1.519829846006d-16
	
	end function hertz_to_hartree
!-----------------------------------------------------
	
	function hartree_to_hertz(x)
	
	implicit none	
	double precision :: hartree_to_hertz
	double precision :: x
	
	hartree_to_hertz = x / 1.519829846006d-16
	
	end function hartree_to_hertz	
!-----------------------------------------------------	

	function joule_to_hartree(x)
	
	implicit none	
	double precision :: joule_to_hartree
	double precision :: x
	
	joule_to_hartree = x * 2.29371269d17
	
	end function joule_to_hartree
!-----------------------------------------------------

	function wavenumber_to_hartree(x)
	
	implicit none	
	double precision :: wavenumber_to_hartree
	double precision :: x
	
	wavenumber_to_hartree = x / 219474.63067d0
	
	end function wavenumber_to_hartree
	
!-----------------------------------------------------

	function hartree_to_wavenumber(x)
	
	implicit none	
	double precision :: hartree_to_wavenumber
	double precision :: x
	
	hartree_to_wavenumber = x * 219474.63067d0
	
	end function hartree_to_wavenumber
	
!-----------------------------------------------------

!			DISTANCE

!-----------------------------------------------------	

	function meter_to_au(x)
	
	implicit none	
	double precision :: meter_to_au
	double precision :: x
	
	meter_to_au = 1.889726133921252d10 * x
	
	end function meter_to_au
!-----------------------------------------------------
	
	function au_to_meter(x)
	
	implicit none	
	double precision :: au_to_meter
	double precision :: x
	
	au_to_meter = x / 1.889726133921252d10
	
	end function au_to_meter
!-----------------------------------------------------
	
	function angstrom_to_au(x)
	
	implicit none	
	double precision :: angstrom_to_au
	double precision :: x
	
	angstrom_to_au = 1.889726133921252 * x
	
	end function angstrom_to_au
!-----------------------------------------------------	
	
	function au_to_angstrom(x)
	
	implicit none	
	double precision :: au_to_angstrom
	double precision :: x
	
	au_to_angstrom = x / 1.889726133921252
	
	end function au_to_angstrom
!-----------------------------------------------------
	
	function nanometer_to_au(x)
	
	implicit none	
	double precision :: nanometer_to_au
	double precision :: x
	
	nanometer_to_au = 18.89726133921252 * x
	
	end function nanometer_to_au	
!-----------------------------------------------------
	
	function au_to_nanometer(x)
	
	implicit none	
	double precision :: au_to_nanometer
	double precision :: x
	
	au_to_nanometer = x / 18.89726133921252
	
	end function au_to_nanometer			
!-----------------------------------------------------
	
	function inversenanometer_to_au(x)
	
	implicit none	
	double precision :: inversenanometer_to_au
	double precision :: x
	
	inversenanometer_to_au = x / 18.89726133921252
	
	end function inversenanometer_to_au	
!-----------------------------------------------------
	
	function au_to_inversenanometer(x)
	
	implicit none	
	double precision :: au_to_inversenanometer
	double precision :: x
	
	au_to_inversenanometer = x * 18.89726133921252
	
	end function au_to_inversenanometer			
!-----------------------------------------------------

!			TIME

!-----------------------------------------------------
	
	function second_to_au(x)
	
	implicit none	
	double precision :: second_to_au
	double precision :: x
	
	second_to_au = 4.134137337414122d16 * x
	
	end function second_to_au
!-----------------------------------------------------
	
	function au_to_second(x)
	
	implicit none	
	double precision :: au_to_second
	double precision :: x
	
	au_to_second = x / 4.134137337414122d16
	
	end function au_to_second	
!-----------------------------------------------------	

!			MASS

!-----------------------------------------------------	

	function amu_to_au(x)
	
	implicit none	
	double precision :: amu_to_au
	double precision :: x
	
	amu_to_au = x * 1822.888479031408
	
	end function amu_to_au		
!-----------------------------------------------------	

!			MAGNETIC FIELD

!-----------------------------------------------------	

	function gauss_to_au(x)
	
	implicit none	
	double precision :: gauss_to_au
	double precision :: x
	
	gauss_to_au = x * 4.254382547308656d-10
	
	end function gauss_to_au	
!-----------------------------------------------------	

!			POWER

!-----------------------------------------------------	

	function watt_to_au(x)
	
	implicit none	
	double precision :: watt_to_au
	double precision :: x
	
	watt_to_au = x * 5.548225828243671
	
	end function watt_to_au		
	
!-----------------------------------------------------	

!-----------------------------------------------------	

!		DIPOLE MOMENT

!-----------------------------------------------------	

	function debye_to_au(x)
	
	implicit none	
	double precision :: debye_to_au
	double precision :: x
	
	debye_to_au = x * 0.393430165d0 
	
	end function debye_to_au		

!-----------------------------------------------------	

	function au_to_debye(x)
	
	implicit none	
	double precision :: au_to_debye
	double precision :: x
	
	au_to_debye = x / 0.393430165d0 
	
	end function au_to_debye
	
!-----------------------------------------------------	
	

	END MODULE unit_conversion_library
	

	
	
	
	
	
