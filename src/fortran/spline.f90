	subroutine cubicsplines(knots,n,y0prime,ynprime,moments)

!	Generates the moments of the cubic splines used to generate interpolation polynomials between two subsequent points
!	in a given table of of N coordinate pairs. We use the condition at the boundaries that the first derivatives are known.
!	See J. Stoer 'Introduction to Numerical Analysis', section 2.4.

	implicit none
	
!	Local variables
	Integer :: n,i
	Double precision:: knots(n,2),upperdiag(n-1),lowerdiag(n-1),maindiag(n),rhsvector(n),moments(n)
	Double precision:: x0,y0,xn,yn,h1,hn,up0,lowN,b,bn,y0prime,ynprime

!	External subroutine variables
	integer :: INFO
	
	h1 = knots(2,1) - knots(1,1)
	hn = knots(n,1) - knots(n-1,1)

	b = 6.d0/h1*((knots(2,2)-knots(1,2))/h1 - y0prime)
	bn = 6.d0/hn*(ynprime - (knots(n,2)-knots(n-1,2))/hn)

	do i =1,n ! generate arrays of dimension n and dimension n-1 
	 if(i.eq.1) then
	  rhsvector(i) = b
	  upperdiag(i) = 1.d0
	 else if (i.eq.n) then
	  rhsvector(i) = bn
!	  lowerdiag(i) = 1.d0
	 else if ((i.ge.2).AND.(i.le.n-1)) then
	  rhsvector(i) = (6.d0/(knots(i+1,1)-knots(i-1,1)))*((knots(i+1,2)-knots(i,2))/(knots(i+1,1)-knots(i,1)) &
	  						    -(knots(i,2)-knots(i-1,2))/(knots(i,1)-knots(i-1,1)))
	  upperdiag(i) = (knots(i+1,1)-knots(i,1))/(knots(i+1,1)-knots(i-1,1))
	  lowerdiag(i) = 1.d0 - upperdiag(i)							    
	 end if 
	end do
	
	do i=1,n
	maindiag(i)=2.d0
	end do
		
	call dgtsv(n,1,lowerdiag,maindiag,upperdiag,rhsvector,n,INFO)

	if (info.eq.0) then 
	moments = rhsvector
	else
	print*,'--- Error in the calculation of moments'
	end if
	
	end subroutine cubicsplines
!========================================================
	
	subroutine splineinterpol(knots,moments,n,x,y)

!	Evaluate the spline for the interval to which x belongs and returns the interpolated value.	
!	See J. Stoer 'Introduction to Numerical Analysis', section 2.4.	
	
	implicit none
	
!	local variables	
	integer:: j,n,i
	double precision::x,xj,y,a,b,c,d
	double precision::moments(n),knots(n,2)
	
!	external subroutine variables
	double precision::array(n)
	
	do i= 1,n
	array(i)= knots(i,1)
	end do
	
	call locate(Array,n,x,j)
		
	a = knots(j,2)
	b = (knots(j+1,2)-knots(j,2))/(knots(j+1,1)-knots(j,1))-(2.d0*moments(j)+moments(j+1))*(knots(j+1,1)-knots(j,1))/6.d0
	c = moments(j)/2.d0
	d = (moments(j+1)-moments(j))/(6.d0*(knots(j+1,1)-knots(j,1)))
	
	xj = knots(j,1)
	y = a + b*(x-xj)+c*(x-xj)**2+d*(x-xj)**3
		
	end subroutine splineinterpol

!========================================================
		
	subroutine locate(grid,n,x,j)

	!Array grid(n) must be in monotonically ascending order. 
	!Subroutine gives the lower index of the smallest interval that includes x.
	! The search is done by bisection.
	
	implicit none
	integer :: n,j,jm, jup, jlow,i
	double precision:: x
	double precision :: grid(n)
	
	jup = n
	jlow = 1

	if(x.gt.(grid(n)).or.(x.lt.grid(1))) then
	WRITE(*,'(A)') 'locate - error - value to interpolate outside range of data'
	WRITE(*,'(A,F10.5)') 'x = ', real(x)
	STOP
	end if
	
	If (DABS(grid(1)-x).lt.1d-9) then
	 j = 1
	 
	else if (DABS(grid(n)-x).lt.1d-9) then
	 j = n-1
	 
	else
	
	 do 
	 
	 if(jup-jlow.gt.1) then
	 jm = (jup + jlow)/2
	
	   if((x.gt.grid(jlow)).and.(x.lt.grid(jm))) then
	     jup = jm
	   end if
	   if((x.gt.grid(jm)).and.(x.lt.grid(jup))) then  
	     jlow = jm
	   end if
	   if(DABS(grid(jm)-x).lt.1d-9) then
	   j=jm
	   exit
	   end if
	
	 else 
	 j = jlow
	 exit
	
	 end if
	 
	 end do
	 
	end if
	
	end subroutine locate
	


