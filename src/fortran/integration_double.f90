subroutine cubint ( ftab, xtab, ntab, ia, ib, result, error )
!
!***********************************************************************
!
!! CUBINT approximates an integral using cubic interpolation of data.
!
!
!  Discussion:
!
!    The integral to be approximated is
! 
!      INTEGRAL (XTAB(IB) to XTAB(IA)) F(X) DX
!
!    The routine estimates the error in integration.
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!    P E Gill and G F Miller
!    An algorithm for the integration of unequally spaced data,
!    Comput J, Number 15, 1972, pages 80-83.
!
!  Modified:
!
!    30 October 2000
!    19 April 2013 'REAL converted to REAL*8'
!  Parameters:
!
!    Input, real*8 FTAB(NTAB), contains the tabulated function
!    values, FTAB(I) = F(XTAB(I)).
!
!    Input, real*8 XTAB(NTAB), contains the points at which the
!    function was tabulated.  XTAB should contain distinct
!    values, given in ascending order.
!
!    Input, integer NTAB, the number of tabulated points.
!    NTAB must be at least 4.
!
!    Input, integer IA, the entry of XTAB at which integration
!    is to begin.  IA must be no less than 1 and no greater
!    than NTAB.
!
!    Input, integer IB, the entry of XTAB at which integration
!    is to end.  IB must be no less than 1 and no greater than
!    NTAB.
!
!    Output, real*8 RESULT, the approximate value of the
!    integral from XTAB(IA) to XTAB(IB) of the function.
!
!    Output, real*8 ERROR, an estimate of the error in
!    integration.
!
  implicit none
!
  integer ntab
!
  real*8 c
  real*8 d1
  real*8 d2
  real*8 d3
  real*8 error
  real*8 ftab(ntab)
  real*8 h1
  real*8 h2
  real*8 h3
  real*8 h4
  integer i
  integer ia
  integer ib
  integer ind
  integer it
  integer j
  integer k
  real*8 r1
  real*8 r2
  real*8 r3
  real*8 r4
  real*8 result
  real*8 s
  real*8 term
  real*8 xtab(ntab)
!
  result = 0.0D+00
  error = 0.0D+00
   
  if ( ia == ib ) then
    return
  end if
 
  if ( ntab < 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i6)' ) '  NTAB must be at least 4, but input NTAB = ',ntab
    stop
  end if
 
  if ( ia < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i6)' ) '  IA must be at least 1, but input IA = ',ia
    stop
  end if
 
  if ( ia > ntab ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i6)' ) '  IA must be <= NTAB, but input IA=',ia
    stop
  end if
 
  if ( ib < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i6)' ) '  IB must be at least 1, but input IB = ',ib
    stop
  end if
 
  if ( ib > ntab ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i6)' ) '  IB must be <= NTAB, but input IB=',ib
    stop
  end if
!
!  Temporarily switch IA and IB, and store minus sign in IND
!  so that, while integration is carried out from low X's
!  to high ones, the sense of the integral is preserved.
!
  if ( ia > ib ) then
    ind = -1
    it = ib
    ib = ia
    ia = it
  else
    ind = 1
  end if
 
  s = 0.0D+00
  c = 0.0D+00
  r4 = 0.0D+00
  j = ntab-2
  if ( ia < ntab-1 .or. ntab == 4 ) then
    j=max(3,ia)
  end if

  k = 4
  if ( ib > 2 .or. ntab == 4 ) then
    k=min(ntab,ib+2)-1
  end if
 
  do i = j, k
 
    if ( i <= j ) then
 
      h2 = xtab(j-1)-xtab(j-2)
      d3 = (ftab(j-1)-ftab(j-2)) / h2
      h3 = xtab(j)-xtab(j-1)
      d1 = (ftab(j)-ftab(j-1)) / h3
      h1 = h2+h3
      d2 = (d1-d3)/h1
      h4 = xtab(j+1)-xtab(j)
      r1 = (ftab(j+1)-ftab(j)) / h4
      r2 = (r1-d1) / (h4+h3)
      h1 = h1+h4
      r3 = (r2-d2) / h1
 
      if ( ia <= 1 ) then
        result = h2 * (ftab(1)+h2*(0.5D0*d3-h2*(d2/6.0D0-(h2+h3+h3)*r3/12.D0)))
        s = -h2**3 * (h2*(3.0D0*h2+5.0D0*h4)+10.0D0*h3*h1)/60.0d0
      end if
 
    else
 
      h4 = xtab(i+1)-xtab(i)
      r1 = (ftab(i+1)-ftab(i))/h4
      r4 = h4+h3
      r2 = (r1-d1)/r4
      r4 = r4+h2
      r3 = (r2-d2)/r4
      r4 = (r3-d3)/(r4+h1)
 
    end if
 
    if ( i > ia .and. i <= ib ) then
 
      term = h3*((ftab(i)+ftab(i-1))*0.5d0-h3*h3*(d2+r2+(h2-h4)*r3) / 12.0d0 )
      result = result+term
      c = h3**3*(2.0d+00 *h3*h3+5.d0*(h3*(h4+h2) + 2.0d0 * h2 * h4 ) ) / 120.0d+00
      error = error+(c+s)*r4
 
      if ( i /= j ) then
        s = c
      else
        s = s+c+c
      end if
 
    else
 
      error = error+r4*s
 
    end if
 
    if ( i >= k ) then
 
      if ( ib >= ntab ) then
        term = h4*(ftab(ntab) - h4*(0.5d0*r1+h4*(r2/6.0d0 +(h3+h3+h4)*r3/12.d0)))
        result = result + term
        error = error - h4**3 * r4 * &
          ( h4 * ( 3.0d0 * h4 + 5.0d0 * h2 ) &
          + 10.0d0 * h3 * ( h2 + h3 + h4 ) ) / 60.0d+00
      end if
 
      if ( ib >= ntab-1 ) error=error+s*r4
    else
      h1 = h2
      h2 = h3
      h3 = h4
      d1 = r1
      d2 = r2
      d3 = r3
    end if
  end do
!
!  Restore original values of IA and IB, reverse signs
!  of RESULT and ERROR, to account for integration
!  that proceeded from high X to low X.
!
  if ( ind /= 1 ) then
    it = ib
    ib = ia
    ia = it
    result = -result
    error = -error
  end if
 
  return
end
