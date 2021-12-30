! ADKL369j.  UNIFIED APPROACH FOR EXACT CALCULATION OF ANGULAR MOMENTUM   
!1   COUPLING AND RECOUPLING COEFFICIENTS.  L. WEI.                    
!REF. IN COMP. PHYS. COMMUN. 120 (1999) 222                           
!The program consists of the following files:!

!drive.f		driver program
!369j.f	

!3j.in		3 files of test data
!6j.in
!9j.in	

       real*8 function threejsymbol(j1,j2,j3,m1,m2,m3)
c
c  This program calculates the exact 3j symbols using two kinds of number
c  representation: prime number representation for prefactor, and 32768-
c  base number representation for summation. The results are tabulated 
c  when the logical variable "pprint" is set to be "true".
c
      parameter (length = 100, ndim = 301, base = 32768.0)
      real j1,j2,j3,m1,m2,m3 
      logical pprint
      integer n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12
      integer i1,i2, sign, numpwr
      integer iac(length), prod(length), sum(length)
      integer init1(length), init2(length), init3(length)
      integer prime(ndim), power(ndim), ppower(ndim)
      integer power1(ndim), power2(ndim)
      integer ppower1(ndim), ppower2(ndim)
      real*8 pprod, ssum, factor, cuthi, cutlo
      common /mc/ iac, sign
      common /pr/ prime
	
      pprint = .false.
	
c  Set the limits for the maximum and the minimum in the real*8 numbers:
c 
      cuthi = 1.0d300
      cutlo = 1.0d-300
c
c  Check validity of input ?
c
      if ( ((m1+m2) .ne. -m3)  .or.
     &     (abs(m1) .gt. j1)  .or.
     &     (abs(m2) .gt. j2)  .or.
     &     (abs(j1-j2) .gt. j3)  .or.
     &     ((j1+j2) .lt. j3) ) then
         
	threejsymbol =0.d0
	 
      else
c
c  Define some integer constants:
c
      n1 = j1+j2-j3
      n2 = j1-j2+j3
      n3 = -j1+j2+j3
      n4 = j1-m1
      n5 = j1+m1
      n6 = j2+m2
      n7 = j3+m3
      n8 = -j2+j3+m1
      n9 = -j1+j3-m2
      n10 = 2*j1
      n11 = 2*j2
      n12 = 2*j3
c
c  DO THE CALCULATION IN THE 32768-BASE NUMBER REPRESENTATION
c
c  Determine the limits of summation:
c
      i1 = max(0,-n8,-n9)
      i2 = min(n1,n4,n6)
c
c  Initialization:
c
      call binom(n1,i1,init1)
      call binom(n2,n4-i1,init2)
      call binom(n3,n6-i1,init3)
      call load(init1)
      call mul(init2)
      call mul(init3)
      if (mod(i1,2) .ne. 0) sign = 0
      call store(prod)
      call store(sum)
c
c  Summation:
c
      do 10 k = i1+1, i2
        call load(prod)
        call muls(n1+1-k)
        call divs(k,l1)
        call muls(n4+1-k)
        call divs(n8+k,l2)
        call muls(n6+1-k)
        call divs(n9+k,l3)
        sign = 1 - sign
        call store(prod)
        call add_sub(sum)
        call store(sum)
   10 continue
c
c  Check whether the sum is zero:
c
      if (iac(1) .eq. 2 .and. iac(2) .eq. 0) then
        threejsymbol = 0.d0
        return
      end if
c
c  DO THE CALCULATION IN THE PRIME NUMBER REPRESENTATION
c      
c  Find the first ndim-1 prime numbers:
c
      call find_prime()
c
c  Decompose the sum into the prime number representation as 
c  completely as possible within the given dimension:
c
      call decomp(power, k)
c
c  Calculate the square of prefactor and its multiplication with the square 
c  of the part in the sum which has been converted to the prime factors:
c
      call binom1(n10, n2, power1, k1)
      call binom1(n10, n5, power2, k2)
      call mul_div(-1, power1, k1, power2, k2, ppower1, kk1)
      call binom1(n11, n1, power1, k1)
      call binom1(n11, n6, power2, k2)
      call mul_div(1, ppower1, kk1, power1, k1, ppower2, kk2)
      call mul_div(-1, ppower2, kk2, power2, k2, ppower1, kk1)
      call binom1(n12, n3, power1, k1)
      call binom1(n12, n7, power2, k2)
      call mul_div(1, ppower1, kk1, power1, k1, ppower2, kk2)
      call mul_div(-1, ppower2, kk2, power2, k2, ppower1, kk1)
      call delta2(j1, j2, j3, ppower2, kk2)
      call mul_div(-1, ppower1, kk1, ppower2, kk2, ppower, kk)
      call mul_div(1, ppower, kk, power, k, ppower1, kk1)
      call mul_div(1, ppower1, kk1, power, k, ppower, kk)
c
      if (mod(int(j1-j2-m3),2) .ne. 0) sign = 1 - sign
c
c  Tabulation of 3j symbols where the part of prime numbers has been squared:
c
      if (pprint) then
        write(6,*) "TABULATION: ", " sign ", sign
        write(6,*) "32768-base number part whose number of digits is",
     &            iac(1)-1," ,"
        write(6,*) (iac(i), i = iac(1), 2, -1)
        write(6,*) "prime part whose number is",kk-1," ,"
        write(6,*) (ppower(i), i = 2, kk)
      end if
c
c  Exact decimal value of 3j symbol:
c
      numpwr = 0      
      ssum = 0.0
      do 20 i = 2, iac(1)
        factor = iac(i)
        j = 1
        do while (j .le. i-2-numpwr)
          factor = factor*base
          if (factor .le. cuthi) then
            j = j + 1
          else
            factor = factor/base
            do 25 jj = 1, (i-2-numpwr-j+1)
              ssum = ssum/base
   25       continue    
            numpwr = numpwr + (i-2-numpwr-j+1)
            j = i-2-numpwr+1
          end if
        end do
        ssum = ssum + factor
        if (ssum .gt. cuthi) then
          ssum = ssum/base
          numpwr = numpwr + 1
        end if
   20 continue    
c
      pprod = 1.0
      do 30 i = 2, kk
        ieven = abs(ppower(i))/2
        irem = mod(abs(ppower(i)), 2)
        if (ppower(i) .gt. 0) then
          do 33 j = 1, ieven
            pprod = pprod*dble(prime(i))
            if (pprod .gt. cuthi) then
              pprod = pprod/base
              numpwr = numpwr + 1
            end if
   33     continue
          if (irem .eq. 1) then
            pprod = pprod*dsqrt(dble(prime(i)))
            if (pprod .gt. cuthi) then
              pprod = pprod/base
              numpwr = numpwr + 1
            end if
          end if
        else if (ppower(i) .lt. 0) then
          do 36 j = 1, ieven
            pprod = pprod/dble(prime(i))
            if (pprod .lt. cutlo) then
              pprod = pprod*base
              numpwr = numpwr - 1
            end if
   36     continue 
          if (irem .eq. 1) then
            pprod = pprod/dsqrt(dble(prime(i)))
            if (pprod .lt. cutlo) then
              pprod = pprod*base
              numpwr = numpwr - 1
            end if
          end if
        end if
   30 continue    
c
      threejsymbol = pprod*ssum
      if (numpwr .gt. 0) then
        do 40 i = 1, numpwr
          threejsymbol = threejsymbol*base
   40   continue
      else if (numpwr .lt. 0) then    
        do 45 i = 1, -numpwr
          threejsymbol = threejsymbol/base
   45   continue
      end if
      if (sign .eq. 0)  threejsymbol = - threejsymbol
      if (threejsymbol.eq. 0.0) write(6,*) "The 3j symbol underflows !"
c
      return
      
      end if
      end
c
c---------------------------------------------------------------------
c
      real*8 function sixj(a,b,c,d,e,f)
c
c  This function calculates the exact 6j symbols using two kinds of number
c  representation: prime number representation for prefactor, and 32768-
c  base number representation for summation. The results are tabulated 
c  when the logical variable "pprint" is set to be "true".
c
      parameter (length = 100, ndim = 301, base = 32768.0)    
      real a,b,c,d,e,f
      logical pprint 
      integer abc,aef,bdf,cde,abde,acdf,bcef
      integer i1,i2, k, sign, numpwr
      integer iac(length), prod(length), sum(length)
      integer init1(length),init2(length),init3(length),init4(length)
      integer prime(ndim),power(ndim)
      integer power1(ndim), power2(ndim), power3(ndim)
      integer ppower1(ndim), ppower2(ndim)
      real*8 ffactor, pprod, ssum, cuthi, cutlo
      common /mc/ iac, sign
      common /pr/ prime
c
      pprint = .false.
c  Set the limits of the maximum and the minimum in the real*8 numbers:
c
      cuthi = 1.0d300
      cutlo = 1.0d-300
c
c  Check the validity of input ?
c
      if (a .lt. abs(b-c) .or. a .gt. b+c
     &  .or. a .lt. abs(e-f) .or. a .gt. e+f
     &  .or. b .lt. abs(d-f) .or. b .gt. d+f
     &  .or. c .lt. abs(d-e) .or. c .gt. d+e) then
            
      sixj = 0.d0
      
      else
      
      
	
c  DO THE SUMMATION IN THE 32768-Base NUMBER REPRESENTATION
c  
c  Define some integer constants:
c
      abc = a+b+c
      aef = a+e+f
      bdf = b+d+f
      cde = c+d+e
      abde = a+b+d+e+1
      acdf = a+c+d+f+1
      bcef = b+c+e+f+1
c
c  Determine the limits of summation:
c
      i1 = max(abc,aef,bdf,cde)
      i2 = min(abde-1,acdf-1,bcef-1)
c
c  Initialization:
c
      call binom(i1+1,abc+1,init1)
      call binom(int(a+b-c),i1-cde,init2)
      call binom(int(a-b+c),i1-bdf,init3)
      call binom(int(-a+b+c),i1-aef,init4)
      call load(init1)
      call mul(init2)
      call mul(init3)
      call mul(init4)
      if (mod(i1,2) .ne. 0)  sign = 0
      call store(prod)
      call store(sum)
c
c  Summation:
c
      do 10 k = i1+1, i2
        call load(prod)
        call muls(k+1)
        call divs(k-abc,l1)
        call muls(abde-k)
        call divs(k-cde,l2)
        call muls(acdf-k)
        call divs(k-bdf,l3)
        call muls(bcef-k)
        call divs(k-aef,l4)
        sign = 1 - sign
        call store(prod)
        call add_sub(sum)
        call store(sum)
   10 continue
c
c  Check whether the sum is zero:
c
      if (iac(1) .eq. 2 .and. iac(2) .eq. 0) then
        sixj = 0.0
        return
      end if
	
c  DO THE CALCULATION IN THE PRIME NUMBER REPRESENTATION
c      
c  Find the first ndim-1 prime numbers:
c
      call find_prime()
c
c  Decompose the sum into the prime number representation as 
c   completely as possible within the given dimension:
c
      call decomp(power, k)
       
c
c  Calculate the square of prefactor and its multiplication with the square 
c  of the part in the sum which has been converted to the prime numbers:
c
      call delta2(a,e,f,power1,k1)
      call delta2(b,d,f,power2,k2)
      call mul_div(1,power1,k1,power2,k2,power3,k3)
      call mul_div(-1,power,k,power3,k3,ppower1,kk1)
      call delta2(a,b,c,power1,k1)
      call delta2(c,d,e,power2,k2)
      call mul_div(-1,power1,k1,power2,k2,power3,k3)
      call mul_div(1,power,k,power3,k3,ppower2,kk2)
      call mul_div(1,ppower1,kk1,ppower2,kk2,power,k)
	
c
c  Tabulation of 6j symbols where the part of prime numbers has been squared:
c
      if (pprint) then
        write(6,*) "TABULATION: ",  " sign ", sign
        write(6,*) "32768-base number part whose number of digits is",
     &            iac(1)-1," ,"
        write(6,*) (iac(i), i = iac(1), 2, -1)
        write(6,*) "prime part whose number is",k-1," ,"
        write(6,*) (power(i), i = 2, k)
      end if
c
c  Exact decimal value of 6j symbol:
c
      numpwr = 0
      ssum = 0.0
      do 20 i = 2, iac(1)
        ffactor = iac(i)
        j = 1
        do while (j .le. i-2-numpwr)
          ffactor = ffactor*base
          if (ffactor .le. cuthi) then
            j = j + 1
          else
            ffactor = ffactor/base
            do 25 jj = 1, i-2-numpwr-j+1
              ssum = ssum/base
   25       continue
            numpwr = numpwr + (i-2-numpwr-j+1)
            j = i-2-numpwr+1    
          end if
        end do   
        ssum = ssum + ffactor
        if (ssum .gt. cuthi) then
          ssum = ssum/base
          numpwr = numpwr + 1
        end if
   20 continue    
c
      pprod = 1.0
      do 30 i = 2, k
        ieven = abs(power(i))/2
        irem = mod(abs(power(i)), 2)
        if (power(i) .gt. 0) then
          do 33 j = 1, ieven
            pprod = pprod*dble(prime(i))
            if (pprod .gt. cuthi) then
              pprod = pprod/base
              numpwr = numpwr + 1
            end if
   33     continue
          if (irem .eq. 1) then
            pprod = pprod*dsqrt(dble(prime(i)))
            if (pprod .gt. cuthi) then
              pprod = pprod/base
              numpwr = numpwr + 1
            end if
          end if                    
        else if (power(i) .lt. 0) then
          do 36 j = 1, ieven
            pprod = pprod/dble(prime(i))
            if (pprod .lt. cutlo) then
              pprod = pprod*base
              numpwr = numpwr - 1
            end if
   36     continue        
          if (irem .eq. 1) then
            pprod = pprod/dsqrt(dble(prime(i)))
            if (pprod .lt. cutlo) then
              pprod = pprod*base
              numpwr = numpwr - 1
            end if
          end if
        end if
   30 continue    
c
      sixj = pprod*ssum
      if (numpwr .gt. 0) then
        do 40 i = 1, numpwr
          sixj = sixj*base
   40   continue    
      else if (numpwr .lt. 0) then
        do 45 i = 1, -numpwr
          sixj = sixj/base
   45   continue
      end if    
      if (sign .eq. 0) sixj = - sixj
      if (sixj .eq. 0.0) write(6,*) "The 6j symbol underflows !"
c
      return
      
      end if
      end
c
c------------------------------------------------------------------------------
c      
      real*8 function ninej(a,b,c,d,e,f,g,h,j)
c
c  This function calculates the exact 9j symbols using two kinds of number 
c  representation: prime number representation for prefactor, and 32768-
c  base number representation for summation. The results are tabulated 
c  when the logical variable "pprint" is set to be "true".
c
      parameter (length = 100, ndim = 301, base = 32768.0)
      real a,b,c,d,e,f,g,h,j,i1,i2,kj
      logical pprint
      integer sign, iac(length), sum(length), numpwr
      integer init1(length), init2(length), init3(length)
      integer prime(ndim), power(ndim), ppower(ndim)
      integer power1(ndim), power2(ndim), power3(ndim)
      integer ppower1(ndim), ppower2(ndim), ppower3(ndim)
      real*8 ffactor, prod, ssum, cuthi, cutlo
      common /mc/ iac, sign
      common /pr/ prime
            
      pprint = .false.
c  Set the limits of the maximum and the minimum in the real*8 numbers:
c
      cuthi = 1.0d300
      cutlo = 1.0d-300
c
c  Check validity of input ?
c
      if ( a .gt. (b+c) .or. a .lt. abs(b-c)
     &  .or. a .gt. (d+g) .or. a .lt. abs(d-g)
     &  .or. d .gt. (e+f) .or. d .lt. abs(e-f)
     &  .or. g .gt. (h+j) .or. g .lt. abs(h-j)
     &  .or. b .gt. (e+h) .or. b .lt. abs(e-h)
     &  .or. c .gt. (f+j) .or. c .lt. abs(f-j) ) then
        
	ninej = 0.d0
	
      else
c
c  DO THE SUMMATION IN THE 32768-BASE NUMBER REPRESENTATION 
c
c  Determine the limits of summation:
c
      i1 = amax1(abs(h-d),abs(b-f),abs(a-j))
      i2 = amin1(h+d,b+f,a+j)
c
c  Initialization of sum as zero and of kj:
c
      sum(1) = 2
      sum(2) = 0
      sum(3) = 1
      kj = i1
c
      do while (kj .ge. i1 .and. kj .le. i2)    
        call factor(a,b,c,f,j,kj, init1)
        call factor(f,d,e,h,b,kj, init2)
        call factor(h,j,g,a,d,kj, init3)
        call load(init1)
        call mul(init2)
        call mul(init3)
        call muls(int(2*kj+1))
        if (mod(int(2*kj), 2) .ne. 0)  sign = 1 - sign
        call add_sub(sum)
        call store(sum)
        kj = kj + 1.0
      end do
c
c  Check whether the sum is zero:
c
      if (iac(1) .eq. 2 .and. iac(2) .eq. 0) then
        ninej = 0.0
        return
      end if
c
c  DO THE CALCULATION IN THE PRIME NUMBER REPRESENTATION
c      
c  Find the first ndim-1 prime numbers:
c
      call find_prime()
c
c  Decompose the sum into the prime number representation as 
c   completely as possible within the given dimension:
c
      call decomp(power, k)
c
c  Calculate the square of prefactor and its multiplication with the square 
c  of the part in the sum which has been converted to the prime factors:
c
      call delta2(a,b,c,power1,k1)
      call delta2(d,e,f,power2,k2)
      call delta2(g,h,j,power3,k3)
      call mul_div(1,power1,k1,power2,k2,ppower1,kk1)
      call mul_div(1,ppower1,kk1,power3,k3,ppower2,kk2)
      call mul_div(-1,power,k,ppower2,kk2,ppower,kk)
      call delta2(a,d,g,power1,k1)
      call delta2(b,e,h,power2,k2)
      call delta2(c,f,j,power3,k3)
      call mul_div(1,power1,k1,power2,k2,ppower1,kk1)
      call mul_div(1,ppower1,kk1,power3,k3,ppower2,kk2)
      call mul_div(-1,power,k,ppower2,kk2,ppower3,kk3)
      call mul_div(1,ppower,kk,ppower3,kk3,power,k)
c
c  Tabulation of 9j symbols where the prime numbers have been squared:
c
      if (pprint) then
        write(6,*) "TABULATION: ", " sign ", sign
        write(6,*) "32768-base number part whose number of digits is",
     &            iac(1)-1," ,"
        write(6,*) (iac(i), i = iac(1), 2, -1)
        write(6,*) "prime part whose number is",k-1," ,"
        write(6,*) (power(i), i = 2, k)
      end if
c
c  Exact decimal value of 9j symbol:
c
      numpwr = 0
      ssum = 0.0
      do 10 i = 2, iac(1)
        ffactor = iac(i)
        i1 = 1
        do while (i1 .le. i-2-numpwr)
          ffactor = ffactor*base
          if (ffactor .le. cuthi) then
            i1 = i1 + 1
          else
            ffactor = ffactor/base
            do 15 jj = 1, (i-2-numpwr-i1+1)
              ssum = ssum/base
   15       continue
            numpwr = numpwr + (i-2-numpwr-i1+1)
            i1 = i-2-numpwr+1    
          end if
        end do    
        ssum = ssum + ffactor
        if (ssum .gt. cuthi) then
          ssum = ssum/base
          numpwr = numpwr + 1
        end if
   10 continue    
c
      prod = 1.0
      do 20 i = 2, k
        ieven = abs(power(i))/2
        irem = mod(abs(power(i)), 2)
        if (power(i) .gt. 0) then
          do 23 i1 = 1, ieven
            prod = prod*dble(prime(i))
            if (prod .gt. cuthi) then
              prod = prod/base
              numpwr = numpwr + 1
            end if
   23     continue
          if (irem .eq. 1) then
            prod = prod*dsqrt(dble(prime(i)))
            if (prod .gt. cuthi) then
              prod = prod/base
              numpwr = numpwr + 1
            end if
          end if
        else if (power(i) .lt. 0) then    
          do 26 i1 = 1, ieven
            prod = prod/dble(prime(i))
            if (prod .lt. cutlo) then
              prod = prod*base
              numpwr = numpwr - 1
            end if
   26     continue
          if (irem .eq. 1) then
            prod = prod/dsqrt(dble(prime(i)))
            if (prod .lt. cutlo) then
              prod = prod*base
              numpwr = numpwr - 1
            end if
          end if
        end if
   20 continue    
c
      ninej = prod*ssum
      if (numpwr .gt. 0) then
        do 30 i = 1, numpwr
          ninej = ninej*base
   30   continue
      else if (numpwr .lt. 0) then
        do 35 i = 1, -numpwr
          ninej = ninej/base
   35   continue   
      end if 
      if (sign .eq. 0) ninej = - ninej
      if (ninej .eq. 0.0) write(6,*) "The 9j symbol underflows !"
c
      return
      
      end if
      end
c
c------------------------------------------------------------------------------
c
      subroutine factor(a,b,c,d,e,f, sum)
c
c  This subroutine calculates the factors in summation steps of 9j symbols 
c  using 32768-base number representation.
c
      parameter (length = 100)    
      real a,b,c,d,e,f
      integer abc,aef,bdf,cde,abde,acdf,bcef
      integer i1,i2, sign
      integer iac(length), prod(length), sum(length)
      integer init1(length),init2(length),init3(length),init4(length)
      common /mc/ iac, sign
c
c  Define some integer constants:
c
      abc = a+b+c
      aef = a+e+f
      bdf = b+d+f
      cde = c+d+e
      abde = a+b+d+e+1
      acdf = a+c+d+f+1
      bcef = b+c+e+f+1
c
c  Determine the limits of summation:
c
      i1 = max(abc,aef,bdf,cde)
      i2 = min(abde-1,acdf-1,bcef-1)
c
c  Initialization:
c
      sign = 1
      call binom(i1+1,aef+1,init1)
      call binom(int(a+e-f),i1-bdf,init2)
      call binom(int(a-e+f),i1-cde,init3)
      call binom(int(-a+e+f),i1-abc,init4)
      call load(init1)
      call mul(init2)
      call mul(init3)
      call mul(init4)
      if (mod(i1,2) .ne. 0) sign = 0
      call store(prod)
      call store(sum)
c
c  Summation:
c
      do 10 k = i1+1, i2
        call load(prod)
        call muls(k+1)
        call divs(k-aef,l1)
        call muls(abde-k)
        call divs(k-bdf,l2)
        call muls(acdf-k)
        call divs(k-cde,l3)
        call muls(bcef-k)
        call divs(k-abc,l4)
        sign = 1 - sign
        call store(prod)
        call add_sub(sum)
        call store(sum)
   10 continue
c
      return
      end
c
c------------------------------------------------------------------------------
c
      subroutine load(ix)
c
c  This subroutine loads the array 'ix' to the array 'iac' which are in 
c  common storage. We use the first element of the array 'iac' to represent 
c  its length. The sign of 'iac' is positive when 'sign' is equal to 1, and
c  negative when it is zero.
c
      parameter (length = 100)
      integer sign, iac(length), ix(length)
      common /mc/ iac, sign
c
      do 10 i = 1, ix(1)
        iac(i) = ix(i)
   10 continue
      sign = ix(ix(1)+1)
      return
      end     
c
c---------------------------------------------------------------------------
c
      subroutine store(ix)
c
c  This subroutine stores the array 'iac' to the array 'ix'. Note that the
c  element ix(ix(1)+1) represents the sign of 'ix'.
c
      parameter (length = 100)
      integer sign, iac(length), ix(length)
      common /mc/ iac, sign
c
      do 10 i = 1, iac(1)
        ix(i) = iac(i)
   10 continue
      ix(iac(1)+1) = sign
      return
      end     
c
c---------------------------------------------------------------------------
c
      subroutine add_sub(ix)
c
c  This subroutine does the algebraic summation of arrays 'iac' and 'ix' and 
c  then load the resultant as 'iac'.
c
      parameter (length = 100)
      integer sign, iac(length), ix(length), a, c
      common /mc/ iac, sign
c
      k = min(iac(1), ix(1))
c
c  ADD ix WITH iac:
c
      if (sign .eq. ix(ix(1)+1)) then
        c = 0
        do 10 i = 2, k
          a = iac(i) + ix(i) + c
          if (a .gt. 32767) then
            iac(i) = mod(a, 32768)
            c = 1
          else
            iac(i) = a
            c = 0
          end if
   10   continue
        if (iac(1) .gt. ix(1)) then
          do 20 i = k+1, iac(1)
            a = iac(i) + c
            if (a .eq. 32768) then
              iac(i) = 0
              c = 1
            else
              iac(i) = a
              c = 0
            end if
   20     continue
        else if (ix(1) .gt. iac(1)) then
          do 30 i = k+1, ix(1)
            a = ix(i) + c
            if (a .eq. 32768) then
              iac(i) = 0
              c = 1
            else
              iac(i) = a
              c = 0
            end if
   30     continue
        end if
        iac(1) = i - 1
        if (c .eq. 1) then
          if (i .gt. length) then
            write(*,*) "increase the dimension of array: iac !"
            stop
          end if
          iac(i) = 1
          iac(1) = i
        end if
      end if
c
c  SUBTRACT ix FROM iac:
c
      if (sign .ne. ix(ix(1)+1)) then
        if (iac(1) .gt. ix(1)) then
          do 40 i = 2, k
            if (iac(i) .ge. ix(i)) then
              iac(i) = iac(i) - ix(i)
            else
              iac(i) = (iac(i) - ix(i)) + 32768
              j = i + 1
              do while (iac(j) .eq. 0)
                iac(j) = 32767
                j = j + 1
              end do
              iac(j) = iac(j) - 1
            end if
   40     continue
          do 50 i = k+1, iac(1)
            iac(i) = iac(i)
   50     continue
        else if (ix(1) .gt. iac(1)) then
          iac(iac(1)+1) = 0
          do 60 i = 2, k
            if (ix(i) .ge. iac(i)) then
              iac(i) = ix(i) - iac(i)
            else
              iac(i) = (ix(i) - iac(i)) + 32768
              j = i + 1
              do while (iac(j) .eq. 32767)
                iac(j) = 0
                j = j + 1
              end do
              iac(j) = iac(j) + 1
            end if
   60     continue
          if (iac(iac(1)+1) .eq. 0) then
            do 70 i = k+1, ix(1)
              iac(i) = ix(i)
   70       continue
          else 
            j = iac(1)+1
            do while (ix(j) .eq. 0)
              iac(j) = 32767
              j = j + 1
            end do
            iac(j) = ix(j) - 1
            do 80 i = j+1, ix(1)
              iac(i) = ix(i)
   80       continue
          end if
          sign = 1 - sign
        else
          do while (iac(k) .eq. ix(k) .and. k .gt. 1)
            iac(k) = 0
            k = k - 1
          end do
          if (k .eq. 1) then
            iac(1) = 2
            iac(2) = 0
            sign = 1
            return
          else
            iac(1) = k
            ix(k+1) = ix(ix(1)+1)
            ix(1) = k
            if (iac(k) .gt. ix(k)) then
              do 90 i = 2, k
                if (iac(i) .ge. ix(i)) then
                  iac(i) = iac(i) - ix(i)
                else
                  iac(i) = (iac(i) - ix(i)) + 32768
                  j = i + 1
                  do while (iac(j) .eq. 0)
                    iac(i) = 32767
                    j = j + 1
                  end do
                  iac(j) = iac(j) - 1
                end if
   90         continue    
            else if (ix(k) .gt. iac(k)) then
              do 100 i = 2, k
                if (ix(i) .ge. iac(i)) then
                  iac(i) = ix(i) - iac(i)
                else
                  iac(i) = (ix(i) - iac(i)) + 32768
                  j = i + 1
                  do while (iac(j) .eq. 32767)
                    iac(j) = 0
                    j = j + 1
                  end do
                  iac(j) = iac(j) + 1
                end if
  100         continue
              sign = 1 - sign
            end if
          end if
        end if
        m = max(iac(1), ix(1))
        if (iac(m) .ne. 0) then
          iac(1) = m
        else
          do while (iac(m) .eq. 0 .and. m .gt. 2)
            m = m - 1
          end do
          iac(1) = m
        end if
      end if      
c         
      return
      end     
c
c--------------------------------------------------------------------------
c
      subroutine muls(n)
c
c  This subroutine multiplies 'iac' by an integer n (<32768) and load it as 
c  'iac'.
c
      parameter (length = 100)
      integer sign, iac(length), a, c ,n
      common /mc/ iac, sign
c
      c = 0
      do 10 i = 2, iac(1)
        a = iac(i)*n + c
        if (a .gt. 32767) then
          iac(i) = mod(a, 32768)
          c = a/32768
        else
          iac(i) = a
          c = 0
        end if
   10 continue
      if (c. ne. 0) then
        if (i .gt. length) then
          write(*,*) "iac overflow !"
          stop
        end if
        iac(i) = c
        iac(1) = i
      end if
      return
      end     
c
c--------------------------------------------------------------------------
c
      subroutine mul(ix)
c  
c  This subroutine multiplies 'iac' by array 'ix' and then load the resultant 
c  as 'iac'.
c
      parameter (length = 100)
      integer sign, ssign
      integer iac(length), ix(length), ix1(length), ix2(length)
      common /mc/ iac, sign
c
c  Initialization of zero:
c
      ix2(1) = 2
      ix2(2) = 0
      ix2(3) = 1
c
      ssign = sign
      sign = 1
      call store(ix1)
      do 10 i = 2, ix(1)
        call load(ix1)
        n = ix(i)
        call muls(n)
        do 20 j = iac(1)+i-2, 2, -1
          if (j .ge. i) then
            iac(j) = iac(j-i+2)
          else
            iac(j) = 0
          end if
   20   continue
        iac(1) = iac(1) + i - 2
        call add_sub(ix2)
        call store(ix2)      
   10 continue
      call load(ix2)
      if (ssign .ne. ix(ix(1)+1)) then
        sign = 0
      else
        sign = 1
      end if
      return
      end     
c
c-------------------------------------------------------------------------
c
      subroutine divs(n, m)
c
c  This subroutine divides 'iac' by an integer n and load the resultant as 
c  'iac' and remainder as an integer m.
c
      parameter (length = 100)
      integer sign, n, m, iac(length), a, c
      common /mc/ iac, sign
c
      c = 0
      do 10 i = iac(1), 2, -1
        a = iac(i) + c*32768
        if (n .gt. a) then
          iac(i) = 0
          c = a
        else
          iac(i) = a/n
          c = mod(a,n)
        end if
   10 continue
      m = c
      do while (iac(iac(1)) .eq. 0 .and. iac(1) .gt. 2)
        iac(1) = iac(1) - 1
      end do
      return
      end     
c
c--------------------------------------------------------------------------
c
      subroutine binom(n,m,ix)
c
c  This subroutine calculates binomial in the 32768-base number representation,
c  and the resultant is array 'ix'.
c
      parameter (length = 100)
      integer sign, n, m, iac(length), ix(length)
      common /mc/ iac, sign
c
      if (n .lt. m) then
        write(*,*) "The binomial coefficient does not exist !"
        stop
      end if
c
c  Initialization of one:
c 
      ix(1) = 2
      ix(2) = 1
      ix(3) = 1
c
      k = min(m, n-m)
      call load(ix)
      do 10 i = 1, k
        call muls(n+1-i)
        call divs(i, l)
   10 continue
      call store(ix)
      return
      end     
c
c-------------------------------------------------------------------------
c
      subroutine find_prime()
c
c  This subroutine finds the first ndim-1 prime numbers and store
c  them into the array "prime" with common storage. Here we assume 
c  that the first prime number is '1' for convenience.
c
      parameter (ndim = 301)
      integer prime(ndim)
      common /pr/ prime
c
      i = 1
      prime(i) = 1
      k = 1
c
      do while (i. le. ndim)
        k = k + 1
        j = 0
        i1 = 2
        do while (i1 .le. k-1 .and. j .eq. 0)
          if ((float(k)/i1 - k/i1) .eq. 0.0) j = j + 1
          i1 = i1 + 1
        end do
        if (j .eq. 0) then
          i = i + 1
          prime(i) = k
        end if
      end do
      return
      end 
c
c------------------------------------------------------------------------------
c
      subroutine decomp(power, k)
c
c  This subroutine decomposes the non-zero 'iac' into the product of prime 
c  numbers as completely as possible within the dimension given and store
c  them into array "power" with dimension k.
c
      parameter (length = 100, ndim = 301)
      integer sign, iac(length), ix(length)
      integer power(ndim), prime(ndim)
      common /mc/ iac, sign
      common /pr/ prime
c
       k = 1
      power(k) = 1
      if (iac(1) .eq. 2 .and. iac(2) .eq. 1) return
c
      n = 1
      call store(ix)
      do while (ix(1) .ge. 2 .and. n .lt. ndim)
        n = n + 1
        m = 0
        j = 0
        do while (j .eq. 0) 
          call load(ix)
          call divs(prime(n), j)
          if (j .eq. 0) then 
            m = m + 1
            call store(ix)
          end if
        end do
        power(n) = m 
        if (m .gt. 0) k = n
        if (ix(1) .eq. 2 .and. ix(2) .eq. 1) then
          call load(ix)
          return
        end if
      end do
      call load(ix)
      return
      end     
c
c-----------------------------------------------------------------------------
c
      subroutine decomp1(num, power, k)
c
c  This subroutine decomposes the integer 'num' into the product of k 
c  prime numbers, each with power(i) equal to m. Also, suppose this 
c  decomposition is complete.
c
      parameter (ndim = 301)
      integer num, num1, power(ndim), prime(ndim)
      common /pr/ prime
c
      k = 1
      power(1) = 1
      num1 = num
      if (num1 .eq. 1) return
c
      i = 1
      do while (num1 .gt. 1 .and. i .lt. ndim)
        i = i + 1
        m = 0
        do while ( (dble(num1)/prime(i)-num1/prime(i)) .eq. 0.0 )
          m = m + 1
          num1 = num1/prime(i)
        end do
        power(i) = m
        if (m .gt. 0)  k = i
      end do
      if (num1 .gt. prime(ndim)) then
        write(*,*) "increase ndim !"
        stop
      end if
      return
      end
c
c------------------------------------------------------------------------------
c
      subroutine mul_div(p_s, power1, k1, power2, k2, power, k)
c
c  This subroutine does the multiplication (p_s=1) or division (p_s=-1) of two
c  arrays of prime numbers and the resultant is stored into the array 'power'
c  with dimension k.
c
      integer power1(k1), power2(k2), power(k), p_s
c
      power(1) = 1
      k = 1
c
      m = min(k1, k2)
      if (p_s .eq. 1) then
        do 10 i = 2, m
          power(i) = power1(i) + power2(i)
   10   continue
      else if (p_s .eq. -1) then
        do 20 i = 2, m
          power(i) = power1(i) - power2(i)
   20   continue
      else
        write(*,*) "Please input 1 or -1 for p_s !"
        stop 
      end if
c
      k = m
      do while (power(k) .eq. 0)
        k = k - 1
      end do
c
      if (k1 .gt. k2) then
        do 30 i = m+1, k1
          power(i) = power1(i)
   30   continue
        k = k1
      else if (k2 .gt. k1 .and. p_s .eq. 1) then
        do 40 i = m+1, k2
          power(i) = power2(i)
   40   continue
        k = k2
      else if (k2 .gt. k1 .and. p_s .eq. -1) then
        do 50 i = m+1, k2
          power(i) = -power2(i)
   50   continue
        k = k2
      end if
      return
      end
c
c------------------------------------------------------------------------------
c
      subroutine binom1(n, m, power, k)
c
c  This subroutine calculates the binomial in the prime number representation 
c  and the resultant is array 'power' with dimension k.
c
      parameter (ndim = 301)
      integer power(ndim), power1(ndim), power2(ndim), ppower(ndim)
c
      if (n .lt. m) then
         write (*,*) "The binomial coeficient does not exist !"
         stop
      end if
c
      k = 1
      power(k) = 1
c
      if (n .eq. m .or. m .eq. 0 .or. n .eq. 0) return     
c
      j = min(m, n-m)
      do 10 i = 1, j
        call decomp1(n+1-i, power1, k1)
        call decomp1(i, power2, k2)
        call mul_div(1, power, k, power1, k1, ppower, kk)
        call mul_div(-1, ppower, kk, power2, k2, power, k)
   10 continue     
      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine delta2(a, b, c, power, k)
c
c  This subroutine calculates the reverse of square of 'delta' defined in the 
c  angular momentum coupling theory in the prime number representation and
c  the resultant is stored into array 'power' with dimension k.
c
      parameter (ndim = 301)
      real a, b, c
      integer power(ndim),ppower(ndim)
      integer power1(ndim),power2(ndim),power3(ndim)
c
      if (a*b*c .eq. 0) then
         if (a .eq. 0) then
            call decomp1(int(2*b+1), power, k)
         else
            call decomp1(int(2*a+1), power, k)
         end if
      else
         call decomp1(int(a+b+c+1), power1, k1)
         call binom1(int(a+b+c), int(2*a), power2, k2)
         call binom1(int(2*a), int(a+b-c), power3, k3)
         call mul_div(1, power1, k1, power2, k2, ppower, kk)
         call mul_div(1, ppower, kk, power3, k3, power, k)
      end if
      return
      end
