Module Math

! Routines
!   - lu_system
!   - norme
!   - inverse (inverseV,inverseM)
!   - det (det33_r,det33_d)
!   - cross product
!   - NewtonRaphson
!     - lnsrch
!   - K-number 

	interface Lu_system
		MODULE PROCEDURE Lu_system_r,Lu_system_d
	end interface
	interface inverse
		MODULE PROCEDURE inverseV_r,inverseV_d,inverseM_r,inverseM_d
	end interface
	interface det
		MODULE PROCEDURE det33_r,det33_d
	end interface
	interface norme
		MODULE PROCEDURE norme_r,norme_d
	end interface
	interface cross_product
		MODULE PROCEDURE cross_product_r,cross_product_d
	end interface
	interface NewtonRaphson
		MODULE PROCEDURE NewtonRaphson_r,NewtonRaphson_d
	end interface
	interface Knumber
	    MODULE PROCEDURE KnumberV_r,KnumberV_d
	end interface


CONTAINS

  subroutine Lu_system_d(A,B,ier,Toler,det)

  Use NumericalRecipes, ONLY: ludcmp,lubksb

  ! External variables

  real(8), dimension(:,:)              :: A
  real(8), dimension(:)                :: B
  integer                              :: ier
  real(8), optional                    :: Toler
  real(8), optional                    :: det

  ! Internal variables

  integer, dimension(size(B))          :: indx  
  integer                              :: i
  integer                              :: jer
  real(8)                              :: d
  real(8)                              :: factor

    ier = 0

factor = 1.d0
!    factor = dsqrt(dot_product(reshape(A,(/size(A)/)),reshape(A,(/size(A)/))))
	A = A / factor

    call ludcmp(A,indx,d,jer)
	if(jer > 0) then
      ier = 2
      B = 0.0
      if(present(det)) det = 0.d0
      return
	endif
    d = 1.d0
    do i=1,size(B)
      d = d*A(i,i)
    enddo
    if(present(det)) det = d

    if(present(Toler)) then
      if(dabs(d) < Toler) then
        ier = 1
        B = 0.d0
        return
	  endif
    else
      if(dabs(d) < 1.d-44) then
        ier = 1
        B = 0.d0
        return
	  endif
    endif

	B = B / factor

    call lubksb(A,indx,B)

  end subroutine Lu_system_d

  subroutine Lu_system_r(A,B,ier,Toler,det)

  Use NumericalRecipes, ONLY: ludcmp,lubksb

  ! External variables

  real(4), dimension(:,:)              :: A
  real(4), dimension(:)                :: B
  integer                              :: ier
  real(4), optional                    :: Toler
  real(4), optional                    :: det

  ! Internal variables

  integer, dimension(size(B))          :: indx  
  integer                              :: i
  integer                              :: jer
  real(4)                              :: d
  real(4)                              :: factor

    ier = 0

factor = 1.
!    factor = sqrt(dot_product(reshape(A,(/size(A)/)),reshape(A,(/size(A)/))))
	A = A / factor

    call ludcmp(A,indx,d,jer)
	if(jer > 0) then
      ier = 2
      B = 0.0
      if(present(det)) det = 0.
      return
	endif

    d = 1.0
    do i=1,size(B)
      d = d*A(i,i)
    enddo
    if(present(det)) det = d

    if(present(Toler)) then
      if(abs(d) < Toler) then
        ier = 1
        B = 0.0
        return
	  endif
    else
      if(abs(d) < 1.e-44) then
        ier = 1
        B = 0.0
        return
	  endif
    endif

	B = B / factor

    call lubksb(A,indx,B)

  end subroutine Lu_system_r

  function det33_r(a) RESULT(x)

  ! External variables

    real(4), dimension(3,3), intent(in) :: a

  ! Internal variables

    real(4)                             :: x


    x = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) + &
        a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) + &
        a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

  end function det33_r

  function det33_d(a) RESULT(x)

  ! External variables

  real(8), dimension(3,3), intent(in) :: a

  ! Internal variables

  real(8)                             :: x


    x = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) + &
        a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) + &
        a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

  end function det33_d


  function cross_product_r(a,b) RESULT(c)

  ! External variables

  real(4), dimension(3), intent(in)  :: a
  real(4), dimension(3), intent(in)  :: b

  ! Internal variables

  real(4), dimension(3)              :: c


    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  end function cross_product_r

  function cross_product_d(a,b) RESULT(c)

  ! External variables

  real(8), dimension(3), intent(in)  :: a
  real(8), dimension(3), intent(in)  :: b

  ! Internal variables

  real(8), dimension(3)              :: c


    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  end function cross_product_d

  function Norme_r(x,wgt) RESULT(xnor)

  ! External variables

  real(4), dimension(:)                :: x
  real(4), dimension(:), optional      :: wgt
  real(4)                              :: xnor

  ! Internal variables

  integer                              :: i,n,adr
  real(4), parameter                   :: TOLER = 1.e-21


    xnor = 0.0

    if(present(wgt)) then
      n = size(x)/size(wgt)
      adr = 0
      do i=1,size(wgt)
        xnor = xnor + sum(x(adr+1:adr+n)**2) * wgt(i)
  	    adr = adr + n
      enddo
    else
      xnor = sum(x**2)
    endif

	if(xnor > TOLER) then
      xnor = sqrt(xnor)
	else
	  xnor = 0.
	endif
  
  end function Norme_r

  function Norme_d(x,wgt) RESULT(xnor)

  ! External variables

  real(8), dimension(:)                :: x
  real(8), dimension(:), optional      :: wgt
  real(8)                              :: xnor

  ! Internal variables

  integer                              :: i,n,adr
  real(4), parameter                   :: TOLER = 1.d-44


    xnor = 0.d0

    if(present(wgt)) then
      n = size(x)/size(wgt)
      adr = 0
      do i=1,size(wgt)
        xnor = xnor + sum(x(adr+1:adr+n)**2) * wgt(i)
  	    adr = adr + n
      enddo
    else
      xnor = sum(x**2)
    endif

	if(xnor > TOLER) then
      xnor = dsqrt(xnor)
	else
	  xnor = 0.d0
	endif  

  end function Norme_d

  
  function inverseV_r(A,det,ier) RESULT(AINV)
  
  ! inverse of A. Use LU decomposition

  Use NumericalRecipes, ONLY: ludcmp,lubksb

  ! External variables

  real(4), dimension(:)                         :: A
  real(4), optional                             :: det
  integer, optional                             :: ier
  real(4), dimension(size(A))                   :: AINV

  ! Internal variables

  integer                                       :: i,j,n,adr
  integer, dimension(size(A))                   :: indx  ! only use sqrt(size(A)) elements
  integer                                       :: jer
  real(4)                                       :: d
  real(4), parameter                            :: TOLER = 1.e-44
  real(4), dimension(:,:), allocatable          :: AUX
  real(4)                                       :: factor


    if(present(ier)) ier = 0

    if(present(det) .or. present(ier)) then
	  if(all(A == 0.0)) then
        if(present(det)) det = 0.0
        if(present(ier)) ier = 2
        AINV = 0.0
        return
	  endif
	endif

    n = sqrt(float(size(A))) 
    factor = maxval(abs(A))

    allocate(AUX(n,n))
    AUX = reshape(A,(/n,n/)) / factor

    AINV = 0.0

    call ludcmp(AUX,indx(1:n),d,jer)
	if(jer > 0 .and. present(ier)) then
      ier = 2
      AINV = 0.0
      return
	endif

    d = 1.0
    do i=1,n
      d = d*AUX(i,i)
    enddo
	d=d*factor
    if(present(det)) det=d

    if(abs(d) < TOLER) then
      if(present(ier)) ier = 1
      AINV = 0.d0
	  deallocate(AUX)
      return
    endif

    adr = 0
    do j=1,n
      AINV(adr+j) = 1.0
      call lubksb(AUX,indx(1:n),AINV(adr+1:adr+n))
      adr = adr + n
    enddo

    AINV = AINV / factor
    deallocate(AUX)

  end function inverseV_r

  function inverseV_d(A,det,ier) RESULT(AINV)
  
  ! inverse of A. Use LU decomposition

  Use NumericalRecipes, ONLY: ludcmp,lubksb

  ! External variables

  real(8), dimension(:)                         :: A
  real(8), optional                             :: det
  integer, optional                             :: ier
  real(8), dimension(size(A))                   :: AINV

  ! Internal variables

  integer                                       :: i,j,n,adr
  integer, dimension(size(A))                   :: indx  ! only use sqrt(size(A)) elements
  integer                                       :: jer
  real(8)                                       :: d
  real(8), parameter                            :: TOLER = 1.d-44
  real(8), dimension(:,:), allocatable          :: AUX
  real(8)                                       :: factor


    if(present(ier)) ier = 0

    if(present(det) .or. present(ier)) then
	  if(all(A == 0.d0)) then
        if(present(det)) det = 0.d0
        if(present(ier)) ier = 2
        AINV = 0.d0
        return
	  endif
	endif

    n = sqrt(float(size(A))) 
    factor = maxval(dabs(A))

    allocate(AUX(n,n))
    AUX = reshape(A,(/n,n/)) / factor

    AINV = 0.d0

    call ludcmp(AUX,indx(1:n),d,jer)
	if(jer > 0 .and. present(ier)) then
      ier = 2
      AINV = 0.0
      return
	endif

    d = 1.d0
    do i=1,n
      d = d*AUX(i,i)
    enddo
	d=d*factor
    if(present(det)) det=d

    if(dabs(d) < TOLER) then
      if(present(ier)) ier = 1
      AINV = 0.d0
	  deallocate(AUX)
      return
    endif

    adr = 0
    do j=1,n
      AINV(adr+j) = 1.d0
      call lubksb(AUX,indx(1:n),AINV(adr+1:adr+n))
      adr = adr + n
    enddo

    AINV = AINV / factor
    deallocate(AUX)

  end function inverseV_d

  function inverseM_r(A,det,ier) RESULT(AINV)
  
  ! inverse of A. Use LU decomposition

  Use NumericalRecipes, ONLY: ludcmp,lubksb

  ! External variables
  
  real(4), dimension(:,:)                         :: A
  real(4), optional                               :: det
  integer, optional                               :: ier
  real(4), dimension(size(A,1),size(A,1))         :: AINV

  ! Internal variables
  
  integer                                         :: i,j,n
  integer, dimension(size(A,1))                   :: indx  
  integer                                         :: jer
  real(4)                                         :: d
  real(4), parameter                              :: TOLER = 1.e-44
  real(4), dimension(:,:), allocatable            :: AUX
  real(4)                                         :: factor


    if(present(ier)) ier = 0
 
    if(present(det) .or. present(ier)) then
	  if(all(A == 0.0)) then
        if(present(det)) det = 0.0
        if(present(ier)) ier = 2
        AINV = 0.0
        return
	  endif
	endif

    n = size(A,1)
    factor = maxval(abs(A))

    allocate(AUX(n,n))
    AUX = A / factor

    AINV = 0.0

    call ludcmp(AUX,indx(1:n),d,jer)
	if(jer > 0 .and. present(ier)) then
      ier = 2
      AINV = 0.0
      return
	endif
    d = 1.d0
    do i=1,n
      d = d*AUX(i,i)
    enddo
	d=d*factor
    if(present(det)) det=d

    if(abs(d) < TOLER) then
      if(present(ier)) ier = 1
      AINV = 0.0
	  deallocate(AUX)
      return
    endif

    do j=1,n
      AINV(j,j) = 1.0
      call lubksb(AUX,indx(1:n),AINV(:,j))
    enddo

    AINV = AINV / factor
    deallocate(AUX)

  end function inverseM_r

  function inverseM_d(A,det,ier) RESULT(AINV)
  
  ! inverse of A. Use LU decomposition

  Use NumericalRecipes, ONLY: ludcmp,lubksb

  ! External variables
  
  real(8), dimension(:,:)                         :: A
  real(8), optional                               :: det
  integer, optional                               :: ier
  real(8), dimension(size(A,1),size(A,1))         :: AINV

  ! Internal variables
  
  integer                                         :: i,j,n
  integer, dimension(size(A,1))                   :: indx  
  integer                                         :: jer
  real(8)                                         :: d
  real(8), parameter                              :: TOLER = 1.d-44
  real(8), dimension(:,:), allocatable            :: AUX
  real(8)                                         :: factor


    if(present(ier)) ier = 0
 
    if(present(det) .or. present(ier)) then
	  if(all(A == 0.d0)) then
        if(present(det)) det = 0.d0
        if(present(ier)) ier = 2
        AINV = 0.d0
        return
	  endif
	endif

    n = size(A,1)
    factor = maxval(dabs(A))

    allocate(AUX(n,n))
    AUX = A / factor

    AINV = 0.d0

    call ludcmp(AUX,indx(1:n),d,jer)
	if(jer > 0 .and. present(ier)) then
      ier = 2
      AINV = 0.0
      return
	endif
    d = 1.d0
    do i=1,n
      d = d*AUX(i,i)
    enddo
	d=d*factor
    if(present(det)) det=d

    if(dabs(d) < TOLER) then
      if(present(ier)) ier = 1
      AINV = 0.d0
	  deallocate(AUX)
      return
    endif

    do j=1,n
      AINV(j,j) = 1.d0
      call lubksb(AUX,indx(1:n),AINV(:,j))
    enddo

    AINV = AINV / factor
    deallocate(AUX)

  end function inverseM_d

  function newt(func,x,fx,tol,niter) RESULT(ier)

  Use NumericalRecipes, ONLY: ludcmp,lubksb

  ! External variables

  interface 
           function func(x,dfdx) RESULT(fx)
                   real(8), dimension(:)             :: x
                   real(8), dimension(:,:), optional :: dfdx
                   real(8), dimension(size(x))       :: fx
           end function
  end interface

  real(8), dimension(:)                              :: x
  real(8), dimension(:)                              :: fx
  real(8), intent(in)                                :: TOL
  integer, intent(in)                                :: niter
  integer                                            :: ier

  ! Internal variables

  real(8), dimension(size(x))                        :: xold
  real(8), dimension(size(x),size(x))                :: dfdx
  real(8), parameter                                 :: TOLF   = 1.d-4
  real(8), parameter                                 :: TOLMIN = epsilon(x)
  real(8), parameter                                 :: STPMX  = 100.d0
  real(8)                                            :: f,fold,stpmax
  real(8), dimension(size(x))                        :: dx,g
  logical                                            :: check

  integer                                            :: iersys
  integer                                            :: k


	fx = func(x)
	f = 0.5d0*dot_product(fx,fx)

	if(maxval(dabs(fx)) < 0.01d0*TOLF) then
	  ier = 1
	  return
	endif
	stpmax = STPMX*max(0.5*dsqrt(dot_product(x,x)),real(size(x)))

    ier = 0
	do k=1,niter
      fx = func(x,dfdx)
	  g = matmul(fx,dfdx)

	  xold = x
	  fold = f

	  dx = -fx
      call lu_system(dfdx,dx,iersys)
      
	  call lnsrch(func,xold,fold,g,dx,x,f,stpmax,check)
      
	  if(maxval(dabs(fx)) < TOLF) then
		ier = 1
		return
	  endif

	  if(check) then
		check=(maxval(dabs(g)*max(dabs(x),1.0d0)/max(f,0.5d0*size(x))) < TOLMIN)
		ier = 2
	    return
	  endif

	  if(maxval(dabs(x-xold)/max(abs(x),1.0d0)) < TOL) return
	enddo

	ier = 99

  CONTAINS
	
	subroutine lnsrch(func,xold,fold,g,dx,x,f,stpmax,check)

  ! External variables

    interface 
             function func(x,dfdx) RESULT(fx)
                     real(8), dimension(:)             :: x
                     real(8), dimension(:,:), optional :: dfdx
                     real(8), dimension(size(x))       :: fx
             end function
    end interface
	real(8), dimension(:), intent(in)                  :: xold
	real(8), intent(in)                                :: fold
	real(8), dimension(:), intent(in)                  :: g
	real(8), dimension(:), intent(inout)               :: dx
	real(8), dimension(:), intent(out)                 :: x
	real(8), intent(out)                               :: f
	real(8), intent(in)                                :: stpmax
	logical, intent(out)                               :: check

  ! Internal variables

	real(8), parameter                                 :: ALF = 1.0e-04
    real(4)                                            :: slope,alamin,alam,alam2
	real(4)                                            :: a,b,rhs1,rhs2
	real(4)                                            :: pabs,tmplam
    real(4)                                            :: disc,f2,fold2
		
	  check=.false.
	
	  pabs = dsqrt(dot_product(dx,dx))
	  if(pabs > stpmax) dx = dx*stpmax/pabs

	  slope  = dot_product(g,dx)
	  alamin = TOLX/maxval(dabs(dx)/max(dabs(xold),1.d0))
	  alam   = 1.0

	  do
		x = xold + alam*dx
		f = 0.5d0*dot_product(func(x),func(x))
		if(alam < alamin) then
		  x = xold
		  check = .true.
		  return
		elseif(f <= fold + ALF*alam*slope) then
	      return
		else
		  if(alam == 1.0) then
			tmplam = -slope/(2.d0*(f-fold-slope))
		  else
			rhs1 = f - fold - alam*slope
			rhs2 = f2 - fold2 - alam2*slope

			a = (rhs1/alam**2 - rhs2/alam2**2)/(alam - alam2)
			b = (-alam2*rhs1/alam**2 + alam*rhs2/alam2**2)/(alam - alam2)
			
			if(a == 0.0) then
			  tmplam = -slope/(2.d0*b)
			else
			  disc = b*b - 3.d0*a*slope
			  if(disc < 0.0) then
			    write(*,'(" roundoff problem in lnsrch ")')
				pause
			  endif
			  tmplam = (-b + sqrt(disc))/(3.d0*a)
			endif
			if(tmplam > 0.5*alam) tmplam = 0.5d0*alam
	      endif
		endif

		alam2 = alam
		f2    = f
		fold2 = fold
		alam = max(tmplam,0.1*alam)
	  end do

	end subroutine lnsrch
  end function newt

  function NewtonRaphson_r(func,x,fx,eps,niter) RESULT(ier)

  ! External variables

  interface 
           function func(x,dfdx) RESULT(fx)
                   real(4), dimension(:)             :: x
                   real(4), dimension(:,:), optional :: dfdx
                   real(4), dimension(size(x))       :: fx
           end function
  end interface

  real(4), dimension(:), intent(inout)               :: x
  real(4), dimension(:), optional, intent(out)       :: fx
  real(4), intent(in)                                :: eps
  integer, intent(in)                                :: niter
  integer                                            :: ier

  ! Internal variables

  integer                                            :: iersys,ierlnsrch
  integer                                            :: k
  real(4), dimension(size(x))                        :: f
  real(4), dimension(size(x))                        :: dx
  real(4), dimension(size(x))                        :: xold
  real(4), dimension(size(x),size(x))                :: dfdx
  real(4)                                            :: err


    ier = 0

    do k=1,niter
      f = func(x,dfdx)
      err = sqrt(dot_product(f,f))
	  xold = x

   	  dx = -f
      call lu_system(dfdx,dx,iersys)
	  if(iersys /= 0) then
	    ier = 99
		return
	  endif

      ierlnsrch = lnsrch(xold,dx,err,x,func)
	  if(ierlnsrch /= 0) then
	    ier = 9
		return
	  endif

	  if(sum(abs(dx)) <= eps) then
	    if(present(fx)) fx = func(x)
	    return
	  endif
  	enddo

	ier = 1


  CONTAINS

    function lnsrch(xold,dx,err0,x,func) RESULT(ierlnsrch)

    ! External variables

      real(4), dimension(:), intent(in)                  :: xold
      real(4), dimension(:), intent(in)                  :: dx
      real(4), intent(in)                                :: err0
      real(4), dimension(:), intent(out)                 :: x
      interface 
               function func(x,dfdx) RESULT(fx)
                     real(4), dimension(:)             :: x
                     real(4), dimension(:,:), optional :: dfdx
                     real(4), dimension(size(x))       :: fx
               end function
      end interface
	  integer                                            :: ierlnsrch

    ! Internal variables
 
      integer                                            :: k
      real(4)                                            :: err
      real(4), dimension(4), parameter                   :: ALPHA = (/1.e+00, &
                                                                      1.e-01, &
                                                                      1.e-02, &
                                                                      1.e-03/)
      real(4), dimension(size(x))                        :: fx

      
	  ierlnsrch = 0

      k= 0 
	  do
	    k = k + 1
	    x = xold + alpha(k)*dx
	    fx = func(x)
	
        err = sqrt(dot_product(fx,fx))
        if(err < err0) return

	    if(k == size(alpha)) exit
      enddo

      x = xold
	  ierlnsrch = 1

    end function lnsrch
  end function NewtonRaphson_r


  function NewtonRaphson_d(func,x,fx,eps,niter) RESULT(ier)

  ! External variables

  interface 
           function func(x,dfdx) RESULT(fx)
                   real(8), dimension(:)             :: x
                   real(8), dimension(:,:), optional :: dfdx
                   real(8), dimension(size(x))       :: fx
           end function
  end interface

  real(8), dimension(:), intent(inout)               :: x
  real(8), dimension(:), optional, intent(out)       :: fx
  real(8), intent(in)                                :: eps
  integer, intent(in)                                :: niter
  integer                                            :: ier

  ! Internal variables

  integer                                            :: iersys,ierlnsrch
  integer                                            :: k
  real(8), dimension(size(x))                        :: f
  real(8), dimension(size(x))                        :: dx
  real(8), dimension(size(x))                        :: xold
  real(8), dimension(size(x),size(x))                :: dfdx
  real(8)                                            :: err


    ier = 0

    do k=1,niter
      f = func(x,dfdx)
      err = dsqrt(dot_product(f,f))
	  xold = x

   	  dx = -f
      call lu_system(dfdx,dx,iersys)
	  if(iersys /= 0) then
	    ier = 99
		return
	  endif

      ierlnsrch = lnsrch(xold,dx,err,x,func)
	  if(ierlnsrch /= 0) then
	    ier = 9
		return
	  endif

	  if(sum(dabs(dx)) <= eps) then
	    if(present(fx)) fx = func(x)
	    return
	  endif
  	enddo

	ier = 1


  CONTAINS

    function lnsrch(xold,dx,err0,x,func) RESULT(ierlnsrch)

    ! External variables

      real(8), dimension(:), intent(in)                  :: xold
      real(8), dimension(:), intent(in)                  :: dx
      real(8), intent(in)                                :: err0
      real(8), dimension(:), intent(out)                 :: x
      interface 
               function func(x,dfdx) RESULT(fx)
                     real(8), dimension(:)             :: x
                     real(8), dimension(:,:), optional :: dfdx
                     real(8), dimension(size(x))       :: fx
               end function
      end interface
	  integer                                            :: ierlnsrch

    ! Internal variables
 
      integer                                            :: k
      real(8)                                            :: err
      real(8), dimension(4), parameter                   :: ALPHA = (/1.d+00, &
                                                                      1.d-01, &
                                                                      1.d-02, &
                                                                      1.d-03/)
      real(8), dimension(size(x))                        :: fx

      
	  ierlnsrch = 0

      k= 0 
	  do
	    k = k + 1
	    x = xold + alpha(k)*dx
	    fx = func(x)
	
        err = dsqrt(dot_product(fx,fx))
        if(err < err0) return

	    if(k == size(alpha)) exit
      enddo

      x = xold
	  ierlnsrch = 1

    end function lnsrch
  end function NewtonRaphson_d


  function KnumberV_r(A,ier) result(k)

  ! External variables

  integer                     :: ier
  real(4), dimension(:)       :: A
  real(4)                     :: k

  ! Internal variables

  real(4), dimension(size(A)) :: Ainv
  real(4), dimension(1)       :: val


    Ainv = inverseV_r(A,ier=ier)
	if(ier /= 0) then
	  k = -1.0
	  return
	endif

	val = maxval(abs(A))*maxval(abs(Ainv))

    k = val(1)

  end function

  function KnumberV_d(A,ier) result(k)
  ! External variables

  integer                     :: ier
  real(8), dimension(:)       :: A
  real(8)                     :: k

  ! Internal variables

  real(8), dimension(size(A)) :: Ainv
  real(8), dimension(1)       :: val


    Ainv = inverseV_d(A,ier=ier)
	if(ier /= 0) then
	  k = -1.d0
	  return
	endif

	val = maxval(dabs(A))*maxval(dabs(Ainv))

    k = val(1)

  end function

end Module Math
