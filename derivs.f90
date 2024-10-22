subroutine cdiff_poly(nr, y, x, dy)
  use constants
  implicit none
  integer, parameter :: dp = selected_real_kind(15)
  integer :: i, nr
  real(dp) :: a, b
  real(dp), dimension(nr), intent(in) :: x, y
  real(dp), dimension(nr), intent(out) :: dy

  dy(1) = (y(2)-y(1))/(x(2)-x(1))
  dy(nr) = (y(nr)-y(nr-1))/(x(nr)-x(nr-1))

  do i = 2,nr-1
!     a = ( (y(i-1) - y(i))*(x(i)-x(i+1)) - (y(i) - y(i+1)) * (x(i-1) - x(i)) ) / &
!          ( (( x(i-1)**2. - x(i)**2. ) * (x(i)-x(i+1)) ) - ( (x(i)**2. - x(i+1)**2.) * ( x(!i-1) - x(i) ) ) )
     a = ( (y(i-1) - y(i))*(x(i)-x(i+1)) - (y(i) - y(i+1)) * (x(i-1) - x(i)) ) / &
          ( (( x(i-1)**2. - x(i)**2. ) * (x(i)-x(i+1)) ) - ( (x(i)**2. - x(i+1)**2.) * ( x(i-1) - x(i) ) ) )
     b = ( y(i-1) - y(i) - a * (x(i-1)**2. - x(i)**2. ) )/ (x(i-1) - x(i)) 
     dy(i) = 2. * a * x(i) + b
!     dy(i)=(y(i-1)-y(i+1))/(x(i-1)-x(i+1))
!     if(i.eq.2)then
!        dy(1) = 2.* a * x(1) + b
!        dy(1) = (y(2)-y(1))/(x(2)-x(1))
!     endif
!     if(i.eq.nr-1)then
!        dy(nr) = 2.* a * x(nr) + b
!        dy(nr)=(y(nr-1)-y(nr))/(x(nr-1)-x(nr))
!     endif
  enddo
  
  RETURN
end subroutine cdiff_poly

subroutine cdiff_poly_second(nr, y, x, dy2)
  use constants
  implicit none
  integer, parameter :: dp = selected_real_kind(15)
  integer :: i, nr
  real(dp) :: a, b
  real(dp), dimension(nr), intent(in) :: x, y
  real(dp), dimension(nr), intent(out) :: dy2

  dy2(1)=0
  dy2(nr)=0
  
  do i = 2,nr-1
     a = ( (y(i-1) - y(i))*(x(i)-x(i+1)) - (y(i) - y(i+1)) * (x(i-1) - x(i)) ) / &
          ( (( x(i-1)**2. - x(i)**2. ) * (x(i)-x(i+1)) ) - ( (x(i)**2. - x(i+1)**2.) * ( x(i-1) - x(i) ) ) )
     b = ( y(i-1) - y(i) - a * (x(i-1)**2. - x(i)**2. ) )/ (x(i-1) - x(i)) 
!     dy(i) = 2. * a * x(i) + b
     dy2(i) = 2. * a
!     dy(i)=(y(i-1)-y(i+1))/(x(i-1)-x(i+1))
!     if(i.eq.2)then
!        dy2(1) = 2.* a! * x(1) + b
!        dy(1) = (y(2)-y(1))/(x(2)-x(1))
!     endif
!     if(i.eq.nr-1)then
!        dy2(nr) = 2.* a! * x(nr) + b
!        dy(nr)=(y(nr-1)-y(nr))/(x(nr-1)-x(nr))
!     endif
  enddo
  
  RETURN
end subroutine cdiff_poly_second

subroutine MidpointD1(n,y,x,dy)
!first derivative on an irregular grid
  IMPLICIT NONE
  integer, parameter :: dp = selected_real_kind(15)
  integer n,i
  real(dp) :: x(n),y(n),dy(n),dx1,dx2

  dy(1)=(y(2)-y(1))/(x(2)-x(1))
  dy(n)=(y(n)-y(n-1))/(x(n)-x(n-1))
  do i=2,n-1
     dx1=(x(i)-x(i-1))
     dx2=(x(i+1)-x(i))
     dy(i)=(y(i)-y(i-1))/dx1+dx1*((y(i+1)-y(i))/dx2-(y(i)-y(i-1))/dx1)/(dx1+dx2)
  enddo

!  print *, '1st', dy(111-2:111+2)
!  read*
  
  return
end subroutine MidpointD1

subroutine MidpointD2(n,y,x,dy2)
  !	second derivative on an irregular grid
  IMPLICIT NONE
  integer, parameter :: dp = selected_real_kind(15)
  integer n,i
  real(dp) ::  x(n),y(n),dy2(n),dy11,dy12,x1,x2
  
  do i=2,n-1
     dy11=(y(i)-y(i-1))/(x(i)-x(i-1))
     dy12=(y(i+1)-y(i))/(x(i+1)-x(i))
     x1=(x(i)+x(i-1))/2.
     x2=(x(i)+x(i+1))/2.
!     if ((i.lt.118+2).and.(i.gt.118-2)) then
!        print *, 'i', i
!        print *, 'dy11 - dy12', dy12 - dy11
!        print *, 'dy12 - dy11', dy11 - dy12
!        print *, 'x2 - x1', x2-x1
!        print *, 'x1 - x2', x1-x2
!        read *
!     endif
!     dy2(i)=(dy11-dy12)/(x1-x2)!/2.
!     dy2(i) = ((2.*y(i-1)/((x(i)-x(i-1))*(x(i+1)-x(i-1)))-&
!          (2.*y(i)/((x(i+1)-x(i))*(x(i)-x(i-1))))+&
!          (2.*y(i+1)/((x(i+1)-x(i))*(x(i+1)-x(i-1))))))
  enddo
!  dy2(1)=dy2(2)
!  dy2(n)=dy2(n-1)
  dy2(1)=0
  dy2(n)=0
  
!  print *, '2nd', dy2(111-2:111+2)
!  read*
  
  return
end subroutine MidpointD2

subroutine rugd2(n,y,x,dy2)
  !	second derivative on an irregular grid
  IMPLICIT NONE
  integer, parameter :: dp = selected_real_kind(15)
  integer n,i
  real(dp) ::  x(n),y(n),dy2(n),hplus,hmin,yplus,ymin
  
  do i=2,n-1
     hplus = x(i+1) - x(i)
     hmin = x(i) - x(i-1)
     yplus = y(i+1)
     ymin = y(i-1)

     dy2(i)=((hmin*yplus) - ((hplus+hmin)*y(i)) + (hplus*ymin))/((1./2.)*&
          (((hplus*hmin)*(hplus+hmin))))

!     if((i.gt.116.).and.(i.lt.120))then
!        print *, 'i', i
!        print *, 'h+', hplus
!        print *, 'h-', hmin
!        print *, 'y+', yplus
!        print *, 'y-', ymin
!        print *, 'num', ((hmin*yplus) - ((hplus+hmin)*y(i)) + (hplus*ymin))
!        print *, 'den', ((1./2.)*(((hplus*hmin)*(hplus+hmin))))
!        print *, 'dy2', dy2(i)
!        read *
!    endif

  enddo
  dy2(1)=0
  dy2(n)=0
  
  return
end subroutine rugd2
	




 
