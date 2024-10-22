subroutine interp1d(mr, x1, y1, x2, y2)
  implicit none
  integer, parameter :: dp = selected_real_kind(15)
  integer :: i, mr
  real(dp), dimension(mr), intent(in) :: x1, y1
  real(dp), dimension(mr-1), intent(in) :: x2
  real(dp), dimension(mr-1), intent(out) :: y2
  real(dp) :: frac
 
  do i = 1, mr-1
     frac = (y1(i+1) - y1(i)) / (x1(i+1) - x1(i))
!     print *, 'frac', frac
!     print *, y1(i) + (x2(i) - x1(i))
!     print *, y1(i) + ( (x2(i) - x1(i)) * frac )     
     y2(i) = (y1(i) + ( (x2(i) - x1(i)) * frac ))    
  enddo
  
  RETURN
end subroutine interp1d

!make 
