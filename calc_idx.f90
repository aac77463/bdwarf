subroutine calc_idx(mr, ra, rpl, idx)

  use constants
  implicit none
  integer, parameter :: dp = selected_real_kind(15)
  integer :: i, idx, mr
  real(dp), dimension(mr) :: diff, ra
  real(dp) :: idxpl
  real(dp), intent(inout) :: rpl
  
  do i = 1, mr
     diff(i) = abs((ra(i)) - (rpl*au))
  enddo
!  print *, 'dist', abs((ra - (rpl*au)))

  !Find the planet location index
  idx = int(minloc(diff(1:mr), dim=1))
  idxpl = ra(idx)
!  print *, 'pl index', idx

end subroutine calc_idx
