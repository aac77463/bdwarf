subroutine remove_dust(mr, ra, rac, du, duc, pebs, pebsc, tcrit, dtmin, pebint)
  use constants
  implicit none
  integer, parameter :: dp = selected_real_kind(15)
  integer :: mr, i
  real(dp), dimension(mr), intent(in) :: ra, rac
  real(dp), dimension(mr), intent(inout) :: du, pebs, duc, pebsc
  real(dp), dimension(mr), intent(inout) :: tcrit
  real(dp), dimension(mr), intent(out) :: pebint
  real(dp), dimension(mr) :: tots, pebintc
  real(dp) :: dtmin, flim
  flim = 0.75 !The amount of material that we want in pebbles

  !Using an exponential such that as dtmin --> inf
  !Then 100% of dust is converted to pebbles
  do i = 1, mr
     pebint(i) = (1. - exp((-dtmin)/tcrit(i)))
  enddo
  call interp1d(mr, ra, pebint, rac, pebintc)

  do i = 1, mr-1
     pebsc(i) = pebsc(i) + ((duc(i)) * pebintc(i))
     duc(i) = duc(i) - ((duc(i)) * pebintc(i))
     if(duc(i).lt.1.e-300)then
        duc(i) = 1.e-300
     endif
  enddo

  tots = pebsc / (duc + pebsc)
  
!!!  !If the amount of pebbles is larger than flim we
!!!  !convert (1-flim) back into dust.
  do i = 1, mr-1
     if(tots(i).gt.flim)then
        tots(i) = duc(i) + pebsc(i)
        duc(i) = (tots(i))*(1.-flim)
        pebsc(i) = (tots(i))*flim
     endif
  enddo

end subroutine remove_dust
