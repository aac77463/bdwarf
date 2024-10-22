subroutine remove(mr, quant, ra, rpl, macc_q, dtmin, mpl, mstar)
  use constants
  implicit none
  integer, parameter :: dp = selected_real_kind(15)
  integer :: mr, i, ixinrh, ixoutrh
  real(dp) :: dtmin, mpl, mstar, r_h, min, mf, rpl
  real(dp), dimension(mr) :: ra, quant
  real(dp), dimension(mr) :: ra2, macc_q
!  rpl = 1.
  
!  !Find the grind for the 1/2 interpolation
  do i = 2, mr
     ra2(i) = ( ra(i) + ra(i-1) ) / 2.  
  enddo
  !
  !  !Finding how big the hill radius is
  r_h = (rpl * au) * ( mpl / (3. * mstar) )**(1./3.)
!
  ixoutrh=minloc(abs((ra-((rpl * au)+r_h))/((rpl * au)+r_h)), dim=1)
  ixinrh=minloc(abs((ra-((rpl * au)-r_h))/((rpl * au)-r_h)), dim=1)
!  ixoutrh=minloc(abs((((rpl * au)+r_h)-ra)/((rpl * au)+r_h)), dim=1)
!  ixinrh=minloc(abs((ra-((rpl * au)-r_h))/((rpl * au)-r_h)), dim=1)
!
  min = 0.
  do i = ixinrh, ixoutrh
     min = min + 2. * pi * quant(i) * ra(i) * (ra2(i+1) - ra2(i))
  enddo

  do i = ixinrh, ixoutrh
     quant(i) = ( 2. * pi * ( ra2(i+1) - ra2(i) ) * quant(i) * ra(i) -&
          (macc_q(i) * dtmin) ) / (2 * pi * (ra2(i+1) - ra2(i)) * ra(i) )
  enddo

!  call interp1d(mr, ra, quant, rac, quantc)
  
  mf = 0.
  do i = ixinrh, ixoutrh
     mf = mf + (2. * pi * quant(i) * ra(i) * (ra2(i+1) - ra2(i)))
  enddo

!  print *, macc_q(ixinrh:ixoutrh)/m_earth
!  print *, (min - mf)/m_earth
!  print *, mpl
  mpl = mpl + (min - mf)

end subroutine remove
