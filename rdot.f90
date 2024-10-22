subroutine migration(mr, ra, sigg, tm, sh, omg, m_star, beta, zeta, mpl, rpl, rmig, dtmin, tinit, ty, idx)

!!$ This subroutine calculates the radial change in [cm/s] for a planet
!!$ undergoing Type 1 migration. I use the formulation in Johansen+19.
!!$ We calculate the total amount of migration. If that is larger than
!!$ the separation of the planet and the i-1 point in the grid, then the
!!$ planet migrates to the i-1 point, which now becomes i. This also resets
!!$ the counter for the migration and we can begin again

  use constants
  implicit none
  integer, parameter :: dp = selected_real_kind(15)
  integer :: mr, idx
  real(dp), dimension(mr), intent(in) :: ra, sigg, tm, sh, omg
  real(dp), dimension(mr) :: beta, zeta, rdot
  real(dp), intent(inout) :: rmig, rpl
  real(dp), intent(in) :: dtmin, tinit, m_star, mpl, ty
  real(dp) :: kmig, rcount, f
  !f=6.e-2
  f=1.

  call calc_idx(mr, ra, rpl, idx)
  
  !Calculate the gradients for the temperature and the gas surface density
  call cdiff_poly(mr, log(sigg(1:mr-1)), log(ra(1:mr-1)), beta)
  call cdiff_poly(mr, log(tm(1:mr-1)), log(ra(1:mr-1)), zeta)
  
  !Calculating the kmig quantity for migration
  !Based on D'Angelo Lubow 3D simulations
  !kmig = 2. * (1.36 + (0.62 * (-1.*beta(idx))) + (0.43 * (-1.*zeta(idx))))
  !Based on not using isothermal
  !kmig = (2.7 + 1.1*beta(idx))
  
  !Calculate the rate of migration [cm/s]
  rdot = (mpl/m_star) * ((sigg * ra**2.)/m_star) * (sh/ra)**(-2.) * (omg * ra) * f

  rcount = rdot(idx) * dtmin

  if(idx.lt.10.)then
  !if(rpl.lt.(0.01*au))then
     rmig = 0.
  !   rcount=0.
  endif
  if(ty.gt.tinit) then
     rmig = rmig + rcount
     if((-1.*rmig).gt.(ra(idx)-ra(idx-1)))then
        rpl = rpl + (rmig/au)
        idx = idx-1
        rmig = 0.
     endif
     if((rmig).gt.(ra(idx)-ra(idx-1)))then
        rpl = rpl - (rmig/au)
        idx = idx+1
        rmig =0.
     endif
  endif

end subroutine migration
