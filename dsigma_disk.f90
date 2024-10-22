subroutine dsigma_disk(mr, sigg, siggc, vs, ra, logsigg, rhs, m_pl, sh, ty, tinit, m_star, &
     gfin, gfout, gfi, rpl, u_g)
!  This subroutine calculates the maximum amount by which the
!  logarithm of the gas surface density can change.
  use constants
  implicit none
  integer, parameter :: dp = selected_real_kind(15)
  integer :: mr, i
  real(dp) :: m_pl, ty, tinit, m_star, gfin, gfout, rpl
  real(dp), dimension(mr), intent(in) :: ra
  real(dp), dimension(mr) :: lograd, radsqr, sqrrad
  real(dp), dimension(mr), intent(in) :: sigg, vs, sh, siggc
  real(dp), dimension(mr), intent(out) :: logsigg
  real(dp), dimension(mr), intent(out) :: gfi
  real(dp), dimension(mr) :: logsinusq, logsinu
  real(dp), dimension(mr) :: sinusq, sinu, sitorad
  real(dp), dimension(mr) :: rhs, tor, u_g
  real(dp), dimension(mr) :: dlogsinusq, dlogsinu, d2logsinusq, dsitorad, dgasv
  
  !log of radial grid, and other quantities for derivatives
  lograd =  log(ra)!
  radsqr = ra**2.!
  sqrrad = ra**(1./2.)!

  !Quantities for calculating right hand term (dlogSigma/dt)
  sinu = sigg * vs!
  logsinu = log(sinu)!
  sinusq = sigg * vs * sqrrad!
  logsinusq = log(sinusq)!
  logsigg = log(sigg)!

  !!!!!!!!!!!!FIRST DERIVATIVES!!!!!!!!!!!!!!
  call cdiff_poly(mr, logsinu, lograd, dlogsinu)!
  call cdiff_poly(mr, logsinusq, lograd, dlogsinusq)!
  call cdiff_poly_second(mr, logsinusq, lograd, d2logsinusq)!

  call torque_vel_calc(mr, ra, sh, m_star, m_pl, tor, rpl)

  !Calculating the third tidal term of the gas surface density to differentiate
  sitorad = tor 
  call cdiff_poly(mr, sitorad, lograd, dsitorad)

!  tinit = 5.e6

  !If time in years is greater than initial planet time
  !We take the gas term with torque
  !Otherwise, disk evolves normally
  if (ty.gt.tinit) then
     rhs = ( (3. * vs ) / radsqr ) * ( (dlogsinu * dlogsinusq) + (d2logsinusq) - &
          ( ( dsitorad ) / ( 3. * vs * sigg) ) )
  else
     rhs = ((3. * vs ) / radsqr) * ( (dlogsinu * dlogsinusq) + (d2logsinusq) )!
  endif

  !Calculating the derivative for the gas velocity
  !(3 / gas * sqr(r)) * d/dr (gas * nu * sqr(r))
  call cdiff_poly(mr, sigg * vs * sqrrad, ra, dgasv)
  
  !Calculating the gas fluxes at the interfaces
  if (ty.lt.tinit) then
     gfi = -(3. /  (sigg * sqrrad)) * dgasv  * sigg
  else
!     gfi = ((-3. /  (sigg * sqrrad)) * dgasv + tor)  * sigg
     gfi = (((-3. /  (sigg * sqrrad)) * dgasv) + tor)  * sigg
  endif
!  print *, gfi

  if(ty.lt.tinit)then
     u_g = -(3. /  (sigg * sqrrad)) * dgasv
  else
!     u_g = ((-3. /  (sigg * sqrrad)) * dgasv) + tor
     u_g = ((-3. /  (sigg * sqrrad)) * dgasv) + tor
  endif
!  print *, 'ug', u_g(1:5)
  
  do i=2, mr-1
     if(u_g(i).lt.0.)then
        gfi(i) = u_g(i) * siggc(i)
     else
        gfi(i) = u_g(i) * siggc(i-1)
     endif
  enddo

  ! Fluxes !
  !<<<----- Inner boundary
  !A positive Flux through the inner boundary is mass lost 
  gfin = gfi(1)

  !----->> Outer boundary
  !A negative flux through the outer boundary is mass lost
  gfout = gfi(mr)

end subroutine dsigma_disk
