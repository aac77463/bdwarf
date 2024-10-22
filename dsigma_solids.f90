subroutine dsigma_solids(mr, ra, ss, sh, quant, quantc, dlogp, pr, omg, rhsquant, st, sigg, nu, logquant, dlogqusigg, &
     dlogqunu, dlogsinusq, dlogssquom, d2logqusigg, d2logsinusq, d2logp, qafin, qafout, qdfin, &
     qdfout, ty, tinit, m_star, m_pl, tor_vel, afiq, dfiq, rpl)
  use constants
  implicit none
  integer, parameter :: dp = selected_real_kind(15)
  integer :: mr, cr, i
  real(dp), dimension(mr), intent(in) :: ss, dlogp, omg, st, sigg, nu, ra, quant, pr, sh, quantc
  real(dp), dimension(mr), intent(out) :: rhsquant, logquant, dlogqusigg, dlogqunu, dlogsinusq, dlogssquom
  real(dp), dimension(mr), intent(out) :: d2logqusigg, d2logsinusq, d2logp, tor_vel
  real(dp), dimension(mr) :: st1, st2, sssq, radsqr, sqrrad, logqusigg, logqunu, logsinusq, logssquom
  real(dp), dimension(mr) :: qutoal, dqutoal, dqusigg, dsinusq, sinusq, lograd, logp, u_r
  real(dp), dimension(mr), intent(out) :: afiq, dfiq
  real(dp) :: qafin, qafout, qdfin, qdfout, ty, tinit, m_star, m_pl, rpl
  
  logp=log(pr)
  st1 = (1./(1. + st**2.))!initial stokes factor
  st2 = (st/(1. + st**2.))!secondary stokes factor
  radsqr = ra**2.!radius squared
  sssq = ss**2. !Sound speed squared
  sqrrad = sqrt(ra)!square root of the radius
  sinusq = sigg * nu * sqrrad!gas * viscosity * sqr(rad) for flux calc
  lograd = log(ra)!log of the radius

  logquant = log(quant)! log of the quantity, be it dust or pebbles
  logqunu = log(quant * nu * st1)!log of the quantity (pebs/dust) and viscosity * stokes1
  logsinusq = log(sigg * nu * sqrrad)!log of the gas viscosity and sqr(rad)
  logssquom = log(st2 * ss**2. * quant / omg)!log of scale height squared quantity and omega * stokes2
  logqusigg = log(quant/sigg)!log of scale height squared quantity and omega

!  call torque_vel_schib(mr, ra, sh, m_star, m_pl, tor_vel, rpl, sigg)
  call torque_vel_calc(mr, ra, sh, m_star, m_pl, tor_vel, rpl)

  !Calculating the tidal term
  qutoal = quant * tor_vel * st1 * ra
  !Taking linear derivative in lograd for torque term
  call cdiff_poly(mr, qutoal, lograd, dqutoal)

  !Calculating the different partial derivatives in log space for equation
  call cdiff_poly(mr, logqunu, lograd, dlogqunu)
  call cdiff_poly(mr, logqusigg, lograd, dlogqusigg)
  call cdiff_poly(mr, logsinusq, lograd, dlogsinusq)
  call cdiff_poly(mr, logssquom, lograd, dlogssquom)
  call cdiff_poly(mr, logp, lograd, dlogp)
  
  !calculating the second derivatives in log space
  call cdiff_poly_second(mr, logqusigg, lograd, d2logqusigg)
  call cdiff_poly_second(mr, logsinusq, lograd, d2logsinusq)
  call cdiff_poly_second(mr, logp, lograd, d2logp)

  if(ty.gt.tinit)then
     !If ty is larger than the initial time, we have a planet and
     !our equation now includes the tidal torques
     cr = 0
  else
     !Otherwise we take the normal pebble drift speeds
     cr = 5
  endif
  
  select case (cr)
     !Here we choose whether we include the torque calculation into the
     !gas evolution equation or not
  case(0) !With torque addition
     rhsquant =  ( ( st1 * nu ) / radsqr )  * ( (dlogqunu * ( dlogqusigg + 3.* dlogsinusq) ) +&
          d2logqusigg + 3. * d2logsinusq - ( dqutoal  / (nu * quant * st1) ) ) - &
          ( ( ( st2 * sssq) / ( omg * radsqr )) * ( ( dlogssquom * dlogp) + d2logp) )
  case default !Without torque addition
     rhsquant =  ( ( st1 * nu ) / radsqr )  * ( (dlogqunu * ( dlogqusigg + 3.* dlogsinusq) ) +&
          d2logqusigg + 3. * d2logsinusq ) - &
          ( ( ( st2 * sssq) / ( omg * radsqr )) * ( ( dlogssquom * dlogp) + d2logp) )
  end select
  
  !Calculating derivatives for the boundary fluxes
  call cdiff_poly(mr, quant/sigg, ra, dqusigg)
  call cdiff_poly(mr, sinusq, ra, dsinusq)
  
  select case(cr)
     !Here we select the case for the different advection velocities.
     !This velocity changes depending on the addition of the torque term
     !Which changes the gas velocity used to calculate the flux.
  case(0)
!  case default
     !Fluxes
     !Negative flux at inner boundary means loss of mass
     qafin = quant(1) * ( ( st1(1) * ( ( (-3. / ( sigg(1) * sqrrad(1) ) ) *  dsinusq(1)  )   + &
          tor_vel(1) ) )   + &
          ( st2(1) * ( sssq(1) / ( omg(1) * ra(1) ) ) * dlogp(1) ) )
     qdfin = (-1. * nu(1) * st1(1)) * dqusigg(1) * sigg(1)
     !Positive flux at inner boundary means loss of mass
     qafout = quant(mr) * ( ( st1(mr) * ( (  -3. / (sigg(mr) * sqrrad(mr)) )  * dsinusq(mr) ) +&
           tor_vel(mr) ) + &
          ( st2(mr) * ( sssq(mr) / ( omg(mr) * ra(mr) ) ) * dlogp(mr) ) ) 
     qdfout = (-1. * nu(mr) * st1(mr)) * dqusigg(mr) * sigg(mr)
  !case(0)
  case default
     !Fluxes
     !Negative flux at inner boundary means loss of mass
     qafin = quant(1) * ( ( st1(1) * ( ( ( -3. / ( sigg(1) * sqrrad(1) ) )  * dsinusq(1) ) ) ) + &
          ( st2(1) * ( sssq(1) / ( omg(1) * ra(1) ) ) * dlogp(1) ) ) 
     qdfin = (-1. * nu(1) * st1(1)) * dqusigg(1) * sigg(1)
     !Positive flux at inner boundary means loss of mass
     qafout = quant(mr) * ( ( st1(mr) * ( ( ( -3. * sigg(mr)) / sqrrad(mr) ) * dsinusq(mr) ) ) + &
          ( st2(mr) * ( sssq(mr) / ( omg(mr) * ra(mr) ) ) * dlogp(mr) ) ) 
     qdfout = (-1. * nu(mr) * st1(mr)) * dqusigg(mr) * sigg(mr)
  end select

  if (ty.gt.tinit) then
     u_r = ( ( st1 * ( ( ( ( -3. / (sigg * sqrrad) ) * dsinusq ) ) + tor_vel) ) + ( st2 * ( sssq / ( omg * ra ) ) * dlogp ) )
  else
     u_r = ( ( st1 * ( ( ( -3. / (sigg * sqrrad) ) * dsinusq ) ) ) + ( st2 * ( sssq / ( omg * ra ) ) * dlogp ) )
  endif
!  u_r = ( ( st1 * ( ( ( -3. / (sigg * sqrrad) ) * dsinusq ) ) ) + ( st2 * ( sssq / ( omg * ra ) ) * dlogp ) )

  afiq = quant * u_r!( ( st1 * ( ( ( -3. / (sigg * sqrrad) ) * dsinusq ) ) ) + ( st2 * ( sssq / ( omg * ra ) ) * dlogp ) )

  do i = 2, mr-1
     if(u_r(i).lt.0.)then
        afiq(i) = u_r(i) * quantc(i)
     else
        afiq(i) = u_r(i) * quantc(i-1)
     endif
  enddo
  afiq(1) = quantc(1) * u_r(1)
  afiq(mr) = quantc(mr-1) * u_r(mr)

!  afiq = quant * ( ( st1 * ( ( ( -3. / (sigg * sqrrad) ) * dsinusq ) ) ) + ( st2 * ( sssq / ( omg * ra ) ) * dlogp ) )

  qafin = afiq(1)
  qafout = afiq(mr-1)
  
  !At the interfaces
  !if we are at ra(1)
  
  
  dfiq = (-1. * nu * st1) * dqusigg * sigg


!  print *,'st', st(1:5)
!  print *, 'ur', u_r(1:5)
 ! read *
!  print *, ' '
!  print *, afiq
!  print *, ' '
!  print *, 'dfq', minval(dfiq)
!  print *, ' '
!  print *, dfiq
!  print *, ' '
!  read *
end subroutine dsigma_solids
