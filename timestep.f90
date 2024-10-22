subroutine timestep(mr, sigg, siggc, logsigg, ts, ty, dtmin, m_pl, du, duc, logdu, &
     logpebs, pebs, pebsc, tcrit, tinit, dafin, dafout, ddfin, ddfout, pafin, pafout, &
     pdfin, pdfout, ra, rac, dmin, dmout, pmin, pmout, vd, gfin, gfout, gmin, gmout, gfi, afip, dfip, afid, &
     dfid, vs, sh, m_star, rpl)

!!$  This subroutine calculates the minimum timestep that can be taken
!!$  away from the gas surface density, assuming that we only take the gas
!!$  surface density into account.
  use constants
  implicit none
  integer :: i, mr
  integer, parameter :: dp = selected_real_kind(15)
  real(dp), dimension(mr), intent(inout) :: logsigg, logdu, logpebs
  real(dp), dimension(mr), intent(inout) :: sigg, du, pebs, siggc, duc, pebsc
  real(dp), dimension(mr), intent(in) :: tcrit, ra, vd, rac, gfi, afip, dfip, afid, dfid, vs, sh
  real(dp), intent(out) :: ts, ty
  real(dp) :: thresh, dtdisk, m_pl, th_plan, ddmin, dgmin, dtdust, dtpebs, dpmin, tinit
  real(dp) :: threshpeb, dtpebdu, tcmin, threshdust, threshra, m_star, rpl
  real(dp), intent(in) :: dafin, dafout, ddfin, ddfout, pafin, pafout, pdfin, pdfout, gfin, gfout
  real(dp), intent(out) :: dtmin, dmin, dmout, pmin, pmout, gmin, gmout
  real(dp), dimension(mr) :: dr, dsigg, dpebs, ddust, tor, dgasv, gasv

  !Setting the timestep quantities to large number so that we
  !can replace them with other minimums as we iterate
  dgmin = 10E30
  ddmin = 10e30
  dpmin = 10e30
  dtmin = 10e30
  dtpebs = 10e33

  !Percentages (in log) by which all of the quantities can change
  thresh = 5.e-2
  threshpeb = 5.e-2
  threshdust = 5.e-2
  th_plan = 1.e-2

  !A second possible timestep limit
  !divide the dr/vd
  do i = 1,mr-1
     dr(i) = ra(i+1) - ra(i)
  enddo
  dr(mr) = dr(mr-1)
  !Make sure that the timestep limit is the minimum
  !of threshpeb or dr/vd
  threshra = minval(abs(dr/vd)) * 0.2


  !Finding minimum timestep for the critical time of pebble growth
  !Chosen to be half of the minimum value of the critical time
  !Find the smallest value of the critical time
  tcmin = minval(tcrit)
  dtpebdu = -1. * tcmin * log(0.5)

!!! Simple prescription for Forward Center Time Space !!!
  do i = 1, mr-1
     dsigg(i) = 2. * ((gfi(i)*ra(i)) - (gfi(i+1)*ra(i+1))) / (ra(i+1)**2. - ra(i)**2.)
     dpebs(i) = 2. * (((afip(i) + dfip(i))*ra(i)) - ((afip(i+1) + dfip(i+1))*ra(i+1)) )/ (ra(i+1)**2. - ra(i)**2.)
     ddust(i) = 2. * (((afid(i) + dfid(i))*ra(i)) - ((afid(i+1) + dfid(i+1))*ra(i+1)) )/ (ra(i+1)**2. - ra(i)**2.)
  enddo
!  print *, 'ddust', ddust
!  print *, 'dpebs', minval(dpebs), minloc(dpebs)
!  print *, 'afip', minval(afip)
!  print *, 'dfip', minval(dfip)
!  print *, 'afid', minval(afid)
!  print *, 'dfid', minval(dfid)

  !Finding minimum of the gas
  do i = 1, mr-1
     dtdisk = thresh * siggc(i) / abs(dsigg(i))
     if ((dtdisk.lt.dgmin)) then
        dgmin = dtdisk
     endif
  enddo

  !Finding minimum of the dust
  do i = 1, mr-1
     dtdust = threshdust * duc(i) / abs(ddust(i))
     if ((dtdust.lt.ddmin)) then
        ddmin = dtdust
     endif
  enddo

  !Finding minimum of the pebbles
  do i = 1, mr-1
     dtpebs = threshpeb * pebsc(i) / abs(dpebs(i))
     if ((dtpebs.lt.dpmin)) then
        dpmin = dtpebs
     endif
  enddo

  if (ty.lt.tinit) then
     dtmin = min(dpmin, dgmin, ddmin)
  else
     dtmin = min(dpmin, dgmin, ddmin)
  endif

  !Computing a CFL timestep condition
  !Calculate first the movement of the gas surface density
  call torque_vel_calc(mr, ra, sh, m_star, m_pl, tor, rpl)
  call cdiff_poly(mr, sigg * vs * sqrt(ra), ra, dgasv)
  gasv = ( -3. / (sigg * sqrt(ra))) * dgasv
  if (ty.lt.tinit) then
     gasv = gasv
  else
     gasv = gasv + tor
  endif
  if (dtmin.lt.1.e-2*yrtos) then
     dtmin = (0.5*minval(dr))/(maxval(abs(gasv)))
  endif

  
  siggc = siggc + (dsigg * dtmin)
  duc = duc + (ddust * dtmin)
  pebsc = pebsc + (dpebs * dtmin)

  !Updating the grid
  call interp1d(mr-1, rac, siggc, ra(2:mr-1), sigg(2:mr-1))
  sigg(1) = siggc(1)
  sigg(mr) = siggc(mr-1)
  call interp1d(mr-1, rac, duc, ra(2:mr-1), du(2:mr-1))
  du(1) = duc(1)
  du(mr) = duc(mr-1)

  call interp1d(mr-1, rac, pebsc, ra(2:mr-1), pebs(2:mr-1))
  pebs(1) = pebsc(1)
  pebs(mr) = pebsc(mr-1)

  !Calculate mass lost through inner boundary!!!!!
     pmin = pmin + ( pafin + pdfin ) * dtmin * ra(1) * 2. * pi / m_earth
     dmin = dmin + ( dafin + ddfin ) * dtmin * 2. * pi * ra(1) / m_earth
     gmin = gmin + ((gfin * dtmin * 2. * pi * ra(1)) / m_earth)

     if((pafout+pdfout).gt.0.)then
        pmout = pmout + ( pafout + pdfout ) * dtmin*2.*pi*ra(mr)/m_earth
     else
        pmout = pmout
     endif

     if((dafout+ddfout).gt.0.)then
        dmout = dmout + ( dafout + ddfout ) * dtmin*2.*pi*ra(mr)/m_earth
     else
        dmout = dmout
     endif

     if((gfout).gt.0.)then
        gmout = gmout + ( gfout * dtmin*2.*pi*ra(mr)/m_earth)
     else
        gmout = gmout
     endif

  !Advance time in seconds and years
  ts = ts + dtmin
  ty = ty + (dtmin/yrtos)

end subroutine timestep
