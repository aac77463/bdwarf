subroutine calculate_quantities(mr3, al, ra, rac, sigg, siggc, tmc, ss, sh, du, duc, pr, vs, omg, m_star, &
     stokes, vdrift, lnpr, pebs, f2, stokes0, rhog, m_pl, idx_loc, ducg)

!!$  This subroutine calculates all of the quantities that we need from the
!!$  gas surface density. These include, midplane temperature, sound speed,
!!$  scale height, dust surface density, pressure, kinematic viscosity, and
!!$  keplerian frequency. All quantities are calculated in [cgs] units. The
!!$  temperature equation includes kinematic viscosity and stellar irradiation.
!!$  Also includes the mass accretion of dust and pebbles onto the star.
  use constants
  implicit none
  integer, parameter :: dp = selected_real_kind(15)
  integer :: ir, mr3, i, idx_loc
  real(dp), dimension(mr3) :: vhw, logp, logra, f, ti, tv, dlogp, vdd, tm, f2c
  real(dp), dimension(mr3), intent(in) :: ra, stokes, pebs, stokes0, rac
  real(dp), dimension(mr3), intent(inout) :: sigg, du, siggc, duc
  real(dp), dimension(mr3), intent(out) :: ss, sh, pr, vs, vdrift, lnpr, omg, f2, rhog, tmc, ducg
  real(dp) :: gran, lumst, m_star, tcap, al, mdotgi, mdotpi, mdotdi, mdotgo, mdotpo, mdotdo, q, m_pl, tot
  !Setting the sublimation temperature
  tcap = 1500.

  !Calculating the GAS to DUST ratio for the temperature
  f = sigg / du

  gran = 0.05
  lumst = lumsol

  !Keplerian frequency
  omg = (((g * m_star) / ra **3.))**(1./2.)

  !Viscous temperature
!  print*, ((27. * sigg * du * kappa * al * k_b * omg))
!  print *,  ((128. * stebol * mu * m_p))**(1./3.)
  tv = ((27. * sigg * du * kappa * al * k_b * omg) / (128. * stebol * mu * m_p))**(1./3.)

  !Define quantities here to cell centers (sigg, du, omg)

  !Irradiation Temperature
  ti = ((gran * lumst) / (4. * pi * stebol * ra**2.))**(1./4.)
!  print *, ti
  !Miplane temperature
  tm = ((tv)**4. +  (ti)**4.)**(1./4.)
  !tm = ti
  !Interpolating the midplane temperature at the interfaces to cell centers
  call interp1d(mr3, ra, tm, rac, tmc)
 
  !Calculating the amount of dust converted to gas due to sublimation
  f2 =  (27. * sigg**2. * kappa * al * k_b * omg) /  (128. * stebol * mu * m_p ) /  ((tcap**4. - ti**4.)**(3./4.))
  !Interpolating the amount of dust conversion factor to the cell centers
  call interp1d(mr3, ra, f2, rac, f2c)
  
!!!!!!!Check total mass of dust lost to dust sublimation
!!!!! do everything on centers
  !This would only happens from 1,mr3-1
  !  print *, 'before', duc(1:5)
  
  do ir=1,mr3-1
     if (tmc(ir).gt.tcap) then
        ducg(ir) = siggc(ir) / f2c(ir)
        tot = duc(ir) + siggc(ir)
        duc(ir) =  siggc(ir)/f2c(ir)
        siggc(ir) = tot - duc(ir)
        tmc(ir) = tcap
     endif
  enddo
  !!!Output the total amount of dust lost to sublimation!!!!
  
  !sound speed, scale height, viscosity, midplane density, and midplane pressure
  ss = ((k_b * tm) / (mu * m_p))**(1./2.)
  sh = ss / omg
  vs = al * ss * sh
  rhog = sigg / (sqrt(2.*pi)*sh)
  pr = rhog * ss**2.

  !Radial gas pressure gradient
  logra = log(ra)
  logp = log(pr)
  call cdiff_poly(mr3, logp, logra, dlogp)
  lnpr = dlogp

  !Calculating headwid velocity and drift velocity
  vhw = lnpr * (ss**2. / (omg * ra))
  vdrift = vhw * (stokes / (1. + stokes**2.))
  vdd = vhw * (stokes0 / (1. + stokes0**2.))
  !Mass lost through the inner boundary
  ![g/cm**2] * [cm/s] * [cm]
  mdotgi = -2. * pi * sigg(1) * ra(1) * vs(1)
  mdotpi = -2. * pi * pebs(1) *  ra(1) * vdrift(1)
  mdotdi = -2. * pi * du(1) *  ra(1) * vdd(1)
  mdotgo = -2. * pi * sigg(mr3) * ra(mr3) * vs(mr3)
  mdotpo = -2. * pi * pebs(mr3) *  ra(mr3) * vdrift(mr3)
  mdotdo = -2. * pi * du(mr3) *  ra(mr3) * vdd(mr3)

end subroutine calculate_quantities

subroutine calculate_quantities_stand(mr3, al, ra, rac, sigg, siggc, tmc, ss, sh, du, duc, pr, vs, omg, m_star, &
     stokes, vdrift, lnpr, pebs, f2, stokes0, rhog, m_pl, idx_loc, ducg)

!!$  This subroutine calculates all of the quantities that we need from the
!!$  gas surface density. These include, midplane temperature, sound speed,
!!$  scale height, dust surface density, pressure, kinematic viscosity, and
!!$  keplerian frequency. All quantities are calculated in [cgs] units. The
!!$  temperature equation includes kinematic viscosity and stellar irradiation.
!!$  Also includes the mass accretion of dust and pebbles onto the star.
  use constants
  implicit none
  integer, parameter :: dp = selected_real_kind(15)
  integer :: ir, mr3, i, idx_loc
  real(dp), dimension(mr3) :: vhw, logp, logra, f, ti, tv, dlogp, vdd, tm, f2c
  real(dp), dimension(mr3), intent(in) :: ra, stokes, pebs, stokes0, rac
  real(dp), dimension(mr3), intent(inout) :: sigg, du, siggc, duc
  real(dp), dimension(mr3), intent(out) :: ss, sh, pr, vs, vdrift, lnpr, omg, f2, rhog, tmc, ducg
  real(dp) :: gran, lumst, m_star, tcap, al, mdotgi, mdotpi, mdotdi, mdotgo, mdotpo, mdotdo, q, m_pl, tot
  !Setting the sublimation temperature
  tcap = 1500.

  !Calculating the GAS to DUST ratio for the temperature
  f = sigg / du

  gran = 0.05
  lumst = lumsol

  !Keplerian frequency
  omg = (((g * m_star) / ra **3.))**(1./2.)

  !Viscous temperature
  tv = ((27. * sigg * du * kappa * al * k_b * omg) / (128. * stebol * mu * m_p))**(1./3.)
  !Define quantities here to cell centers (sigg, du, omg)

  !Irradiation Temperature
  ti = ((gran * lumst) / (4. * pi * stebol * ra**2.))**(1./4.)
  !Miplane temperature
  tm = ((tv)**4. +  (ti)**4.)**(1./4.)

  !Interpolating the midplane temperature at the interfaces to cell centers
  call interp1d(mr3, ra, tm, rac, tmc)
 
  !Calculating the amount of dust converted to gas due to sublimation
  f2 =  (27. * sigg**2. * kappa * al * k_b * omg) /  (128. * stebol * mu * m_p ) /  ((tcap**4. - ti**4.)**(3./4.))
  !Interpolating the amount of dust conversion factor to the cell centers
  call interp1d(mr3, ra, f2, rac, f2c)
  
  !sound speed, scale height, viscosity, midplane density, and midplane pressure
  ss = ((k_b * tm) / (mu * m_p))**(1./2.)
  sh = ss / omg
  vs = al * ss * sh
  rhog = sigg / (sqrt(2.*pi)*sh)
  pr = rhog * ss**2.

  !Radial gas pressure gradient
  logra = log(ra)
  logp = log(pr)
  call cdiff_poly(mr3, logp, logra, dlogp)
  lnpr = dlogp

  !Calculating headwid velocity and drift velocity
  vhw = lnpr * (ss**2. / (omg * ra))
  vdrift = vhw * (stokes / (1. + stokes**2.))
  vdd = vhw * (stokes0 / (1. + stokes0**2.))
  !Mass lost through the inner boundary
  ![g/cm**2] * [cm/s] * [cm]
  mdotgi = -2. * pi * sigg(1) * ra(1) * vs(1)
  mdotpi = -2. * pi * pebs(1) *  ra(1) * vdrift(1)
  mdotdi = -2. * pi * du(1) *  ra(1) * vdd(1)
  mdotgo = -2. * pi * sigg(mr3) * ra(mr3) * vs(mr3)
  mdotpo = -2. * pi * pebs(mr3) *  ra(mr3) * vdrift(mr3)
  mdotdo = -2. * pi * du(mr3) *  ra(mr3) * vdd(mr3)

end subroutine calculate_quantities_stand

subroutine calculate_quantities_const(mr3, al, ra, rac, sigg, siggc, tmc, ss, sh, du, duc, pr, vs, omg, m_star, &
     stokes, vdrift, lnpr, pebs, f2, stokes0, rhog, m_pl, idx_loc, ducg)

!!$  This subroutine calculates all of the quantities that we need from the
!!$  gas surface density. These include, midplane temperature, sound speed,
!!$  scale height, dust surface density, pressure, kinematic viscosity, and
!!$  keplerian frequency. All quantities are calculated in [cgs] units. The
!!$  temperature equation includes kinematic viscosity and stellar irradiation.
!!$  Also includes the mass accretion of dust and pebbles onto the star.
  use constants
  implicit none
  integer, parameter :: dp = selected_real_kind(15)
  integer :: ir, mr3, i, idx_loc
  real(dp), dimension(mr3) :: vhw, logp, logra, f, ti, tv, dlogp, vdd, tm, f2c
  real(dp), dimension(mr3), intent(in) :: ra, stokes, pebs, stokes0, rac
  real(dp), dimension(mr3), intent(inout) :: sigg, du, siggc, duc
  real(dp), dimension(mr3), intent(out) :: ss, sh, pr, vs, vdrift, lnpr, omg, f2, rhog, tmc, ducg
  real(dp) :: gran, lumst, m_star, tcap, al, mdotgi, mdotpi, mdotdi, mdotgo, mdotpo, mdotdo, q, m_pl, tot
  !Setting the sublimation temperature
  tcap = 1500.

  !Calculating the GAS to DUST ratio for the temperature
  f = sigg / du

  gran = 0.05
  lumst = lumsol

  !Keplerian frequency
  omg = (((g * m_star) / ra **3.))**(1./2.)

  !Viscous temperature
  tv = ((27. * sigg * du * kappa * al * k_b * omg) / (128. * stebol * mu * m_p))**(1./3.)
  !Define quantities here to cell centers (sigg, du, omg)

  !Irradiation Temperature
  ti = ((gran * lumst) / (4. * pi * stebol * ra**2.))**(1./4.)
  !Miplane temperature
  tm = ti

  !Interpolating the midplane temperature at the interfaces to cell centers
  call interp1d(mr3, ra, tm, rac, tmc)
 
  !Calculating the amount of dust converted to gas due to sublimation
  f2 =  (27. * sigg**2. * kappa * al * k_b * omg) /  (128. * stebol * mu * m_p ) /  ((tcap**4. - ti**4.)**(3./4.))
  !Interpolating the amount of dust conversion factor to the cell centers
  call interp1d(mr3, ra, f2, rac, f2c)
  
  !sound speed, scale height, viscosity, midplane density, and midplane pressure
  ss = ((k_b * tm) / (mu * m_p))**(1./2.)
  sh = ss / omg
  vs = al * ss * sh
  rhog = sigg / (sqrt(2.*pi)*sh)
  pr = rhog * ss**2.

  !Radial gas pressure gradient
  logra = log(ra)
  logp = log(pr)
  call cdiff_poly(mr3, logp, logra, dlogp)
  lnpr = dlogp

  !Calculating headwid velocity and drift velocity
  vhw = lnpr * (ss**2. / (omg * ra))
  vdrift = vhw * (stokes / (1. + stokes**2.))
  vdd = vhw * (stokes0 / (1. + stokes0**2.))
  !Mass lost through the inner boundary
  ![g/cm**2] * [cm/s] * [cm]
  mdotgi = -2. * pi * sigg(1) * ra(1) * vs(1)
  mdotpi = -2. * pi * pebs(1) *  ra(1) * vdrift(1)
  mdotdi = -2. * pi * du(1) *  ra(1) * vdd(1)
  mdotgo = -2. * pi * sigg(mr3) * ra(mr3) * vs(mr3)
  mdotpo = -2. * pi * pebs(mr3) *  ra(mr3) * vdrift(mr3)
  mdotdo = -2. * pi * du(mr3) *  ra(mr3) * vdd(mr3)

end subroutine calculate_quantities_const

subroutine calculate_quantities_bdwarf(mr3, al, ra, rac, sigg, siggc, tmc, ss, sh, du, duc, pr, vs, omg, m_star, &
     stokes, vdrift, lnpr, pebs, f2, stokes0, rhog, m_pl, idx_loc, ducg)

!!$  This subroutine calculates all of the quantities that we need from the
!!$  gas surface density. These include, midplane temperature, sound speed,
!!$  scale height, dust surface density, pressure, kinematic viscosity, and
!!$  keplerian frequency. All quantities are calculated in [cgs] units. The
!!$  temperature equation includes kinematic viscosity and stellar irradiation.
!!$  Also includes the mass accretion of dust and pebbles onto the star.
  use constants
  implicit none
  integer, parameter :: dp = selected_real_kind(15)
  integer :: ir, mr3, i, idx_loc
  real(dp), dimension(mr3) :: vhw, logp, logra, f, ti, tv, dlogp, vdd, tm, f2c
  real(dp), dimension(mr3), intent(in) :: ra, stokes, pebs, stokes0, rac
  real(dp), dimension(mr3), intent(inout) :: sigg, du, siggc, duc
  real(dp), dimension(mr3), intent(out) :: ss, sh, pr, vs, vdrift, lnpr, omg, f2, rhog, tmc, ducg
  real(dp) :: gran, lumst, m_star, tcap, al, mdotgi, mdotpi, mdotdi, mdotgo, mdotpo, mdotdo, q, m_pl, tot
  !Setting the sublimation temperature
  tcap = 1500.

  !Calculating the GAS to DUST ratio for the temperature
  f = sigg / du

  gran = 0.05
  lumst = lumsol * 0.03

  !Keplerian frequency
  omg = (((g * m_star) / ra **3.))**(1./2.)

  !Viscous temperature
  tv = ((27. * sigg * du * kappa * al * k_b * omg) / (128. * stebol * mu * m_p))**(1./3.)
  !Define quantities here to cell centers (sigg, du, omg)

  !Irradiation Temperature
  ti = ((gran * lumst) / (4. * pi * stebol * ra**2.))**(1./4.)
  !Miplane temperature
  tm = ti

  !Interpolating the midplane temperature at the interfaces to cell centers
  call interp1d(mr3, ra, tm, rac, tmc)
 
  !Calculating the amount of dust converted to gas due to sublimation
  f2 =  (27. * sigg**2. * kappa * al * k_b * omg) /  (128. * stebol * mu * m_p ) /  ((tcap**4. - ti**4.)**(3./4.))
  !Interpolating the amount of dust conversion factor to the cell centers
  call interp1d(mr3, ra, f2, rac, f2c)
  
  !sound speed, scale height, viscosity, midplane density, and midplane pressure
  ss = ((k_b * tm) / (mu * m_p))**(1./2.)
  sh = ss / omg
  vs = al * ss * sh
  rhog = sigg / (sqrt(2.*pi)*sh)
  pr = rhog * ss**2.

  !Radial gas pressure gradient
  logra = log(ra)
  logp = log(pr)
  call cdiff_poly(mr3, logp, logra, dlogp)
  lnpr = dlogp

  !Calculating headwid velocity and drift velocity
  vhw = lnpr * (ss**2. / (omg * ra))
  vdrift = vhw * (stokes / (1. + stokes**2.))
  vdd = vhw * (stokes0 / (1. + stokes0**2.))
  !Mass lost through the inner boundary
  ![g/cm**2] * [cm/s] * [cm]
  mdotgi = -2. * pi * sigg(1) * ra(1) * vs(1)
  mdotpi = -2. * pi * pebs(1) *  ra(1) * vdrift(1)
  mdotdi = -2. * pi * du(1) *  ra(1) * vdd(1)
  mdotgo = -2. * pi * sigg(mr3) * ra(mr3) * vs(mr3)
  mdotpo = -2. * pi * pebs(mr3) *  ra(mr3) * vdrift(mr3)
  mdotdo = -2. * pi * du(mr3) *  ra(mr3) * vdd(mr3)

end subroutine calculate_quantities_bdwarf
