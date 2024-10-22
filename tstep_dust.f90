subroutine tstep_dust(mr, sigg, du, pebs, omg, tcrit, lnpr, cs, al, ra, adrift, afrag, &
     acrit, agrow, adiff, stokes_crit, stokes0, ts)
  use constants
  implicit none
  integer :: i, mr
  integer, parameter :: dp = selected_real_kind(15)
  real(dp), dimension(mr), intent(in) :: lnpr, cs, ra, sigg, du, omg, pebs
  real(dp), dimension(mr) :: tau_grow
  real(dp), dimension(mr), intent(out) :: tcrit, afrag, adrift, acrit, agrow, adiff, stokes0, stokes_crit
  real(dp), intent(in) :: ts
  real(dp) :: rhos, vfrag, al

  !Fragmentation velocity in [cm/s]
  vfrag = 200.

  !The material density of the material
  rhos = 1.6 ![g/cm^3]

  !Calculate the growth time of the particles from d2g
  tau_grow = sigg / ( abs(du+pebs) * omg )  

  !The different growth limits, see Birnstiel+12
  adrift = (2. * (du + pebs) * (ra * omg)**2.) / (pi * rhos * cs**2.) * abs(lnpr)**(-1.) * .55
  afrag = (2. * sigg * vfrag**2.) / (3. * pi * rhos * al * cs**2.) * .37
  adiff = (4. * sigg * vfrag * omg * ra) / (pi * rhos * cs**2.) * abs(lnpr)**(-1.) * .37
  agrow = a0 * exp( ts / tau_grow )
  
  !Assuring sizes are never smaller than  micron size
  do i=1,mr
     if(.not.afrag(i).gt.a0)then
        afrag(i)=a0
     endif
     if(.not.adrift(i).gt.a0)then
        adrift(i)=a0
     endif
     if(.not.adiff(i).gt.a0)then
        adiff(i)=a0
     endif
     if(.not.agrow(i).gt.a0)then
        agrow(i)=a0
     endif
  enddo

  !Taking the smallest of each growth limit
  do i=1,mr
     acrit(i) = min(min(adrift(i),afrag(i),adiff(i)), agrow(i))
  enddo

  !Calculating the stokes number for the critical size and 0 for the dust
  stokes_crit = (pi * acrit * rhos) / (2. * sigg)
  stokes0 = (pi * a0 * rhos) / (2. * sigg)

  !Calculate the critical time to grow to acrit
  do i = 1, mr
     tcrit(i) = log(acrit(i) / a0) * tau_grow(i)
     if(.not.tcrit(i).ge.0.)then
     endif
  enddo

end subroutine tstep_dust

subroutine tstep_dust_frags(mr, sigg, du, pebs, tm, omg, tcrit, lnpr, cs, al, ra, adrift, afrag, &
     acrit, agrow, adiff, stokes_crit, stokes0, ts, p)
  use constants
  implicit none
  integer :: i, mr
  integer, parameter :: dp = selected_real_kind(15)
  real(dp), dimension(mr), intent(in) :: lnpr, cs, ra, sigg, du, omg, pebs, tm, p
  real(dp), dimension(mr) :: tau_grow, rhos, vfrag, dpr
  real(dp), dimension(mr), intent(out) :: tcrit, afrag, adrift, acrit, agrow, adiff, stokes0, stokes_crit
  real(dp), intent(in) :: ts
  real(dp) :: al

  !Calculating the dP factor for the azimuthal fragmentation velocity
  call cdiff_poly(mr, p, ra, dpr)
  
  !Fragmentation velocity in [cm/s]
  vfrag = 100.

  !The material density of the material
  rhos = 1.6 ![g/cm^3]

  do i = 1, mr
     if(tm(i).lt.273.)then
        vfrag(i) = 100.
        rhos(i) = 1.6
     endif
  enddo

  !Calculate the growth time of the particles from d2g
  tau_grow = sigg / ( abs(du+pebs) * omg )  

  !The different growth limits, see Birnstiel+12
  adrift = (2. * (du + pebs) * (ra * omg)**2.) / (pi * rhos * cs**2.) * abs(lnpr)**(-1.) * .55
  !Fragmentation velocity due to Brownian motion
  afrag = (2. * sigg * vfrag**2.) / (3. * pi * rhos * al * cs**2.) * .37
  adiff = (4. * sigg * vfrag * omg * ra) / (pi * rhos * cs**2.) * abs(lnpr)**(-1.) * .37
  agrow = a0 * exp( ts / tau_grow )
  
  !Assuring sizes are never smaller than  micron size
  do i=1,mr
     if(.not.afrag(i).gt.a0)then
        afrag(i)=a0
     endif
     if(.not.adrift(i).gt.a0)then
        adrift(i)=a0
     endif
     if(.not.adiff(i).gt.a0)then
        adiff(i)=a0
     endif
     if(.not.agrow(i).gt.a0)then
        agrow(i)=a0
     endif
  enddo

  !Taking the smallest of each growth limit
  do i=1,mr
     acrit(i) = min(min(adrift(i),afrag(i),adiff(i)), agrow(i))
  enddo

  !Fragmentation velocity due to azimuthal mixing
!  afrag = abs(-1. * (dp/(rhos*omg)) * (1./(1.+stc**2.)))

  !Calculating the stokes number for the critical size and 0 for the dust
  stokes_crit = (pi * acrit * rhos) / (2. * sigg)
  stokes0 = (pi * a0 * rhos) / (2. * sigg)

  !Calculate the critical time to grow to acrit
  do i = 1, mr
     tcrit(i) = log(acrit(i) / a0) * tau_grow(i)
     if(.not.tcrit(i).ge.0.)then
     endif
  enddo

end subroutine tstep_dust_frags
