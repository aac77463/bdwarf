module constants
implicit none
  integer, parameter :: dq = selected_real_kind(15)
  real(dq), parameter :: pi = 3.14159
  real(dq), parameter :: au = 1.496E13
  real(dq), parameter :: m_earth = 5.9723657E27
  real(dq), parameter :: m_sol = 1.9884754E33
  real(dq), parameter :: epsilon = 1.E-2
  real(dq), parameter :: ism = 1.E-4
  real(dq), parameter :: mu = 2.43
  real(dq), parameter :: m_p = 1.673E-24
  real(dq), parameter :: k_b = 1.381E-16
  real(dq), parameter :: stebol = 5.67037E-5
  real(dq), parameter :: kappa = 508d0
  real(dq), parameter :: lumsol = 3.839E33
  real(dq), parameter :: g = 6.674E-8
  real(dq), parameter :: yrtos = 3.15576E7
  real(dq), parameter :: a0 = 1.e-4
end module constants

program test
  use constants
  implicit none
  real, external :: subroutine
  integer :: i, nth, fp, gp, hp, ep, mth, ip, nr, index_loc, oth, ixinrh, ixoutrh
  integer, parameter :: nr1 = 1000.,  nr2=0., nt = 10000., status=1
  integer, parameter :: dp = selected_real_kind(15)
  real(dp), dimension(nr1+nr2) :: rad, omega, x, macc_p, radc
  real(dp) :: times, timey, pl_mass, tinit, xin, xout, dxref, dxhr
  real(dp), dimension(nr1+nr2) :: sigmag, logsigmag, rthdsd, logrhop, rthdsddu, logpebs, rhspebs, vdrift
  real(dp), dimension(nr1+nr2) :: tmid, cs, h, rhop, p, nu, pebbles, tcrit, logp, adrift, afrag, adiff
  real(dp), dimension(nr1+nr2) :: acrit, agrow, stcrit, stokes0, pebint, f2, rhog
  real(dp), dimension(nr1+nr2) :: dlogqusigg, dlogqunu, dlogsinusq, dlogp, dlogshquom, d2logqusigg, d2logsinusq
  real(dp), dimension(nr1+nr2) :: d2logp, tor_vel, afip, dfip, afid, dfid, gfi
  real(dp), dimension(nr1+nr2) :: sigmagc, rhopc, pebblesc, gasv, ducg, beta, zeta, rmig
  real(dp) :: alpha, mdisk, gamma, rcrit, r_planet, m_planet, macc_d, macc_g
  real(dp) :: star_coeff, m_star, tstep, r_hill, mdotgi, mdotpi, mdotdi, mdotgo, mdotpo, mdotdo
  real(dp) :: dmin, dmout, pmin, pmout, gmin, gmout, gfin, gfout
  real(dp) :: dafin, dafout, ddfin, ddfout, pafin, pafout, pdfin, pdfout
  logical :: keep_running

  !open(11, file='parameters.txt')
  !read(11, *) alpha, tinit, r_planet
  !close(11)

  fp = 20
  gp = 30
  hp = 40
  ep = 50
  ip = 60
  !planet insertion time
  tinit = 1.e3
  oth = 500!0
  nth = 10000!0
  mth = 10000
  !alpha = 10**(-3.5)
  alpha = 1.e-4
  gamma = 0.8
  rcrit = 10.0
  !initial planet mass in [m_earth]
  m_planet = 1.e-3 !scale down mass to disk size
  !initial planet location in [au]
  r_planet = 1.
  star_coeff = 0.1
  m_star = star_coeff * m_sol
  mdisk = m_star * 0.1
  !Radial extent of grid refinement
  dxref = .5
  !Minimum radial resolution in Hill radii
  dxhr = 0.1
  !Innner and outer disk radius
  xin = 0.005  !xin=0.15 AU
  xout = 100. !xout=30 AU
!  stokes0=1.e-8

  !Opening the file that will be written to
  open(unit=ep, file='./data/bdwarf/grid.dat')
  open(unit=fp, file='./data/bdwarf/2darrs.dat')
  open(unit=gp, file='./data/bdwarf/1darrs.dat')
  open(unit=hp, file='./data/bdwarf/params.dat')

!  open(unit=ep, file='./grid.dat')
!  open(unit=fp, file='./2darrs.dat')
!  open(unit=gp, file='./1darrs.dat')
!  open(unit=hp, file='./params.dat')
  
  i = 1
  keep_running=.true.
  do while(keep_running)
     if (i.eq.1) then
        call sigma_init(mdisk, gamma, rcrit, nr1, rad, radc, sigmag, sigmagc, m_planet, &
             pl_mass, index_loc, r_planet, rhop, rhopc, pebbles, pebblesc, times, x, nr2, nr, &
             acrit, xin, xout, dxref, dxhr, m_star, &
             dmin, dmout, pmin, pmout, gmin, gmout)
        write (hp, *) alpha, mdisk, gamma, rcrit, nr, index_loc, tinit, dxref, dxhr
        close(hp)
        call calculate_quantities_bdwarf(nr, alpha, rad, radc, sigmag, sigmagc, tmid, cs, h, rhop, rhopc, p, &
             nu, omega, m_star, stcrit, vdrift, logp, pebbles, &
             f2, stokes0, rhog, pl_mass, index_loc, ducg)
        call tstep_dust_frags(nr, sigmag, rhop, pebbles, tmid, omega, tcrit, logp, cs, alpha, rad, &
             adrift, afrag, acrit, agrow, adiff, stcrit, stokes0, times, p)
        write (ep, *) rad, x, omega
        close(ep)
        call write_2d(fp, nr, sigmag, rhop, pebbles, tcrit, vdrift, tmid, logp, &
             adrift, afrag, acrit, agrow, adiff, stcrit, pebint, dlogqusigg, dlogqunu, &
             dlogsinusq, dlogp, dlogshquom, d2logqusigg, d2logsinusq, d2logp, nu, h, p, tor_vel, cs, &
             sigmagc, rhopc, pebblesc, afid, dfid, afip, dfip, stokes0, gasv, ducg)
        print *, 'planet mass', pl_mass
        write (gp, *) times, 0, pl_mass, tstep, r_hill, ixinrh, ixoutrh, mdotgi, mdotpi, mdotdi, mdotgo, &
             mdotpo, mdotdo, dmin, dmout, pmin, pmout, dafin, dafout, ddfin, ddfout, pafin, pafout, pdfin, &
             pdfout, gfin, gfout, gmin, gmout, r_planet, index_loc
     elseif (i.eq.2) then
        call calculate_quantities_bdwarf(nr, alpha, rad, radc, sigmag, sigmagc, tmid, cs, h, rhop, rhopc, p, &
             nu, omega, m_star, stcrit, vdrift, logp, pebbles, &
             f2, stokes0, rhog, pl_mass, index_loc, ducg)
        call dsigma_disk(nr, sigmag, sigmagc, nu, rad, logsigmag, rthdsd, pl_mass, h, timey, tinit, &
             m_star, gfin, gfout, gfi, r_planet, gasv)
        call dsigma_solids(nr, rad, cs, h, rhop, rhopc, dlogp, p, omega, rthdsddu, stokes0, sigmag, nu, logrhop, dlogqusigg, &
             dlogqunu, dlogsinusq, dlogshquom, d2logqusigg, d2logsinusq, d2logp, dafin, dafout, ddfin, &
             ddfout, timey, tinit, m_star, pl_mass, tor_vel, afid, dfid, r_planet)
        call dsigma_solids(nr, rad, cs, h, pebbles, pebblesc, dlogp, p, omega, rhspebs, stcrit, sigmag, nu, logpebs, dlogqusigg, &
             dlogqunu, dlogsinusq, dlogshquom, d2logqusigg, d2logsinusq, d2logp, pafin, pafout, pdfin, &
             pdfout, timey, tinit, m_star, pl_mass, tor_vel, afip, dfip, r_planet)
        call timestep(nr, sigmag, sigmagc, logsigmag, times, timey, tstep, pl_mass, &
             rhop, rhopc, logrhop, logpebs, pebbles, pebblesc, tcrit, tinit, &
             dafin, dafout, ddfin, ddfout, pafin, pafout, pdfin, pdfout, rad, radc, dmin, dmout, pmin, &
             pmout, vdrift, gfin, gfout, gmin, gmout, gfi, afip, dfip, afid, dfid, nu, h, m_star, &
             r_planet)
       elseif (timey.gt.tinit) then
        call calculate_quantities_bdwarf(nr, alpha, rad, radc, sigmag, sigmagc, tmid, cs, h, rhop, rhopc, p, &
             nu, omega, m_star, stcrit, vdrift, logp, pebbles, &
             f2, stokes0, rhog, pl_mass, index_loc, ducg)
        call tstep_dust_frags(nr, sigmag, rhop, pebbles, tmid, omega, tcrit, logp, cs, alpha, rad, &
             adrift, afrag, acrit, agrow, adiff, stcrit, stokes0, times, p)
        call remove_dust(nr, rad, radc, rhop, rhopc, pebbles, pebblesc, tcrit, times, pebint)
        call dsigma_disk(nr, sigmag, sigmagc, nu, rad, logsigmag, rthdsd, pl_mass, h, timey, tinit, &
             m_star, gfin, gfout, gfi, r_planet, gasv)
        call dsigma_solids(nr, rad, cs, h, rhop, rhopc, dlogp, p, omega, rthdsddu, stokes0, sigmag, nu, logrhop, dlogqusigg, &
             dlogqunu, dlogsinusq, dlogshquom, d2logqusigg, d2logsinusq, d2logp, dafin, dafout, ddfin, &
             ddfout, timey, tinit, m_star, pl_mass, tor_vel, afid, dfid, r_planet)
        call dsigma_solids(nr, rad, cs, h, pebbles, pebblesc, dlogp, p, omega, rhspebs, stcrit, sigmag, nu, logpebs, dlogqusigg, &
             dlogqunu, dlogsinusq, dlogshquom, d2logqusigg, d2logsinusq, d2logp, pafin, pafout, pdfin, &
             pdfout, timey, tinit, m_star, pl_mass, tor_vel, afip, dfip, r_planet)
        call timestep(nr, sigmag, sigmagc, logsigmag, times, timey, tstep, pl_mass, &
             rhop, rhopc, logrhop, logpebs, pebbles, pebblesc, tcrit, tinit, &
             dafin, dafout, ddfin, ddfout, pafin, pafout, pdfin, pdfout, rad, radc, dmin, dmout, pmin, &
             pmout, vdrift, gfin, gfout, gmin, gmout, gfi, afip, dfip, afid, dfid, nu, h, m_star, &
             r_planet)
        call macc_rates(nr, rad, pl_mass, r_planet, omega, macc_g, macc_d, macc_p, &
             r_hill, pebbles, stcrit, m_star, ixinrh, ixoutrh)
        call remove(nr, pebbles, rad, r_planet, macc_p, tstep, pl_mass, m_star)
        call migration(nr, rad, sigmag, tmid, h, omega, m_star, beta, zeta, pl_mass, r_planet, rmig, tstep, tinit, timey, index_loc)
     else
        call calculate_quantities_bdwarf(nr, alpha, rad, radc, sigmag, sigmagc, tmid, cs, h, rhop, rhopc, p, &
             nu, omega, m_star, stcrit, vdrift, logp, pebbles, &
             f2, stokes0, rhog, pl_mass, index_loc, ducg)
        call tstep_dust_frags(nr, sigmag, rhop, pebbles, tmid, omega, tcrit, logp, cs, alpha, rad, &
             adrift, afrag, acrit, agrow, adiff, stcrit, stokes0, times, p)
        call remove_dust(nr, rad, radc, rhop, rhopc, pebbles, pebblesc, tcrit, times, pebint)
        call dsigma_disk(nr, sigmag, sigmagc, nu, rad, logsigmag, rthdsd, pl_mass, h, timey, tinit, &
             m_star, gfin, gfout, gfi, r_planet, gasv)
        call dsigma_solids(nr, rad, cs, h, rhop, rhopc, dlogp, p, omega, rthdsddu, stokes0, sigmag, nu, logrhop, dlogqusigg, &
             dlogqunu, dlogsinusq, dlogshquom, d2logqusigg, d2logsinusq, d2logp, dafin, dafout, ddfin, &
             ddfout, timey, tinit, m_star, pl_mass, tor_vel, afid, dfid, r_planet)
        call dsigma_solids(nr, rad, cs, h, pebbles, pebblesc, dlogp, p, omega, rhspebs, stcrit, sigmag, nu, logpebs, dlogqusigg, &
             dlogqunu, dlogsinusq, dlogshquom, d2logqusigg, d2logsinusq, d2logp, pafin, pafout, pdfin, &
             pdfout, timey, tinit, m_star, pl_mass, tor_vel, afip, dfip, r_planet)
        call timestep(nr, sigmag, sigmagc, logsigmag, times, timey, tstep, pl_mass, &
             rhop, rhopc, logrhop, logpebs, pebbles, pebblesc, tcrit, tinit, &
             dafin, dafout, ddfin, ddfout, pafin, pafout, pdfin, pdfout, rad, radc, dmin, dmout, pmin, &
             pmout, vdrift, gfin, gfout, gmin, gmout, gfi, afip, dfip, afid, dfid, nu, h, m_star, &
             r_planet)
     endif
     if((times/yrtos).lt.5.e4)then
        if (mod(i, oth)==0) then
           call write_2d(fp, nr, sigmag, rhop, pebbles, tcrit, vdrift, tmid, logp, &
                adrift, afrag, acrit, agrow, adiff, stcrit, pebint, dlogqusigg, dlogqunu, &
                dlogsinusq, dlogp, dlogshquom, d2logqusigg, d2logsinusq, d2logp, nu, h, p, tor_vel, cs, &
                sigmagc, rhopc, pebblesc, afid, dfid, afip, dfip, stokes0, gasv, ducg)
           write(gp, *) times, timey, pl_mass, tstep, r_hill, ixinrh, ixoutrh, mdotgi, mdotpi, mdotdi, &
                mdotgo, mdotpo, mdotdo, dmin, dmout, pmin, pmout, dafin, dafout, ddfin, ddfout, pafin, pafout, &
                pdfin, pdfout, gfin, gfout, gmin, gmout, r_planet, index_loc
        endif
     else
        if (mod(i, nth)==0) then
!           print *, macc_p
           call write_2d(fp, nr, sigmag, rhop, pebbles, tcrit, vdrift, tmid, logp, &
                adrift, afrag, acrit, agrow, adiff, stcrit, pebint, dlogqusigg, dlogqunu, &
                dlogsinusq, dlogp, dlogshquom, d2logqusigg, d2logsinusq, d2logp, nu, h, p, tor_vel, cs, &
                sigmagc, rhopc, pebblesc, afid, dfid, afip, dfip, stokes0, gasv, ducg)
           write(gp, *) times, timey, pl_mass, tstep, r_hill, ixinrh, ixoutrh, mdotgi, mdotpi, mdotdi, &
                mdotgo, mdotpo, mdotdo, dmin, dmout, pmin, pmout, dafin, dafout, ddfin, ddfout, pafin, pafout, &
                pdfin, pdfout, gfin, gfout, gmin, gmout, r_planet, index_loc
        endif
     endif
     if (mod(i, mth)==0) then
        print *, 'it', i, 'ty', timey, 'pl', pl_mass/m_earth
        print *, 'dt', tstep/yrtos, r_hill/au, ixinrh, ixoutrh
        print *, 'pebs max', maxval(pebbles)
        print *, 'pdfin, pafin', pdfin, pafin
        print *, 'pdfout, pafout', pdfout, pafout
        print *, 'planet location', r_planet
        flush(gp)
        flush(fp)
     endif
     if ((times/yrtos).gt.3.e6) then
        keep_running=.false.
        print *, 'model exited with message: time in years has reached disk dissipation'
        call EXIT(STATUS)
     endif
     i = i + 1
  enddo

  close(hp)
end program test
