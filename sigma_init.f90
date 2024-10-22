subroutine sigma_init(md, ga, rc, mr, ra, rac, sigg, siggc, mplinit, m_pl, &
     idx_loc, rpl, du, duc, pebs, pebsc, ts, x, mr2, mr3, &
     acrit, xin, xout, dxref, dxhr, m_star, &
     dmin, dmout, pmin, pmout, gmin, gmout)
!!$  This subroutine calculates the first iteration of the gas surface density.
!!$  The inputs are the disk mass (md) in [g]. The gamma (ga) parameter in
!!$  a real format. The critical radius (rc) in [au]. It also takes in the
!!$  radial (mr) dimension of the grid. Lastly, the subroutine outputs the
!!$  radial (ra) grid in [cm] and gas surface density (sigg) in [g/cm^2]
  use constants
  implicit none

  integer :: mr, mr2, mr3, i
  integer, intent(out) :: idx_loc
  integer, parameter :: dp = selected_real_kind(15)
  real(dp), dimension(mr+mr2), intent(out) :: x, ra, acrit
  real(dp), dimension(mr+mr2), intent(out) :: sigg, du, pebs, siggc, duc, pebsc, rac
  real(dp), dimension(mr+mr2) :: gauss
  real(dp), dimension(mr+mr2) :: x2
  real(dp) :: rcau, rpl, dxtot, dxmin, xc
  real(dp), intent(in) :: md, ga, rc, mplinit, dxref, dxhr, xin, xout, m_star
  real(dp), intent(out) :: m_pl, ts
  real(dp), dimension(mr+mr2) :: diff
  real(dp) :: dmin, dmout, pmin, pmout, gmin, gmout

  !Changing the planet from earth masses to [g]
  m_pl = m_earth * mplinit

  !Making an equidistant logarithmic grid in [AU] COARSE GRID
  do i = 1, mr
     x(i) = log(xin) + (i-1) * ( (log(xout) - log(xin)) / (mr-1))
  enddo

  !Making the refined equidistant grid [AU]
!  dxtot = dxref * rpl
!  dxmin = (rpl * (m_earth / (3. * m_star))**(1./3.)) * dxhr
  x = exp(x)
!  call gridrefine(x, rpl, dxtot, dxmin, mr, mr2, mr3)

  !Skipping the creating of the refined grid
  !We have to make sure that every mr is now mr3
  mr3 = mr
  
  !Now the the refined grid is made
  !We continue by making the cell centers
  do i = 1, mr3-1
     x2(i) = sqrt(x(i+1)*x(i))
  enddo

  do i = 1, mr3
     diff(i) = abs((x(i)) - rpl)
  enddo

  !Find the planet location index
  idx_loc = int(minloc(diff(1:mr3), dim=1))
  xc = x(idx_loc)

  !Converting grids to cm
  ra = x * au
  rac = x2 * au

  !Converting rcrit to [cm]
  rcau = rc*au

  !Initializing the first gas surface density
  sigg = (2. - ga) * (md / (2. * pi * rcau**2.)) * &
       (ra / rcau) ** (-ga) * &
       exp(-(ra / rcau)**(2. - ga))

  call interp1d(mr3, ra, sigg, rac, siggc)
  
  !Initializing the dust density
  do i = 1, mr3
!     du(i) = sigg(i) *  1.e-30
     du(i) = sigg(i) * epsilon
  enddo
  call interp1d(mr3, ra, du, rac, duc)

  !Initializing the pebble density floor
  do i = 1, mr3
     pebs(i) = sigg(i) * 1.e-30
!     pebs(i) = sigg(i) * epsilon * 0.5
  enddo
  call interp1d(mr3, ra, pebs, rac, pebsc)

  !Initializing time in [s]
  ts = 0.

  !Initializing the critical size of pebbles
  acrit = 1.e10

  dmin = 0.
  dmout = 0.
  pmin = 0.
  pmout = 0.
  gmin = 0.
  gmout = 0.

end subroutine sigma_init
