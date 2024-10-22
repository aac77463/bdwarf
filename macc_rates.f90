 subroutine macc_rates(mr, ra, m_pl, rpl, omg, macc_g, macc_d, macc_p, &
     r_h, pebs, stokes, m_star, ixinrh, ixoutrh)

  !!!!IMPORTANT!!!!!!!
  ! I need to read in the scale height of the gas and the stokes numbers of the particles
  ! in the disk in order to calculate the 3D accretion regime.

  !!$  This subroutine calculates the size of the hill sphere around the planet
  !!$  at the indicated mass. It also calculated the mass accretion rate of the
  !!$  gas and refractory components.
  use constants
  implicit none
  integer :: mr, i
  integer, intent(out) :: ixinrh, ixoutrh
  integer, parameter :: dp = selected_real_kind(15)
  real(dp) :: macc_g, macc_d, tot, wtz, r_h, m_star, m_pl
  real(dp), intent(in) :: rpl!,m_pl
  real(dp), dimension(mr), intent(in) :: ra, omg, stokes, pebs
  real(dp), dimension(mr) :: ra2, macc_p
  macc_g = 0.
  macc_d = 0.


  !Calculating the grid r1/2 for calculating interpolations
  do i = 2, mr
     ra2(i) = ( ra(i) + ra(i-1) ) / 2.
  enddo

  !Finding how big the hill radius is
  r_h = (rpl * au) * ( m_pl / (3. * m_star) )**(1./3.)


  !Calculating the pebble scale height
  !h_peb = h * sqrt(al/stokes)
  !cond_hpeb = (pi * (stokes/0.1)**(1./3.)*r_h) / (2. * sqrt(2. * pi))
  
  
  !Choosing the inner and outer points from where to accrete material
  !They depends on how big the hill radius becomes, and are updated as
  !The plane grows in mass
  ixoutrh=minloc(abs((((rpl * au)+r_h)-ra)/((rpl * au)+r_h)), dim=1)
  ixinrh=minloc(abs((ra-((rpl * au)-r_h))/((rpl * au)-r_h)), dim=1)

  !Tot is a weighting factor to make sure more material is accreted from
  !the grid points near the center of the planet, and less from those
  !at the edges
  tot = 0.
  macc_p = 0.
  do i = ixinrh+1, ixoutrh-1
     wtz = sqrt(max((r_h**2. - (ra(i) - (rpl * au))**2.),0.))
!     if (h_peb(i).lt.cond_hpeb(i)) then
!        macc_p(i) = 2. * ( ra2(i+1) - ra2(i) ) *  omg(i) * stokes(i)**(2./3.) * pebs(i) * wtz *&
!             (cond_hpeb(i))
!     else
!     macc_p(i) = 2. * ( ra2(i+1) - ra2(i) ) *  omg(i) * stokes(i)**(2./3.) * pebs(i) * wtz
!     endif
     macc_p(i) = 2. * ( ra2(i+1) - ra2(i) ) *  omg(i) * (stokes(i)/0.1)**(2./3.) * pebs(i) * wtz
     tot = tot + ( ra2(i+1) - ra2(i) ) * wtz
  enddo

  !!! Checking the points inner to outer points
  do i=ixoutrh,ixoutrh
     if(((r_h+(rpl * au)).gt.(ra2(i))))then
        wtz = ( ((rpl * au)+r_h) + ra2(i) )/2.
        wtz = sqrt(max((r_h**2. - (wtz - (rpl * au))**2.),0.))
        macc_p(i) = 2. * (((rpl * au)+r_h) - ra2(i)) * omg(i) *&
             stokes(i)**(2./3.) * pebs(i) * wtz
        tot = tot + ( ((rpl * au)+r_h) - ra2(i) ) * wtz
     else
        macc_p(i) = 2. * ( ra2(i+1) - ra2(i) ) *  omg(i) * stokes(i)**(2./3.) * pebs(i)
        tot = tot + (ra2(i+1) - ra2(i))
        print *, 'taking else'
     endif
  enddo

  do i=ixinrh,ixinrh
     if((((rpl * au)-r_h).lt.(ra2(i+1))))then
        wtz = ( ((rpl * au)+r_h) + ra2(i+1) )/2.
        wtz = sqrt(max((r_h**2. - (wtz - (rpl * au))**2.),0.))
        macc_p(i) = 2. * (ra2(i+1) - ((rpl * au) - r_h)) * omg(i) *&
             stokes(i)**(2./3.) * pebs(i) * wtz
        tot = tot + ( ra2(i+1) - ((rpl * au) - r_h) ) * wtz
     else
        macc_p(i) = 2. * ( ra2(i+1) - ra2(i) ) *  omg(i) * stokes(i)**(2./3.) * pebs(i)
        tot = tot + (ra2(i+1) - ra2(i))
        print *, 'taking else 2'
     endif
  enddo
  
  !Dividing the mass accretion rate using the previous integral

!  print *, 'macc_p', macc_p(ixinrh:ixoutrh)
!  print *, 'tot', tot
!  print *, 'r hill sq', r_h**2.
  macc_p = (macc_p / tot) * r_h**2.
!  print *, macc_p


end subroutine macc_rates

 subroutine macc_rates_3d(mr, ra, m_pl, rpl, omg, macc_g, macc_d, macc_p, &
     r_h, pebs, stokes, m_star, ixinrh, ixoutrh)

  !!!!IMPORTANT!!!!!!!
  ! I need to read in the scale height of the gas and the stokes numbers of the particles
  ! in the disk in order to calculate the 3D accretion regime.

  !!$  This subroutine calculates the size of the hill sphere around the planet
  !!$  at the indicated mass. It also calculated the mass accretion rate of the
  !!$  gas and refractory components.
  use constants
  implicit none
  integer :: mr, i
  integer, intent(out) :: ixinrh, ixoutrh
  integer, parameter :: dp = selected_real_kind(15)
  real(dp) :: macc_g, macc_d, tot, wtz, r_h, m_star, m_pl
  real(dp), intent(in) :: rpl!,m_pl
  real(dp), dimension(mr), intent(in) :: ra, omg, stokes, pebs
  real(dp), dimension(mr) :: ra2, macc_p
  macc_g = 0.
  macc_d = 0.

  !Calculating the pebble scale height
!  hpeb = h * sqrt(al/stokes)

  
  
  !Calculating the grid r1/2 for calculating interpolations
  do i = 2, mr
     ra2(i) = ( ra(i) + ra(i-1) ) / 2.
  enddo

  !Finding how big the hill radius is
  r_h = (rpl * au) * ( m_pl / (3. * m_star) )**(1./3.)

  !Choosing the inner and outer points from where to accrete material
  !They depends on how big the hill radius becomes, and are updated as
  !The plane grows in mass
  ixoutrh=minloc(abs((((rpl * au)+r_h)-ra)/((rpl * au)+r_h)), dim=1)
  ixinrh=minloc(abs((ra-((rpl * au)-r_h))/((rpl * au)-r_h)), dim=1)

  !Tot is a weighting factor to make sure more material is accreted from
  !the grid points near the center of the planet, and less from those
  !at the edges
  tot = 0.
  macc_p = 0.
  do i = ixinrh+1, ixoutrh-1
     wtz = sqrt(max((r_h**2. - (ra(i) - (rpl * au))**2.),0.))
     macc_p(i) = 2. * ( ra2(i+1) - ra2(i) ) *  omg(i) * stokes(i)**(2./3.) * pebs(i) * wtz
     tot = tot + ( ra2(i+1) - ra2(i) ) * wtz
  enddo

  !!! Checking the points inner to outer points
  do i=ixoutrh,ixoutrh
     if(((r_h+(rpl * au)).gt.(ra2(i))))then
        wtz = ( ((rpl * au)+r_h) + ra2(i) )/2.
        wtz = sqrt(max((r_h**2. - (wtz - (rpl * au))**2.),0.))
        macc_p(i) = 2. * (((rpl * au)+r_h) - ra2(i)) * omg(i) *&
             stokes(i)**(2./3.) * pebs(i) * wtz
        tot = tot + ( ((rpl * au)+r_h) - ra2(i) ) * wtz
     else
        macc_p(i) = 2. * ( ra2(i+1) - ra2(i) ) *  omg(i) * stokes(i)**(2./3.) * pebs(i)
        tot = tot + (ra2(i+1) - ra2(i))
        print *, 'taking else'
     endif
  enddo

  do i=ixinrh,ixinrh
     if((((rpl * au)-r_h).lt.(ra2(i+1))))then
        wtz = ( ((rpl * au)+r_h) + ra2(i+1) )/2.
        wtz = sqrt(max((r_h**2. - (wtz - (rpl * au))**2.),0.))
        macc_p(i) = 2. * (ra2(i+1) - ((rpl * au) - r_h)) * omg(i) *&
             stokes(i)**(2./3.) * pebs(i) * wtz
        tot = tot + ( ra2(i+1) - ((rpl * au) - r_h) ) * wtz
     else
        macc_p(i) = 2. * ( ra2(i+1) - ra2(i) ) *  omg(i) * stokes(i)**(2./3.) * pebs(i)
        tot = tot + (ra2(i+1) - ra2(i))
        print *, 'taking else 2'
     endif
  enddo
  
  !Dividing the mass accretion rate using the previous integral

  macc_p = (macc_p / tot) * r_h**2.
!  print *, macc_p


end subroutine macc_rates_3d
