subroutine torque_vel_calc(mr, ra, sh, m_star, m_pl, tor_vel, rpl)
  use constants
  implicit none
  integer, parameter :: dp = selected_real_kind(15)
  integer :: mr, i
  real(dp), dimension(mr), intent(in) :: ra, sh
  real(dp), dimension(mr), intent(out) :: tor_vel
  real(dp), dimension(mr) :: signpm, rnew, rhill
  real(dp), intent(in) :: m_star, m_pl, rpl
  real(dp) :: q, tf
  tf = 0.8
  
  q = m_pl / m_star

  rhill = (rpl * au) * (q/3.)**(1./3.)
!  rnew = max((rhill),(sh))
  
  do i = 1, mr
     !Drift velocity is (-). (-) means tor_velque velocity is moving inward w/ vdrift
     !(+) means tor_velque velocity is moving outward against vdrift
     !     if (ra(i).gt.(rpl * au)+rhill(i)) then
     if (ra(i).gt.(rpl * au)) then
        tor_vel(i) = sqrt(ra(i) * g * m_star) * tf * q**2. * (ra(i)/max(abs(ra(i) - (rpl * au)), sh(i)))**4. / ra(i)
!        tor_vel(i) = sqrt(ra(i) * g * m_star) * tf * q**2. * (ra(i)/max(rnew(i), sh(i)))**4. / ra(i)
        signpm(i) = 1.
!       elseif (ra(i).eq.(rpl * au)) then
     !elseif(ra(i).lt.(rpl * au)-rhill(i))then
     elseif (ra(i).eq.(rpl * au)) then
        !tor_vel(i) = sqrt(ra(i) * g * m_star) * tf * q**2. * (ra(i)/max(rnew(i), sh(i)))**4. / ra(i)
        tor_vel(i) = 0.
        signpm(i) = 1.
     else
        tor_vel(i) =  sqrt(ra(i) * g * m_star) * tf * q**2. * (ra(i)/max(abs(ra(i) - (rpl * au)), sh(i)))**4. / ra(i)
        !tor_vel(i) = sqrt(ra(i) * g * m_star) * tf * q**2. * (ra(i)/max(rnew(i), sh(i)))**4. / ra(i)
!        tor_vel(i) = 0.
        signpm(i) = -1.
     endif
  enddo

  do i = 1, mr
     if (abs(ra(i) - (rpl * au)).lt.sh(i)) then
        signpm(i) = ( ra(i) - (rpl * au) ) / sh(i)
     endif
  enddo

  !Doing the approximation of the torque from Dangelo & Lubow 2010 where the near the planet location
  !the torque is equal to zero to prevent an early gap
  do i = 1, mr
!     print *, abs((ra(i) - (rpl*au))) / rhill(i), m_pl/m_earth
     !     if (((abs(ra(i)-(rpl*au)))/rhill(i).gt.0.5).and.((abs(ra(i)-(rpl*au)))/rhill(i).lt.1.5)) then
     if((ra(i).gt.(rpl*au-(3.*rhill(i)))).and.(ra(i).lt.(rpl*au+(3.*rhill(i)))))then
        tor_vel(i) = 0.
!        print *, tor_vel(i)
!        read *
     endif
  enddo


  tor_vel = tor_vel * signpm

end subroutine torque_vel_calc
