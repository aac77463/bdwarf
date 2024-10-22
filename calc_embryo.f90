subroutine calc_embryo(mr, idx, du, gas, ra, h, m_plt)
  use constants
  implicit none
  integer, parameter :: dp = selected_real_kind(15)
  integer :: mr, idx
  real(dp), dimension(mr), intent(in) :: ra, gas, du, h
  real(dp), dimension(mr) :: h_gas, z
  real(dp), intent(out) :: m_plt
  z = du/gas
  h_gas = h/ra
  print *, 'gas scale height', h_gas

  m_plt = m_earth * 1.2e-6 * (z(idx)/0.02)**(0.5) * (gas(idx)/1700.)**(3./2.) * (h_gas(idx)/0.033)**(3./2.) *&
       (ra(idx)/(3.*au))**(9./8.)

  print*, 'initial planet mass [Earth]', m_plt/m_earth
  read *
end subroutine calc_embryo
