subroutine write_2d(fp, mr3, sigg, du, pebs, tcrit, vdrift, tmid, logp, adrift, afrag, acrit, agrow, &
     adiff, stcrit, pebint, dlogqusigg, dlogqunu, dlogsinusq, dlogp, dlogshquom, d2logqusigg, d2logsinusq, &
     d2logp, nu, sh, p, tor_vel, cs, siggc, duc, pebsc, afid, dfid, afip, dfip, st0, gasv, ducg)
  implicit none
  integer, parameter :: dp = selected_real_kind(15)
  integer, intent(in) :: fp, mr3
  real(dp), dimension(mr3) :: sigg, du, pebs, tcrit, vdrift, tmid, logp, adrift, afrag, acrit, agrow
  real(dp), dimension(mr3) :: adiff, stcrit, pebint, dlogqusigg, dlogqunu, dlogsinusq, dlogp, dlogshquom
  real(dp), dimension(mr3) :: d2logqusigg, d2logsinusq, d2logp, nu, sh, p, tor_vel, cs, siggc, duc, pebsc
  real(dp), dimension(mr3) :: afid, dfid, afip, dfip, st0, gasv, ducg
  write(fp, *) sigg, du, pebs, tcrit, vdrift, tmid, logp, adrift, afrag, acrit, agrow, &
       adiff, stcrit, pebint, dlogqusigg, dlogqunu, dlogsinusq, dlogp, dlogshquom, &
       d2logqusigg, d2logsinusq, d2logp, nu, sh, p, tor_vel, cs, siggc, duc, pebsc, &
       afid, dfid, afip, dfip, st0, gasv, ducg
end subroutine write_2d
