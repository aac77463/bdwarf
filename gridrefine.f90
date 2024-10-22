subroutine gridrefine(x,xc,dxtot,dxmin,nr,nr2,nr3)
! this subroutine creates a refined grid around the location xc
! grid runs from xc-dxtot to xc+dxtot
! the grid has a resolution of dxmin at the location xc
! The subroutine assumes the grid x(1:nr) is filled with a larger global grid
! the number of points used for the refined region is nr2
! the final total number of points is output as nr3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! So for a general use in a disk:
! x(1:nr) is filled with the low resolution grid of the entire disk
! xc is the location of the planet
! dxtot is typically set to 0.25*xc
! dxmin is typically set to 0.1*Rhill
! nr and nr2 can be chosen
! nr3 is set by the subroutine (and will typically be smaller than nr+nr2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IMPLICIT NONE
  integer, parameter :: dq = selected_real_kind(15)
  integer nr,nr2,nr3
  real(dq) ::  x(nr+nr2),x1,x2,xc,dxtot,dxmin,xx(nr2)
  real(dq) :: p
  integer i,i1,i2

  x1=xc-dxtot
  x2=xc+dxtot

!      find first index
  do i=1,nr
     if(x(i).gt.x1) exit
  enddo
  i1=i
  !       find second index
  do i=i1,nr
     if(x(i).gt.x2) exit
  enddo
  i2=i

  nr3=i1+nr-i2+nr2+1 !removing the point at the planet location
  if(nr3.gt.nr+nr2) then
     print*,'nr3 > nr+nr2 : exit'
     stop
  endif

  p=2d0
  p=-log(abs(dxmin/dxtot))/log(real(nr2/2))
  if(p.lt.1d0) p=1d0
  print*,p
  do i=1,nr2/2
     xx(i)=xc+dxtot*(real(i)/real(nr2/2))**p
  enddo
  do i=1,nr2/2
     xx(i+nr2/2)=xc-dxtot*(real(i)/real(nr2/2))**p
  enddo
  call sort(xx,nr2)
    
  !       move array
  x(i1)=xc !removing the point at the planet location
  x(i1+1:i1+nr-i2+1)=x(i2:nr)
!  x(i1:i1+nr-i2)=x(i2:nr)
  x(i1+nr-i2+2:nr3)=xx(1:nr2)
  !       sort the new array
  call sort(x,nr3)
  
  !       check for duplicates
1 continue
  do i=1,nr3-1
     if(x(i).eq.x(i+1)) then
        x(i:nr3-1)=x(i+1:nr3)
        nr3=nr3-1
        goto 1
     endif
  enddo
  
  return
end subroutine gridrefine
	
	
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
	
SUBROUTINE sort(arr,n)
  INTEGER n,M,NSTACK
    integer, parameter :: dq = selected_real_kind(15)
    REAL(dq) :: arr(n)
    PARAMETER (M=7,NSTACK=50)
    INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
    REAL(dq) ::  a,temp
    if(n.le.1) return
    jstack=0
    l=1
    ir=n
1   if(ir-l.lt.M)then
       do j=l+1,ir
          a=arr(j)
          do i=j-1,l,-1
             if(arr(i).le.a)goto 2
             arr(i+1)=arr(i)
          enddo
          i=l-1
2         arr(i+1)=a
       enddo
       if(jstack.eq.0)return
       ir=istack(jstack)
       l=istack(jstack-1)
       jstack=jstack-2
    else
       k=(l+ir)/2
       temp=arr(k)
       arr(k)=arr(l+1)
       arr(l+1)=temp
       if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
       endif
       if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
       endif
       if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
       endif
       i=l+1
       j=ir
       a=arr(l+1)
3      continue
       i=i+1
       if(arr(i).lt.a)goto 3
4      continue
       j=j-1
       if(arr(j).gt.a)goto 4
       if(j.lt.i)goto 5
       temp=arr(i)
       arr(i)=arr(j)
       arr(j)=temp
       goto 3
5      arr(l+1)=arr(j)
       arr(j)=a
       jstack=jstack+2
       if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
       if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
       else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
       endif
    endif
    goto 1
  END SUBROUTINE sort


	
