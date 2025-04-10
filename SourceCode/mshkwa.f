      subroutine mshkwa(hkw,f,k,z,n)
      implicit none
c
      integer n
      double precision f,k,z
      double complex hkw(2,2)
c
      include 'msglobal.h'
c
      double complex cx,cem,cch,csh
c
      cx=kp(n)*dcmplx(2.d0*z,0.d0)
      if(z.gt.0.d0)then
        cem=cdexp(-cx)
        cch=(0.5d0,0.d0)*((1.d0,0.d0)+cem)
        csh=(0.5d0,0.d0)*((1.d0,0.d0)-cem)
      else
        cem=cdexp(cx)
        cch=(0.5d0,0.d0)*((1.d0,0.d0)+cem)
        csh=-(0.5d0,0.d0)*((1.d0,0.d0)-cem)
      endif
      hkw(1,1)=cch
      hkw(1,2)=-kp(n)*csh/acc(n)
      hkw(2,1)=-acc(n)*csh/kp(n)
      hkw(2,2)=cch
c
      return
      end