      subroutine mssource(ros,vps,vss)
      implicit none
c
      double precision ros,vps,vss
c
      include 'msglobal.h'
c
      double precision pi2
c
      pi2=8.d0*datan(1.d0)
c
c     istp = 1
c     explosion source in water (m11=m22=m33=1)
c
      sfct=-1.d0/(pi2*ros*vps*vps)
c
      return
      end
