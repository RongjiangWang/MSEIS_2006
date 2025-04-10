      subroutine msbsj(dk)
      implicit none
      double precision dk
c
      include 'msglobal.h'
c
      integer i,ir,ik
      double precision k,x
      double precision bessj0,bessj1,bessj
c
      do ik=1,nbsjmax
        k=dble(ik)*dk
        do ir=1,nr
          x=k*r(ir)
          bsj(ik,0,ir)=bessj0(x)
          bsj(ik,1,ir)=bessj1(x)
          if(x.gt.2.d0)then
            bsj(ik,2,ir)=bsj(ik,1,ir)*2.d0/x-bsj(ik,0,ir)
          else
            bsj(ik,2,ir)=bessj(2,x)
          endif
          if(x.gt.3.d0)then
            bsj(ik,3,ir)=bsj(ik,2,ir)*4.d0/x-bsj(ik,1,ir)
          else
            bsj(ik,3,ir)=bessj(3,x)
          endif
          bsj(ik,-1,ir)=-bsj(ik,1,ir)
        enddo
      enddo
c
      return
      end