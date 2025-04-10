      subroutine mswaveno(f,k)
      implicit none
c
      double precision f,k
c
      include 'msglobal.h'
c
      integer n
      double complex ck,ck2
c
      ck=dcmplx(k,0.d0)
      ck2=dcmplx(k*k,0.d0)
c
      kpair=cdsqrt((ck+comega/cvpair)*(ck-comega/cvpair))
c
      do n=1,n0
        kp(n)=cdsqrt((ck+comega/cvp(n))*(ck-comega/cvp(n)))
        wb(n)=(2.d0,0.d0)*cmu(n)*ck
        if(vs(n).gt.0.d0)then
          ks(n)=cdsqrt((ck+comega/cvs(n))*(ck-comega/cvs(n)))
          wa(n)=cmu(n)*(ck2+ks(n)*ks(n))
        endif
      enddo
c
      return
      end
