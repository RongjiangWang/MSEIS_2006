      subroutine msqmodel(f)
      implicit none
c
c     calculate q based on the constant q model
c
c     f = frequency
c
      double precision f
c
      include 'msglobal.h'
c
      integer n
      double precision alfa
      double complex cmp,cms
c
      double precision pi2
      data pi2/6.28318530717959d0/
c
      comega=dcmplx(pi2*f,pi2*fi)
      cvpair=dcmplx(vpair,0.d0)
      accair=dcmplx(roair,0.d0)*comega**2
c
      alfa=0.5d0*f/dabs(fi)
c
      do n=1,n0
        if(qp(n).le.0.d0)then
          cmp=(1.d0,0.d0)
        else if(vs(n).le.0.d0)then
c
c         for water
c
          cmp=dcmplx(1.d0,0.5d0/(qp(n)*f**qahz))
        else
          cmp=dcmplx(1.d0,dmin1(0.5d0/qp(n),alfa))
        endif
        if(qs(n).le.0.d0.or.vs(n).le.0.d0)then
          cms=(1.d0,0.d0)
        else
          cms=dcmplx(1.d0,dmin1(0.5d0/qs(n),alfa))
        endif
        cvp(n)=dcmplx(vp(n),0.d0)*cmp
        cvs(n)=dcmplx(vs(n),0.d0)*cms
        cmu(n)=dcmplx(ro(n)*vs(n)**2,0.d0)*cms**2
        cla(n)=dcmplx(ro(n)*vp(n)**2,0.d0)*cmp**2-(2.d0,0.d0)*cmu(n)
        acc(n)=dcmplx(ro(n),0.d0)*comega**2
      enddo
c
      return
      end
