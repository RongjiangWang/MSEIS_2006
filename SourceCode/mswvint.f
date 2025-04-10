      subroutine mswvint(srate)
      implicit none
c
      double precision srate
c
      include 'msglobal.h'
c
      integer n,l,lf,i,i1,i2,ir,ik,nk1,nk2,id,jd
      integer mj(3)
      double precision f,fcut,kmax,dk
      double precision pi,pi2,thickness,slwmax,fac,x
      double precision k(0:ndmax),ph(3),kcut(4),kcut1(4),kcut2(4)
      double complex ck,ck2,c2dk,cdk2,cfac
      double complex cy(4),y(3,-ndmax:ndmax,0:ndmax),cm2(3)
      double precision taper
c
      double complex c2
      data c2/(2.d0,0.d0)/
c
      pi=4.d0*datan(1.d0)
      pi2=2.d0*pi
c
      if(hydroseis)then
        i1=3
        i2=3
      else
        i1=1
        i2=2
      endif
      mj(1)=0
      mj(2)=1
      mj(3)=0
      ph(1)=0.25d0*pi
      ph(2)=0.75d0*pi
      ph(3)=0.25d0*pi
      cm2(1)=(0.d0,0.d0)
      cm2(2)=(1.d0,0.d0)
      cm2(3)=(0.d0,0.d0)
c
      fcut=dble(nf)*df
c
      thickness=0.d0
      do l=1,lp-1
        thickness=thickness+hp(l)
      enddo
c
      dk=pi/dmax1(srate*r(nr),3.d0*thickness)
      c2dk=dcmplx(2.d0*dk,0.d0)
      cdk2=dcmplx(dk*dk,0.d0)
c
      zrs2=(zs-zr)**2
      if(nd.eq.0)then
        do ir=1,nr
          geospr(ir)=1.d0
        enddo
      else
        do ir=1,nr
          geospr(ir)=1.d0/(zrs2+r(ir)*r(ir))**nd
        enddo
      endif
      call msbsj(dk)
c
      kcut1(1)=0.d0
      kcut1(2)=0.d0
      kcut1(3)=0.d0
      kcut1(4)=0.d0
      do i=1,4
        kcut2(i)=pi2*fcut*slw(i)
      enddo
c
      slwmax=1.d0/vp(nno(ls))
      kmax=pi2*fcut*slwmax
      if(autoslwcut)then
        kcut2(1)=0.d0
        kcut2(2)=0.d0
        kcut2(3)=kmax
        kcut2(4)=1.1d0*kmax
      else
        kcut2(1)=dmin1(kcut2(1),0.9d0*kmax)
        kcut2(2)=dmin1(kcut2(2),kmax)
        kcut2(3)=dmin1(kcut2(3),kmax)
        kcut2(4)=dmin1(kcut2(4),1.1d0*kmax)
      endif
c
      do lf=1,nf
        do i=i1,i2
          do ir=1,nr
            grns(lf,i,ir)=(0.d0,0.d0)
          enddo
        enddo
      enddo
c
      write(*,'(a,2(f10.7,a))')' Min./max. slowness at f_cut: ',
     &     1000.d0*kcut2(1)/(pi2*fcut),' / ',
     &     1000.d0*kcut2(4)/(pi2*fcut),' s/km'
c
      do lf=2,nf
        f=dble(lf-1)*df
        call msqmodel(f)
c
        do i=1,4
          kcut(i)=kcut1(i)+dble(lf-1)*(kcut2(i)-kcut1(i))/dble(nf)
        enddo
c
        nk1=max0(1,1+idint(kcut(1)/dk)-nd)
        nk2=nd+idint(kcut(4)/dk)
c
        do jd=0,nd
          do id=-nd,nd
            do i=i1,i2
              y(i,id,jd)=(0.d0,0.d0)
            enddo
          enddo
        enddo
c
        do ik=nk1,nk2
          k(0)=dble(ik)*dk
          ck=dcmplx(k(0),0.d0)
          call mskern(cy,f,k(0))
          y(1,nd,0)=cy(1)
          y(2,nd,0)=-cy(3)
          y(3,nd,0)=-cy(2)
          do jd=1,nd
            k(jd)=k(jd-1)-dk
            if(k(jd).gt.0.d0)then
              id=nd-jd
              ck=dcmplx(k(jd),0.d0)
              ck2=ck*ck
              do i=i1,i2
                y(i,id,jd)=y(i,id,jd-1)*(zrs2+cm2(i)/ck2)
     &             -(y(i,id+1,jd-1)-y(i,id-1,jd-1))/ck/c2dk
     &             -(y(i,id+1,jd-1)-c2*y(i,id,jd-1)+y(i,id-1,jd-1))/cdk2
              enddo
            endif
          enddo
          if(ik.ge.nk1+2*nd.and.ik.le.nk2-nd)then
            cfac=dcmplx(k(nd)*dk
     &          *taper(k(nd),kcut(1),kcut(2),kcut(3),kcut(4)))
            do i=i1,i2
              y(i,0,nd)=y(i,0,nd)*cfac
            enddo
            if(ik-nd.le.nbsjmax)then
              do ir=1,nr
                do i=i1,i2
                  grns(lf,i,ir)=grns(lf,i,ir)+y(i,0,nd)
     &                *dcmplx(bsj(ik-nd,mj(i),ir)*geospr(ir),0.d0)
                enddo
              enddo
            else
              do ir=1,nr
                x=k(nd)*r(ir)
                fac=dsqrt(2.d0/(pi*x))
                do i=i1,i2
                  grns(lf,i,ir)=grns(lf,i,ir)+y(i,0,nd)
     &                *dcmplx(fac*dcos(x-ph(i))*geospr(ir),0.d0)
                enddo
              enddo
            endif
          endif
          do jd=0,nd-1
            do id=-(nd-jd),(nd-jd)-1
              do i=i1,i2
                y(i,id,jd)=y(i,id+1,jd)
              enddo
            enddo
          enddo
        enddo
c
        write(*,'(i6,a,E13.6,a,i7)')lf,'.',f,
     &      'Hz: slowness samples = ',1+nk2-nk1-nd
      enddo
c
      return
      end