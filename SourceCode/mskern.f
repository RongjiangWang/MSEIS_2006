      subroutine mskern(y,f,k)
      implicit none
c
c     calculation of response to p-sv source
c     y(4): solution vector (complex)
c     f: frequency
c     k: wave number
c
      double precision f,k
      double complex y(4)
c
      include 'msglobal.h'
c
c     work space
c
      integer i,j,l,n,lup,llw,key
      double complex ck,ch0,pwave,swave,delta
      double complex y0(4,2),c0(4,2),c1(4,2),b(2)
      double complex cinc(4,6),hkw(2,2)
      double complex y1(4,2),yup(4,2),ylw(4,2),orth(2,2)
      double complex coef(2,2)
c
      double precision eps,pi2
      double complex c2
      data eps,pi2/1.0d-12,6.28318530717959d0/
	data c2/(2.d0,0.d0)/
c
      ck=dcmplx(k,0.d0)
      call mswaveno(f,k)
      do i=1,4
        y(i)=(0.d0,0.d0)
      enddo
c
      lup=1
c
c===============================================================================
c
c     matrix propagation from surface to source
c
c     determination of starting upper sublayer
c
      do j=1,2
        do i=1,4
          yup(i,j)=(0.d0,0.d0)
        enddo
      enddo
c
      yup(1,1)=(1.d0,0.d0)
      yup(2,1)=-accair/kpair
c
      if(lup.eq.lzr)call cmemcpy(yup,y0,8)
c
      do l=lup+1,ls
        ch0=dcmplx(hp(l-1),0.d0)
        n=nno(l-1)
c
c       determination of propagation matrix
c
        pwave=cdexp(-kp(n)*ch0)
c
        if(l.gt.lzr)then
c
c         additional normalization to avoid overflow
c
          do i=1,2
            y0(i,1)=y0(i,1)*pwave
          enddo
        endif
        call mshkwa(hkw,f,k,hp(l-1),n)
        call caxcb(hkw,yup(1,1),2,2,1,y1(1,1))
        yup(1,1)=y1(1,1)
        yup(2,1)=y1(2,1)
c
        if(l.eq.lzr)call cmemcpy(yup,y0,8)
      enddo
c
c===============================================================================
c
c     matrix propagation from half-space to source
c
      if(lp.ge.lob)then
        n=nno(lp)
        ylw(1,1)=-kp(n)
        ylw(2,1)=wa(n)
        ylw(3,1)=ck
        ylw(4,1)=-wb(n)*kp(n)
c
        ylw(1,2)=ck
        ylw(2,2)=-wb(n)*ks(n)
        ylw(3,2)=-ks(n)
        ylw(4,2)=wa(n)
c
        if(lp.gt.ls.and.lp.eq.lzr)call cmemcpy(ylw,y0,8)
c
        do l=lp-1,lob,-1
          ch0=dcmplx(hp(l),0.d0)
          n=nno(l)
c
c         determination of propagation matrix
c
          call msve2am(n,ck,ylw,c0,2)
          pwave=cdexp(-kp(n)*ch0)
          swave=cdexp(-ks(n)*ch0)
c
c         orthonormalization of the p-sv modes
c
          delta=(1.d0,0.d0)/(c0(4,2)*c0(2,1)-c0(2,2)*c0(4,1))
          orth(1,1)=c0(4,2)*delta
          orth(1,2)=-c0(2,2)*delta
          orth(2,1)=-c0(4,1)*delta
          orth(2,2)=c0(2,1)*delta
          call caxcb(c0,orth,4,2,2,c1)
c
          if(l.lt.lzr)then
c
c           additional normalization to avoid overflow
c
            do i=1,2
              orth(i,1)=orth(i,1)*pwave
              orth(i,2)=orth(i,2)*swave
            enddo
            call caxcb(y0,orth,4,2,2,y1)
            call cmemcpy(y1,y0,8)
          endif
c
          c1(1,1)=c1(1,1)*pwave*pwave
          c1(2,1)=(1.d0,0.d0)
          c1(3,1)=c1(3,1)*pwave*swave
          c1(4,1)=(0.d0,0.d0)
c
          c1(1,2)=c1(1,2)*swave*pwave
          c1(2,2)=(0.d0,0.d0)
          c1(3,2)=c1(3,2)*swave*swave
          c1(4,2)=(1.d0,0.d0)
c
          call msam2ve(n,ck,ylw,c1,2)
          if(l.gt.ls.and.l.eq.lzr)call cmemcpy(ylw,y0,8)
        enddo
c
        b(1)=ylw(4,2)
        b(2)=-ylw(4,1)
        do i=1,4
          ylw(i,1)=b(1)*ylw(i,1)+b(2)*ylw(i,2)
          ylw(i,2)=(0.d0,0.d0)
        enddo
        if(lob.eq.lzr)then
          call cmemcpy(ylw,y0,8)
        else if(lob.lt.lzr)then
          do i=1,4
            y0(i,1)=b(1)*y0(i,1)+b(2)*y0(i,2)
            y0(i,2)=(0.d0,0.d0)
          enddo
        endif
        ylw(3,1)=(0.d0,0.d0)
        llw=lob
      else
        n=nno(lp)
	  ylw(1,1)=kp(n)
	  ylw(2,1)=acc(n)
        ylw(3,1)=(0.d0,0.d0)
        ylw(4,1)=(0.d0,0.d0)
        do i=1,4
          ylw(i,2)=(0.d0,0.d0)
        enddo
        if(lp.gt.ls.and.lp.eq.lzr)call cmemcpy(ylw,y0,8)
        llw=lp
      endif
c
      do l=llw-1,ls,-1
        ch0=dcmplx(hp(l),0.d0)
        n=nno(l)
        pwave=cdexp(-kp(n)*ch0)
        if(l.lt.lzr)then
c
c         additional normalization to avoid overflow
c
          do i=1,4
            y0(i,1)=y0(i,1)*pwave
          enddo
        endif
        call mshkwa(hkw,f,k,-hp(l),n)
        call caxcb(hkw,ylw(1,1),2,2,1,y1(1,1))
        ylw(1,1)=y1(1,1)
        ylw(2,1)=y1(2,1)
c
        if(l.gt.ls.and.l.eq.lzr)call cmemcpy(ylw(1,1),y0(1,1),2)
      enddo
c
c===============================================================================
c     source function
c===============================================================================
c
      b(1)=dcmplx(sfct,0.d0)
      b(2)=(0.d0,0.d0)
      do i=1,2
        coef(i,1)=yup(i,1)
        coef(i,2)=-ylw(i,1)
      enddo
      key=0
      call cdgemp(coef,b,2,1,0.d0,key)
      if(key.eq.0)then
        print *,'warning in mskern: anormal exit from cdgemp!'
        return
      endif
      if(lzr.le.ls)then
        do i=1,4
          y(i)=b(1)*y0(i,1)
        enddo
      else
        do i=1,4
          y(i)=b(2)*y0(i,1)
        enddo
      endif
      return
      end
