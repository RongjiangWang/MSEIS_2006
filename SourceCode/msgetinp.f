      subroutine msgetinp(unit,srate)
      implicit none
      integer unit
      double precision srate
c
      include 'msglobal.h'
c
c     work space
c
      integer i,ir,j,l,l1,lcut,n,ierr
      integer irlast,irnow,is,ns
      integer ieqdis,flen0,nup,nlw
      double precision z,pi,taunorm,rnow,rlast,tnow,tlast
      double precision v00,depth,hpmin
      double precision s1,s2,s,ds,smin,shead,twinmin
      double precision r1,r2,dr,dm,t1,rdis
      double precision suppress,ros,vps,vss,fcut
      double precision resolut(3),t0(nrmax)
      character*80 outfile0,comments*180
c
c     source parameters
c     =================
c
      pi=4.d0*datan(1.d0)
      call getdata(unit,comments)
      read(comments,*)zs
      zs=dmax1(0.d0,km2m*zs)
c
c     receiver parameters
c     ===================
c
      call getdata(unit,comments)
      read(comments,*)zr
      zr=km2m*zr
      call getdata(unit,comments)
      read(comments,*)ieqdis
      call getdata(unit,comments)
      read(comments,*)nr
      if(nr.gt.nrmax)then
        stop 'Error in msgetinp: nr > nrmax!'  
      endif
      if(ieqdis.eq.1.and.nr.gt.1)then
        read(unit,*)r1,r2
        if(nr.eq.1)then
          dr=0.d0
        else
          dr=(r2-r1)/dble(nr-1)
        endif
        do i=1,nr
          r(i)=km2m*(r1+dr*dble(i-1))
        enddo
      else
        read(unit,*)(r(i),i=1,nr)
        do i=1,nr
          r(i)=km2m*r(i)
        enddo
      endif
      do ir=2,nr
        if(r(ir).lt.r(ir-1))then
          stop ' Error in msgetinp: distance array wrong ordered!'
        endif
      enddo
      if(r(1).le.0.d0)then
        stop ' Error in msgetinp: distance <= 0!'
      endif
c
      call getdata(unit,comments)
      read(comments,*)tstart,twindow,nt
      if(twindow.le.0.d0.or.nt.le.0)then
        stop 'Error in msgetinp: time window or sampling no <= 0!'  
      endif
c
      call getdata(unit,comments)
      read(comments,*)v0
      v0=km2m*v0
c
c     wavenumber integration parameters
c     =================================
c
      call getdata(unit,comments)
      read(comments,*)nd
      if(nd.lt.0.or.nd.gt.ndmax)then
        stop 'Error in msgetinp: wrong select of integration algorithm!'
      endif
      call getdata(unit,comments)
      read(comments,*)(slw(j),j=1,4)
      do j=1,4
        slw(j)=slw(j)/km2m
      enddo
      if(slw(1).le.0.d0.or.slw(2).le.0.d0.or.
     +   slw(3).le.0.d0.or.slw(4).le.0.d0.or.
     +   slw(2).lt.slw(1).or.slw(3).lt.slw(2).or.
     +   slw(4).lt.slw(3))then
        autoslwcut=.true.
      else
        autoslwcut=.false.
      endif
      call getdata(unit,comments)
      read(comments,*)srate
      if(srate.lt.1.d0)srate=1.d0
c
      call getdata(unit,comments)
      read(comments,*)suppress
      if(suppress.le.0.d0.or.suppress.ge.1.d0)then
        suppress=dexp(-1.d0)
        print *,'warning in qsmain: aliasing suppression'
        print *,'factor is replaced by the default value of 1/e.'
      endif
      fi=dlog(suppress)/(2.d0*pi*twindow)
c
c     wavelet parameters
c     ==================
c
      call getdata(unit,comments)
      read(comments,*)taunorm,wdeg
      if(wdeg.lt.0.or.wdeg.gt.2)then
        stop ' Error in qsmain: wrong wavelet selection!'
      else if(wdeg.eq.0)then
        call getdata(unit,comments)
        read(comments,*)nn0
        read(unit,*)(wv0(i),i=1,nn0)
      endif
c
c     seimometer parameters
c     =====================
c
      call getdata(unit,comments)
      read(comments,*)asm
      call getdata(unit,comments)
      read(comments,*)nroot
      read(unit,*)(root(i),i=1,nroot)
      call getdata(unit,comments)
      read(comments,*)npole
      read(unit,*)(pole(i),i=1,npole)
c
c     output files
c     ============
c
      varbtxt='U'
      call getdata(unit,comments)
      read(comments,*)outfile0
      do flen0=80,1,-1
        if(outfile0(flen0:flen0).ne.' ')goto 100
      enddo
100   continue
      outfile(1)=outfile0(1:flen0)//'.tz'
      outfile(2)=outfile0(1:flen0)//'.tr'
      outfile(3)=outfile0(1:flen0)//'.tp'
      do i=1,3
        flen(i)=flen0+3
      enddo
c
c     global model parameters
c     =======================
c
      call getdata(unit,comments)
      read(comments,*)(resolut(i),i=1,3)
      do i=1,3
        if(resolut(i).le.0.d0)resolut(i)=0.1d0
        resolut(i)=1.d-02*resolut(i)
      enddo
      call getdata(unit,comments)
      read(comments,*)l
      if(l.gt.lmax)then
        stop ' Error in msgetinp: to large number of layers!'
      endif
c
c     multilayered model parameters
c     =============================
c
      do i=1,l
        call getdata(unit,comments)
        read(comments,*)j,h(i),vp(i),vs(i),ro(i),qp(i),qs(i)
c
c       input units:    -,km,  km/s, km/s, g/cm^3,-,-
c
        h(i)=km2m*h(i)
        vp(i)=km2m*vp(i)
        vs(i)=km2m*vs(i)
        ro(i)=km2m*ro(i)
      enddo
      if(vs(1).gt.0.d0)then
        stop ' Error in msgetinp: no ocean layer at top!'
      endif
c
c     end of inputs
c     =============
c
      dt=twindow/dble(nt-1)
      nf=1
200   nf=2*nf
      if(nf.lt.nt)goto 200
      nf=nf/2
      if(nf.gt.nfmax)then
        print *,'Error in msgetinp: time sampling no > ',2*nfmax,'!'
        stop
      endif
      df=1.d0/(dble(2*nf)*dt)
      fcut=0.5d0/dt
      tau=taunorm*dt
      if(taunorm.le.0.d0)tau=2.d0*dt
c
      comptxt(1)='z'
      comptxt(2)='r'
      comptxt(3)='p'
      do j=1,nr
        i=j/1000
        rcvtxt(j)(1:1)=char(ichar('0')+i)
        i=mod(j,1000)/100
        rcvtxt(j)(2:2)=char(ichar('0')+i)
        i=mod(j,100)/10
        rcvtxt(j)(3:3)=char(ichar('0')+i)
        i=mod(j,10)
        rcvtxt(j)(4:4)=char(ichar('0')+i)
      enddo
c
c     determine upper und lower parameter values of each layer
c
      l0=1
      z1(l0)=0.d0
      do i=2,l
        if(h(i).gt.h(i-1))then
          z1(l0)=h(i-1)
          vp1(l0)=vp(i-1)
          vs1(l0)=vs(i-1)
          ro1(l0)=ro(i-1)
          qp1(l0)=qp(i-1)
          qs1(l0)=qs(i-1)
c
          z2(l0)=h(i)
          vp2(l0)=vp(i)
          vs2(l0)=vs(i)
          ro2(l0)=ro(i)
          qp2(l0)=qp(i)
          qs2(l0)=qs(i)
          l0=l0+1
        else
          z1(l0)=h(i)
          vp1(l0)=vp(i)
          vs1(l0)=vs(i)
          ro1(l0)=ro(i)
          qp1(l0)=qp(i)
          qs1(l0)=qs(i)
        endif
      enddo
      z1(l0)=h(l)
      vp1(l0)=vp(l)
      vs1(l0)=vs(l)
      ro1(l0)=ro(l)
      qp1(l0)=qp(l)
      qs1(l0)=qs(l)
c
c     construction of sublayers at the cutoff frequency
c
      call mssublay(resolut,fcut)
      write(*,*)' The layered model:'
      write(*,'(7a)')'    no ','  z(km)  ',
     &               '  vp(km/s) ','  vs(km/s) ',' ro(g/cm^3)',
     &               '     qp   ','     qs'
      depth=0.d0
      do i=1,n0
        write(*,1000)i,depth/km2m,vp(i)/km2m,
     &               vs(i)/km2m,ro(i)/km2m,qp(i),qs(i)
        depth=depth+h(i)
        if(i.lt.n0)then
          write(*,1000)i,depth/km2m,vp(i)/km2m,
     &               vs(i)/km2m,ro(i)/km2m,qp(i),qs(i)
        endif
      enddo
c
      call mslayer(ierr)
      n=nno(ls)
      ros=ro(n)
      vps=vp(n)
      vss=vs(n)
      call mssource(ros,vps,vss)
c
c     find ocean bottom
c
      lob=lp+1
      do l=1,lp
        n=nno(l)
        if(vs(n).gt.0.d0)then
          lob=l
          goto 300
        endif
      enddo
300   continue
      if(ls.ge.lob)then
        stop ' Error in msgetinp: source not in ocean layer!'
      endif
c
c     for marine seismic application (seismograms only for receiver in solid medium)
c
      if(lzr.lt.lob)then
        hydroseis=.true.
        fsel(1)=0
        fsel(2)=0
        fsel(3)=1
      else
        hydroseis=.false.
        fsel(1)=1
        fsel(2)=1
        fsel(3)=0
      endif
c
      if(v0.gt.0.d0)then
        v00=1.d0/v0
      else
        v00=0.d0
      endif
      do ir=1,nr
        t0(ir)=tstart+r(ir)*v00
      enddo
c
      hpmin=hp(min0(ls,lzr))
      smin=1.d0/vp(nno(min0(ls,lzr)))
      do l=min0(ls,lzr),max0(ls,lzr)-1
        if(smin.gt.1.d0/vp(nno(l)))then
          smin=1.d0/vp(nno(l))
          hpmin=hp(l)
        endif
      enddo
c
c     compare direct p, reflected p and head wave phase
c
      lcut=0
      do l=max0(ls,lzr),lp
c
        twinmin=twindow
c
c       1. direct or reflected p wave
c
        irlast=0
        rlast=0.d0
        tlast=0.d0
        do l1=min0(ls,lzr),max0(ls,lzr)-1
          tlast=tlast+hp(l1)/vp(nno(l1))
        enddo
        do l1=max0(ls,lzr),l-1
          tlast=tlast+2.d0*hp(l1)/vp(nno(l1))
        enddo
        s1=0.d0
        s2=smin/dsqrt(1.d0+0.5d0*(hpmin/r(nr))**2)
        ns=2*nr+10
        ds=s2/dble(ns)
        do is=1,ns
          s=s1+dble(is)*ds
          rnow=0.d0
          tnow=0.d0
          do l1=min0(ls,lzr),max0(ls,lzr)-1
            rdis=hp(l1)*s/dsqrt(1.d0/vp(nno(l1))**2-s**2)
            rnow=rnow+rdis
            tnow=tnow+dsqrt(rdis**2+hp(l1)**2)/vp(nno(l1))
          enddo
          do l1=max0(ls,lzr),l-1
            rdis=hp(l1)*s/dsqrt(1.d0/vp(nno(l1))**2-s**2)
            rnow=rnow+2.d0*rdis
            tnow=tnow+2.d0*dsqrt(rdis**2+hp(l1)**2)/vp(nno(l1))
          enddo
          irnow=irlast
          do ir=irlast+1,nr
            if(r(ir).gt.rnow)then
              goto 400
            else if(r(ir).ge.rlast)then
              t1=(tlast*(rnow-r(ir))+tnow*(r(ir)-rlast))
     &          /(rnow-rlast)
              irnow=ir
              twinmin=dmin1(twinmin,t1-t0(ir))
            endif
          enddo
400       rlast=rnow
          tlast=tnow
          irlast=irnow
        enddo
c
c       2. head wave
c
        shead=1.d0/vp(nno(l))
        if(smin.gt.shead)then
          rlast=0.d0
          tlast=0.d0
          do l1=min0(ls,lzr),max0(ls,lzr)-1
            rdis=hp(l1)*shead/dsqrt(1.d0/vp(nno(l1))**2-shead**2)
            rlast=rlast+rdis
            tlast=tlast+dsqrt(rdis**2+hp(l1)**2)/vp(nno(l1))
          enddo
          do l1=max0(ls,lzr),l-1
            rdis=hp(l1)*shead/dsqrt(1.d0/vp(nno(l1))**2-shead**2)
            rlast=rlast+2.d0*rdis
            tlast=tlast+2.d0*dsqrt(rdis**2+hp(l1)**2)/vp(nno(l1))
          enddo
          irlast=1
          do ir=1,nr
            if(r(ir).lt.rlast)irlast=ir+1
          enddo
          do ir=irlast,nr
            t1=tlast+shead*(r(ir)-rlast)
            twinmin=dmin1(twinmin,t1-t0(ir))
          enddo
          smin=shead
          hpmin=hp(l)
        endif
c
        if(twinmin.lt.twindow)lcut=l
      enddo
      if(lcut.lt.1)then
        stop ' time window too small!'
      else if(lcut.lt.lp)then
        lp=lcut
        hp(lp)=0.d0
        n0=nno(lp)
        write(*,'(a,i3)')' actually used number of layers: ',n0
      endif
1000  format(i5,f12.2,3f11.4,2f9.1)
c
      return
      end
