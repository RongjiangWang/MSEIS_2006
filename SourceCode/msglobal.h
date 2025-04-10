c     GLOBAL INDEX PARAMETERS FOR DEFINING ARRAYS
c     ===========================================
c     nzmax: max. interface index;
c     lmax: max. number of total homogeneous layers (lmax <= nzmax-2);
c     nrmax: max. number of traces;
c     nfmax: max. number of frequency samples
c     nbsjmax: max. number of high-precission Bessel function samples
c
      integer nzmax,lmax,nrmax,nfmax,nbsjmax,ndmax
      parameter(lmax=100)
      parameter(nzmax=lmax+2)
      parameter(nrmax=3001,nfmax=2048)
      parameter(nbsjmax=10000)
      parameter(ndmax=4)
c
c     INDEX PARAMETERS FOR SEISMOMETER CHARACTERISTICS
c     ================================================
c     (max. number of roots and poles)
c
      integer nrootmax,npolemax
      parameter(nrootmax=10,npolemax=10)
c
c     LENGTH UNIT
c     ===========
c
      double precision km2m
      parameter(km2m=1.0d+03)
c
c     ATMOSPHERIC PARAMETERS
c     ======================
      double precision roair,vpair
      parameter(roair=0.1300d+01,vpair=0.3318d+03)
c
c     water quality factor: Q(f)=Q(1 Hz)*f^qahz
c     derived from J.F. Boehme (2001)
c
      double precision qahz
      parameter(qahz=-0.389d0)
c
      double complex accair,cvpair,kpair,comega
      common /airpara/ accair,cvpair,kpair,comega
c
c     zr: receiver depth
c     lzr: sublayer no of receiver
c     lob: sublayer no of ocean bottom
c
      integer lzr,lob
      double precision zr
      logical hydroseis
      common /dreceiver/ zr
      common /ireceiver/ lzr,lob
      common /lreceiver/ hydroseis
c
      integer nr
      double precision r(nrmax)
      common /distance/ r,nr
c
      integer lp,nno(nzmax)
      double precision hp(nzmax)
      common /dsublayer/ hp
      common /isublayer/ lp,nno
c
c     original model parameters
c
      integer l0
      double precision z1(lmax),z2(lmax),ro1(lmax),ro2(lmax)
      double precision vp1(lmax),vp2(lmax),vs1(lmax),vs2(lmax)
      double precision qp1(lmax),qp2(lmax),qs1(lmax),qs2(lmax)
      common /dmodel0/ z1,z2,ro1,ro2,vp1,vp2,vs1,vs2,qp1,qp2,qs1,qs2
      common /imodel0/ l0
c       
c     layered model parameter:
c     n0: number of homogeneous layers
c
      integer n0
      double precision h(lmax),ro(lmax),vp(lmax),vs(lmax)
      double precision qp(lmax),qs(lmax)
      common /dmodel/ h,ro,vp,vs,qp,qs
      common /imodel/ n0
c
      double complex acc(lmax),kp(lmax),ks(lmax),cla(lmax),cmu(lmax)
      double complex cvp(lmax),cvs(lmax),wa(lmax),wb(lmax)
      common /dpara/ acc,kp,ks,cla,cmu,cvp,cvs,wa,wb
c
c     source parameters
c
      integer ls
      double precision zs,sfct
      common /dsource/ zs,sfct
      common /isource/ ls
c
c     slowness cut-offs
c
      double precision slw(4)
      logical autoslwcut
      common /dslwcutoffs/ slw
      common /lslwcutoffs/ autoslwcut
c
c     path filtering
c
      integer nd
      common /dtransform/ nd
c
c     table of J_n(x), n = -1, 0, 1, 2, 3
c
      double precision bsj(nbsjmax,-1:3,nrmax)
      double precision zrs2,geospr(nrmax)
      common /dbessels/ bsj,zrs2,geospr
c
      integer nt,nf
      double precision dt,df,fi
      common /dsampling/ dt,df,fi
      common /isampling/ nt,nf
c
      double precision tstart,twindow,tau,v0
      common /dtparas/ tstart,twindow,tau,v0
c
      integer nnmax,nn0,iexist,wdeg
      parameter(nnmax=1024)
      double precision wv0(nnmax)
      common /dwavelets/ wv0
      common /iwavelets/ nn0,iexist,wdeg
c
c     seismometer filtering
c
      integer nroot,npole
      double precision asm
      double complex root(nrootmax),pole(npolemax) 
      common /dseismometer/ root,pole,asm
      common /iseismometer/ nroot,npole
c
c     green's functions
c
      double complex grns(nfmax,3,nrmax)
      common /dgrnfcts/ grns
c
c     title text
c
      character*1 comptxt(3),varbtxt
      character*4 rcvtxt(nrmax)
      common /ctitle/ comptxt,varbtxt,rcvtxt
c
c     input and output data files
c
      character*80 inputfile
      common /cinputddata/ inputfile
      integer ssel(3),fsel(3),flen(3)
      character*83 outfile(3)
      common /ioutsel/ ssel,fsel,flen
      common /foutdata/ outfile
