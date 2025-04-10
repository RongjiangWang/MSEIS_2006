      program mseis
      implicit none
c
      include 'msglobal.h'
c
c     work space
c
      integer i,runtime
      double precision pi,srate,tuser
      integer time
      character*15 passwd,passwd0
      data passwd0/'sanyuye27022008'/
c
c     read input file file
c
      print *,'######################################################'
      print *,'#                                                    #'
      print *,'#               Welcome to the program               #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#       M   M    SSSS    EEEEE    III     SSSS       #'
      print *,'#       MM MM   S        E         I     S           #'
      print *,'#       M M M    SSS     EEEE      I      SSS        #'
      print *,'#       M   M       S    E         I         S       #'
      print *,'#       M   M   SSSS     EEEEE    III    SSSS        #'
      print *,'#                                                    #'
      print *,'#                  (Version 2006)                    #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#                      by                            #'
      print *,'#                 Rongjiang Wang                     #'
      print *,'#              (wang@gfz-potsdam.de)                 #'
      print *,'#                                                    #'
      print *,'#           GeoForschungsZentrum Potsdam             #'
      print *,'#             Last modified: Feb 2008                #'
      print *,'######################################################'
      print *,'                          '
c      tuser=dble(time()/24/3600)-38.d0*365.25d0-150.d0
c      if(tuser.gt.180.d0)then
c        write(*,'(a,$)')' Password: '
c        read(*,'(a)')passwd
c        if(passwd.ne.passwd0)then
c          stop ' Program interupted: wrong pass word.'
c        endif
c      endif
      write(*,'(a,$)')' the input data file is '
      read(*,'(a)')inputfile
      runtime=time()
c
      pi=4.d0*datan(1.d0)
c
      open(10,file=inputfile,status='old')
      call msgetinp(10,srate)
      close(10)
c
      call mswvint(srate)
      iexist=0
      do i=1,3
        if(fsel(i).eq.1)call msfftinv(i)
      enddo
c
      runtime=time()-runtime
      write(*,'(a)')' #############################################'
      write(*,'(a)')' #                                           #'
      write(*,'(a)')' #      End of computations with mseis06     #'
      write(*,'(a)')' #                                           #'
      write(*,'(a,i10,a)')' #       Run time: ',runtime,
     +                                           ' sec            #'
      write(*,'(a)')' #############################################'
1001  format(2i7,E12.4,a)
1002  format(i4,a,E12.4,a,$)
1003  format(E12.5,$)
1004  format(2E12.4,$)
1005  format(2E12.4)
 500  stop
      end
