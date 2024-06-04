!------------------------------------------------------------
!This module contains the objective functions for optimization.
!Users should modifiy this module according to their objective functions.
!XX is a vector of nparm control parameters
!nrow is the row location of each control parameter in the nominal Input.in file
!ncol is the column location of each control parameter in the nominal Input.in file
!ImpactT.in is the true input file for the ImpactT code after the control parameter update. 
!The user should use their simulation code's input file name and run their executable code.
!------------------------------------------------------------
!Author: Ji Qiang, LBNL
!
!
      module objfuncs

      contains
        !test objective functions
        !//Run a simulation through accelerator.
        subroutine objfunc(nrow,ncol,XX,objval,nparm,nobj)
        implicit none
        !include 'mpif.h'
!---------------------------------------------------
!for parallel optimization
        integer, intent(in) :: nparm,nobj
        real*8, dimension(nparm), intent(in) :: XX
        integer, dimension(nparm), intent(in) :: nrow,ncol
        real*8, dimension(nobj), intent(out) :: objval
!---------------------------------------------------
        integer :: i
        character*220 :: comst
        character*250 :: comst2
        character :: xbc
        integer :: j,j1,j2,j2b,iline,nlinemax,ii
        integer :: nchars,code,j3,j4
        integer :: j5,j6,j7,j8,j9,j50,j4b,j8b,j7b
        character(30) :: tmpstr
        character(30) :: tmpstr2
        real*8,dimension(8) :: fvalue
        real*8,dimension(10) :: objtmp1
        real*8 :: sigz,eng,rr,xtmp,xxtmp
        integer :: imod,ip

        !nominal input file name (Input.in) and 
        !simulation input file name (impactT.in).
        !The content (control parameters) in the nominal input file
        !will be modified by the control parameters and output into
        !the true input file for simulation run
        open(unit=1,file='Input.in',status='old')
        open(unit=2,file='ImpactT.in',status='unknown')

!-----------------
! prepare ImpactT simulation inputs using the control parameters
!
        nlinemax = 1000000
        iline = 0
        ii = 0
        do i = 1, nlinemax
          READ(1, "(A)", ADVANCE="no", SIZE=nchars, IOSTAT=code,end=111) comst
 
          iline = iline + 1
         do ip = 1, nparm
          if(i.eq.nrow(ip)) then
            j1 = 1
            j2 = 0
            j2b = 0
            j3 = 0
            j4 = 0
            j7 = 0
            j4b = 0
            j7b = 0
            do j = 1, nchars
              if(comst(j:j).ne." ") then
                j2 = j2 + 1
                xbc = comst(j:j)
                j2b = 0
              else
                j2b = j2b + 1
                j2 = 0
                if(j2b.eq.1) j1 = j1+1
              endif

              if(j1.lt.ncol(ip)) then
                  j3 = j3 + 1
              else if(j1.gt.ncol(ip)) then
                  j4 = j4 + 1
              else
                  j7 = j7 + 1
              endif
            enddo
            write(tmpstr,*)XX(ip)
            j8 = len(tmpstr)
            j8b = 0
            j4b = 0

            j5 = j3+j8+j4+j8b+j4b
            do j = 1, j3
              comst2(j:j) = comst(j:j)
            enddo
            do j = j3+1, j3+j8
              j9 = j - j3
              comst2(j:j) = tmpstr(j9:j9)
            enddo
            do j = j3+j8+1,j3+j8+j4
              j6 = j-j8+j7
              comst2(j:j) = comst(j6:j6)
            enddo
            comst(1:j5) = comst2(1:j5)
            nchars = j5
          else
            imod = 0
          endif
         enddo
         write(2,"(A)")comst(1:nchars)

        enddo
111       continue
        close(1)
        close(2)

        call flush(2)

!run the simulation
        call system("./ImpactTexe > bb")

        objtmp1 = 0.0d0
        !this loop is a dummy delay loop to ensure that fort.* files have been
        !produced.
        do i = 1, 10000
          xxtmp = (i-1)*0.00005d0
          xtmp = cos((i-1)*sqrt(xxtmp))
        enddo

        open(unit=12,file='fort.24',status='old')
        do i = 1, nlinemax
          read(12,*,end=444)fvalue(1:8)
        enddo
444       continue
        close(12)
        objtmp1(1) = fvalue(8) !x rms projected emittance

        open(unit=12,file='fort.25',status='old')
        do i = 1, nlinemax
          read(12,*,end=445)fvalue(1:8)
        enddo
445       continue
        close(12)
        objtmp1(2) = fvalue(8) !y rms emittance

        open(unit=14,file='fort.26',status='old')
        do i = 1, nlinemax
          read(14,*,end=446)fvalue(1:7)
        enddo
446       continue
        close(14)
        sigz = fvalue(3) !z rms bench length

        open(unit=15,file='fort.18',status='old')
        do i = 1, nlinemax
          read(15,*,end=447)fvalue(1:7)
        enddo
447       continue
        close(15)
        eng = fvalue(4) !beam kinetic energy (MeV)

        if(eng.gt.5.0) then
          objval(1) = (objtmp1(1)+objtmp1(2))/2e-6 !mm-mrad
          objval(2) = sigz/1e-3 !in mm
        else
          call random_number(rr)
          objval(1) = 4000.0d0 + rr
          call random_number(rr)
          objval(2) = 5000.0d0 + rr
        endif
        
        end subroutine objfunc

      end module objfuncs
