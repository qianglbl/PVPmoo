!.......................................................................
! PVPMOO: A Parallel Variable Population Multi-objective Optimizer
! with an Adaptive Unified Differential Evolution Algorithm or a Genetic Algorithm: v.1.0
!.......................................................................
!****************************
!
!*** Copyright Notice ***
!
!Parallel Variable Population Multi-Objective Optimizer (pvpmoo)
!Copyright (c) 2024, The Regents of the University of California,
!through Lawrence Berkeley National Laboratory (subject to receipt of
!any required approvals from the U.S. Dept. of Energy). All rights reserved.
!
!If you have questions about your rights to use or distribute this software,
!please contact Berkeley Lab's Intellectual Property Office at
!IPO@lbl.gov.
!
!NOTICE.  This Software was developed under funding from the U.S. Department
!of Energy and the U.S. Government consequently retains certain rights.  As
!such, the U.S. Government has been granted for itself and others acting on
!its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
!Software to reproduce, distribute copies to the public, prepare derivative
!works, and perform publicly and display publicly, and to permit others to do so.
!
!****************************
! This program was developed by Ji Qiang (jqiang@lbl.gov) of the Lawrence Berkeley National Lab.

!.........................................................................
! Some features include:
! 1) The population size varies from generation to generation.
! 2) The population is uniformly distributed to a number of parallel processors
! 3) The objective function evalution is evaluated by an external serial program.
! Ref: J. Qiang, "A parallel variable population multi-objective optimizer for 
! accelerator beam dynamics optimization," Nuclear Inst. and Methods in Physics Research, 
! A 1054 (2023) 168402.
! The optimization input parameters are in pvpmoo.in.
! The objective function should be modified according to user's application program.
!.........................................................................
!
      program main
      use parallel
      use objfuncs 
      implicit none
      include 'mpif.h'
      real*8, allocatable, dimension(:) :: XCmin,XCmax
      integer, allocatable, dimension(:) :: nrow,ncol
      real*8, dimension(2,4) :: frange
      real*8, dimension(2) :: crange
      real*8, dimension(5) :: rdrange
      integer :: flagchild,ierr
      real*8 :: etam,etae
      real*8 :: t0,t1,dt,dtmax
      !ndim : initial number of the control parameters
      integer ::  ndim=20 
      !nobj : initial number of the objective functions
      integer  ::  nobj=2
      integer :: itermax=5000
      integer :: npini=100
      integer  :: npmin=20,npmax=80,npextmax=500
      integer :: i
!----------------------------------------------
      integer*8 :: iseed
      integer :: flagndom,flaginit,flaginitsamp
      real*8 :: rancheck,probde
 
      call init_parallel()

      open(11,file="pvpmoo.in",status="old")
      !ndim - # of control parameters; nobj - # of objective functions; itermax - maximum # of generations
      read(11,*)ndim,nobj,itermax
      !npini - initial population size; npmin - mini. pop size; npmax - maxi. pop. size; 
      !npextmax - maxi. external storage size
      read(11,*)npini,npmin,npmax,npextmax
      !flagndom - switch of nomdominated sorting (1 for 1st obj. sorted first,  2 for no-order sorting, 3 for
      !last obj. sorted first, otherwise 1st or last obj. randomly sorted.)
      !flaginit - switch for generating initial parents (1 for random sampling, 5 for read-in with re-evalution
      !of objective function values, otherise, just read-in.)
      !flaginitsamp - sample initial population using a random sampling (1) or a quasi-random sampling (2)
      read(11,*)flagndom,flaginit,flaginitsamp
      !iseed - random number seed
      read(11,*)iseed
      !mini. and maxi. ranges of 5 hyper-parameters used in uDE
      read(11,*)frange(1,1),frange(1,2),frange(1,3),frange(1,4),crange(1)
      read(11,*)frange(2,1),frange(2,2),frange(2,3),frange(2,4),crange(2)
      !probablity to keep these 5 hyper-parameters from generation to generation
      read(11,*)rdrange(1:5)
      !flagchild - swith for uDE (1) or GA (2) or with probde probability of uDE 
      !probde - probability of using uDE for next generation children production
      !etam - parameters used in GA
      !etae - parameters used in GA
      read(11,*)flagchild,probde,etam,etae
      
      allocate(nrow(ndim))
      allocate(ncol(ndim))
      allocate(XCmin(ndim))
      allocate(XCmax(ndim))
      !nrow - row line # of the control parameters in the input.in file
      read(11,*)nrow(1:ndim)
      !ncol - column variable # of the control parameters in the input.in file
      read(11,*)ncol(1:ndim)
      !lower bounds of the control parameters
      read(11,*)XCmin(1:ndim)
      !upper bounds of the control parameters
      read(11,*)XCmax(1:ndim)

      close(11)

      if(idproc.eq.0) then
        print*,"-------------------------------------------------------------------------"
        print*,"PVPmoo: Parallel Variable Population Multi-Objective Opitimizer: v.1.0"
        print*,"-------------------------------------------------------------------------"
      endif

      t0 = MPI_WTIME()

      call PVPmoo(objfunc, ndim, nobj, XCmin, XCmax,npini,&
             itermax,mygroupid,idproc,ngroup,ngproc,iseed,&
             npextmax,npmin,npmax,flagndom,flaginit,flaginitsamp,&
             frange,crange,flagchild,etam,etae,nrow,ncol,rdrange,probde)

      t1 = MPI_WTIME()
      dt = t1 - t0
      call MPI_ALLREDUCE(dt,dtmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
              MPI_COMM_WORLD,ierr)

      print*,"time:",dtmax

      deallocate(nrow)
      deallocate(ncol)
      deallocate(XCmin)
      deallocate(XCmax)

      call end_parallel()
 
      end program main

subroutine PVPmoo(objfunc,ndim,nobj,XCmin,XCmax,npini,&
           itermax,myid,myidpr,Ngroupfn,NprocPergrp,iseed,&
           npextmax,npmin,npmax,flagndom,flaginit,flaginitsamp,&
           frange,crange,flagchild,etam,etae,nrow,ncol,rdrange,probde)
!.........................................................................
!                objfunc : The user provided file for evlauting the objective function.
!      subroutine objfunc(nrow(dim),ncol(dim),pop(ndim,np),val(nobj,np),ndim,nobj)
!                      where "pop" is the real control parameter vector.(input)
!                            "val" is the objective functions value.(output)
!             ndim : ndimension of the real control parameters.
!             nobj   : ndimension of the objective functions.
!      XCmin(ndim) : The lower bound of the real control parameters.
!      XCmax(ndim) : The upper bound of the real control parameters.
!                 np : Population size.
!          itermax : The maximum number of iteration.
!         myid       : processor group ID.
!         myidpr     : processor ID in the whole processor group.
!         Ngroupfn   : number of processor groups in the simulation.
!        NprocPergrp : number of processors per group.
!          iseed     : seed for random number generator.

     implicit none
     include 'mpif.h'
     integer, intent(in) :: npini,ndim,nobj,itermax,  &
                        npmin,npmax,npextmax
     integer, intent(in) :: myid,myidpr,Ngroupfn,nprocpergrp,&
                            flagchild
     integer*8, intent(inout) :: iseed
     real*8, dimension(ndim), intent(in) :: XCmin, XCmax
     real*8, intent(in) :: etam,etae,probde
     integer, dimension(ndim), intent(in) :: nrow,ncol
     integer  :: nfevaltot
     !pop, val current parent generation
     real*8, dimension(ndim,npmax) :: pop, bm, mui, mpo,   &
                     pop2, rand, pop1,poptmp1
     real*8, dimension(nobj,npmax) :: val,val2,val1,objtmp1
     !temporary large storage 
     real*8, dimension(ndim,npextmax+npmax) :: poptmp,poptmp2
     !external large storage
     real*8, dimension(ndim,npextmax) :: pop2np
     real*8, dimension(nobj,npextmax) :: val2np
     integer :: i, ibest, iter, j, ig 
     real*8 :: xminj,xmaxj,h
     real*8, dimension(npextmax+npmax) :: tmpval
     real*8, dimension(nobj,npextmax+npmax) :: objvals,objtmp,valtmp2
     real*8, dimension(ndim) :: rand_C1
     integer :: np,nplc,i1,inp,ind1,indd,npcount,ndg,npop,nhalf,itmp,ih
     external  objfunc
     real*8 :: ran22,ran222
     integer :: ierr,ii,npext,npcountmp
     integer*8 :: iseed2
     real*8 :: rtmp
     integer :: iseed0,ilow0,ihigh0,ind2
     integer :: taus,ig1,npopkeep
     logical flag(2)
     real*8 quasi(ndim)
     integer :: ntmp
     integer :: flagndom,flaginit,flaginitsamp
     real*8, dimension(nobj) :: valmin,valmax
     integer :: nrank,k,nkeep
     real*8, dimension(npmax) :: f1np,f2np,f3np,f4np,c1np 
     real*8 :: f1,f2,f3,f4,c1,rdp1,rdp2,rdp3,rdp4,rdp5
     real*8, dimension(2,4) :: frange
     real*8, dimension(2) :: crange
     real*8, dimension(5) :: rdrange
     integer :: ip,nfeval
     real*8 :: rr, tmp1

!!-----Initialize a population --------------------------------------------!!

     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! pop contains the global total population array (size np).
     np = npini
     pop=0.0d0
     nplc = np/Ngroupfn
     np = nplc*Ngroupfn 
     !iseed < 0
     iseed2 = iseed - myid*nplc*5-10*myid
     !iseed2 = -iseed - myid*1.0d0*nplc*5-1.0d0*myid

     if(flaginit.eq.1) then
       !random sampling
       ! only a local fraction of population is given sampled value
       if(flaginitsamp.eq.1) then
        do i=1,nplc
           ig = myid*nplc + i
           do ii = 1, ndim
             rand_C1(ii) = ran222(iseed2)
           enddo
           do j = 1, ndim
             pop(j,ig)=XCmin(j)+rand_C1(j)*(XCmax(j)-XCmin(j))
           enddo
         end do
       !quasi-random sampling
       !initialize the initial population using a quiet start sampling.
       else if(flaginitsamp.eq.2) then
         call inihalt( flag, ndim, np, quasi)
         ig = myid*nplc
         ig1 = (myid+1)*nplc
         do i = 1, np
          call gohalt(ndim,quasi)
          if((i>ig).and.(i<=ig1)) then
           do j = 1, ndim
             pop(j,i)=XCmin(j)+quasi(j)*(XCmax(j)-XCmin(j))
           enddo
          endif
         enddo
       else
         print*,"no initial sampling!!!!"
         stop
       endif
!!--------------------------------------------------------------------------!!

!!------Evaluate fitness functions and find the initial best member-----------------!!
       val=0.0d0
       i1 = 1 + myid*nplc
       inp = myid*nplc+nplc
! Begin counting evaluation calls.
       nfeval=0
! loop through local population size in a global population array.
       do i=i1,inp
          call objfunc(nrow,ncol,pop(:,i),val(:,i),ndim,nobj)
          nfeval=nfeval+1
       enddo
  
       pop2 = 0.0d0
       call MPI_ALLREDUCE(pop,pop2,np*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,&
              MPI_COMM_WORLD,ierr)
! now the pop contains all global population and is the same among
! all processors.
       pop = pop2/NprocPergrp

! sum up local objective values into global popul. obj. values.
       !this causes the objective function value * NprocPergrp
       val2=0.0d0
       call MPI_ALLREDUCE(val,val2,np*nobj,MPI_DOUBLE_PRECISION,MPI_SUM,&
              MPI_COMM_WORLD,ierr) 
       val = val2/NprocPergrp

       !output initial population of control parameters (fort.1) and 
       !corresponding objective function values (fort.2).
       if(myidpr.eq.0) then
         open(1,file="pop.out",status="unknown",position="append")
         open(2,file="obj.out",status="unknown",position="append")
         do i = 1, np
           write(1,222)pop(1:ndim,i)
         enddo
         call flush(1)
         do i = 1, np
           write(2,333)val(1:nobj,i)
         enddo
         call flush(2)
         close(1)
         close(2)
         call system("cp pop.out obj.out ../")
       endif

!redefine the parent popluation size from the initial sampling
!using non-dominate sorting

       pop2 = pop
       val2 = val
       ntmp = npextmax+npmax

       !find the initial non-dominated group
       !do non-domianted sorting and find the next generation population
       npop = np
       npcount = 0 !counter of how many population
       ndg = 0 !counter of how many non-dominated selection call
       ind1 = 0 !counter of the rank 1 non-dominated solution
       ind2 = 0 !counter of the rank 2 non-dominated solution
       !only 1st non-dominate sorted groups or so are used as parents
       !do non-dominate sorting and save the sorted results in pop2np and val2np
       do while(ndg<2 .or. npcount<npmin)
          if(flagndom.eq.1) then
            call ndomgroup(pop,ndim,val,nobj,npop,pop2,val2,indd)
          else if(flagndom.eq.2) then
            call ndomgroup2(pop,ndim,val,nobj,npop,pop2,val2,indd)
          else if(flagndom.eq.3) then
            call ndomgroup3(pop,ndim,val,nobj,npop,pop2,val2,indd)
          else
            rtmp = ran22(iseed)
            if(rtmp.lt.0.5d0) then
              call ndomgroup(pop,ndim,val,nobj,npop,pop2,val2,indd)
            else
              call ndomgroup3(pop,ndim,val,nobj,npop,pop2,val2,indd)
            endif
          endif

          ndg = ndg + 1
          do i = npcount + 1, npcount+indd
            ii = i - npcount
            pop2np(:,i) = pop2(:,ii)
            val2np(:,i) = val2(:,ii)
          enddo
          npcount = npcount + indd
          if(ndg.eq.1) ind1 = indd
          if(ndg.eq.2) ind2 = indd
 
          do i = indd+1,npop
            ii = i - indd
            pop(:,ii) = pop2(:,i)
            val(:,ii) = val2(:,i)
          enddo
          !reduce the sorting range
          npop = npop-indd
          !!set the penalty number
          !do i = npop-indd+1,npop
          !  val(:,i) = 1.0e5+10*i
          !enddo
       enddo

       !define the new parent population (has to be multiple of Ngroupfn).
       nplc = npcount/Ngroupfn
       np = Ngroupfn*nplc
       do i = 1, np
            pop(:,i) = pop2np(:,i)
            val(:,i) = val2np(:,i)
       enddo

       npext = np
  
       !the best value will be randomly chosen from the non-dominated solution
       bm=0.0d0
       !for non-dominated solution group, best is itself
       do i = 1,ind1 
         bm(:,i) = pop2np(:,i)
       enddo
       do i = ind1+1, np
         itmp = ran22(iseed)*ind1+1
         bm(:,i) = pop2np(:,itmp)
       enddo
     !!--------------------------------------------------------------------------!!
     else
       !initialize initial population from read-in
       !for restart, one needs: Nptot,npext,ind1
       !popXC: contain 1st Nptot of the external archive
       !pop2np:external archive contains non-dom. sorted pop.
       !if one wants to increase the number of population size,
       !one set set np = npext
       open(12,file="pop.in",status="old")
       read(12,*)np,npext,ind1
       do i = 1, npext
          read(12,*)pop2np(1:ndim,i)
       enddo
       close(12)
       open(11,file="obj.in",status="old")
       do i = 1, npext
          read(11,*)val2np(1:nobj,i)
       enddo
       close(11)
       nplc = np/Ngroupfn
       npcount = nplc*Ngroupfn
       ind2 = 0
       np = npcount
       npext = npcount
       ind1 = min(np,ind1)

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       !the following is needed when we do a restart simulation but
       !changing the resolution of the simulations
       !the objective function values from the read-in will not be right
       if(flaginit.eq.5) then !update the obj. values using the new input contrl. param.
          val=0.0
          val2=0.0
          nfeval=0
          i1 = 1 + myid*nplc
          inp = myid*nplc+nplc
          do i=i1,inp
             call objfunc(nrow,ncol,pop2np(:,i),val(:,i),ndim,nobj)
             nfeval=nfeval + 1
          enddo

          ! sum up local objective values into global popul. obj. values.
          !this causes the objective function value * NprocPergrp
          call MPI_ALLREDUCE(val,val2,npmax*nobj,MPI_DOUBLE_PRECISION,MPI_SUM,&
              MPI_COMM_WORLD,ierr)
          val2np = val2/NprocPergrp

          !output initial population of control parameters (fort.1) and
          !corresponding objective function values (fort.2).
          if(myidpr.eq.0) then
            open(1,file="pop.out",status="unknown",position="append")
            open(2,file="obj.out",status="unknown",position="append")
            do i = 1, np
              write(1,222)pop2np(1:ndim,i)
            enddo
            do i = 1, np
              write(2,333)val2np(1:nobj,i)
            enddo
            call flush(1)
            call flush(2)
            close(1)
            close(2)
            call system("cp pop.out obj.out ../")
          endif
       endif
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     endif
!!--------------------------------------------------------------------------!!

     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     rdp1 = rdrange(1)
     rdp2 = rdrange(2)
     rdp3 = rdrange(3)
     rdp4 = rdrange(4)
     rdp5 = rdrange(5)
     ! generate initial individual control parameter from random range.
     call genFCnp(f1np,f2np,f3np,f4np,c1np,frange,crange,iseed,npmax)
     iter=0 
!!--------------------------------------------------------------------------!!
!!------Perform evolutionary iteration----------------------!! 
!!--------------------------------------------------------------------------!!
     nfevaltot = 0
     call MPI_ALLREDUCE(nfeval,nfevaltot,1,MPI_INTEGER,MPI_SUM,&
              MPI_COMM_WORLD,ierr)
     if(myidpr.eq.0) print*,"it: ",iter,npext,ind1,nfevaltot
     do while (iter < itermax)
        pop2=pop
        val2 = val

        !choose the 5 control parameters
        do ip = 1,np
          rr = ran22(iseed)
          if(rr.lt.rdp1) then !generate the new parameters.
            tmp1 = ran22(iseed)
            f1 = frange(1,1)+tmp1*(frange(2,1)-frange(1,1))
          else
            f1 = f1np(ip)
          endif
          rr = ran22(iseed)
          if(rr.lt.rdp2) then !generate the new parameters.
            tmp1 = ran22(iseed)
            f2 = frange(1,2)+tmp1*(frange(2,2)-frange(1,2))
          else
            f2 = f2np(ip)
          endif
          rr = ran22(iseed)
          if(rr.lt.rdp3) then !generate the new parameters.
            tmp1 = ran22(iseed)
            f3 = frange(1,3)+tmp1*(frange(2,3)-frange(1,3))
          else
            f3 = f3np(ip)
          endif
          rr = ran22(iseed)
          if(rr.lt.rdp4) then !generate the new parameters.
            tmp1 = ran22(iseed)
            f4 = frange(1,4)+tmp1*(frange(2,4)-frange(1,4))
          else
            f4 = f4np(ip)
          endif
          rr = ran22(iseed)
          if(rr.lt.rdp5) then !generate the new parameters.
            tmp1 = ran22(iseed)
            c1 = crange(1)+tmp1*(crange(2)-crange(1))
          else
            c1 = c1np(ip)
          endif
 
          f1np(ip) = f1
          f2np(ip) = f2
          f3np(ip) = f3
          f4np(ip) = f4
          c1np(ip) = c1
        enddo

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !produce np children using np parents following the DE algorithm with immediate replacement
        if(flagchild.eq.1) then
          call prodchildrenDE(np,nobj,ndim,nprocpergrp,nplc,myid,Ngroupfn,&
                   iseed,iseed2,nfeval,f1np,f2np,f3np,f4np,c1np,&
                   XCmin,XCmax,bm,pop2,val2,pop1,val1,objfunc,npmax,nrow,ncol) 
        !produce np children using np parents following the GA algorithm.
        else if(flagchild.eq.2) then
          call prodchildrenGA(np,nobj,ndim,nprocpergrp,nplc,myid,Ngroupfn,&
                   iseed,iseed2,nfeval,etam,etae,&
                   XCmin,XCmax,pop2,val2,pop1,val1,objfunc,npmax,nrow,ncol) 
        !random selection of DE and GA with "probde" probability of DE.
        else
          rr = ran22(iseed)
          if(rr.lt.probde) then !generate the new parameters.
            call prodchildrenDE(np,nobj,ndim,nprocpergrp,nplc,myid,Ngroupfn,&
                   iseed,iseed2,nfeval,f1np,f2np,f3np,f4np,c1np,&
                   XCmin,XCmax,bm,pop2,val2,pop1,val1,objfunc,npmax,nrow,ncol) 
          else
            call prodchildrenGA(np,nobj,ndim,nprocpergrp,nplc,myid,Ngroupfn,&
                   iseed,iseed2,nfeval,etam,etae,&
                   XCmin,XCmax,pop2,val2,pop1,val1,objfunc,npmax,nrow,ncol) 
          endif
        endif

        !combine the new trial solution pop1, val1 with the parent solution in external storage
        !pop2np,val2np, into a np+npext array for non-dominated sorting
        do i = 1, np
          pop(:,i) = pop1(:,i)
          objvals(:,i) = val1(:,i)
        enddo
        do i = np+1, np+npext
          pop(:,i) = pop2np(:,i-np)
          objvals(:,i) = val2np(:,i-np)
        enddo
        !do non-domianted sorting and find the next generation population
        !use 3 temporary storage: pop,poptmp,poptmp2
        !here, poptmp2 contains the solution after sorting
        !here, the non-dominated sorting is done using the np children +
        !npext external stored parents
        npop = np + npext
        npopkeep = npop
        npcount = 0 !counter of how many population
        ndg = 0 !counter of how many non-dominated selection call
        ind1 = 0
        ind2 = 0
        nrank = 1
        !do while(ndg<1 .or. npcount<npmin)
        do while(ndg<nrank .or. npcount<npmin)
        !do a complete sorting
        !do while(npop.gt.1)
          if(flagndom.eq.1) then
            call ndomgroup(pop,ndim,objvals,nobj,npop,poptmp,objtmp,indd)
          else if(flagndom.eq.2) then
            call ndomgroup2(pop,ndim,objvals,nobj,npop,poptmp,objtmp,indd)
          else if(flagndom.eq.3) then
            call ndomgroup3(pop,ndim,objvals,nobj,npop,poptmp,objtmp,indd)
          else
            rtmp = ran22(iseed)
            if(rtmp.lt.0.5d0) then
              call ndomgroup(pop,ndim,val,nobj,npop,pop2,val2,indd)
            else
              call ndomgroup3(pop,ndim,val,nobj,npop,pop2,val2,indd)
            endif
          endif
          ndg = ndg + 1
          !do i = npcount + 1, npcount+indd
          do i = npcount + 1, npopkeep
              ii = i - npcount
              poptmp2(:,i) = poptmp(:,ii)
              valtmp2(:,i) = objtmp(:,ii)
          enddo
          npcount = npcount + indd
          if(ndg.eq.1) ind1 = indd
          if(ndg.eq.2) ind2 = indd

          do i = indd+1,npop
            ii = i - indd
            pop(:,ii) = poptmp(:,i)
            objvals(:,ii) = objtmp(:,i)
          enddo
          npop = npop - indd
          !!set the penalty number
          !do i = npop-indd+1,npop
          !  objvals(:,i) = 1.0e5+10*i
          !enddo
        enddo

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !here min(npmax,npcount) is used as the number of next generation parents
        !to generate np next generation children.
        if(npcount.le.npmax) then
            nplc = npcount/Ngroupfn
        else
            nplc = npmax/Ngroupfn
        endif
        np = Ngroupfn*nplc
        
        !select np best solutions (non-dominated solutions)
        if(ind1.le.npmax) then
          !for non-dominated solution group, best is itself
          do i = 1, ind1
            bm(:,i) = poptmp2(:,i)
          enddo
          !otherwise, random selection from non-dominated group
          do i = ind1+1, np
            itmp = ran22(iseed)*ind1+1
            bm(:,i) = poptmp2(:,itmp)
          enddo
        else
          if(np>ind1/2) then !remove ind1-np from ind1
             call seldshm2(poptmp2,ndim,valtmp2,nobj,ind1,&
                       bm,val,ind1-np,npmax,ntmp)
          else !select np out ind1
            call seldshm3(poptmp2,ndim,valtmp2,nobj,ind1,&
                       bm,val,np,npmax,ntmp)
          endif
        endif

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !select np parents solutions from the npcount sorted solutions.
        if(npcount.le.npmax) then
          do i = 1, np
            pop(:,i) = poptmp2(:,i)
            val(:,i) = valtmp2(:,i)
          enddo
        else
          if(np>npcount/2) then !remove npcount-np from npcount
            call seldshm2(poptmp2,ndim,valtmp2,nobj,npcount,&
                       pop,val,npcount-np,npmax,ntmp)
          else
            call seldshm3(poptmp2,ndim,valtmp2,nobj,npcount,&
                       pop,val,np,npmax,ntmp)
          endif
        endif

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !fill in the external storage (pop2np,val2np) using the non-dominated 
        !sorted parent+children solutions (poptmp2,valtmp2)(min(npcount,npextmax)).
        !If we do only 1 non-dominated sorting, npcount = ind1
        !if ind1< npmin, some dominated solutions will also be put into the storage.
        if(ind1.le.npextmax) then
          if(npcount.le.npextmax) then
            npext = npcount
          else
            npext = npextmax
          endif
          do i = 1, npext
            pop2np(:,i) = poptmp2(:,i)
            val2np(:,i) = valtmp2(:,i)
          enddo
        else
            npext = npextmax
            !select npextmax from ind1 non-dominated solution
            !the following algorithm is try to keep the new non-dominated generation 
            !as diverse as possible
            nkeep = 0
             if(npext>ind1/2) then !remove ind1-npext from ind1
               call seldshm2(poptmp2,ndim,valtmp2,nobj,ind1,&
                      pop2np,val2np,ind1-npext,npextmax,ntmp)
             else
               call seldshm3(poptmp2,ndim,valtmp2,nobj,ind1,&
                      pop2np,val2np,npext,npextmax,ntmp)
             endif
        endif

        iter=iter+1
!        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!output population of control parameters (fort.1) and 
!corresponding objective function values (fort.2) in each generation.
        call MPI_ALLREDUCE(nfeval,nfevaltot,1,MPI_INTEGER,MPI_SUM,&
              MPI_COMM_WORLD,ierr)
        if(myidpr.eq.0) print*,"it: ",iter,npext,ind1,nfevaltot

        if(myidpr.eq.0) then
          open(1,file="pop.out",status="unknown",position="append")
          open(2,file="obj.out",status="unknown",position="append")
          write(1,*)iter,npext,ind1
          do i = 1, npext
            write(1,222)pop2np(1:ndim,i)
          enddo
          write(2,710)iter,npext,ind1,npopkeep,nfevaltot
          do i = 1, npext
            write(2,333)val2np(1:nobj,i)
          enddo
          call flush(1)
          call flush(2)
          close(1)
          close(2)
          call system("cp pop.out obj.out ../")
        endif
 
     end do

!!------end the evolutionary computation------------------------------!!
222 format(1x, 50(1x,e20.10))
333 format(1x, 10(e18.9))
710 format(1x, 5I8)

end subroutine PVPmoo

function randperm(num,iseed)
  implicit none
  integer*8,intent(inout) :: iseed
  integer, intent(in) :: num
  integer :: number, i, j, k
  integer, dimension(num) :: randperm
  real*8, dimension(num) :: rand2
  intrinsic random_number
  real*8 :: ran22
  !call random_number(rand2)
  do i = 1, num
    rand2(i) = ran22(iseed)
  enddo
  do i=1,num
     number=1
     do j=1,num
        if (rand2(i) > rand2(j)) then
	       number=number+1
        end if
     end do
     do k=1,i-1
        if (rand2(i) <= rand2(k) .and. rand2(i) >= rand2(k)) then
	       number=number+1
        end if
     end do
     randperm(i)=number
  end do
  return
end function randperm

      FUNCTION ran22(idum)
      integer*8 :: idum
      real*8 ran22
      real*8 :: rtmp
        call random_number(rtmp)
        ran22 = rtmp
      END Function ran22

      FUNCTION ran222(idum)
      INTEGER*8 idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 ran222,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.d0/IM1,IMM1=IM1-1,&
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
      NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-7,RNMX=1.d0-EPS)
      INTEGER*8 idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran222=min(AM*iy,RNMX)
      return
      END Function ran222

      !find a group of non-dominated solution (poptmp,objtmp,ind)
      !from the input solution (pop,objvals,npop)
      subroutine ndomgroup(pop,ndim,objvals,nobj,npop,poptmp,objtmp,ind)
      implicit none
      integer :: ndim,npop,nobj
      real*8, dimension(ndim,npop) :: pop,poptmp
      real*8, dimension(nobj,npop) :: objvals,objtmp
      integer :: i,j,k,ind,icc,ii
      integer, dimension(npop) :: iwksp
      real*8, dimension(npop) :: wksp,ra
      logical, dimension(npop) :: ndom
      real*8 :: eps

      eps = 1.0d-12

      !rearrage the pop and objective function value
      !using the ascending sorting of the 1st objective row
      ra(:) = objvals(1,1:npop)
      call indexx(npop,ra,iwksp)
      wksp = ra
      do i = 1, npop
          ra(i) = wksp(iwksp(i))
      enddo
      objvals(1,1:npop) = ra
      do j = 2, nobj
        wksp = objvals(j,1:npop)
        do i = 1, npop
          ra(i) = wksp(iwksp(i))
        enddo
        objvals(j,:) = ra
      enddo 
      do j = 1, ndim
        wksp = pop(j,1:npop)
        do i = 1, npop
          ra(i) = wksp(iwksp(i))
        enddo
        pop(j,1:npop) = ra
      enddo

      ind = 1
      !1st one after sorting is always a nondominated solution
      poptmp = 0.0d0
      objtmp = 0.0d0
!      iwksp = 0
      ndom = .True.
      poptmp(:,ind) = pop(:,ind) 
      objtmp(:,ind) = objvals(:,ind)
!      iwksp(ind) = 1
      ndom(ind) = .False.
      do i = 2, npop
        do j = 2, nobj
          icc = 0
          do k = 1, i-1 !?can the operation of this loop be constant
            !print*,"i,j,k:",i,j,k,objvals(j,i),objvals(j,k)
            !if(iwksp(k).eq.0) then
            if(ndom(k)) then
              icc = icc + 1
            else if(objvals(j,i).lt.objvals(j,k)) then
            !else if((-objvals(j,i)+objvals(j,k)).gt.eps) then
              !i dominates k in objective function j 
              icc = icc + 1
            else
              exit
            endif
          enddo
          !print*,"icc: ",icc,i-1
          !i is not dominated by previous solutions.
          if(icc.eq.(i-1)) then !find a non-dominated solution
            !i dominates 1 -> i-1 in objective function j
            ind = ind + 1
            !print*,"ind: ",ind,icc,i
            !copy the non-dominated solutions into poptmp,objtmp
            poptmp(:,ind) = pop(:,i)
            objtmp(:,ind) = objvals(:,i)
            !iwksp(i) = 1
            ndom(i) = .False.
            exit
          endif
        enddo
      enddo
      
      ii = 0
      do i = 1, npop
        !if(iwksp(i).eq.0) then
        if(ndom(i)) then
          ii = ii + 1
          !print*,"ii: ",i,ii,ind+ii
          poptmp(:,ind+ii) = pop(:,i)
          objtmp(:,ind+ii) = objvals(:,i)
        endif
      enddo

      end subroutine ndomgroup

      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      REAL*8 arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL*8 a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END

!produce np children using np parents following the uDE algorithm with 
!immediate replacement.
!"iseed" is the global random seed. "iseed2" is the local random seed.
     subroutine prodchildrenDE(np,nobj,ndim,nprocpergrp,nplc,myid,Ngroupfn,iseed,&
                iseed2,nfeval,F1,F2,F3,F4,CR,XCmin,XCmax,bm,popold,&
                val,ui,val1,objfunc,npmax,nrow,ncol)
     implicit none
     include 'mpif.h'
     integer, intent(in) :: np,ndim,nobj,npmax
     integer, intent(in) :: myid,nplc,nprocpergrp,Ngroupfn
     integer*8, intent(inout) :: iseed,iseed2
     integer, intent(inout) :: nfeval
     integer, dimension(ndim), intent(in) :: nrow,ncol
     real*8, dimension(ndim,np) :: bm, mui, mpo,   &
                              popold,ui
     real*8, dimension(nobj,np) :: val,val1,val2
     real*8, dimension(nobj+ndim,np) :: tmp1,tmp2
     real*8, dimension(ndim) :: XCmin,XCmax
     real*8, dimension(npmax) :: F1,F2,F3,F4
     real*8, dimension(npmax) :: CR
     integer, dimension(np) :: rot, a1, a2, a3, a4, a5, rt
     external objfunc
     real*8 :: ran22,ran222
     integer :: ip,ii,i,j,flagdom,ierr,nsend,ndim1,ndimobj
     real*8 :: rtmp,rnd
     integer :: jj,itmp,ll
     integer :: irand
     real*8 :: tmp0
     integer :: npexplore
     real*8 :: fracexplore,rr

     nsend = (ndim+nobj)*np
     ndim1 = ndim+1
     ndimobj = ndim+nobj

     !here, np > 5
     if(np<5) then
       print*,"np is too small: ",np
       stop
     endif
     do i = 1, np
11     itmp = ran22(iseed)*np+1
       if(itmp.eq.i) goto 11
       a1(i) = itmp

12     itmp = ran22(iseed)*np+1
       if(itmp.eq.i .or. itmp.eq.a1(i)) goto 12
       a2(i) = itmp

13     itmp = ran22(iseed)*np+1
       if(itmp.eq.i .or. itmp.eq.a1(i) .or. itmp.eq.a2(i)) goto 13
       a3(i) = itmp

14     itmp = ran22(iseed)*np+1
       if(itmp.eq.i .or. itmp.eq.a1(i) .or. itmp.eq.a2(i) &
          .or. itmp.eq.a3(i)) goto 14
       a4(i) = itmp

15     itmp = ran22(iseed)*np+1
       if(itmp.eq.i .or. itmp.eq.a1(i) .or. itmp.eq.a2(i) &
          .or. itmp.eq.a3(i) .or. itmp.eq.a4(i)) goto 15
       a5(i) = itmp
     enddo

!!------for Crossover operation-------------------------------------------------!!
        mui=0.0d0
        mpo=0.0d0
        do j = 1, np
          tmp0 = ran22(iseed)
          do i = 1, ndim
            rr = ran22(iseed)
            irand = ndim*tmp0+1
            if((rr < CR(j)) .or. (i.eq.irand)) then
              mui(i,j) = 1.0d0
            else
              mpo(i,j) = 1.0d0
            endif
          enddo
        enddo

!---- select a mutation strategy-----------------------------------------!
        ui = 0.0d0
        val1 = 0.0d0
        do i=1,nplc 
           ip = myid*nplc + i

           ui(:,ip)=popold(:,ip)+F1(ip)*(bm(:,ip)-popold(:,ip)) &
                +F2(ip)*(popold(:,a5(ip))-popold(:,ip)) &
                +F3(ip)*(popold(:,a1(ip))-popold(:,a2(ip))) &
                +F4(ip)*(popold(:,a3(ip))-popold(:,a4(ip)))

!!------   do Crossover operation-------------------------------------------------!!
     	   ui(:,ip)=popold(:,ip)*mpo(:,ip)+ui(:,ip)*mui(:,ip)

! check control parameter range
           do ii = 1, ndim
               if(ui(ii,ip).lt.XCmin(ii)) then
                 rtmp = ran222(iseed2)
                 ui(ii,ip) = XCmin(ii) + rtmp*(XCmax(ii)-XCmin(ii))
               else if(ui(ii,ip).gt.XCmax(ii)) then
                 rtmp = ran222(iseed2)
                 ui(ii,ip) = XCmin(ii) + rtmp*(XCmax(ii)-XCmin(ii))
               else
               endif
           enddo

!!--------------------------------------------------------------------------!!
!!------Evaluate new objective functions -----------------!!
           call objfunc(nrow,ncol,ui(:,ip),val1(:,ip),ndim,nobj)
           nfeval=nfeval+1
!check val1(ip) against val(ip)
           flagdom = 1
           do j = 1, nobj
             if(val1(j,ip).gt.val(j,ip)) flagdom=0 
           enddo
           if(flagdom.eq.1) then !val1 dominate val,immediate replacement
             popold(:,ip) = ui(:,ip)
             val(:,ip) = val1(:,ip)
           endif
!use allreduce to update the whole popold with the new non-dominated solution.
           tmp1 = 0.0d0
           tmp2 = 0.0d0
           tmp1(1:ndim,ip) = popold(1:ndim,ip)
           tmp1(ndim1:ndimobj,ip) = val(1:nobj,ip)
           call MPI_ALLREDUCE(tmp1,tmp2,nsend,MPI_DOUBLE_PRECISION,MPI_SUM,&
              MPI_COMM_WORLD,ierr)
           tmp2 = tmp2/NprocPergrp
           do j = 1, Ngroupfn
             ii = (j-1)*nplc+i
             popold(1:ndim,ii) = tmp2(1:ndim,ii)
             val(1:nobj,ii) = tmp2(ndim1:ndimobj,ii)
           enddo
        enddo

!!--------------------------------------------------------------------------!!
        popold = 0.0d0
        val2 = 0.0d0
        tmp1 = 0.0d0
        tmp2 = 0.0d0
        do i = 1, nplc
          ip = myid*nplc + i
          tmp1(1:ndim,ip) = ui(1:ndim,ip)
          tmp1(ndim1:ndimobj,ip) = val1(1:nobj,ip)
        enddo
        call MPI_ALLREDUCE(tmp1,tmp2,nsend,MPI_DOUBLE_PRECISION,MPI_SUM,&
              MPI_COMM_WORLD,ierr)
        tmp1 = tmp2/NprocPergrp
        do i = 1, np
          ui(1:ndim,i) = tmp1(1:ndim,i)
          val1(1:nobj,i) = tmp1(ndim1:ndimobj,i)
        enddo

     end subroutine prodchildrenDE

!produce np children using np parents following the real code GA algorithm 
!"iseed" is the global random seed. "iseed2" is the local random seed.
     subroutine prodchildrenGA(np,nobj,ndim,nprocpergrp,nplc,myid,Ngroupfn,iseed,&
                iseed2,nfeval,etam,etae,XCmin,XCmax,popold,&
                val,ui,val1,objfunc,npmax,nrow,ncol)
     implicit none
     include 'mpif.h'
     integer, intent(in) :: np,ndim,nobj,npmax
     integer, intent(in) :: myid,nplc,nprocpergrp,Ngroupfn
     integer*8, intent(inout) :: iseed,iseed2
     integer, intent(inout) :: nfeval
     integer, dimension(ndim), intent(in) :: nrow,ncol
     real*8, dimension(ndim,np) :: popold,ui
     real*8, dimension(nobj,np) :: val,val1
     real*8, dimension(nobj+ndim,np) :: tmp1,tmp2
     real*8, dimension(ndim) :: XCmin,XCmax
     real*8, intent(in) :: etam,etae
     external objfunc
     real*8 :: ran22,ran222
     integer :: ip,ii,i,j,flagdom,ierr,nsend,ndim1,ndimobj
     integer :: jj,itmp,ll
     integer :: irand
     real*8 :: delta1,delta2,delta3,delta4,delta,betak,etapm,etape,eta
     real*8 :: bb, bb2,rr,deltak,rtmp
     integer :: i1,i2,nr,nhalf,nleft

     nsend = (ndim+nobj)*np
     ndim1 = ndim+1
     ndimobj = ndim+nobj

! etam is usually > 5, a typical value can be 5
! etae > 1 & <5, a typical value can be 1
!
!!--------------------------------------------------------------------------!!
        ui = 0.0d0
        val1 = 0.0d0


        etapm = 1.0d0/(etam+1.0d0)
        etape = 1.0d0/(etae+1.0d0)

!!------   do crossover operation--------------------------------------!!
        !nplc has to be even
        nhalf = nplc/2
        nleft = nplc-2*nhalf
        do i=1, nhalf
          i1 = myid*nplc+2*i-1
          i2 = myid*nplc+2*i
          do ii = 1, ndim
            delta1 = XCmax(ii) - popold(ii,i1)
            delta2 = XCmax(ii) - popold(ii,i2)
            delta3 = popold(ii,i1)-XCmin(ii)
            delta4 = popold(ii,i2)-XCmin(ii)
            delta = min(delta1,delta2,delta3,delta4) 
            bb = 1.0d0+2*delta/abs(popold(ii,i1)-popold(ii,i2))
            bb2 = 2.0d0-1.0d0/bb**(1+etae)
            rr = ran222(iseed2)
            if(rr<=0.5d0) then
                 betak = (2*rr)**etape
            else
                 !betak = (bb2*(1.0d0-rr))**(-etape)
                 betak = (2*(1.0d0-rr))**(-etape)
            endif
            !ui(ii,i1) = 0.5d0*(popold(ii,i1)+popold(ii,i2))- &
            !      0.5d0*betak*abs(popold(ii,i1)-popold(ii,i2)) 
            !ui(ii,i2) = 0.5d0*(popold(ii,i1)+popold(ii,i2))+ &
            !      0.5d0*betak*abs(popold(ii,i1)-popold(ii,i2)) 
            ui(ii,i1) = 0.5d0*(popold(ii,i1)+popold(ii,i2))+ &
                  0.5d0*betak*(popold(ii,i2)-popold(ii,i1)) 
            ui(ii,i2) = 0.5d0*(popold(ii,i1)+popold(ii,i2))- &
                  0.5d0*betak*(popold(ii,i2)-popold(ii,i1)) 
          enddo
        enddo
        !for the left single solution if nplc is odd
        do i=1, nleft
          i1 = myid*nplc+2*nhalf+i
          !random parent solution from the other 1:2*nhalf
          nr = 2*nhalf*ran222(iseed2)
          if(nr.eq.0) nr = 1
          i2 = myid*nplc+nr
          do ii = 1, ndim
            delta1 = XCmax(ii) - popold(ii,i1)
            delta2 = XCmax(ii) - popold(ii,i2)
            delta3 = popold(ii,i1)-XCmin(ii)
            delta4 = popold(ii,i2)-XCmin(ii)
            delta = min(delta1,delta2,delta3,delta4)
            bb = 1.0d0+2*delta/abs(popold(ii,i1)-popold(ii,i2))
            bb2 = 2.0d0-1.0d0/bb**(1+etae)
            rr = ran222(iseed2)
            if(rr<0.5d0) then
                 betak = (2*rr)**etape
            else
                 !betak = (bb2*(1.0d0-rr))**(-etape)
                 betak = (2*(1.0d0-rr))**(-etape)
            endif
            !ui(ii,i1) = 0.5d0*(popold(ii,i1)+popold(ii,i2))- &
            !      0.5d0*betak*abs(popold(ii,i1)-popold(ii,i2))
            ui(ii,i1) = 0.5d0*(popold(ii,i1)+popold(ii,i2))+ &
                  0.5d0*betak*(popold(ii,i2)-popold(ii,i1))
          enddo
        enddo

!!------   do mutation operation--------------------------------------!!
        do i=1,nplc 
           ip = myid*nplc + i
           !polynomial mutataion
           do ii = 1, ndim
              delta1 = XCmax(ii) - ui(ii,ip)
              delta2 = ui(ii,ip)-XCmin(ii)
              !delta = min(delta1,delta2) 
              delta = XCmax(ii)-XCmin(ii)
              rr = ran222(iseed2)
              if(rr<0.5d0) then
                 deltak = (2*rr)**etapm-1.0d0
              else
                 deltak = 1.0d0-(2*(1.0d0-rr))**etapm
              endif
              ui(ii,ip) = ui(ii,ip) + deltak*delta
           enddo
        enddo

        ! check constraint range
        do i=1,nplc 
           ip = myid*nplc + i
           do ii = 1, ndim
               if(ui(ii,ip).lt.XCmin(ii)) then
                 rtmp = ran222(iseed2)
                 ui(ii,ip) = XCmin(ii) + rtmp*(XCmax(ii)-XCmin(ii))
               else if(ui(ii,ip).gt.XCmax(ii)) then
                 rtmp = ran222(iseed2)
                 ui(ii,ip) = XCmin(ii) + rtmp*(XCmax(ii)-XCmin(ii))
               else
               endif
           enddo
        enddo

        do i=1,nplc 
           ip = myid*nplc + i
           call objfunc(nrow,ncol,ui(:,ip),val1(:,ip),ndim,nobj)
           nfeval=nfeval+1
        enddo

!!--------------------------------------------------------------------------!!
        popold = 0.0d0
        tmp1 = 0.0d0
        tmp2 = 0.0d0
        do i = 1, nplc
          ip = myid*nplc + i
          tmp1(1:ndim,ip) = ui(1:ndim,ip)
          tmp1(ndim1:ndimobj,ip) = val1(1:nobj,ip)
        enddo
        call MPI_ALLREDUCE(tmp1,tmp2,nsend,MPI_DOUBLE_PRECISION,MPI_SUM,&
              MPI_COMM_WORLD,ierr)
        tmp1 = tmp2/NprocPergrp
        do i = 1, np
          ui(1:ndim,i) = tmp1(1:ndim,i)
          val1(1:nobj,i) = tmp1(ndim1:ndimobj,i)
        enddo

     end subroutine prodchildrenGA

      !find a group of non-dominated solution (poptmp,objtmp,ind)
      !from the input solution (pop,objvals,npop)
      subroutine ndomgroup3(pop,ndim,objvals,nobj,npop,poptmp,objtmp,ind)
      implicit none
      integer :: ndim,npop,nobj
      real*8, dimension(ndim,npop) :: pop,poptmp
      real*8, dimension(nobj,npop) :: objvals,objtmp
      integer :: i,j,k,ind,icc,ii
      integer, dimension(npop) :: iwksp
      real*8, dimension(npop) :: wksp,ra
      logical, dimension(npop) :: ndom
      real*8 :: eps

      eps = 1.0d-12

      !rearrage the pop and objective function value
      !using the ascending sorting of the last objective row
      ra(:) = objvals(nobj,1:npop)
      call indexx(npop,ra,iwksp)
      wksp = ra
      do i = 1, npop
          ra(i) = wksp(iwksp(i))
      enddo
      objvals(nobj,1:npop) = ra
      do j = 1, nobj-1
        wksp = objvals(j,1:npop)
        do i = 1, npop
          ra(i) = wksp(iwksp(i))
        enddo
        objvals(j,1:npop) = ra
      enddo 
      do j = 1, ndim
        wksp = pop(j,1:npop)
        do i = 1, npop
          ra(i) = wksp(iwksp(i))
        enddo
        pop(j,1:npop) = ra
      enddo

      ind = 1
      !1st one after sorting is always a nondominated solution
      poptmp = 0.0d0
      objtmp = 0.0d0
!      iwksp = 0
      ndom = .True.
      poptmp(:,ind) = pop(:,ind) 
      objtmp(:,ind) = objvals(:,ind)
!      iwksp(ind) = 1
      ndom(ind) = .False.
      do i = 2, npop
        do j = 1, nobj-1
          icc = 0
          do k = 1, i-1
            !print*,"i,j,k:",i,j,k,objvals(j,i),objvals(j,k)
            !if(iwksp(k).eq.0) then
            if(ndom(k)) then
              icc = icc + 1
            else if(objvals(j,i).lt.objvals(j,k)) then
            !else if((-objvals(j,i)+objvals(j,k)).gt.eps) then
              !i dominates k in objective function j 
              icc = icc + 1
            else
              exit
            endif
          enddo
          !print*,"icc: ",icc,i-1
          if(icc.eq.(i-1)) then !find a non-dominated solution
            !i dominates 1 -> i-1 in objective function j
            ind = ind + 1
            !print*,"ind: ",ind,icc,i
            !copy the non-dominated solutions into poptmp,objtmp
            poptmp(:,ind) = pop(:,i)
            objtmp(:,ind) = objvals(:,i)
            !iwksp(i) = 1
            ndom(i) = .False.
            exit
          endif
        enddo
      enddo
      
      ii = 0
      do i = 1, npop
        !if(iwksp(i).eq.0) then
        if(ndom(i)) then
          ii = ii + 1
          !print*,"ii: ",i,ii,ind+ii
          poptmp(:,ind+ii) = pop(:,i)
          objtmp(:,ind+ii) = objvals(:,i)
        endif
      enddo

      end subroutine ndomgroup3

      !find a group of non-dominated solution (poptmp,objtmp,ind)
      !from the input solution (pop,objvals,npop)
      subroutine ndomgroup2(pop,ndim,objvals,nobj,npop,poptmp,objtmp,ind)
      implicit none
      integer :: ndim,npop,nobj
      real*8, dimension(ndim,npop) :: pop,poptmp
      real*8, dimension(nobj,npop) :: objvals,objtmp
      integer :: i,j,k,ind,icc,ii,idom
      integer, dimension(npop) :: iwksp
      real*8, dimension(npop) :: wksp,ra
      logical, dimension(npop) :: ndom
      logical :: tmpdom,tmpdom2
      real*8 :: eps

      eps = 1.0d-12

      !rearrage the pop and objective function value
      !using the ascending sorting of the 1st objective row
!      ra(:) = objvals(1,:)
!      call indexx(npop,ra,iwksp)
!      wksp = ra
!      do i = 1, npop
!          ra(i) = wksp(iwksp(i))
!      enddo
!      objvals(1,:) = ra
!      do j = 2, nobj
!        wksp = objvals(j,:)
!        do i = 1, npop
!          ra(i) = wksp(iwksp(i))
!        enddo
!        objvals(j,:) = ra
!      enddo 
!      do j = 1, ndim
!        wksp = pop(j,:)
!        do i = 1, npop
!          ra(i) = wksp(iwksp(i))
!        enddo
!        pop(j,:) = ra
!      enddo

      ind = 1
      poptmp = 0.0d0
      objtmp = 0.0d0
      ndom = .True.
      poptmp(:,ind) = pop(:,ind) 
      objtmp(:,ind) = objvals(:,ind)
      ndom(ind) = .False.
      do i = 2, npop
        tmpdom2 = .True. !assume initially i is nondominated with respect to i-1
        do k = 1, i-1
          if(ndom(k)) then
          else
            tmpdom = .False.
            idom = 0
            do j = 1, nobj
              if(objvals(j,i).le.objvals(j,k)) then
                idom = idom + 1
              endif
              if(objvals(j,i).lt.objvals(j,k)) then
              !if((-objvals(j,i)+objvals(j,k)).gt.eps) then
                tmpdom=.True.
              endif
            enddo
            if((idom.eq.nobj).and.tmpdom) then !i dominates k
               ndom(k) = .True.
            else if(idom.eq.nobj) then
               tmpdom=.True.
            endif
            if(tmpdom) then
            else
              tmpdom2 = .False.
            endif
          endif
        enddo
        if(tmpdom2) then
          ndom(i) = .False.
        endif
      enddo
      ind = 0
      do i = 1, npop
        if(ndom(i)) then
        else
          ind = ind + 1
          !copy the non-dominated solutions into poptmp,objtmp
          poptmp(:,ind) = pop(:,i)
          objtmp(:,ind) = objvals(:,i)
        endif
      enddo
      
      ii = 0
      do i = 1, npop
        if(ndom(i)) then
          ii = ii + 1
          !print*,"ii: ",i,ii,ind+ii
          poptmp(:,ind+ii) = pop(:,i)
          objtmp(:,ind+ii) = objvals(:,i)
        endif
      enddo

      end subroutine ndomgroup2

subroutine inihalt ( flag, dimen, atmost, quasi )

!*****************************************************************************80
!
!! INIHALT initializes the Halton quasirandom number generator.
!
!  Discussion:
!
!    INIHALT first checks whether the user-supplied dimension DIMEN of
!    the quasirandom vectors is acceptable (between 2 and 50).
!
!    INIHALT then calculates a tolerance parameter E to make the program work
!    correctly in finite precision arithmetic and a parameter DELTA
!    to check that E works.  If the test is not passed, then ATMOST
!    is too big relative to the machine precision.
!
!    Otherwise, INIHALT computes and returns the first vector QUASI.
!    For the following values of QUASI, it is necessary to call GOHALT.
!
!  Modified:
!
!    4 June 2024
!
!  Parameters:
!
!    Output, logical FLAG(2), error flags.
!    FLAG(1) is FALSE if the input value of DIMEN is unacceptable.
!    FLAG(2) is FALSE if the input value of ATMOST is unacceptable.
!
!    Input, integer DIMEN, the spatial dimension.  DIMEN should
!    satisfy: 2 <= DIMEN <= 50.
!
!    Input, integer ATMOST, the maximum number of quasirandom
!    vectors to be computed.
!
!    Output, double precision QUASI(DIMEN), the first element of
!    the Halton sequence.
!
!  Local Parameters:
!
!    Local, double precision PRIME(50), the first 50 primes.
!
!  Global Parameters:
!
!    Stored in common block /HALTON/:
!
!    Global, double precision E, a tolerance.
!
!    Global, double precision PRIME_INV(50), the reciprocals of the
!    first 50 primes.
!
  implicit none

  integer, intent(in) :: dimen
  integer, intent(in) :: atmost
  double precision, intent(out) :: quasi(dimen)
  integer, parameter :: dim_max = 50

  double precision delta
  double precision e
  logical flag(2)
  integer prime(dim_max)
  double precision prime_inv(dim_max)
  double precision small

  common /halton/ e, prime_inv
  save /halton/
!
!  Check DIMEN.
!
  flag(1) = .true.
  flag(2) = .true.

  if ( dimen < 2 .or. dim_max < dimen ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INIHALT - Fatal error!'
    write ( *, '(a)' ) '  The spatial dimension should satisfy:'
    write ( *, '(a,i6)' ) '    2 <= dimen <= ', dim_max
    write ( *, '(a,i6)' ) '  But this input value is dimen = ', dimen
    flag(1) = .false.
    return
  end if
!
!  Set the primes.
!
  prime(1:dim_max) = (/ &
      2,   3,   5,   7,  11,  13,  17,  19,  23,  29, &
     31,  37,  41,  43,  47,  53,  59,  61,  67,  71, &
     73,  79,  83,  89,  97, 101, 103, 107, 109, 113, &
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173, &
    179, 181, 191, 193, 197, 199, 211, 223, 227, 229  /)
!
!  Compute the tolerance and make the check.
!
  small = epsilon ( small )

  e = 0.9D+00 * ( 1.0D+00 / ( dble ( atmost * prime(dimen) ) ) - 10.0D+00 * small )

  delta = 100.0D+00 * small * dble ( atmost + 1 ) * log10 ( dble ( atmost ) )

  if ( 0.09D+00 * ( e - 10.0D+00 * small ) < delta ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INIHALT - Fatal error!'
    write ( *, '(a)' ) '  The value of ATMOST is too great.'
    flag(2) = .false.
    return
  end if
!
!  Set the inverse primes.
!
  prime_inv(1:dim_max) = 1.0D+00 / dble ( prime(1:dim_max) )
!
!  Compute the first vector.
!
  quasi(1:dimen) = prime_inv(1:dimen)

end subroutine inihalt

subroutine gohalt (dimen,quasi )

!*****************************************************************************80
!
!! GOHALT generates a new quasirandom Halton vector with each call.
!
!  Discussion:
!
!    The routine adapts key ideas from Halton and Smith.
!
!    The routine INIHALT must be called once before using
!    this routine.
!
!  Parameters:
!
!    Input/output, real*8 QUASI(DIMEN), on input, the previous
!    quasirandom vector; on output, the next quasirandom vector.
!    On the first call, the input value should be the output
!    value given by INIHALT.
!
!  Global Parameters:
!
!    In labelled common block /HALTON/:
!
!    Global, double precision E, a tolerance.
!
!    Global, double precision PRIME_INV(50), the reciprocal of
!    the first 50 prime numbers.
!
  implicit none

  integer, intent(in) :: dimen
  double precision, intent(inout) :: quasi(dimen)
  integer, parameter :: dim_max = 50

  double precision e
  double precision f
  double precision g
  double precision h
  integer i
  double precision prime_inv(dim_max)
  double precision t

  common /halton/ e, prime_inv
  save /halton/
!
!  Generate QUASI one component at a time, using radix 1/PRIME(K) for 
!  component K.
!
  do i = 1, dimen

    t = prime_inv(i)
    f = 1.0D+00 - quasi(i)
    g = 1.0D+00
    h = prime_inv(i)

    do

      if ( e <= f - h ) then
        exit
      end if
!
!  This checks whether Q + H > 1 - E.
!
      g = h
      h = h * t
!
!  If this is the (K-1)-st time this statement is reached, check whether
!  QUASI(I) + R**(-K) > 1-E.
!
    end do
!
!  For the appropriate I (depending on how many times the loop above
!  is executed), add H**(I+1) + H**(I) - 1
!  to the old QUASI(I) to get the next QUASI(I).
!
    quasi(i) = g + h - f

  end do

end subroutine gohalt

      subroutine distharm(pop,ndim,objvals,nobj,npop,dist)
      implicit none
      integer :: ndim,npop,nobj
      real*8, dimension(ndim,npop) :: pop
      real*8, dimension(nobj,npop) :: objvals
      integer :: i,j,k
      integer, dimension(npop) :: iwksp
      real*8, dimension(npop) :: wksp,ra
      logical, dimension(npop) :: ndom
      real*8, dimension(npop,npop) :: ds
      real*8, dimension(npop) :: dist
      real*8 :: tol,ds2,x1,xmin,xmin2
      integer :: imin,imin2

      tol = 1.0d-15

      ds = 0.0d0
      do i = 1, npop
        do j = 1,npop
          if(j.ne.i) then
            ds2 = 0.0d0
            do k = 1, nobj
              ds2 = ds2 + (objvals(k,i)-objvals(k,j))**2 
            enddo
            ds(i,j) = sqrt(ds2)
          endif
        enddo
      enddo
      do j = 1, npop
        x1 = 1.0d10
        do i = 1, npop
          if(x1.ge.ds(i,j).and.(i.ne.j)) then
            x1 = ds(i,j)
            imin = i
          endif
        enddo
        xmin = x1
!        print*,"xmin: ",j,xmin,imin
        x1 = 1.0d10
        do i = 1, npop
          if((x1.ge.ds(i,j)) .and. (i.ne.imin) .and. (i.ne.j) ) then
            x1 = ds(i,j)
            imin2 = i
          endif
        enddo
        xmin2 = x1
!        print*,"xmin2: ",j,xmin2,imin2
        if((xmin.gt.tol)) then
          dist(j) = 2.0d0/(1.0d0/xmin+1.0d0/xmin2)
        else if(xmin2.gt.tol) then
          dist(j) = xmin
        else
          dist(j) = (xmin+xmin2)/2
        endif
      enddo

!      do i = 1, npop
!        print*,"dist: ",dist(i)
!      enddo
 
      ra = dist
      call indexx(npop,ra,iwksp)

      wksp = ra
      do i = 1, npop
          ra(i) = wksp(iwksp(i))
      enddo
      dist = ra
      do j = 1, nobj
        wksp = objvals(j,1:npop)
        do i = 1, npop
          ra(i) = wksp(iwksp(i))
        enddo
        objvals(j,1:npop) = ra
      enddo 
      do j = 1, ndim
        wksp = pop(j,1:npop)
        do i = 1, npop
          ra(i) = wksp(iwksp(i))
        enddo
        pop(j,1:npop) = ra
      enddo

      end subroutine distharm

      !remove nrmv solutions out of npop solutions using the
      !harmonic distance as a criteria.
      !modified the original scheme proposed by Huang et al. to account for
      !the effects that two same solutions could be in the pool.
      !nm2 is the maximum size of the input solutions.
      !nm1 is the maximum size of the output solutions.
      subroutine seldshm2(pop,ndim,objvals,nobj,npop,popout,objout,nrmv,nm1,nm2)
      implicit none
      integer :: ndim,npop,nobj,nrmv,nm2,nm1
      real*8, dimension(ndim,nm2) :: pop
      real*8, dimension(nobj,nm2) :: objvals
      real*8, dimension(ndim,nm1) :: popout
      real*8, dimension(nobj,nm1) :: objout
      integer :: i,j,k
      real*8, dimension(npop,npop) :: ds
      integer, dimension(npop,npop) :: ijds
      integer, dimension(npop) :: iwksp
      real*8, dimension(npop) :: dist,ra,wksp
      real*8 :: ds2,dmin,dtmp,tol0
      integer :: imin,imin2,idshm,ii,kk,kh
      logical, dimension(npop) :: ipop

      tol0 = -1.0d-10
      ds = 0.0d0
      do i = 1, npop
        do j = 1,npop
          if(j.ne.i) then
            ds2 = 0.0d0
            do k = 1, nobj
              ds2 = ds2 + (objvals(k,i)-objvals(k,j))**2 
            enddo
            ds(i,j) = sqrt(ds2)
          else
            ds(i,j) = -2.0d0
          endif
        enddo
      enddo

      do i = 1, npop
        ra = ds(:,i)
        call indexx(npop,ra,iwksp)
        ijds(:,i) = iwksp
      enddo

      kh = 2

      ipop = .true.
      do ii = 1, nrmv
        dmin = 1.0d11
        idshm = 100000
        do i = 1, npop
          if(ipop(i)) then
            iwksp = ijds(:,i)
            wksp = ds(:,i)
            do j = 1, npop
              ra(j) = wksp(iwksp(j))
            enddo
            kk = 0
            dtmp = 0.0d0
            do j = 1, npop
              if(ra(j).gt.tol0 .and. kk.lt.kh) then     
                kk = kk + 1
                if(ra(j).eq.0.0d0) then
                  ra(j) = 1.0d-10
                endif 
                dtmp = dtmp + 1.0d0/ra(j) 
              endif
            enddo
            if(dtmp.gt.0.0d0) then
              dist(i) = kh*1.0d0/dtmp
            else
              dist(i) = 1.0d10 
            endif
          else
            dist(i) = 1.0d10
          endif
          if(dmin.gt.dist(i)) then
             dmin = dist(i)
             idshm = i
          endif
        enddo
        ds(idshm,:) = -1.0d0
        ds(:,idshm) = -1.0d0
        !remove the one with the smallest harm. dist.
        ipop(idshm) = .false.
      enddo

      i = 0
      do j = 1, npop
          if(ipop(j)) then
            i = i + 1
            popout(:,i) = pop(:,j)
            objout(:,i) = objvals(:,j)
          endif
      enddo

!          write(28,*)npop,i
!          do j = 1, i
!            write(28,129)objout(1,j),objout(2,j)
!          enddo
!          call flush(28)
! 
!          write(29,*)npop,i
!          do j = 1, npop
!            write(29,129)objvals(1,j),objvals(2,j)
!          enddo
! 
!          call flush(29)

129       format(3(1x,e16.7))
      end subroutine seldshm2

      !select nout solutions out of npop solutions using the
      !harmonic distance as a criteria.
      !selection can miss some good but clustered solutions. 
      !their h-distance can be small due to cluster. however, one of them
      !could have a large-distance if the other clustered solutions are removed.
      !the advantage is that it can be much faster than removing the solutions.
      !modified the original scheme proposed by Huang et al. to account for
      !the effects that two same solutions could be in the pool.
      !nm2 is the maximum size of the input solutions.
      !nm1 is the maximum size of the output solutions.
      subroutine seldshm3(pop,ndim,objvals,nobj,npop,popout,objout,nout,nm1,nm2)
      implicit none
      integer :: ndim,npop,nobj,nout,nm2,nm1
      real*8, dimension(ndim,nm2) :: pop
      real*8, dimension(nobj,nm2) :: objvals
      real*8, dimension(ndim,nm1) :: popout
      real*8, dimension(nobj,nm1) :: objout
      integer :: i,j,k
      real*8, dimension(npop,npop) :: ds
      integer, dimension(npop,npop) :: ijds
      integer, dimension(npop) :: iwksp
      real*8, dimension(npop) :: dist,ra,wksp
      real*8 :: ds2,dmin,dtmp,tol0
      integer :: imin,imin2,idshm,ii,kk,kh,imloc,icount,iflag,nedgept
      logical, dimension(npop) :: ipop
      integer, dimension(nobj) :: iminst
   
      ipop = .true.
      !select the minimum edge values
      do i = 1, nobj
        dmin = 1.0d10
        imin = 10000
        do j = 1, npop
          if(dmin.ge.objvals(i,j)) then
            dmin = objvals(i,j)
            imin = j
          endif
        enddo
        ipop(imin) = .false.
!        print*,"imin: ",imin,npop,objvals(i,imin),i,nm1,nm2
        iminst(i) = imin
      enddo

      !find the number of edge points to be included in the selection
      icount = 0
      do i = 1, nobj
        do j = 1, nobj
          if(iminst(i).eq.iminst(j)) then
            icount = icount + 1
          endif
        enddo
      enddo
      nedgept = nobj-(icount - nobj)

      tol0 = -1.0d-10
      ds = 0.0d0
      do i = 1, npop
        do j = 1,npop
          if(j.ne.i) then
            ds2 = 0.0d0
            do k = 1, nobj
              ds2 = ds2 + (objvals(k,i)-objvals(k,j))**2 
            enddo
            ds(i,j) = sqrt(ds2)
          else
            ds(i,j) = -2.0d0
          endif
        enddo
      enddo

      do i = 1, npop
        ra = ds(:,i)
        call indexx(npop,ra,iwksp)
        ijds(:,i) = iwksp
      enddo

      kh = 2

      do i = 1, npop
            iwksp = ijds(:,i)
            wksp = ds(:,i)
            do j = 1, npop
              ra(j) = wksp(iwksp(j))
            enddo
            kk = 0
            dtmp = 0.0d0
            do j = 1, npop
              if(ra(j).gt.tol0 .and. kk.lt.kh) then     
                kk = kk + 1
                if(ra(j).eq.0.0d0) then
                  ra(j) = 1.0d-10
                endif 
                dtmp = dtmp + 1.0d0/ra(j) 
              endif
            enddo
            if(dtmp.gt.0.0d0) then
              dist(i) = kh*1.0d0/dtmp
            else
              dist(i) = -1.0d10 
            endif
      enddo

      ra = dist
      call indexx(npop,ra,iwksp)

      !reduce the nout to be selected solution be icount in order
      !to include the minimum boundary solutions.
      !select nout-nedgept largest h-averaged distance.
      icount = 0
      do i = npop, 1, -1
        if(ipop(iwksp(i))) then
          ipop(iwksp(i)) = .false.
          icount = icount + 1
          if(icount.eq.nout-nedgept) then
            exit
          endif
        endif
      enddo

      i = 0
      do j = 1, npop
          if(ipop(j)) then
          else
            i = i + 1
            popout(:,i) = pop(:,j)
            objout(:,i) = objvals(:,j)
          endif
      enddo

!      print*,"i: ",i,icount,nedgept
!          write(28,*)npop,i
!          do j = 1, i
!            write(28,129)objout(1,j),objout(2,j)
!          enddo
!          call flush(28)
 
!          write(29,*)npop,i
!          do j = 1, npop
!            write(29,129)objvals(1,j),objvals(2,j)
!          enddo
 
!          call flush(29)

129       format(3(1x,e16.7))
      end subroutine seldshm3

!generate the new set of parameters. 
subroutine genFCnp(f1,f2,f3,f4,c1,frange,crange,iseed,nplc)
implicit none
integer, intent(in) :: nplc
integer*8, intent(inout) :: iseed
real*8, intent(out), dimension(nplc) :: f1,f2,f3,f4,c1
real*8, dimension(2,4), intent(in) :: frange
real*8, dimension(2), intent(in) :: crange
real*8 :: tmp1
real*8 :: ran22
integer :: i
 
  do i = 1, nplc
    tmp1 = ran22(iseed)
    f1(i) = frange(1,1)+tmp1*(frange(2,1)-frange(1,1))
    tmp1 = ran22(iseed)
    f2(i) = frange(1,2)+tmp1*(frange(2,2)-frange(1,2))
    tmp1 = ran22(iseed)
    f3(i) = frange(1,3)+tmp1*(frange(2,3)-frange(1,3))
    tmp1 = ran22(iseed)
    f4(i) = frange(1,4)+tmp1*(frange(2,4)-frange(1,4))
    tmp1 = ran22(iseed)
    c1(i) = crange(1)+tmp1*(crange(2)-crange(1))
  enddo
 
end subroutine genFCnp
