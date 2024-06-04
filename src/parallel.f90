      module parallel
!---------------------------------------------------
!for parallel parameter scan
        integer, parameter :: ngproc = 1
        integer :: ngroup,mygroupid
        integer :: orig_group,new_group,gworld,lworld
        integer :: nproc,idproc,mpicommwd
!---------------------------------------------------
      contains

!------------------------------------------------------------------
!This function sets up the preparation for the parallel optimization
!
        subroutine init_parallel()
        implicit none
        include 'mpif.h'
        integer :: myid,myidx,myidy,ierr,Flagbc,npyhalf
!---------------------------------------------------
!for parallel parameter scan
	integer, dimension(ngproc) :: ranks
        character*4 name1
        character*5 name2
        character*6 name3
        character*7 name4
        integer :: i,j,k,l,m,n,lcid


        !start up MPI.
        call MPI_INIT(ierr)

!---------------------------------------------------
!for parallel parameter scan
        lworld = MPI_COMM_WORLD
        call MPI_COMM_RANK(lworld,idproc,ierr)
        call MPI_COMM_SIZE(lworld,nproc,ierr)

        ngroup = nproc/ngproc
        if(mod(nproc,ngproc).ne.0) then
          print*,"wrong number of processors!"
          print*,"number of processors has to be multiple of the processor per group: ",ngproc
          stop
        endif
 
        !this is the group id for each processor id
        mygroupid=idproc/ngproc

        ! array "ranks" lists the global PE# of each PE in mygroupnum:
        do n=1,ngproc
          ranks(n)=mygroupid*ngproc + n-1
        enddo

        call MPI_COMM_GROUP(lworld,orig_group,ierr)
        call MPI_GROUP_INCL(orig_group,ngproc,ranks,new_group,ierr)
        call MPI_COMM_CREATE(lworld,new_group,gworld,ierr)

        call MPI_COMM_RANK(gworld,lcid,ierr)

        mpicommwd = gworld
        call MPI_BARRIER(lworld,ierr)
        if(mygroupid.eq.0) print*,"finish initalization.....",lcid


        name1 = 'tmpx'
        name2 = 'tmpxx'
        name3 = 'tmpxxx'
        name4 = 'tmpxxxx'
        if(mygroupid.lt.10) then
            name1(4:4) = char(mygroupid+48)
            call system("rm -r "//name1)
            call system("mkdir "//name1)
!            call system("cp rf* input.in opt.in partcl.data ImpactTexe* "//name1)
            call system("cp * "//name1)
            call chdir(name1)
          else if((mygroupid.ge.10).and.(mygroupid.lt.100)) then
            i = mygroupid/10
            j = mygroupid - 10*i
            name2(4:4) = char(i+48)
            name2(5:5) = char(j+48)
            call system("mkdir "//name2)
            call system("cp * "//name2)
            call chdir(name2)
          else if((mygroupid.ge.100).and.(mygroupid.lt.1000)) then
            i = mygroupid/100
            j = mygroupid - 100*i
            k = j/10
            l = j - 10*k
            name3(4:4) = char(i+48)
            name3(5:5) = char(k+48)
            name3(6:6) = char(l+48)
            call system("mkdir "//name3)
            call system("cp * "//name3)
            call chdir(name3)
          else
            i = mygroupid/1000
            j = mygroupid - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(4:4) = char(i+48)
            name4(5:5) = char(k+48)
            name4(6:6) = char(m+48)
            name4(7:7) = char(n+48)
            call system("mkdir "//name4)
            call system("cp * "//name4)
            call chdir(name4)
          endif

        end subroutine init_parallel

        subroutine end_parallel()
        implicit none
        include 'mpif.h'
        integer :: ierr

        call system("cd .. ")
        call MPI_Finalize(ierr)

        end subroutine end_parallel

      end module parallel
