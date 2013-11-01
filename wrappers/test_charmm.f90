module charmm_io
contains
  subroutine charmm_cor_read(n,x,q,size,filename,numex,natex,nat,atype,rscale,gscale,fgscale)
    implicit none
    logical qext
    integer i, nat, im, in, n, natex_size, vdw_size
    real(8) size, sizex, sizey, sizez
    integer, allocatable, dimension(:) :: numex, natex, atype
    real(8), allocatable, dimension(:) :: x, q, rscale, gscale, fgscale
    character(len=100) lin
    character(len=*) filename

    open(unit=1,file=filename,status='old')
    read(1,'(a100)')lin
    read(lin(9:100),*)sizex,sizey,sizez
    size=sizex
    do while (lin(1:1) == '*')
       read(1,'(a100)')lin
    enddo
    qext=index(lin,'EXT') > 0
    if(qext)then
       read(lin,'(i10)')n
    else
       read(lin,'(i5)')n
    endif
    allocate( x(3*n), q(n) )
    do i=1,n
       if(qext)then
          read(1,'(2i10,20x,3f20.10,20x,f20.10)') &
               in,im,x(i*3-2),x(i*3-1),x(i*3),q(i)
       else
          read(1,'(2i5,10x,3f10.5,10x,f10.5)') &
               in,im,x(i*3-2),x(i*3-1),x(i*3),q(i)
       endif
    enddo
    read(1,'(a100)')lin
    allocate(numex(n))
    read(1,'(10i10)')numex(1:n)
    read(1,'(a100)')lin
    read(lin(6:100),*)natex_size
    if(natex_size /= sum(numex(1:n))) then
       print*,'Something is wrong with the input file'
       print*,filename
       stop
    endif
    allocate(natex(natex_size))
    read(1,'(10i10)')natex(1:natex_size)
    read(1,'(a100)')lin
    read(lin(16:100),*)nat
    vdw_size=nat**2
    allocate(rscale(vdw_size),gscale(vdw_size),fgscale(vdw_size))
    read(1,'(5f20.10)')rscale(1:vdw_size)
    read(1,'(a100)')lin
    read(1,'(5f20.10)')gscale(1:vdw_size)
    read(1,'(a100)')lin
    read(1,'(5f20.10)')fgscale(1:vdw_size)
    allocate(atype(n))
    read(1,'(a100)')lin
    read(1,'(20i5)')atype(1:n)
    return
  end subroutine charmm_cor_read
end module charmm_io

subroutine split_range(ista, iend, isplit, nsplit)
  implicit none
  integer ista, iend, isplit, nsplit
  integer increment, remainder, size
  size = iend - ista + 1
  increment = size / nsplit
  remainder = mod(size, nsplit)
  ista = ista + isplit * increment + min(isplit,remainder)
  iend = ista + increment - 1
  if(remainder.gt.isplit) iend = iend + 1
  return
end subroutine split_range

program main
  use charmm_io
  implicit none
  include 'mpif.h'
  character(128) filename
  integer d, i, ierr, images, ista, iend, istat, ksize, lnam, mpirank, mpisize
  integer nat, nglobal, prange, ix, iy, iz, j
  real(8) alpha, sigma, cuton, cutoff, average, norm, pcycle, pi, ccelec
  real(8) coef, dx, dy, dz, fx, fy ,fz, pp, R2, R3inv, Rinv
  real(8) accDif, accDifGlob
  real(8) accNrm, accNrmGlob
  real(8) accNrm2, accNrmGlob2
  real(8) potDif, potDifGlob
  real(8) potNrmGlob2
  real(8) potSum, potSumGlob
  real(8) potSum2, potSumGlob2
  integer, dimension (128) :: iseed
  real(8), dimension (3) :: dipole = (/0, 0, 0/)
  real(8), dimension (3) :: xperiodic
  integer, allocatable, dimension(:) :: icpumap, numex, natex, atype
  real(8), allocatable, dimension(:) :: x, q, p, f, x2, q2, p2, f2
  real(8), allocatable, dimension(:) :: rscale, gscale, fgscale
  parameter(pi=3.14159265358979312d0, ccelec=332.0716d0)

  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, mpisize, ierr)
  call mpi_comm_rank(mpi_comm_world, mpirank, ierr)
  nglobal = 1000
  images = 3
  ksize = 11
  pcycle = 20 * pi
  sigma = .25 / pi
  cuton = 9.5;
  cutoff = 10.0;
  alpha = 10 / pcycle
  nat = 16
  if (command_argument_count() > 0) then
     call get_command_argument(1,filename,lnam,istat)
     call charmm_cor_read(nglobal,x,q,pcycle,filename,numex,natex,nat,atype,rscale,gscale,fgscale)
     allocate( p(nglobal),  f(3*nglobal), icpumap(nglobal) )
     allocate( x2(3*nglobal), q2(nglobal), p2(nglobal), f2(3*nglobal) )
     alpha = 10 / pcycle
  else
     allocate( x(3*nglobal),  q(nglobal),  p(nglobal),  f(3*nglobal), icpumap(nglobal) )
     allocate( x2(3*nglobal), q2(nglobal), p2(nglobal), f2(3*nglobal) )
     allocate( numex(nglobal), natex(nglobal), atype(nglobal) )
     allocate( rscale(nat*nat), gscale(nat*nat), fgscale(nat*nat) )
#if 1
     do i = 1, 128
        iseed(i) = 0
     end do
     call random_seed(put=iseed)
     call random_number(x)
     call random_number(q)
     average = 0
     do i = 1, nglobal
        x(3*i-2) = x(3*i-2) * pcycle - pcycle / 2
        x(3*i-1) = x(3*i-1) * pcycle - pcycle / 2
        x(3*i-0) = x(3*i-0) * pcycle - pcycle / 2
        p(i) = 0
        f(3*i-2) = 0
        f(3*i-1) = 0
        f(3*i-0) = 0
        icpumap(i) = 0
        average = average + q(i)
     end do
     average = average / nglobal
     do i = 1, nglobal
        q(i) = q(i) - average
     end do
#else
     write (filename, '("source", i4.4, ".dat")') 0
     open(10, file=filename)
     do i = 1, nglobal
        read(10,*) x(3*i-2), x(3*i-1), x(3*i-0), q(i)
        p(i) = 0
        f(3*i-2) = 0
        f(3*i-1) = 0
        f(3*i-0) = 0
        icpumap(i) = 0
     end do
#endif
     call random_number(rscale)
     call random_number(gscale)
     do i = 1, nglobal
        numex(i) = 1
        if(mod(i,2).eq.1)then
           natex(i) = i+1
        else
           natex(i) = i-1
        endif
        atype(i) = rand() * nat
     enddo
     do i = 1, nat
        do j = 1, nat
           if( i.ne.j ) then
              gscale((i-1)*nat+j) = sqrt(gscale((i-1)*nat+i)*gscale((j-1)*nat+j))
              rscale((i-1)*nat+j) = (sqrt(rscale((i-1)*nat+i)) + sqrt(rscale((j-1)*nat+j))) * 0.5;
              rscale((i-1)*nat+j) = rscale((i-1)*nat+j) * rscale((i-1)*nat+j);
           endif
        enddo
     enddo
     do i = 1, nat*nat
        fgscale(i) = gscale(i)
     enddo
  end if
  ista = 1
  iend = nglobal
  call split_range(ista,iend,mpirank,mpisize)
  do i = ista, iend
     icpumap(i) = 1
  end do
  call fmm_init(images)
  call fmm_partition(nglobal, icpumap, x, q, pcycle)
  !call fmm_coulomb(nglobal, icpumap, x, q, p, f, pcycle)
  call fmm_vanderwaals(nglobal, icpumap, atype, x, p, f, cuton, cutoff,&
       pcycle, nat, rscale, gscale, fgscale)
  do i = 1, nglobal
     x2(3*i-2) = x(3*i-2)
     x2(3*i-1) = x(3*i-1)
     x2(3*i-0) = x(3*i-0)
     q2(i) = q(i)
     p2(i) = 0
     f2(3*i-2) = 0
     f2(3*i-1) = 0
     f2(3*i-0) = 0
  end do
#if 1
  !call ewald_coulomb(nglobal, icpumap, x2, q2, p2, f2, ksize, alpha, sigma, cutoff, pcycle)
#else
  call direct_coulomb(nglobal, icpumap, x2, q2, p2, f2, pcycle)
#endif
  call direct_vanderwaals(nglobal, icpumap, atype, x2, p2, f2, cuton, cutoff,&
       pcycle, nat, rscale, gscale, fgscale)
  !call coulomb_exclusion(nglobal, icpumap, x, q, p, f, pcycle, numex, natex)
  !call coulomb_exclusion(nglobal, icpumap, x2, q2, p2, f2, pcycle, numex, natex)
  !call vanderwaals_exclusion(nglobal, icpumap, atype, x, p, f, cuton, cutoff,&
  !     pcycle, nat, rscale, gscale, fgscale, numex, natex)
  !call vanderwaals_exclusion(nglobal, icpumap, atype, x2, p2, f2, cuton, cutoff,&
  !     pcycle, nat, rscale, gscale, fgscale, numex, natex)
  potSum = 0
  potSum2 = 0
  accDif = 0
  accNrm = 0
  accNrm2 = 0
  do i = 1, nglobal
     if (icpumap(i).eq.1) then
        potSum  = potSum  + p(i) * q(i)
        potSum2 = potSum2 + p2(i) * q(i)
        accDif  = accDif  + (f(3*i-2) - f2(3*i-2)) * (f(3*i-2) - f2(3*i-2))&
             + (f(3*i-1) - f2(3*i-1)) * (f(3*i-1) - f2(3*i-1))&
             + (f(3*i-0) - f2(3*i-0)) * (f(3*i-0) - f2(3*i-0))
        accNrm  = accNrm  + f(3*i-2) * f(3*i-2) + f(3*i-1) * f(3*i-1) + f(3*i-0) * f(3*i-0)
        accNrm2  = accNrm2  + f2(3*i-2) * f2(3*i-2) + f2(3*i-1) * f2(3*i-1) + f2(3*i-0) * f2(3*i-0)
     end if
  end do
  potSumGlob = 0
  potSumGlob2 = 0
  accDifGlob = 0
  accNrmGlob = 0
  accNrmGlob2 = 0
  call mpi_reduce(potSum,  potSumGlob,  1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(potSum2, potSumGlob2, 1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(accDif,  accDifGlob,  1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(accNrm,  accNrmGlob,  1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(accNrm2, accNrmGlob2, 1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr)
  potDifGlob = (potSumGlob - potSumGlob2) * (potSumGlob - potSumGlob2)
  potNrmGlob2 = potSumGlob2 * potSumGlob2
  if (mpirank.eq.0) then
     print"(a)",'--- FMM vs. direct ---------------'
     print"(a,f9.7)",'Rel. L2 Error (pot)  : ', sqrt(potDifGlob/potNrmGlob2)
     print"(a,f9.7)",'Rel. L2 Error (acc)  : ', sqrt(accDifGlob/accNrmGlob2)
     print"(a,f12.4)",'Energy (FMM)         : ', ccelec*potSumGlob/2.0
     print"(a,f12.4)",'Energy (Ewald)       : ', ccelec*potSumGlob2/2.0
     print"(a,f12.4)",'GRMS (FMM)           : ', ccelec*sqrt(accNrmGlob/3.0/nglobal)/2.0
     print"(a,f12.4)",'GRMS (Ewald)         : ', ccelec*sqrt(accNrmGlob2/3.0/nglobal)/2.0
  end if

  deallocate( x, q, p, f, icpumap, x2, q2, p2, f2 )
  call mpi_finalize(ierr)
end program main
