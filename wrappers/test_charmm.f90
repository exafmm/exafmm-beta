module charmm_io

  real(8), allocatable, dimension(:) :: rscale,gscale,fgscale
  integer idist

contains

  subroutine charmm_cor_read(n,x,q,size,nam,numex,natex)
    implicit none
    integer i, im, in, n
    integer, allocatable, dimension(:) :: numex, natex
    real(8) size, sizex, sizey, sizez
    real(8), allocatable, dimension(:) :: x, q
    character(len=100) lin
    character(len=*) nam
    logical qext
    integer vdw_size,natex_size

    open(unit=1,file=nam,status='old')
    read(1,'(a100)')lin
    read(lin(9:100),*)sizex,sizey,sizez
    size=sizex
    do while (lin(1:1) == '*')
       read(1,'(a100)')lin
    enddo
    write(*,*)'lin=',lin(1:30)
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
       write(*,*)'Something is wrong with the input file'
       write(*,*)nam
       stop
    endif
    allocate(natex(natex_size))
    read(1,'(10i10)')natex(1:natex_size)
    read(1,'(a100)')lin
    read(lin(16:100),*)idist
    vdw_size=idist**2
    allocate(rscale(vdw_size),gscale(vdw_size),fgscale(vdw_size))
    read(1,'(5f20.10)')rscale(1:vdw_size)
    read(1,'(a100)')lin
    read(1,'(5f20.10)')gscale(1:vdw_size)
    read(1,'(a100)')lin
    read(1,'(5f20.10)')fgscale(1:vdw_size)

    return
  end subroutine charmm_cor_read
end module charmm_io

subroutine split_range(ista, iend, isplit, nsplit)
  implicit real*8(a-h,o-z)
  integer size, remainder
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
  integer d, i, iend, ierr, images, ista, istat, ksize, lnam, mpirank, mpisize
  integer nglobal, prange
  real(8) accDif, accDifGlob, accNrm, accNrm2, accNrmGlob, accNrmGlob2
  real(8) alpha, average, ccelec, cutoff, norm, pcycle, pi
  real(8) potDif, potDifGlob, potNrmGlob2, potSum, potSum2, potSumGlob, potSumGlob2, sigma
  integer, dimension (128) :: iseed
  real(8), dimension (3) :: dipole = (/0, 0, 0/)
  real(8), dimension (3) :: xperiodic
  integer, allocatable, dimension(:) :: icpumap, numex, natex
  real(8), allocatable, dimension(:) :: x, q, p, f, x2, q2, p2, f2
  parameter(pi = 3.14159265358979312d0, ccelec=332.0716d0)

  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, mpisize, ierr)
  call mpi_comm_rank(mpi_comm_world, mpirank, ierr)
  nglobal = 1000
  images = 3
  ksize = 11
  pcycle = 2 * pi
  sigma = .25 / pi
  cutoff = 10.
  alpha = 10 / pcycle
  if (command_argument_count() > 0) then
     call get_command_argument(1,filename,lnam,istat)
     call charmm_cor_read(nglobal,x,q,pcycle,filename,numex,natex)
     allocate( p(nglobal),  f(3*nglobal), icpumap(nglobal) )
     allocate( x2(3*nglobal), q2(nglobal), p2(nglobal), f2(3*nglobal) )
     alpha = 10 / pcycle
  else
     allocate( x(3*nglobal),  q(nglobal),  p(nglobal),  f(3*nglobal), icpumap(nglobal) )
     allocate( x2(3*nglobal), q2(nglobal), p2(nglobal), f2(3*nglobal) )
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
  end if
  ista = 1
  iend = nglobal
  call split_range(ista,iend,mpirank,mpisize)
  do i = ista, iend
     icpumap(i) = 1
  end do
  call fmm_init(images)
  call fmm_partition(nglobal, icpumap, x, q, pcycle)
  call fmm(nglobal, icpumap, x, q, p, f, pcycle)
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
  call ewald(nglobal, icpumap, x2, q2, p2, f2, ksize, alpha, sigma, cutoff, pcycle)
#else
  prange = 0
  do i = 0, images-1
     prange = prange + 3**i
  end do
  do i = 1, nglobal
     if (icpumap(i).eq.1) then
        pp = 0
        fx = 0
        fy = 0
        fz = 0
        do ix = -prange, prange
           do iy = -prange, prange
              do iz = -prange, prange
                 xperiodic(1) = ix * pcycle
                 xperiodic(2) = iy * pcycle
                 xperiodic(3) = iz * pcycle
                 do j = 1, nglobal
                    dx = x(3*i-2) - x2(3*j-2) - xperiodic(1)
                    dy = x(3*i-1) - x2(3*j-1) - xperiodic(2)
                    dz = x(3*i-0) - x2(3*j-0) - xperiodic(3)
                    R2 = dx * dx + dy * dy + dz * dz
                    Rinv = 1 / sqrt(R2)
                    if (R2.eq.0) Rinv = 0
                    R3inv = q2(j) * Rinv * Rinv * Rinv
                    pp = pp + q2(j) * Rinv
                    fx = fx + dx * R3inv
                    fy = fy + dy * R3inv
                    fz = fz + dz * R3inv
                 end do
              end do
           end do
        end do
        p2(i) = p2(i) + pp
        f2(3*i-2) = f2(3*i-2) - fx
        f2(3*i-1) = f2(3*i-1) - fy
        f2(3*i-0) = f2(3*i-0) - fz
     end if
  end do
  do i = 1, nglobal
     do d = 1, 3
        dipole(d) = dipole(d) + x(3*i+d-3) * q(i)
     end do
  end do
  norm = dipole(1) * dipole(1) + dipole(2) * dipole(2) + dipole(3) * dipole(3)
  coef = 4 * pi / (3 * pcycle * pcycle * pcycle)
  do i = 1, nglobal
     p2(i) = p2(i) - coef * norm / nglobal / q(i)
     f2(3*i-2) = f2(3*i-2) - coef * dipole(1)
     f2(3*i-1) = f2(3*i-1) - coef * dipole(2)
     f2(3*i-0) = f2(3*i-0) - coef * dipole(3)
  end do
#endif
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
     print"(a,f10.6)",'GRMS (FMM)           : ', ccelec*sqrt(accNrmGlob/3.0/nglobal)/2.0
     print"(a,f10.6)",'GRMS (Ewald)         : ', ccelec*sqrt(accNrmGlob2/3.0/nglobal)/2.0
  end if

  deallocate( x, q, p, f, icpumap, x2, q2, p2, f2 )
  call mpi_finalize(ierr)
end program main
