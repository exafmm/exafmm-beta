program main
  implicit none
  include 'mpif.h'
  character(len=128) filename
  integer i, j, ni, nj, ierr, mpisize, mpirank
  integer images, ncrit, nworkers, verbose
  integer, dimension (128) :: iseed
  real(8) pi, pcycle, theta, eps2, diff1, norm1, diff2, norm2
  real(8), allocatable, dimension(:) :: xi, yi, vi, xj, yj, vj, vd
  parameter(pi=3.14159265358979312d0)

  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, mpisize, ierr)
  call mpi_comm_rank(mpi_comm_world, mpirank, ierr)
  ni = 1000
  nj = 1000
  pcycle = 2 * pi
  ncrit = 32
  nworkers = 1
  eps2 = 0.000001
  verbose = 0
  allocate( xi(2*ni),  yi(2*ni), vi(2*ni), vd(2*ni) )
  allocate( xj(2*nj),  yj(2*nj), vj(2*nj) )
  do i = 1, 128
     iseed(i) = mpirank
  enddo
  call random_seed(put=iseed)
  call random_number(xi)
  call random_number(yi)
  call random_number(xj)
  call random_number(yj)
  do i = 1, ni
     xi(i) = xi(i) * pcycle - pcycle / 2
     yi(i) = yi(i) * pcycle - pcycle / 2
     vi(i) = 0
  enddo
  do i = 1, nj
     xj(i) = xj(i) * pcycle - pcycle / 2
     yj(i) = yj(i) * pcycle - pcycle / 2
     vj(i) = 1.0 / nj
  enddo
  call fmm_init(eps2, ncrit, nworkers, ni, xi, yi, vi, nj, xj, yj, vj)
  call fmm_partition(ni, xi, yi, vi, nj, xj, yj, vj)
  call fmm_buildtree()
  call fmm_v2b(vi, vj, verbose)
  do i = 1, ni
     vd(i) = 0
  enddo
  call direct(100, xi, yi, vd, nj, xj, yj, vj)

  diff1 = 0
  norm1 = 0
  do i = 1, 100
     diff1 = diff1 + (vi(i) - vd(i)) * (vi(i) - vd(i))
     norm1 = norm1 + vd(i) * vd(i)
  enddo
  diff2 = 0
  norm2 = 0
  call mpi_reduce(diff1, diff2, 1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(norm1, norm2, 1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr)
  if (mpirank == 0) then
     print"(a,f9.6)",'Rel. L2 Error (pot)  : ', sqrt(diff2/norm2)
  endif

  deallocate( xi, yi, vi, vd, xj, yj, vj )
  call fmm_finalize()
  call mpi_finalize(ierr)
end program main
