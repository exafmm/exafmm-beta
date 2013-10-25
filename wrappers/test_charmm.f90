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
      end

      program main
      use iso_c_binding, only: c_ptr
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      parameter(nmax = 1000000, pi = 3.14159265358979312d0)
      integer prange, d
      real(8) norm, norm1, norm2, norm3, norm4
      type(c_ptr) ctx
      integer, dimension (128) :: iseed
      real(8), dimension (3) :: dipole = (/0, 0, 0/)
      real(8), dimension (3) :: xperiodic
      allocatable :: x(:)
      allocatable :: q(:)
      allocatable :: p(:)
      allocatable :: f(:)
      allocatable :: icpumap(:)
      allocatable :: x2(:)
      allocatable :: q2(:)
      allocatable :: p2(:)
      allocatable :: f2(:)

      call mpi_init(ierr)
      call mpi_comm_size(mpi_comm_world, mpisize, ierr);
      call mpi_comm_rank(mpi_comm_world, mpirank, ierr);
      nglobal = 1000;
      images = 2
      ksize = 11
      pcycle = 2 * pi
      alpha = 10 / pcycle
      allocate( x(3*nmax),  q(nmax),  p(nmax),  f(3*nmax), icpumap(nmax) )
      allocate( x2(3*nmax), q2(nmax), p2(nmax), f2(3*nmax) )

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
#if 0
      call fmm_ewald(nglobal, x2, q2, p2, f2, ksize, alpha, pcycle)
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
                  if(R2.eq.0) Rinv = 0
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
      do i = 1, nglobal
        if (icpumap(i).eq.1) then
          potSum  = potSum  + p(i) * q(i)
          potSum2 = potSum2 + p2(i) * q(i)
          accDif  = accDif  + (f(3*i-2) - f2(3*i-2)) * (f(3*i-2) - f2(3*i-2))&
                            + (f(3*i-1) - f2(3*i-1)) * (f(3*i-1) - f2(3*i-1))&
                            + (f(3*i-0) - f2(3*i-0)) * (f(3*i-0) - f2(3*i-0))
          accNrm  = accNrm  + f2(3*i-2) * f2(3*i-2) + f2(3*i-1) * f2(3*i-1) + f2(3*i-0) * f2(3*i-0)
        end if
      end do
      potSumGlob = 0
      potSumGlob2 = 0
      accDifGlob = 0
      accNrmGlob = 0
      call mpi_reduce(potSum,  potSumGlob,  1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr);
      call mpi_reduce(potSum2, potSumGlob2, 1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr);
      call mpi_reduce(accDif,  accDifGlob,  1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr);
      call mpi_reduce(accNrm,  accNrmGlob,  1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr);
      potDifGlob = (potSumGlob - potSumGlob2) * (potSumGlob - potSumGlob2)
      potNrmGlob = potSumGlob2 * potSumGlob2
      if (mpirank.eq.0) then
        print"(a)",'--- FMM vs. direct ---------------'
	print"(a,f9.7)",'Rel. L2 Error (pot)  : ', sqrt(potDifGlob/potNrmGlob)
        print"(a,f9.7)",'Rel. L2 Error (acc)  : ', sqrt(accDifGlob/accNrmGlob)
      end if

      deallocate( x, q, p, f, icpumap, x2, q2, p2, f2 )
      call mpi_finalize(ierr);
      end
