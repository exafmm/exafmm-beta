      subroutine mpi_shift(var, nold, mpisize, mpirank)
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      parameter(nmax = 1000000)
      integer istatus(mpi_status_size)
      real(8),intent(inout) :: var(nmax)
      real(8),allocatable :: buf(:)
      isend = mod(mpirank + 1,           mpisize)
      irecv = mod(mpirank - 1 + mpisize, mpisize)

      call mpi_isend(nold, 1, mpi_integer, irecv, 0, mpi_comm_world, ireqs, ierr)
      call mpi_irecv(nnew, 1, mpi_integer, isend, 0, mpi_comm_world, ireqr, ierr)
      call mpi_wait(ireqs, istatus, ierr)
      call mpi_wait(ireqr, istatus, ierr)

      allocate( buf(nnew) )
      call mpi_isend(var, nold, mpi_real8, irecv, 1, mpi_comm_world, ireqs, ierr)
      call mpi_irecv(buf, nnew, mpi_real8, isend, 1, mpi_comm_world, ireqr, ierr)
      call mpi_wait(ireqs, istatus, ierr)
      call mpi_wait(ireqr, istatus, ierr)
      do i = 1,nnew
        var(i) = buf(i)
      end do
      nold = nnew
      deallocate(buf)
      return
      end

      program test_laplace
      use iso_c_binding, only: c_ptr
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      parameter(nmax = 1000000, pi = 3.14159265358979312d0)
      integer prange, d
      real(8) norm, norm1, norm2, norm3, norm4
      type(c_ptr) ctx
      integer, dimension (128) :: iseed
      real(8), dimension (3) :: dipole = (/0, 0, 0/)
      real(8), dimension (3) :: dipole2
      real(8), dimension (3) :: xperiodic
      allocatable :: x(:)
      allocatable :: q(:)
      allocatable :: p(:)
      allocatable :: f(:)
      allocatable :: x2(:)
      allocatable :: q2(:)
      allocatable :: p2(:)
      allocatable :: f2(:)

      call mpi_init(ierr)
      call mpi_comm_size(mpi_comm_world, mpisize, ierr);
      call mpi_comm_rank(mpi_comm_world, mpirank, ierr);
      ni = 2;
      images = 0
      ksize = 11
      pcycle = 2 * pi
      alpha = 10 / pcycle
      allocate( x(3*nmax),  q(nmax),  p(nmax),  f(3*nmax)  )
      allocate( x2(3*nmax), q2(nmax), p2(nmax), f2(3*nmax) )

      do i = 1, 128
        iseed(i) = mpirank
      end do
      call random_seed(put=iseed)
      call random_number(x)
      call random_number(q)
      average = 0
      do i = 1, ni
        x(3*i-2) = x(3*i-2) * pcycle - pcycle / 2
        x(3*i-1) = x(3*i-1) * pcycle - pcycle / 2
        x(3*i-0) = x(3*i-0) * pcycle - pcycle / 2
        p(i) = 0
        f(3*i-2) = 0
        f(3*i-1) = 0
        f(3*i-0) = 0
        average = average + q(i)
      end do
      average = average / ni
      do i = 1, n
        q(i) = q(i) - average
      end do
      call fmm_init(images)
      call fmm_partition(ni, x, q, pcycle)
      call fmm(ni, x, q, p, f, pcycle)
      do i = 1, ni
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
      call fmm_ewald(ni, x2, q2, p2, f2, ksize, alpha, pcycle)
#else
      prange = 0
      do i = 0, images-1
        prange = prange + 3**i
      end do
      do i = 1, ni
        do d = 1, 3
          dipole(d) = dipole(d) + x(3*i+d-3) * q(i)
        end do
      end do
      call mpi_allreduce(dipole, dipole2, 3, mpi_real8, mpi_sum, mpi_comm_world, ierr)
      norm = dipole2(1) * dipole2(1) + dipole2(2) * dipole2(2) + dipole2(3) * dipole2(3)
      coef = 4 * pi / (3 * pcycle * pcycle * pcycle)
      nj = ni
      if(mpirank.eq.0) print"(a)",'--- MPI direct sum ---------------'
      do irank = 0, mpisize-1
        if (mpirank.eq.0) print"(a,i1,a,i1)",'Direct loop          : ',irank+1,'/',mpisize
        call mpi_shift(x2, 3*nj, mpisize, mpirank)
        call mpi_shift(q2, nj,   mpisize, mpirank)
        do i = 1, ni
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
                do j = 1, nj
                  dx = x(3*i-2) - x2(3*j-2) - xperiodic(1)
                  dy = x(3*i-1) - x2(3*j-1) - xperiodic(2)
                  dz = x(3*i-0) - x2(3*j-0) - xperiodic(3)
                  R2 = dx * dx + dy * dy + dz * dz
                  Rinv = 1 / sqrt(R2)
                  if(irank.eq.mpisize-1.and.i.eq.j) Rinv = 0
                  R3inv = q2(j) * Rinv * Rinv * Rinv
                  pp = pp + q2(j) * Rinv
                  fx = fx + dx * R3inv
                  fy = fy + dy * R3inv
                  fz = fz + dz * R3inv
                end do
              end do
            end do
          end do
          p2(i) = p2(i) + pp - coef * norm / n / q(i)
          f2(3*i-2) = f2(3*i-2) - fx - coef * dipole2(1)
          f2(3*i-1) = f2(3*i-1) - fy - coef * dipole2(2)
          f2(3*i-0) = f2(3*i-0) - fz - coef * dipole2(3)
        end do
      end do
#endif
      diff2 = 0
      norm2 = 0
      diff3 = 0
      norm3 = 0
      diff4 = 0
      norm4 = 0
      v = 0
      v2 = 0
      do i = 1, ni
        v = v + p(i) * q(i)
        v2 = v2 + p2(i) * q(i)
        diff2 = diff2 + (f(3*i-2) - f2(3*i-2)) * (f(3*i-2) - f2(3*i-2))&
                      + (f(3*i-1) - f2(3*i-1)) * (f(3*i-1) - f2(3*i-1))&
                      + (f(3*i-0) - f2(3*i-0)) * (f(3*i-0) - f2(3*i-0))
        norm2 = norm2 + f2(3*i-2) * f2(3*i-2) + f2(3*i-1) * f2(3*i-1) + f2(3*i-0) * f2(3*i-0)
      end do
      diff1 = (v - v2) * (v - v2)
      norm1 = v2 * v2
      call mpi_reduce(diff1, diff3, 1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr);
      call mpi_reduce(norm1, norm3, 1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr);
      call mpi_reduce(diff2, diff4, 1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr);
      call mpi_reduce(norm2, norm4, 1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr);
      if (mpirank.eq.0) then
        print"(a)",'--- FMM vs. direct ---------------'
	print"(a,f10.7)",'Rel. L2 Error (pot)  : ', sqrt(diff3/norm3)
      end if
      if (abs(diff4).gt.0) then
        print"(a,f10.7)",'Rel. L2 Error (acc)  : ', sqrt(diff4/norm4)
      end if

      deallocate( x, q, p, f, x2, q2, p2, f2 )
      call mpi_finalize(ierr);
      end
