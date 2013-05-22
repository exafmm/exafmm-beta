      subroutine mpi_shift(var, n, mpisize, mpirank)
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      parameter(nmax = 1000000)
      integer istatus(mpi_status_size)
      real*8,intent(inout) :: var(nmax)
      real*8,allocatable :: buf(:)
      allocate( buf(n) )
      isend = mod(mpirank + 1,           mpisize)
      irecv = mod(mpirank - 1 + mpisize, mpisize)

      call mpi_isend(var, n, mpi_double, irecv, 1, mpi_comm_world, ireqs, ierr)
      call mpi_irecv(buf, n, mpi_double, isend, 1, mpi_comm_world, ireqr, ierr)
      call mpi_wait(ireqs, istatus, ierr)
      call mpi_wait(ireqr, istatus, ierr)
      do i = 1,n
        var(i) = buf(i)
      end do
      deallocate(buf)
      return
      end

      program test_laplace
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      parameter(nmax = 1000000, pim = 3.14159265358979312d0)
      dimension nn(8)
      allocatable :: xi(:)
      allocatable :: pi(:)
      allocatable :: fi(:)
      allocatable :: pd(:)
      allocatable :: fd(:)
      allocatable :: xj(:)
      allocatable :: qj(:)

      call mpi_init(ierr)
      call mpi_comm_size(mpi_comm_world, mpisize, ierr);
      call mpi_comm_rank(mpi_comm_world, mpirank, ierr);
      n = 1000000 / mpisize;
      allocate( xi(3*n), pi(n), fi(3*n), pd(n), fd(3*n) )
      allocate( xj(3*n), qj(n) )

      call random_seed(put=nn(1:8))
      call random_number(xi)
      call random_number(xj)
      do i = 1, n
        xi(3*i-2) = xi(3*i-2) * 2 * pim - pim
        xi(3*i-1) = xi(3*i-1) * 2 * pim - pim
        xi(3*i-0) = xi(3*i-0) * 2 * pim - pim
        pi(i) = 0
        fi(3*i-2) = 0
        fi(3*i-1) = 0
        fi(3*i-0) = 0
        pd(i) = 0
        fd(3*i-2) = 0
        fd(3*i-1) = 0
        fd(3*i-0) = 0
        xj(3*i-2) = xj(3*i-2) * 2 * pim - pim
        xj(3*i-1) = xj(3*i-1) * 2 * pim - pim
        xj(3*i-0) = xj(3*i-0) * 2 * pim - pim
        qj(i) = 1. / n
      end do
      call fmm(n, xi, pi, fi, n, xj, qj, 0)
      if(mpirank.eq.0) print"(a)",'--- MPI direct sum ---------------'
      do irank = 0, mpisize-1
        if (mpirank.eq.0) print"(a,i1.1,a,i1.1)",'Direct loop          : ',irank+1,'/',mpisize
        call mpi_shift(xj, 3*n, mpisize, mpirank)
        call mpi_shift(qj, n,   mpisize, mpirank)
        do i = 1, 100
          p = 0
          fx = 0
          fy = 0
          fz = 0
          do j = 1, n
            dx = xi(3*i-2) - xj(3*j-2)
            dy = xi(3*i-1) - xj(3*j-1)
            dz = xi(3*i-0) - xj(3*j-0)
            R2 = dx * dx + dy * dy + dz * dz
            Rinv = 1 / sqrt(R2)
            if(irank.eq.mpisize-1.and.i.eq.j) Rinv = 0
            R3inv = qj(j) * Rinv * Rinv * Rinv
            p = p + qj(j) * Rinv
            fx = fx + dx * R3inv
            fy = fy + dy * R3inv
            fz = fz + dz * R3inv
          end do
          pd(i) = pd(i) + p
          fd(3*i-2) = fd(3*i-2) - fx
          fd(3*i-1) = fd(3*i-1) - fy
          fd(3*i-0) = fd(3*i-0) - fz
        end do
      end do
      dif1 = 0
      val1 = 0
      dif2 = 0
      val2 = 0
      dif3 = 0
      val3 = 0
      dif4 = 0
      val4 = 0
      do i = 1, 100
        dif1 = dif1 + (pi(i) - pd(i)) * (pi(i) - pd(i));
        val1 = val1 + pd(i) * pd(i);
        dif2 = dif2 + (fi(3*i-2) - fd(3*i-2)) * (fi(3*i-2) - fd(3*i-2))&
              + (fi(3*i-1) - fd(3*i-1)) * (fi(3*i-1) - fd(3*i-1))&
              + (fi(3*i-0) - fd(3*i-0)) * (fi(3*i-0) - fd(3*i-0));
        val2 = val2 + fd(3*i-2) * fd(3*i-2) + fd(3*i-1) * fd(3*i-1) + fd(3*i-0) * fd(3*i-0);
      end do
      call mpi_reduce(dif1, dif3, 1, mpi_double, mpi_sum, 0, mpi_comm_world, ierr);
      call mpi_reduce(val1, val3, 1, mpi_double, mpi_sum, 0, mpi_comm_world, ierr);
      call mpi_reduce(dif2, dif4, 1, mpi_double, mpi_sum, 0, mpi_comm_world, ierr);
      call mpi_reduce(val2, val4, 1, mpi_double, mpi_sum, 0, mpi_comm_world, ierr);
      if (mpirank.eq.0) then
        print"(a)",'--- FMM vs. direct ---------------'
	print"(a,f9.7)",'Rel. L2 Error (pot)  : ', sqrt(dif3/val3)
      end if
      if (abs(dif3).gt.0) then
        print"(a,f9.7)",'Rel. L2 Error (acc)" : ', sqrt(dif4/val4)
      end if

      deallocate( xi, pi, fi, pd, fd, xj, qj )
      call mpi_finalize(ierr);
      return
      end
