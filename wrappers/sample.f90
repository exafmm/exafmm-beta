module simulation_setting
  integer natom ! 全原子数
  integer nstep_total     ! ステップ数
  integer dump_step       ! 構造を書き出すのは何ステップに１回？
  integer dump_ene_step       ! energyを書き出すのは何ステップに１回？
  
  double precision box_length(3)
  double precision  mass, dt
  double precision  half_dt_mass
end module simulation_setting

module lj_param
  double precision epsi, sigma
  double precision sig6, eps24, eps4
  double precision spring_k, spring_length
  double precision cutoff, cutoff2
end module lj_param

module energy_module
  type energy_var
    double precision pot ! potential
    double precision p_lj, p_sp ! potential energy of Lennard-Jones and spring
    double precision kin ! kinetic
  end type energy_var
end module energy_module

program main
  use simulation_setting
  use energy_module
  use lj_param, only : cutoff
  implicit none
  include 'mpif.h'
  ! global
  type ( energy_var ) ene1, ene2
  ! local
  integer nstep
  double precision buf
  character(80) filename, restart_file
  integer ierror, mpi_all, mpi_rank
  integer , parameter :: nsize = 100000
  integer res_index1(nsize), res_index2(nsize)
  double precision x1(3*nsize), v1(3*nsize), q(nsize), & ! for old style
    f1(3*nsize), fsp1(3*nsize), flj1(3*nsize), p1(nsize)
  double precision x2(3*nsize), v2(3*nsize), & ! for exafmm API
    f2(3*nsize), fsp2(3*nsize), flj2(3*nsize), p2(nsize)
  integer ni, nj

  ! initialize
  call read_parameter(restart_file)

  ene1%pot  = 0.0d0
  ene1%kin  = 0.0d0
  ene1%p_lj = 0.0d0
  ene1%p_sp = 0.0d0
  ene2 = ene1

  x1 = 0.0d0; f1 = 0.0d0 ; v1 = 0.0d0 ; q = 0.0d0
  call read_restartfile(restart_file, natom, nsize, &
    res_index1, x1, v1, q)

  call check_bounds(natom, nsize, res_index1, x1, box_length)

  call MPI_Init(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD , mpi_all  , ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD , mpi_rank , ierror)
  call fmm_init_wrapper(1,1,0.5d0, cutoff, 1)

  if( mpi_rank .eq. 0)then
    res_index2 = res_index1
    x2 = x1
    v2 = v1
    ni = natom
  else
    ni = 0
  endif
  nj = nsize
  call fmm_partition_wrapper(ni, nj , res_index2, x2, q, v2, box_length)
  f2(:) = 0.0d0
  
  do nstep = 1, nstep_total
    ! velocity verlet method start
    ! 1, update velocity
    v1(1:3*natom) = v1(1:3*natom) + half_dt_mass*f1(1:3*natom)
    v2(1:3*ni)    = v2(1:3*ni)    + half_dt_mass*f2(1:3*ni)
    
    ! 2, update position
    x1(1:3*natom) = x1(1:3*natom) + dt*v1(1:3*natom)
    x2(1:3*ni)    = x2(1:3*ni)    + dt*v2(1:3*ni)
    
    if( mod(nstep , 10) .eq. 1)then
      call check_bounds(natom , nsize , res_index1 , x1 , box_length)
      call check_bounds(ni    , nsize , res_index2 , x2 , box_length)
    endif

    ! 3, calc force
    call calc_force_spring(natom, nsize, res_index1, x1, fsp1, p1, ene1%p_sp)
    call calc_force_lj(natom, nj, nsize, res_index1, x1, p1, flj1, ene1%p_lj, .true., box_length)
    f1(1:3*natom) = fsp1(1:3*natom) + flj1(1:3*natom)
    
    ! exafmm api start
    if( mod(nstep, 2) .eq. 1)then
      nj = nsize
      call fmm_partition_wrapper(ni, nj , res_index2, x2, q, v2, box_length)
    endif
    nj = nsize
    call fmm_fmm_wrapper(ni, nj, res_index2, x2, q, p2, flj2, box_length)
    flj2 = 0.0d0
    ! calc_force 
    call calc_force_spring(ni, nsize, res_index2, x2, fsp2, p2, ene2%p_sp)
    call calc_force_lj(ni, nj, nsize, res_index2, x2, p2, flj2, ene2%p_lj, .false., box_length)
    f2(1:3*ni) = fsp2(1:3*ni) + flj2(1:3*ni)
    ! exafmm api end

    ! 4, update velocity
    v1(1:3*natom) = v1(1:3*natom) + half_dt_mass*f1(1:3*natom)
    v2(1:3*ni)    = v2(1:3*ni)    + half_dt_mass*(f2(1:3*ni))
    ! velocity verlet method end

    ! ***** output start *****
    if( mod( nstep-1, dump_ene_step ) .eq. 0 )then
      if( mpi_rank .eq. 0)then
        write(*,*) '# nstep ', nstep
        
        ene1%pot= ene1%p_lj + ene1%p_sp
        ene1%kin = 0.5d0*mass*sum(v1(1:3*natom)*v1(1:3*natom))
        call output_energy_dat(dt, nstep,natom,ene1)
      endif
      
      buf = 0.5d0*mass*sum(v2(1:3*ni)*v2(1:3*ni))
      call MPI_Allreduce( buf, ene2%kin, 1, MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD, ierror)
      buf = ene2%p_sp
      call MPI_Allreduce( buf, ene2%p_sp, 1, MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD, ierror)
      buf = ene2%p_lj
      call MPI_Allreduce( buf, ene2%p_lj, 1, MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD, ierror)
      ene2%pot= ene2%p_lj + ene2%p_sp
      call output_energy_dat_fmm(dt, nstep,natom,ene2, mpi_rank)
    endif

    if( mod(nstep-1, dump_step).eq. 0)then
      call output_xyz_fmm(natom, natom, nsize, x1, box_length, ene1, 10)
      call output_xyz_fmm(ni, nj, nsize, x2, box_length, ene2, mpi_rank)
      
      filename = 'restart_out.d'
      call output_checkpoint(filename, natom, nsize, res_index1, x1, v1, q)
    endif
  enddo
    
  call fmm_finalize()
  call MPI_finalize(ierror)

  close(34)
end program
