subroutine read_parameter(restart_file)
  use simulation_setting
  use lj_param
  implicit none
  ! local
  character (len=255) :: finp ! input file name
  character(80) restart_file

  namelist /input1/ natom, nstep_total, dump_step, dump_ene_step
  namelist /data_file/ restart_file
  namelist /phys_param/ box_length, mass, dt
  namelist /lennard_jones_parameter/ epsi, sigma, spring_k, spring_length, cutoff

  call get_command_argument(1, finp)
  open(99, file=finp, status='old')
  read(99, input1)
  read(99, data_file)
  read(99, phys_param )
  half_dt_mass = (dt/2.0d0)/mass

  read(99,lennard_jones_parameter)
  sig6 = sigma**6
  eps24 = 24.0d0*epsi
  eps4 = 4.0d0*epsi
  spring_k = 100.0d0
  spring_length = 0.171d0
  cutoff2 = (1.0d0/cutoff)**2

  close(99)
  
  ! dump input data
  open(34, file='04_input_copy.dat')

  write(34,"(a30)")    "&input1"
  write(34,"(a30,i10)")"  natom = ", natom
  write(34,"(a30,i10)")"  nstep_total   = ", nstep_total
  write(34,"(a30,i10)")"  dump_step     = ", dump_step
  write(34,"(a30)")     "&end"
  write(34,"(a30)")"&data_file"
  write(34,"(a30,a40)")"  restart_file           = ", restart_file
  write(34,"(a30)")     "&end"
  write(34,"(a30)")""
  write(34,"(a30)")"&phys_param"
  write(34,"(a30,f15.5)")"  box_length(1) = ", box_length(1)
  write(34,"(a30,f15.5)")"  box_length(2) = ", box_length(2)
  write(34,"(a30,f15.5)")"  box_length(3) = ", box_length(3)
  write(34,"(a30,f15.5)")"  mass        = ", mass
  write(34,"(a30,f15.5)")"  dt          = ", dt
  write(34,"(a30)")"&end"
  write(34,"(a30)")""
  write(34,"(a30)")"&lennard_jones_parameter"
  write(34,"(a30,f15.5)")"  epsi          = ", epsi         
  write(34,"(a30,f15.5)")"  sigma         = ", sigma        
  write(34,"(a30,f15.5)")"  spring_k      = ", spring_k     
  write(34,"(a30,f15.5)")"  spring_length = ", spring_length
  write(34,"(a30,f15.5)")"  cutoff        = ", cutoff       
  write(34,"(a30)")"&end"
  write(34,"(a30)")""

end subroutine

subroutine read_restartfile(filename, natom, nsize, &
    res_index, x1, v1, q)
  implicit none
  ! global
  character(80) filename
  integer natom, nsize, res_index(natom)
  double precision x1(3*nsize), v1(3*nsize), q(nsize)
  ! local
  integer ii

  open(20,file=filename)

  do ii = 1, natom
    read(20,*) res_index(ii), &
      x1(3*ii-2), x1(3*ii-1), x1(3*ii), v1(3*ii-2), v1(3*ii-1), v1(3*ii), q(3*ii-2)
  end do
  close(20)
end subroutine

subroutine check_bounds(ni, nsize, res_index, x, box_length)
  implicit none
  include 'mpif.h'
  ! global
  integer ni, nsize, res_index(nsize)
  double precision x(3*nsize), box_length(3)
  ! local
  integer iatom, ii, counter
  double precision dtemp(3)

  do iatom = 1, ni
    if( res_index(iatom) .le. 0) cycle

    counter = 0
    do ii = 1, ni - iatom
      if( res_index(iatom  + ii) .gt. 0)then
        exit
      else
        counter = counter + 1
      endif
    enddo

    dtemp(1:3) = x(3*iatom-2:3*iatom)
    do ii = iatom, iatom + counter
      x(3*ii-2:3*ii) = x(3*ii-2:3*ii) &
        - nint(dtemp(1:3)/box_length(1:3))*box_length(1:3)
    enddo
  enddo
end subroutine

subroutine calc_force_spring(ni, nsize, &
  res_index_fmm, x, f_sp_fmm, p, E_sp)
  use lj_param, only : spring_k, spring_length
  implicit none
  ! global
  integer ni, nsize, res_index_fmm(nsize)
  double precision x(3*nsize), f_sp_fmm(3*nsize), p(nsize)
  double precision E_sp
  ! local
  integer ii, iatom,  iatom1, iatom2, counter
  double precision delta(3), r, f(3), spring_x

  ! intra molecule (spring)
  p = 0.0d0
  f_sp_fmm = 0.0d0
  do iatom = 1, ni
    if( res_index_fmm(iatom) .le. 0) cycle

    counter =0
    do ii = 1, ni - iatom
      if( res_index_fmm(iatom  + ii) .gt. 0)then
        exit
      else
        counter = counter + 1
      endif
    enddo

    if( counter .ne. 0)then
      do iatom1 = iatom, iatom + counter -1 
        do iatom2 = iatom1+1, iatom + counter
          delta(1:3) = x(3*iatom2-2:3*iatom2) - x(3*iatom1-2:3*iatom1)
          r = dsqrt(sum( delta(1:3)*delta(1:3)))
          spring_x = r-spring_length
          f(1:3) = spring_x*(delta(1:3)/r)
          f_sp_fmm(3*iatom1-2:3*iatom1) = f_sp_fmm(3*iatom1-2:3*iatom1) + f(1:3)
          f_sp_fmm(3*iatom2-2:3*iatom2) = f_sp_fmm(3*iatom2-2:3*iatom2) - f(1:3)
          p(iatom1) = p(iatom1) + 0.25d0*spring_k*spring_x*spring_x
          p(iatom2) = p(iatom2) + 0.25d0*spring_k*spring_x*spring_x
        enddo
      enddo
    endif
  enddo
  f_sp_fmm(1:3*ni) = f_sp_fmm(1:3*ni)*spring_k
  E_sp = sum( p(1:ni))
end subroutine

subroutine calc_force_lj(ni, nj, nsize, &
    res_index, x, p, force, E_lj, pbc_flag, box_length)
  use lj_param, only : eps4
  implicit none
  ! global
  integer ni, nj, nsize, res_index(nsize)
  double precision x(3*nsize), p(nsize), force(3*nsize), E_lj
  double precision box_length(3)
  logical pbc_flag
  ! local
  integer ii, im, iatom, counter, iatom1, iatom2
  double precision delta(3), E_each, f(3)
  double precision half_length(3)

  half_length(:) = box_length(:) / 2.0d0

  force(:) = 0.0d0
  p(:) = 0.0d0
  E_lj = 0.0d0
  ! ni - ni
  do iatom = 1, ni
    if( res_index(iatom) .le. 0) cycle

    counter =0
    do ii = 1, ni - iatom
      if( res_index(iatom  + ii) .gt. 0)then
        exit
      else
        counter = counter + 1
      endif
    enddo

    do iatom1 = iatom, iatom + counter
      do iatom2 = iatom+counter +1 , ni
        if( (iatom2 .ge. iatom ).and.(iatom2 .le. iatom+counter))then
          cycle ! intra molecule
        endif
        delta(1:3) = x(3*iatom1-2:3*iatom1) - x(3*iatom2-2:3*iatom2)
        if( pbc_flag )then
          do im = 1, 3
            if(delta(im) .gt. half_length(im))then
              delta(im) = delta(im) - box_length(im)
            elseif( delta(im) .lt. -half_length(im))then
              delta(im) = delta(im) + box_length(im)
            endif
          enddo
        endif
        call calc_force_body(delta, f, E_each )
        force(3*iatom1-2:3*iatom1) = force(3*iatom1-2:3*iatom1) + f(1:3)
        force(3*iatom2-2:3*iatom2) = force(3*iatom2-2:3*iatom2) - f(1:3)
        ! calc potential energy
        p(iatom1) = p(iatom1) + 0.5d0*E_each
        p(iatom2) = p(iatom2) + 0.5d0*E_each
      enddo
    enddo
  enddo

  if( pbc_flag )then
    p =  eps4*p
    E_lj = sum(p(1:ni))
    return
  endif

  ! ni - nj
  do iatom1 = 1, ni
    do iatom2 = ni+1, nj
      delta(1:3) = x(3*iatom1-2:3*iatom1) - x(3*iatom2-2:3*iatom2)
      call calc_force_body(delta, f, E_each )
      force(3*iatom1-2:3*iatom1) = force(3*iatom1-2:3*iatom1) + f(1:3)
      ! calc potential energy
      p(iatom1) = p(iatom1) + 0.5d0*E_each
    enddo
  enddo
  p =  eps4*p
  E_lj = sum(p(1:ni))
end subroutine

subroutine output_energy_dat(dt, nstep,natom, ene)
  use energy_module
  implicit none
  ! global
  double precision dt
  integer nstep, natom
  type ( energy_var ) ene
  ! local
  type ( energy_var ) ene_local

  ene_local = ene
  ene_local%pot       = ene%pot       / dble(natom)
  ene_local%p_lj      = ene%p_lj      / dble(natom)
  ene_local%p_sp      = ene%p_sp      / dble(natom)
  ene_local%kin       = ene%kin       / dble(natom)

  if( nstep .eq. 1 )then
    open(31, file='01_energy.log')
    write(31,"(a50)",advance='no') "#old   1      2          3       4           5 "
    write(31,"(a50)") " 6        (/atom)   "
    write(31,"(a50)",advance='no') "#  time, total, pot_total, pot_lj, pot_spring,"
    write(31,"(a50)") " kin_total"
  endif
  write(31,"(6e20.12)") dt*dble(nstep), ene_local%pot+ene_local%kin,  ene_local%pot, ene_local%p_lj, ene_local%p_sp, ene_local%kin

end subroutine

subroutine output_energy_dat_fmm(dt, nstep,natom, ene, mpi_rank)
  use energy_module
  implicit none
  ! global
  double precision dt
  integer nstep, natom, mpi_rank
  type ( energy_var ) ene
  character(80) filename
  ! local
  type ( energy_var ) ene_local

  ene_local%pot       = ene%pot       / dble(natom)
  ene_local%p_lj      = ene%p_lj      / dble(natom)
  ene_local%p_sp      = ene%p_sp      / dble(natom)
  ene_local%kin       = ene%kin       / dble(natom)

  if( nstep .eq. 1 )then
    write(filename,'("05_mpi_ene_",i2.2, ".dat")') mpi_rank
    open(33, file=filename)
    write(33,"(a50)",advance='no') "#API   1      2          3       4           5 "
    write(33,"(a50)") " 6        (/atom)   "
    write(33,"(a50)",advance='no') "#  time, total, pot_total, pot_lj, pot_spring,"
    write(33,"(a50)") " kin_total"
  endif
  write(33,"(6e20.12)") dt*dble(nstep), ene_local%pot+ene_local%kin,  ene_local%pot, ene_local%p_lj, ene_local%p_sp, ene_local%kin

end subroutine

subroutine output_xyz_fmm(ni, nj, nsize, x, box_length, ene, irank)
  use energy_module
  implicit none
  ! global
  integer ni, nj, nsize, irank
  double precision x(3*nsize), box_length(3)
  type ( energy_var ) ene
  ! local
  integer ii
  double precision half_length(1:3), temp_coord(3)
  double precision , parameter :: factor = 4.0
  character(80) filename

  half_length = factor*box_length / 2.0d0

  write(filename,'("06_test_",i3.3, ".xyz")') irank
  open(35+irank,file=filename)

  write(35+irank,"(i6)") nj+8
  write(35+irank,"(3f15.5)") ene%pot, ene%kin, ene%pot + ene%kin
  write(35+irank,"(a3,3f15.3)") "Xe " ,   half_length(1) ,   half_length(2) ,  half_length(3)
  write(35+irank,"(a3,3f15.3)") "Xe " ,   half_length(1) , - half_length(2) ,  half_length(3)
  write(35+irank,"(a3,3f15.3)") "Xe " , - half_length(1) ,   half_length(2) ,  half_length(3)
  write(35+irank,"(a3,3f15.3)") "Xe " , - half_length(1) , - half_length(2) ,  half_length(3)
  write(35+irank,"(a3,3f15.3)") "Xe " ,   half_length(1) ,   half_length(2) ,- half_length(3)
  write(35+irank,"(a3,3f15.3)") "Xe " ,   half_length(1) , - half_length(2) ,- half_length(3)
  write(35+irank,"(a3,3f15.3)") "Xe " , - half_length(1) ,   half_length(2) ,- half_length(3)
  write(35+irank,"(a3,3f15.3)") "Xe " , - half_length(1) , - half_length(2) ,- half_length(3)
  do ii = 1, ni
    temp_coord(:) = factor*x(3*ii-2:3*ii)
    write(35+irank,"(a3,3f15.3)") "He ", temp_coord(1) , temp_coord(2), temp_coord(3)
  enddo
  do ii = ni+1, nj
    temp_coord(:) = factor*x(3*ii-2:3*ii)
    write(35+irank,"(a3,3f15.3)") "Ar ", temp_coord(1) , temp_coord(2), temp_coord(3)
  enddo
end subroutine

subroutine output_checkpoint(filename, natom, nsize, res_index, &
    x1, v1, q)
  implicit none
  ! global
  character(80) filename
  integer natom, nsize, res_index(nsize)
  double precision x1(3*nsize), v1(3*nsize), q(nsize)
  ! local
  integer ii

  open(32, file= filename)
  do ii = 1, natom
    write(32,"(i8,7f18.5)") res_index(ii), &
      x1(3*ii-2), x1(3*ii-1), x1(3*ii), v1(3*ii-2), v1(3*ii-1), v1(3*ii), q(ii)
  enddo

  close(32)
end subroutine

subroutine calc_force_body(delta, f, E_each)
  use lj_param, only : sig6, eps24, cutoff2
  implicit none
  ! global
  double precision delta(3), f(3), E_each
  ! local
  double precision r2_inv, sig6_r6_inv

  f = 0.0d0
  E_each = 0.0d0
  r2_inv = 1.0d0/(sum( delta(1:3)*delta(1:3)))
  if( r2_inv .le. cutoff2 ) return
  sig6_r6_inv = sig6*r2_inv*r2_inv*r2_inv
  ! potential = 4*epsi*( (sigma/r)**12 - (sigma/r)**6 )
  ! force     = 4*epsi*( 12*sigma_12*(r**13) - 6*sigma_6*(r**7))
  !           = 24*(epsi/r)* ( 2*sigma_12*(r**12) - sigma_6*(r**6))
  !f(1:3) = delta(1:3)*(24*epsi*(1/r**2))*(sig_6/r**6)*(2*sig_6*(1/r**6) - 1.0d0 )
  f(1:3) = delta(1:3)*(eps24*r2_inv)*sig6_r6_inv*( 2.0d0*sig6_r6_inv - 1.0d0)
  E_each = sig6_r6_inv*(sig6_r6_inv - 1.0d0)

end subroutine

