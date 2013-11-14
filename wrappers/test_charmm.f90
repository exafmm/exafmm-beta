module charmm_io
contains
  subroutine charmm_cor_read(n,x,q,size,filename,numex,natex,nat,atype,&
       rscale,gscale,fgscale,nbonds,ntheta,ib,jb,it,jt,kt,rbond,cbond,&
       aangle,cangle,mass,xc,v,xnew,nres,ires)
    implicit none
    logical qext
    integer i, nat, im, in, n, natex_size, vdw_size, idist, resn, nres
    real(8) size, sizex, sizey, sizez
    integer, allocatable, dimension(:) :: numex,natex,atype,ib,jb,it,jt,kt,ires
    real(8), allocatable, dimension(:) :: x, q, rscale, gscale, fgscale
    real(8), allocatable, dimension(:) :: mass,xc,v,xnew
    real(8), allocatable, dimension(:,:) :: rbond,cbond
    real(8), allocatable, dimension(:,:,:) :: aangle,cangle
    integer nbonds,ntheta
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
    allocate( x(3*n), q(n), ires(n)) ! too much for ires!!
    ! read also the residue info from cor file
    ! residue info is enough for pure water system
    ! for protein + water the segments are needed, too!!!
    resn=1
    ires(1:n)=1 ! at which atom the residue starts
    do i=1,n
       if(qext)then
          read(1,'(2i10,20x,3f20.10,20x,f20.10)') &
               in,im,x(i*3-2),x(i*3-1),x(i*3),q(i)
       else
          read(1,'(2i5,10x,3f10.5,10x,f10.5)') &
               in,im,x(i*3-2),x(i*3-1),x(i*3),q(i)
       endif
       if(im/=resn)then
          resn=resn+1
          ires(resn)=i
       endif
    enddo
    nres=resn
    !write(*,*)nres,ires(1:nres)
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
    read(1,'(a100)',end=99)lin
    ! we maybe able to read the bonding info here
    ! it is used for dynamics only
    if(lin(1:5)=='BONDS')then
       read(lin(6:100),*)nbonds
       allocate(ib(nbonds),jb(nbonds))
       read(1,'(10i10)')(ib(i),jb(i),i=1,nbonds)
    endif
    read(1,'(a100)',end=99)lin
    if(lin(1:6)=='ANGLES')then
       read(lin(7:100),*)ntheta
       allocate(it(ntheta),jt(ntheta),kt(ntheta))
       read(1,'(9i10)')(it(i),jt(i),kt(i),i=1,ntheta)
    endif
    ! parameters...
    allocate(rbond(nat,nat),cbond(nat,nat),mass(nat))
    allocate(aangle(nat,nat,nat),cangle(nat,nat,nat))
    read(1,'(a100)',end=99)lin
    read(1,'(5f20.10)')rbond(1:nat,1:nat)
    read(1,'(a100)',end=99)lin
    read(1,'(5f20.10)')cbond(1:nat,1:nat)
    read(1,'(a100)',end=99)lin
    read(1,'(5f20.10)')aangle(1:nat,1:nat,1:nat)
    read(1,'(a100)',end=99)lin
    read(1,'(5f20.10)')cangle(1:nat,1:nat,1:nat)
    read(1,'(a100)',end=99)lin
    read(1,'(5f20.10)')mass(1:nat)
    ! read all vectors from CHARMM restart file:
    ! coordinates from restart are needed for proper continuation!
    ! dual coordinates are nedeed anyway for image centering problem!
    ! equilibrated velocities from CHARMM restart
    ! and the coordinate corrections /xnew/ ??? Are they useful ???
    allocate(xc(3*n),v(3*n),xnew(3*n))
    read(1,'(3d28.18)')(xc(3*i-2),xc(3*i-1),xc(3*i),i=1,n)
    read(1,'(3d28.18)')(v(3*i-2),v(3*i-1),v(3*i),i=1,n)
    read(1,'(3d28.18)')(xnew(3*i-2),xnew(3*i-1),xnew(3*i),i=1,n)
!!!    
99  continue
    return
  end subroutine charmm_cor_read

  subroutine bonded_terms(nglobal,icpumap,nat,atype,x,f,nbonds,ntheta,&
       ib,jb,it,jt,kt,rbond,cbond,aangle,cangle,eb,et)
    implicit none
    integer nglobal,nat
    real(8), allocatable, dimension(:) :: x, f
    integer i,ii,jj,kk,nbonds,ntheta
    real(8), optional :: eb,et
    real(8) dx,dy,dz,r,db,df
    integer, allocatable, dimension(:) :: ib,jb,it,jt,kt,icpumap,atype
    real(8), allocatable, dimension(:,:) :: rbond,cbond
    real(8), allocatable, dimension(:,:,:) :: aangle,cangle
    real(8) dx1,dy1,dz1,dx2,dy2,dz2,r1,r2,r1r,r2r
    real(8) dx1r,dy1r,dz1r,dx2r,dy2r,dz2r,cst,at,da
    real(8) dxi,dyi,dzi,st2r,str,dtx1,dtx2,dty1,dty2,dtz1,dtz2,dt1,dt2
    if (present(eb)) then
       eb = 0.0
       do i = 1, nbonds
          ii=ib(i)
          jj=jb(i)
          if (icpumap(ii) == 0) cycle
          dx=x(3*ii-2)-x(3*jj-2)
          dy=x(3*ii-1)-x(3*jj-1)
          dz=x(3*ii)-x(3*jj)
          r=sqrt(dx*dx+dy*dy+dz*dz)
          db=r-rbond(atype(ii),atype(jj))
          df=cbond(atype(ii),atype(jj))*db
          eb=eb+df*db
          df=2.0*df/r
          dxi=dx*df
          dyi=dy*df
          dzi=dz*df
          f(3*ii-2)=f(3*ii-2)+dxi
          f(3*ii-1)=f(3*ii-1)+dyi
          f(3*ii)=f(3*ii)+dzi
          f(3*jj-2)=f(3*jj-2)-dxi
          f(3*jj-1)=f(3*jj-1)-dyi
          f(3*jj)=f(3*jj)-dzi
       enddo
    endif

    if (present(et)) then
       et = 0.0
       do i = 1, ntheta
          ii=it(i)
          if (icpumap(ii) == 0) cycle
          jj=jt(i)
          kk=kt(i)
          dx1=x(3*ii-2)-x(3*jj-2)
          dy1=x(3*ii-1)-x(3*jj-1)
          dz1=x(3*ii)-x(3*jj)
          dx2=x(3*kk-2)-x(3*jj-2)
          dy2=x(3*kk-1)-x(3*jj-1)
          dz2=x(3*kk)-x(3*jj)
          r1=sqrt(dx1*dx1+dy1*dy1+dz1*dz1)
          r2=sqrt(dx2*dx2+dy2*dy2+dz2*dz2)
          r1r=1.0/r1
          r2r=1.0/r2
          dx1r=dx1*r1r
          dy1r=dy1*r1r
          dz1r=dz1*r1r
          dx2r=dx2*r2r
          dy2r=dy2*r2r
          dz2r=dz2*r2r
          cst=dx1r*dx2r+dy1r*dy2r+dz1r*dz2r
          at=acos(cst)
          da=at-aangle(atype(ii),atype(jj),atype(kk))
          df=da*cangle(atype(ii),atype(jj),atype(kk))
          et=et+df*da
          st2r=1.0/(1.0-cst*cst)
          str=sqrt(st2r)
          df=-2.0*df*str
          dtx1=r1r*(dx2r-cst*dx1r)
          dtx2=r2r*(dx1r-cst*dx2r)
          dty1=r1r*(dy2r-cst*dy1r)
          dty2=r2r*(dy1r-cst*dy2r)
          dtz1=r1r*(dz2r-cst*dz1r)
          dtz2=r2r*(dz1r-cst*dz2r)
          dt1=df*dtx1
          dt2=df*dtx2
          f(3*ii-2)=f(3*ii-2)+dt1
          f(3*kk-2)=f(3*kk-2)+dt2
          f(3*jj-2)=f(3*jj-2)-dt1-dt2
          dt1=df*dty1
          dt2=df*dty2
          f(3*ii-1)=f(3*ii-1)+dt1
          f(3*kk-1)=f(3*kk-1)+dt2
          f(3*jj-1)=f(3*jj-1)-dt1-dt2
          dt1=df*dtz1
          dt2=df*dtz2
          f(3*ii)=f(3*ii)+dt1
          f(3*kk)=f(3*kk)+dt2
          f(3*jj)=f(3*jj)-dt1-dt2
       enddo
    endif

  end subroutine bonded_terms

  subroutine energy(nglobal,nat,nbonds,ntheta,ksize,&
       alpha,sigma,cutoff,cuton,ccelec,pcycle,&
       xc,p,f,q,v,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
       ib,jb,it,jt,kt,atype,icpumap,numex,natex,etot,eb,et,efmm,evdw)

    implicit none
    integer nglobal,nat,nbonds,ntheta,ksize
    real(8),optional :: eb,et,efmm,evdw
    real(8) alpha,sigma,cutoff,cuton,ccelec,etot,pcycle
    real(8), allocatable, dimension(:) :: x,xc,p,f,fl,q,v,gscale,fgscale,rscale
    real(8), allocatable, dimension(:,:) :: rbond,cbond
    real(8), allocatable, dimension(:,:,:) :: aangle,cangle
    integer, allocatable, dimension(:) :: ib,jb,it,jt,kt,atype,icpumap,numex,natex
    integer i,j

    ! zero the force
    f(1:3*nglobal) = 0.0
    allocate(x(3*nglobal))

    ! there is a prblem here: I need coordinates from surrounding region,
    ! not only from icpumap atoms to perform the bond calculations.
    ! And they must be by residue, which they are not

    !x(1:3*nglobal) = xc(1:3*nglobal) !copy coordinates
    call fmm_partition(nglobal, icpumap, xc, q, v, pcycle)

    if (present(eb).and.present(et)) then
       call bonded_terms(nglobal,icpumap,nat,atype,xc,f,nbonds,ntheta,&
            ib,jb,it,jt,kt,rbond,cbond,aangle,cangle,eb,et)
    elseif (present(eb)) then
       call bonded_terms(nglobal,icpumap,nat,atype,xc,f,nbonds,ntheta,&
            ib,jb,it,jt,kt,rbond,cbond,aangle,cangle,eb=eb)
    elseif (present(et)) then
       call bonded_terms(nglobal,icpumap,nat,atype,xc,f,nbonds,ntheta,&
            ib,jb,it,jt,kt,rbond,cbond,aangle,cangle,et=et)
    endif

    if(present(efmm)) then
       allocate(fl(3*nglobal))
       fl(1:3*nglobal)=0.0
       p(1:nglobal)=0.0
       call fmm_coulomb(nglobal, icpumap, xc, q, p, fl, pcycle)
       !call ewald_coulomb(nglobal, icpumap, xc, q, p, fl, ksize, alpha, sigma, cutoff, pcycle)
       call coulomb_exclusion(nglobal, icpumap, xc, q, p, fl, pcycle, numex, natex)
       efmm=0.0
       do i=1, nglobal
          if (icpumap(i) == 0) cycle
          efmm=efmm+p(i)*q(i)
          f(3*i-2)=f(3*i-2)+fl(3*i-2)*q(i)*ccelec
          f(3*i-1)=f(3*i-1)+fl(3*i-1)*q(i)*ccelec
          f(3*i)=f(3*i)+fl(3*i)*q(i)*ccelec
       enddo
       efmm=efmm*ccelec*0.5
       deallocate(fl)
    endif
    if (present(evdw)) then
       allocate(fl(3*nglobal))
       p(1:nglobal)=0.0
       fl(1:3*nglobal)=0.0
       call fmm_vanderwaals(nglobal, icpumap, atype, xc, p, fl, cuton, cutoff,&
            pcycle, nat, rscale, gscale, fgscale)
       call vanderwaals_exclusion(nglobal, icpumap, atype, xc, p, fl, cuton, cutoff,&
            pcycle, nat, rscale, gscale, fgscale, numex, natex)
       evdw=0.0
       do i = 1, nglobal
          if (icpumap(i) == 0) cycle
          evdw=evdw+p(i)
          f(3*i-2)=f(3*i-2)-fl(3*i-2)
          f(3*i-1)=f(3*i-1)-fl(3*i-1)
          f(3*i)=f(3*i)-fl(3*i)
       enddo
       evdw=evdw*0.5
       deallocate(fl)
    endif
    etot=0.0
    if(present(eb)) etot=etot+eb
    if(present(et)) etot=etot+et
    if(present(efmm)) etot=etot+efmm
    if(present(evdw)) etot=etot+evdw
    deallocate(x)
    return
  end subroutine energy

  subroutine force_testing(nglobal,nat,nbonds,ntheta,ksize,&
       alpha,sigma,cutoff,cuton,ccelec,pcycle,&
       xc,p,f,q,v,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
       ib,jb,it,jt,kt,atype,icpumap,numex,natex,etot,eb,et,efmm,evdw)

    implicit none
    integer nglobal,nat,nbonds,ntheta,ksize
    real(8) eb,et,efmm,evdw
    real(8) alpha,sigma,cutoff,cuton,ccelec,etot,pcycle
    real(8), allocatable, dimension(:) :: x,xc,p,f,q,v,gscale,fgscale,rscale
    real(8), allocatable, dimension(:,:) :: rbond,cbond
    real(8), allocatable, dimension(:,:,:) :: aangle,cangle
    integer, allocatable, dimension(:) :: ib,jb,it,jt,kt,atype,icpumap,numex,natex
    real(8), allocatable, dimension(:) :: f_analytic
    real(8) xsave, e0, step, step2, eplus, eminus, nforce
    integer i,j,istart,iend

    call energy(nglobal,nat,nbonds,ntheta,ksize,&
         alpha,sigma,cutoff,cuton,ccelec,pcycle,&
         xc,p,f,q,v,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
         ib,jb,it,jt,kt,atype,icpumap,numex,natex,etot,eb,et,efmm,evdw)

!!!
!!! this routine is not parallelized ??
!!!
!!! And also it doesn't test forces but derivatives!!!
!!!
    e0=etot
    allocate(f_analytic(3*nglobal))
    f_analytic(1:3*nglobal)=f(1:3*nglobal)
    istart=1
    iend=10
    step=0.0001
    step2=1.0/(2.0*step)

    write(*,*)'   I      J             Numeric           Analytic             Difference'
    do i = istart, iend
       do j = 1, 3
          xsave=xc(3*(i-1)+j) ! save it
          xc(3*(i-1)+j)=xc(3*(i-1)+j)-step

          call energy(nglobal,nat,nbonds,ntheta,ksize,&
               alpha,sigma,cutoff,cuton,ccelec,pcycle,&
               xc,p,f,q,v,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
               ib,jb,it,jt,kt,atype,icpumap,numex,natex,eminus,eb,et,efmm,evdw)

          xc(3*(i-1)+j)=xc(3*(i-1)+j)+2.0*step

          call energy(nglobal,nat,nbonds,ntheta,ksize,&
               alpha,sigma,cutoff,cuton,ccelec,pcycle,&
               xc,p,f,q,v,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
               ib,jb,it,jt,kt,atype,icpumap,numex,natex,eplus,eb,et,efmm,evdw)

          nforce=step2*(eplus-eminus)
          xc(3*(i-1)+j)=xsave  ! restore
          write(*,'(2i7,3f20.9)')i,j,nforce,f_analytic(3*(i-1)+j),nforce-f_analytic(3*(i-1)+j)
       end do
    end do

    deallocate(f_analytic)

    return
  end subroutine force_testing


  subroutine print_energy(timstart,nglobal,nat,nbonds,ntheta,ksize,&
       alpha,sigma,cutoff,cuton,ccelec,pcycle,&
       xc,p,f,q,v,mass,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
       ib,jb,it,jt,kt,atype,icpumap,numex,natex,etot,eb,et,efmm,evdw)

    use mpi

    implicit none
    integer nglobal,nat,nbonds,ntheta,ksize
    real(8),optional :: eb,et,efmm,evdw,timstart
    real(8) alpha,sigma,cutoff,cuton,ccelec,etot,pcycle
    real(8), allocatable, dimension(:) :: x,xc,p,f,fl,q,v,mass,gscale,fgscale,rscale
    real(8), allocatable, dimension(:,:) :: rbond,cbond
    real(8), allocatable, dimension(:,:,:) :: aangle,cangle
    integer, allocatable, dimension(:) :: ib,jb,it,jt,kt,atype,icpumap,numex,natex
    integer i,j, ierr, mpirank
    real(8) temp,kboltz,ekinetic,grms,ekineticglob,grmsglob,etotglob

    kboltz=1.987191d-03 !from CHARMM

    call energy(nglobal,nat,nbonds,ntheta,ksize,&
         alpha,sigma,cutoff,cuton,ccelec,pcycle,&
         xc,p,f,q,v,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
         ib,jb,it,jt,kt,atype,icpumap,numex,natex,etot,eb,et,efmm,evdw)

    ! calculate kinetic energy, temperature:
    ekinetic=0.0
    grms=0.0
    do i=1, nglobal
       if(icpumap(i) == 0) cycle
       ekinetic=ekinetic+mass(atype(i))*(v(3*i-2)**2+v(3*i-1)**2+v(3*i)**2)
       grms=grms+f(3*i-2)**2+f(3*i-1)**2+f(3*i)**2
    enddo
    ! Summ up the terms in parallel
    call mpi_reduce(etot, etotGlob,  1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr)
    call mpi_reduce(ekinetic, ekineticGlob,  1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr)
    call mpi_reduce(grms, grmsGlob,  1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr)
    grmsglob=sqrt(grmsglob/3.0/real(nglobal))
    temp=ekineticglob/3.0/nglobal/kboltz
    ekineticglob=ekineticglob/2.0
    call mpi_comm_rank(mpi_comm_world, mpirank, ierr)
    if(mpirank==0) &
         write(*,'(''time:'',f9.3,'' Etotal:'',f14.5,'' Ekin:'',f14.5,'' Epot:'',f14.5,'' T:'',f12.3,'' Grms:'',f12.5)')&
         timstart,etotglob+ekineticglob,ekineticglob,etotglob,temp,grmsglob

    return
  end subroutine print_energy

  subroutine run_dynamics(dynsteps,imcentfrq,printfrq,&
       nglobal,nat,nbonds,ntheta,ksize,&
       alpha,sigma,cutoff,cuton,ccelec,pcycle,&
       xc,p,f,q,v,mass,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
       ib,jb,it,jt,kt,atype,icpumap,numex,natex,etot,eb,et,efmm,evdw,nres,ires)
    use mpi
    implicit none
    integer dynsteps,nglobal,nat,nbonds,ntheta,ksize,imcentfrq,printfrq,nres
    real(8) eb,et,efmm,evdw
    real(8) alpha,sigma,cutoff,cuton,ccelec,etot,pcycle
    real(8), allocatable, dimension(:) :: x,v,mass,xc,p,f,q,gscale,fgscale,rscale
    real(8), allocatable, dimension(:,:) :: rbond,cbond
    real(8), allocatable, dimension(:,:,:) :: aangle,cangle
    integer, allocatable, dimension(:) :: ib,jb,it,jt,kt,atype,icpumap,numex,natex,ires
    real(8), allocatable, dimension(:) :: xnew,vnew,fnew,fac1,fac2
    real(8) xsave, e0, step, step2, eplus, eminus, nforce,tstep,timstart,time,tstep2
    integer i,j,istart,iend,unit,istep,ierr,mpirank
    real(8),parameter :: TIMFAC=4.88882129D-02
    logical leapfrog_maybe

    unit=1
    call mpi_comm_rank(mpi_comm_world, mpirank, ierr)
    if(mpirank==0)open(unit=unit,file='water.pdb',status='new')

    leapfrog_maybe=.false.
!!!
!!! this synchronized seems to be working:
!!! http://en.wikipedia.org/wiki/Leapfrog_integration
!!! also do Bernie's algorithm from CHARMM
!!!
    tstep=0.001/timfac !ps -> akma
    tstep2=tstep**2
    timstart = 100.0 ! first 100ps was equilibration with standard CHARMM
    time = timstart

    allocate(xnew(3*nglobal),vnew(3*nglobal),fnew(3*nglobal)) ! do we need first two ??
    allocate(fac1(nglobal),fac2(nglobal))
    do i=1,nglobal
       if(icpumap(i)==0)cycle
       fac1(i) = tstep2/mass(atype(i))/2.0
       fac2(i) = tstep/mass(atype(i))/2.0
    enddo

    call energy(nglobal,nat,nbonds,ntheta,ksize,&
         alpha,sigma,cutoff,cuton,ccelec,pcycle,&
         xc,p,f,q,v,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
         ib,jb,it,jt,kt,atype,icpumap,numex,natex,etot,eb,et,efmm,evdw)

    do istep = 1, dynsteps

       if (leapfrog_maybe) then
          call energy(nglobal,nat,nbonds,ntheta,ksize,&
               alpha,sigma,cutoff,cuton,ccelec,pcycle,&
               xc,p,f,q,v,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
               ib,jb,it,jt,kt,atype,icpumap,numex,natex,etot,eb,et,efmm,evdw)
          do j = 1, nglobal
             if(icpumap(j)==0)cycle
             vnew(3*j-2) = v(3*j-2)  - f(3*j-2)*tstep/mass(atype(j))
             vnew(3*j-1) = v(3*j-1)  - f(3*j-1)*tstep/mass(atype(j))
             vnew(3*j)   = v(3*j)    - f(3*j)*tstep/mass(atype(j))
             xnew(3*j-2) = xc(3*j-2) + vnew(3*j-2)*tstep
             xnew(3*j-1) = xc(3*j-1) + vnew(3*j-1)*tstep
             xnew(3*j)   = xc(3*j)   + vnew(3*j)*tstep
          enddo
       else
          do j = 1, nglobal
             if(icpumap(j)==0)cycle
             xnew(3*j-2) = xc(3*j-2) + v(3*j-2)*tstep - f(3*j-2)*fac1(j)
             xnew(3*j-1) = xc(3*j-1) + v(3*j-1)*tstep - f(3*j-1)*fac1(j)
             xnew(3*j)   = xc(3*j)   + v(3*j)*tstep   - f(3*j)*fac1(j)
          enddo

          call energy(nglobal,nat,nbonds,ntheta,ksize,&
               alpha,sigma,cutoff,cuton,ccelec,pcycle,&
               xnew,p,fnew,q,v,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
               ib,jb,it,jt,kt,atype,icpumap,numex,natex,etot,eb,et,efmm,evdw)

 !!!         vnew(1:3*nglobal) = 0.0   ! to use mpi_allreduce()
          do j = 1, nglobal
             if(icpumap(j)==0)cycle
             vnew(3*j-2) = v(3*j-2)  - fac2(j)*(f(3*j-2) + fnew(3*j-2))
             vnew(3*j-1) = v(3*j-1)  - fac2(j)*(f(3*j-1) + fnew(3*j-1))
             vnew(3*j)   = v(3*j)    - fac2(j)*(f(3*j)   + fnew(3*j))
          enddo
          do j = 1, nglobal
             if(icpumap(j)==0)cycle
             f(3*j-2)=fnew(3*j-2)
             f(3*j-1)=fnew(3*j-1)
             f(3*j)  =fnew(3*j)
          enddo
       endif

!!! These copies are not really neded :-(
       do j = 1, nglobal
          if(icpumap(j)==1) then
             xc(3*j-2)= xnew(3*j-2)
             xc(3*j-1)= xnew(3*j-1)
             xc(3*j)  = xnew(3*j)
             v(3*j-2) = vnew(3*j-2)
             v(3*j-1) = vnew(3*j-1)
             v(3*j)   = vnew(3*j)
          endif
       enddo

       ! FIXME: WE MUST GET RID OF THIS MPI_ALLREDUCE()!!!!!
       ! velocities communication for parallel: It must be done in FMM library ???
       ! broadcast/globalsum would be OK for small number of processes for testing purposes:
!!!       call mpi_allreduce(vnew, v , 3*nglobal, mpi_real8, mpi_sum, mpi_comm_world, ierr)

       if (mod(istep,imcentfrq) == 0) call image_center(nglobal,xc,nres,ires,pcycle,icpumap)

       time=time+tstep*timfac  ! for printing only
       if (mod(istep,printfrq) == 0) call print_energy(time,nglobal,nat,nbonds,ntheta,ksize,&
            alpha,sigma,cutoff,cuton,ccelec,pcycle,&
            xc,p,f,q,v,mass,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
            ib,jb,it,jt,kt,atype,icpumap,numex,natex,etot,eb,et,efmm,evdw)

       if (mod(istep,printfrq) == 0) call pdb_frame(unit,time,nglobal,xc,nres,ires,icpumap)

    enddo

    deallocate(xnew,vnew,fnew)

  end subroutine run_dynamics

  subroutine pdb_frame(unit,time,nglobal,xc,nres,ires,icpumap)
    use mpi
    implicit none
    integer nglobal,nres,ipt,i,unit,ierr, mpirank
    integer,allocatable,dimension(:) :: ires,icpumap
    real(8),allocatable,dimension(:) :: xc,xl
    real(8) time

    ! Broadcast the coordinates for parallel
    allocate(xl(3*nglobal))
    xl(1:3*nglobal)=0.0
    do i=1,nglobal
       if(icpumap(i)==1) xl(3*i-2:3*i) = xc(3*i-2:3*i)
    enddo
    call mpi_allreduce(xl, xc, 3*nglobal, mpi_real8, mpi_sum, mpi_comm_world, ierr)
    deallocate(xl)

    call mpi_comm_rank(mpi_comm_world, mpirank, ierr)
    if (mpirank /= 0) return

    write(unit,'(''REMARK frame at time:'',f14.4,'' ps'')')time
    ipt=0
    do i = 1, nres
       ipt=ipt+1
       write(unit,'(a6,i5,2x,a3,1x,a3,1x,a1,i4,3x,3f8.3,2f6.2,11x,a1)')&
            'HETATM',ipt,'OH2','HOH','w',i,xc(3*ipt-2),xc(3*ipt-1),xc(3*ipt),1.0,1.0,'O'
       ipt=ipt+1
       write(unit,'(a6,i5,2x,a3,1x,a3,1x,a1,i4,3x,3f8.3,2f6.2,11x,a1)')&
            'HETATM',ipt,'HO1','HOH','w',i,xc(3*ipt-2),xc(3*ipt-1),xc(3*ipt),1.0,1.0,'H'
       ipt=ipt+1
       write(unit,'(a6,i5,2x,a3,1x,a3,1x,a1,i4,3x,3f8.3,2f6.2,11x,a1)')&
            'HETATM',ipt,'HO2','HOH','w',i,xc(3*ipt-2),xc(3*ipt-1),xc(3*ipt),1.0,1.0,'H'
    enddo
    write(unit,'(''END'')')
    return

  end subroutine pdb_frame

  subroutine image_center(nglobal,xc,nres,ires,pcycle,icpumap)
    implicit none
    integer nglobal,nres,flag
    integer, allocatable, dimension(:) :: ires,icpumap
    real(8), allocatable, dimension(:) :: xc
    real(8) pcycle,xmin,xmax,ymin,ymax,zmin,zmax,xcen,ycen,zcen
    integer istart, iend, nsel, i, j

    ! only simple translations, so velocities stay the same! "?"

    xmin=-pcycle/2.0
    ymin=-pcycle/2.0
    zmin=-pcycle/2.0
    xmax=pcycle/2.0
    ymax=pcycle/2.0
    zmax=pcycle/2.0

    residues: do i = 1, nres
       istart=ires(i)
       if(i==nres) then
          iend=nglobal
       else
          iend=ires(i+1)-1
       endif
!!!
!!! for parallel we need to treat them by residues
!!! if one of the atoms is in my processor, I need to move
!!! the whole residue!! No extra communication needed :-)
!!!
       flag=0
       do j = istart, iend
          if(icpumap(j)==1)flag=1
       enddo
       if (flag==0) cycle residues
!!!
       xcen=0.0
       ycen=0.0
       zcen=0.0
       nsel=iend-istart+1
       do j = istart, iend
          xcen=xcen+xc(3*j-2)
          ycen=ycen+xc(3*j-1)
          zcen=zcen+xc(3*j)
       enddo
       xcen=xcen/real(nsel)
       ycen=ycen/real(nsel)
       zcen=zcen/real(nsel)
       if(xcen<xmin)then
          do j = istart,iend
             xc(3*j-2)=xc(3*j-2)+pcycle
          enddo
       endif
       if(xcen>xmax)then
          do j = istart,iend
             xc(3*j-2)=xc(3*j-2)-pcycle
          enddo
       endif
       if(ycen<ymin)then
          do j = istart,iend
             xc(3*j-1)=xc(3*j-1)+pcycle
          enddo
       endif
       if(ycen>ymax)then
          do j = istart,iend
             xc(3*j-1)=xc(3*j-1)-pcycle
          enddo
       endif
       if(zcen<zmin)then
          do j = istart,iend
             xc(3*j)=xc(3*j)+pcycle
          enddo
       endif
       if(zcen>zmax)then
          do j = istart,iend
             xc(3*j)=xc(3*j)-pcycle
          enddo
       endif
    enddo residues

  end subroutine image_center

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
  character(len=128) filename
  integer dynsteps
  integer d, i, ierr, images, ista, iend, istat, ksize, lnam, mpirank, mpisize
  integer nat, nglobal, prange, ix, iy, iz, j, verbose
  real(8) alpha, sigma, cuton, cutoff, average, norm, pcycle, pi, ccelec
  real(8) coef, dx, dy, dz, fx, fy ,fz, pp, R2, R3inv, Rinv, theta
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
  integer, allocatable, dimension(:) :: icpumap,numex,natex,atype,ib,jb,it,jt,kt,ires
  real(8), allocatable, dimension(:) :: x, q, p, f, p2, f2, v, mass,xc,xnew
  real(8), allocatable, dimension(:) :: rscale, gscale, fgscale
  real(8), allocatable, dimension(:,:) :: rbond,cbond
  real(8), allocatable, dimension(:,:,:) :: aangle,cangle
  integer nbonds,ntheta,imcentfrq,printfrq,nres
  real(8) efmm, evdw, etot, eb,et,timstart
  logical test_first,integrate
  parameter(pi=3.14159265358979312d0, ccelec=332.0716d0)

  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, mpisize, ierr)
  call mpi_comm_rank(mpi_comm_world, mpirank, ierr)
  nglobal = 1000
  images = 3
  theta = 0.4
  verbose = 0
  ksize = 11
  pcycle = 10 * pi
  sigma = .25 / pi
  cuton = 9.5
  cutoff = 10.0
  alpha = 10 / pcycle
  nat = 16
  CHARMMIO: if (command_argument_count() > 0) then
     call get_command_argument(1,filename,lnam,istat)
     call charmm_cor_read(nglobal,x,q,pcycle,filename,numex,natex,nat,atype,&
          rscale,gscale,fgscale,nbonds,ntheta,ib,jb,it,jt,kt,rbond,cbond,&
          aangle,cangle,mass,xc,v,xnew,nres,ires)
     allocate( p(nglobal),  f(3*nglobal), icpumap(nglobal) )
     allocate( p2(nglobal), f2(3*nglobal) )
     alpha = 10 / pcycle

  else
     allocate( x(3*nglobal),  q(nglobal),  p(nglobal),  f(3*nglobal), icpumap(nglobal) )
     allocate( p2(nglobal), f2(3*nglobal) )
     allocate( numex(nglobal), natex(nglobal), atype(nglobal) )
     allocate( rscale(nat*nat), gscale(nat*nat), fgscale(nat*nat) )
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
     call random_number(rscale)
     call random_number(gscale)
     do i = 1, nglobal
        numex(i) = 1
        if(mod(i,2).eq.1)then
           natex(i) = i+1
        else
           natex(i) = i-1
        endif
        atype(i) = rand() * nat + 1
     enddo
     do i = 1, nat*nat
        rscale(i) = 1
        gscale(i) = 0.0001
        fgscale(i) = gscale(i)
     enddo
  end if CHARMMIO
  ista = 1
  iend = nglobal
  call split_range(ista,iend,mpirank,mpisize)
  do i = ista, iend
     icpumap(i) = 1
  end do
  call fmm_init(images,theta,verbose)
  call fmm_partition(nglobal, icpumap, x, q, v, pcycle)
  call fmm_coulomb(nglobal, icpumap, x, q, p, f, pcycle)
  do i = 1, nglobal
     p2(i) = 0
     f2(3*i-2) = 0
     f2(3*i-1) = 0
     f2(3*i-0) = 0
  end do
#if 1
  call ewald_coulomb(nglobal, icpumap, x, q, p2, f2, ksize, alpha, sigma, cutoff, pcycle)
#else
  call direct_coulomb(nglobal, icpumap, x, q, p2, f2, pcycle)
#endif
  call coulomb_exclusion(nglobal, icpumap, x, q, p, f, pcycle, numex, natex)
  call coulomb_exclusion(nglobal, icpumap, x, q, p2, f2, pcycle, numex, natex)

  !print for charmm
  !do i = 1, nglobal
  !   write(1001,'(f20.10)')p2(i)
  !enddo

  potSum = 0
  potSum2 = 0
  accDif = 0
  accNrm = 0
  accNrm2 = 0
  do i = 1, nglobal
     if (icpumap(i).eq.1) then
        p(i) = p(i) * q(i)
        f(3*i-2) = f(3*i-2) * q(i)
        f(3*i-1) = f(3*i-1) * q(i)
        f(3*i-0) = f(3*i-0) * q(i)
        p2(i) = p2(i) * q(i)
        f2(3*i-2) = f2(3*i-2) * q(i)
        f2(3*i-1) = f2(3*i-1) * q(i)
        f2(3*i-0) = f2(3*i-0) * q(i)
        potSum  = potSum  + p(i)
        potSum2 = potSum2 + p2(i)
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
     print"(a)",'--- Coulomb FMM vs. Ewald -------'
     print"(a,f9.7)",'Rel. L2 Error (pot)  : ', sqrt(potDifGlob/potNrmGlob2)
     print"(a,f9.7)",'Rel. L2 Error (acc)  : ', sqrt(accDifGlob/accNrmGlob2)
     print"(a,f12.4)",'Energy (FMM)         : ', ccelec*potSumGlob/2.0
     print"(a,f12.4)",'Energy (Ewald)       : ', ccelec*potSumGlob2/2.0
     print"(a,f12.4)",'GRMS (FMM)           : ', ccelec*sqrt(accNrmGlob/3.0/nglobal)
     print"(a,f12.4)",'GRMS (Ewald)         : ', ccelec*sqrt(accNrmGlob2/3.0/nglobal)
  end if

  do i = 1, nglobal
     p(i) = 0
     f(3*i-2) = 0
     f(3*i-1) = 0
     f(3*i-0) = 0
  end do
  call fmm_vanderwaals(nglobal, icpumap, atype, x, p, f, cuton, cutoff,&
       pcycle, nat, rscale, gscale, fgscale)
  do i = 1, nglobal
     p2(i) = 0
     f2(3*i-2) = 0
     f2(3*i-1) = 0
     f2(3*i-0) = 0
  end do
  call direct_vanderwaals(nglobal, icpumap, atype, x, p2, f2, cuton, cutoff,&
       pcycle, nat, rscale, gscale, fgscale)
  call vanderwaals_exclusion(nglobal, icpumap, atype, x, p, f, cuton, cutoff,&
       pcycle, nat, rscale, gscale, fgscale, numex, natex)
  call vanderwaals_exclusion(nglobal, icpumap, atype, x, p2, f2, cuton, cutoff,&
       pcycle, nat, rscale, gscale, fgscale, numex, natex)
  potSum = 0
  potSum2 = 0
  accDif = 0
  accNrm = 0
  accNrm2 = 0
  do i = 1, nglobal
     if (icpumap(i).eq.1) then
        potSum  = potSum  + p(i)
        potSum2 = potSum2 + p2(i)
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
     print"(a)",'--- VdW FMM vs. Direct ----------'
     print"(a,f9.7)",'Rel. L2 Error (pot)  : ', sqrt(potDifGlob/potNrmGlob2)
     print"(a,f9.7)",'Rel. L2 Error (acc)  : ', sqrt(accDifGlob/accNrmGlob2)
     print"(a,f12.4)",'Energy (FMM)         : ', potSumGlob/2.0
     print"(a,f12.4)",'Energy (Direct)      : ', potSumGlob2/2.0
     print"(a,f12.4)",'GRMS (FMM)           : ', sqrt(accNrmGlob/3.0/nglobal)
     print"(a,f12.4)",'GRMS (Direct)        : ', sqrt(accNrmGlob2/3.0/nglobal)
  end if
 
  ! run dynamics if second command line argument specified
  if (command_argument_count() > 1) then
     call get_command_argument(2,filename,lnam,istat)
     read(filename,*)dynsteps
     write(*,*)'will run dynamics for ',dynsteps,' steps'
     ! for pure water systems there is no need for nbadd14() :-)

     test_first=.false.
     integrate=.true.
     printfrq=1
     imcentfrq=10
     timstart=100.0 ! time of restart file (later: get it from there)

     if (test_first) call force_testing(nglobal,nat,nbonds,ntheta,ksize,&
          alpha,sigma,cutoff,cuton,ccelec,pcycle,&
          xc,p,f,q,v,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
          ib,jb,it,jt,kt,atype,icpumap,numex,natex,etot,eb,et,efmm,evdw)

     call print_energy(timstart,nglobal,nat,nbonds,ntheta,ksize,&
          alpha,sigma,cutoff,cuton,ccelec,pcycle,&
          xc,p,f,q,v,mass,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
          ib,jb,it,jt,kt,atype,icpumap,numex,natex,etot,eb,et,efmm,evdw)

     if (integrate) call run_dynamics(dynsteps,imcentfrq,printfrq,&
          nglobal,nat,nbonds,ntheta,ksize,&
          alpha,sigma,cutoff,cuton,ccelec,pcycle,&
          xc,p,f,q,v,mass,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
          ib,jb,it,jt,kt,atype,icpumap,numex,natex,etot,eb,et,efmm,evdw,nres,ires)

  endif

  deallocate( x, q, p, f, icpumap, p2, f2 )
  call mpi_finalize(ierr)
end program main
