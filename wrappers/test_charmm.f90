module charmm_io
contains

  subroutine charmm_cor_write(n,x,q,size,outfile,numex,natex,nat,atype,&
       rscale,gscale,fgscale,nbonds,ntheta,ib,jb,it,jt,kt,rbond,cbond,&
       aangle,cangle,mass,xc,v,time)
    implicit none
    character(len=*) outfile
    character(len=3) atom_name(3)
    integer i,nat,n,natex_size,vdw_size,nbonds,ntheta
    real(8) size,time
    integer,allocatable,dimension(:) :: numex,natex,atype,ib,jb,it,jt,kt
    real(8),allocatable,dimension(:) :: x,q,rscale,gscale,fgscale
    real(8),allocatable,dimension(:) :: mass,xc,v
    real(8),allocatable,dimension(:,:) :: rbond,cbond
    real(8),allocatable,dimension(:,:,:) :: aangle,cangle

    atom_name(1)='OH2'
    atom_name(2)='H1 '
    atom_name(3)='H2 '

    natex_size=sum(numex(1:n))
    vdw_size=nat**2

    open(unit=2,file=outfile,status='unknown')
    write(2,'(a8,3x,3f22.14)')'* SIZE =',size,size,size
    write(2,'(''* This is test_charmm restart file'')')
    write(2,'(''*'')')
    write(2,'(i10,''  EXT'')')n
    do i=1,n
       write(2,'(2i10,''  TIP3      '',a3,5x,3f20.10,''  WAT       '',i5,3x,f20.10)') &
            i,(i-1)/3+1,atom_name(mod(i-1,3)+1),x(i*3-2),x(i*3-1),x(i*3),(i-1)/3+1,q(i)
    enddo
    write(2,'(''NUMEX'')')
    write(2,'(10i10)')numex(1:n)
    write(2,'(''NATEX'',i6)')natex_size
    write(2,'(10i10)')natex(1:natex_size)
    write(2,'(''RSCALE==FRSCALE'',i3)')nat
    write(2,'(5f20.10)')rscale(1:vdw_size)
    write(2,'(''GSCALE'',i3)')nat
    write(2,'(5f20.10)')gscale(1:vdw_size)
    write(2,'(''FGSCALE'',i3)')nat
    write(2,'(5f20.10)')fgscale(1:vdw_size)
    write(2,'(''IACSOR'')')
    write(2,'(20i5)')atype(1:n)
    write(2,'(''BONDS'',i6)')nbonds
    write(2,'(10i10)')(ib(i),jb(i),i=1,nbonds)
    write(2,'(''ANGLES'',i5)')ntheta
    write(2,'(9i10)')(it(i),jt(i),kt(i),i=1,ntheta)
    write(2,'(''RBOND'',i3)')nat
    write(2,'(5f20.10)')rbond(1:nat,1:nat)
    write(2,'(''CBOND'',i3)')nat
    write(2,'(5f20.10)')cbond(1:nat,1:nat)
    write(2,'(''AANGLE'',i3)')nat
    write(2,'(5f20.10)')aangle(1:nat,1:nat,1:nat)
    write(2,'(''CANGLE'',i3)')nat
    write(2,'(5f20.10)')cangle(1:nat,1:nat,1:nat)
    write(2,'(''MASS'',i3)')nat
    write(2,'(5f20.10)')mass(1:nat)
    write(2,'(3d28.18)')(xc(3*i-2),xc(3*i-1),xc(3*i-0),i=1,n)
    write(2,'(3d28.18)')(v(3*i-2),v(3*i-1),v(3*i-0),i=1,n)
    write(2,'(f20.10)')time
    return
  end subroutine charmm_cor_write

  subroutine charmm_cor_read(n,x,q,size,infile,numex,natex,nat,atype,&
       rscale,gscale,fgscale,nbonds,ntheta,ib,jb,it,jt,kt,rbond,cbond,&
       aangle,cangle,mass,xc,v,nres,ires,time)
    implicit none
    logical qext
    character(len=100) lin
    character(len=*) infile
    integer i,nat,im,in,n,natex_size,vdw_size,resn,nres,nbonds,ntheta
    real(8) size,sizex,sizey,sizez,time
    integer,allocatable,dimension(:) :: numex,natex,atype,ib,jb,it,jt,kt,ires
    real(8),allocatable,dimension(:) :: x,q,rscale,gscale,fgscale
    real(8),allocatable,dimension(:) :: mass,xc,v
    real(8),allocatable,dimension(:,:) :: rbond,cbond
    real(8),allocatable,dimension(:,:,:) :: aangle,cangle

    open(unit=1,file=infile,status='old')
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
    allocate( x(3*n),q(n),v(3*n),ires(n) ) ! too much for ires
    ! read also the residue info from cor file
    ! residue info is enough for pure water system
    ! for protein + water the segments are needed, too
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
       v(3*i-2) = 0
       v(3*i-1) = 0
       v(3*i-0) = 0
       if(im /= resn)then
          resn=resn+1
          ires(resn)=i
       endif
    enddo
    nres=resn
    read(1,'(a100)')lin
    allocate(numex(n))
    read(1,'(10i10)')numex(1:n)
    read(1,'(a100)')lin
    read(lin(6:100),*)natex_size
    if(natex_size /= sum(numex(1:n))) then
       print*,'Something is wrong with the input file'
       print*,infile
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
    read(1,'(a100)',end=99)lin ! variables from here are only used only for dynamics
    if(lin(1:5) == 'BONDS')then
       read(lin(6:100),*)nbonds
       allocate(ib(nbonds),jb(nbonds))
       read(1,'(10i10)')(ib(i),jb(i),i=1,nbonds)
    endif
    read(1,'(a100)',end=99)lin
    if(lin(1:6) == 'ANGLES')then
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
    allocate(xc(3*n))
    read(1,'(3d28.18)')(xc(3*i-2),xc(3*i-1),xc(3*i-0),i=1,n)
    read(1,'(3d28.18)')(v(3*i-2),v(3*i-1),v(3*i-0),i=1,n)
    read(1,'(f20.10)')time
99  continue
    return
  end subroutine charmm_cor_read

  subroutine bcast1(nglobal,icpumap,recv)
    use mpi
    implicit none
    integer i,nglobal,ierr
    integer,allocatable,dimension(:) :: icpumap
    real(8),allocatable,dimension(:) :: send,recv
    allocate(send(nglobal))
    send(1:nglobal)=0.0
    do i=1,nglobal
       if(icpumap(i) == 1) send(i) = recv(i)
    enddo
    call mpi_allreduce(send,recv,nglobal,mpi_real8,mpi_sum,mpi_comm_world,ierr)
    deallocate(send)
  end subroutine bcast1

  subroutine bcast3(nglobal,icpumap,recv)
    use mpi
    implicit none
    integer i,nglobal,ierr
    integer,allocatable,dimension(:) :: icpumap
    real(8),allocatable,dimension(:) :: send,recv
    allocate(send(3*nglobal))
    send(1:3*nglobal)=0.0
    do i=1,nglobal
       if(icpumap(i) == 1) send(3*i-2:3*i) = recv(3*i-2:3*i)
    enddo
    call mpi_allreduce(send,recv,3*nglobal,mpi_real8,mpi_sum,mpi_comm_world,ierr)
    deallocate(send)
  end subroutine bcast3

  subroutine bonded_terms(icpumap,atype,x,f,nbonds,ntheta,&
       ib,jb,it,jt,kt,rbond,cbond,aangle,cangle,eb,et)
    implicit none
    integer i,ii,jj,kk,nbonds,ntheta
    real(8) dx,dy,dz,r,db,df,dfr,eb,et
    real(8) dx1,dy1,dz1,dx2,dy2,dz2,r1,r2,r1r,r2r
    real(8) dx1r,dy1r,dz1r,dx2r,dy2r,dz2r,cst,at,da
    real(8) dxi,dyi,dzi,st2r,str,dtx1,dtx2,dty1,dty2,dtz1,dtz2
    integer,allocatable,dimension(:) :: ib,jb,it,jt,kt,icpumap,atype
    real(8),allocatable,dimension(:) :: x,f
    real(8),allocatable,dimension(:,:) :: rbond,cbond
    real(8),allocatable,dimension(:,:,:) :: aangle,cangle
    eb = 0.0
    do i = 1,nbonds
       ii=ib(i)
       jj=jb(i)
       if (icpumap(ii) /= 1 .and. icpumap(jj) /= 1) cycle
       if (icpumap(ii) == 0) print*,"missing i = ",ii
       if (icpumap(jj) == 0) print*,"missing j = ",jj
       dx=x(3*ii-2)-x(3*jj-2)
       dy=x(3*ii-1)-x(3*jj-1)
       dz=x(3*ii-0)-x(3*jj-0)
       r=sqrt(dx*dx+dy*dy+dz*dz)
       db=r-rbond(atype(ii),atype(jj))
       df=cbond(atype(ii),atype(jj))*db
       dfr=2.0*df/r
       dxi=dx*dfr
       dyi=dy*dfr
       dzi=dz*dfr
       if (icpumap(ii) == 1) then
          eb=eb+df*db/2
          f(3*ii-2)=f(3*ii-2)+dxi
          f(3*ii-1)=f(3*ii-1)+dyi
          f(3*ii-0)=f(3*ii-0)+dzi
       endif
       if (icpumap(jj) == 1) then
          eb=eb+df*db/2
          f(3*jj-2)=f(3*jj-2)-dxi
          f(3*jj-1)=f(3*jj-1)-dyi
          f(3*jj-0)=f(3*jj-0)-dzi
       endif
    enddo
    et = 0.0
    do i = 1,ntheta
       ii=it(i)
       jj=jt(i)
       kk=kt(i)
       if (icpumap(ii) /= 1 .and. icpumap(jj) /= 1 .and. icpumap(kk) /= 1) cycle
       if (icpumap(ii) == 0) print*,"missing i = ",ii
       if (icpumap(jj) == 0) print*,"missing j = ",jj
       if (icpumap(kk) == 0) print*,"missing k = ",kk
       dx1=x(3*ii-2)-x(3*jj-2)
       dy1=x(3*ii-1)-x(3*jj-1)
       dz1=x(3*ii-0)-x(3*jj-0)
       dx2=x(3*kk-2)-x(3*jj-2)
       dy2=x(3*kk-1)-x(3*jj-1)
       dz2=x(3*kk-0)-x(3*jj-0)
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
       st2r=1.0/(1.0-cst*cst)
       str=sqrt(st2r)
       dfr=-2.0*df*str
       dtx1=dfr*r1r*(dx2r-cst*dx1r)
       dtx2=dfr*r2r*(dx1r-cst*dx2r)
       dty1=dfr*r1r*(dy2r-cst*dy1r)
       dty2=dfr*r2r*(dy1r-cst*dy2r)
       dtz1=dfr*r1r*(dz2r-cst*dz1r)
       dtz2=dfr*r2r*(dz1r-cst*dz2r)
       if (icpumap(ii) == 1) then
          et=et+df*da/3
          f(3*ii-2)=f(3*ii-2)+dtx1
          f(3*ii-1)=f(3*ii-1)+dty1
          f(3*ii-0)=f(3*ii-0)+dtz1
       endif
       if (icpumap(jj) == 1) then
          et=et+df*da/3
          f(3*jj-2)=f(3*jj-2)-dtx1-dtx2
          f(3*jj-1)=f(3*jj-1)-dty1-dty2
          f(3*jj-0)=f(3*jj-0)-dtz1-dtz2
       endif
       if (icpumap(kk) == 1) then
          et=et+df*da/3
          f(3*kk-2)=f(3*kk-2)+dtx2
          f(3*kk-1)=f(3*kk-1)+dty2
          f(3*kk-0)=f(3*kk-0)+dtz2
       endif
    enddo
  end subroutine bonded_terms

  subroutine verify(nglobal,icpumap,p,p2,f,f2,&
    pl2err,fl2err,enerf,enere,grmsf,grmse)
    use mpi
    implicit none
    integer nglobal,i,ierr
    real(8) potSum,potSum2,potDifGlob,potSumGlob,potSumGlob2,potNrmGlob2
    real(8) accDif,accNrm,accNrm2,accDifGlob,accNrmGlob,accNrmGlob2
    real(8) pl2err,fl2err,enerf,enere,grmsf,grmse
    integer,allocatable,dimension(:) :: icpumap
    real(8),allocatable,dimension(:) :: p,p2,f,f2
    potSum = 0
    potSum2 = 0
    accDif = 0
    accNrm = 0
    accNrm2 = 0
    do i = 1,nglobal
       if (icpumap(i) == 1) then
          p(i) = p(i)
          f(3*i-2) = f(3*i-2)
          f(3*i-1) = f(3*i-1)
          f(3*i-0) = f(3*i-0)
          p2(i) = p2(i)
          f2(3*i-2) = f2(3*i-2)
          f2(3*i-1) = f2(3*i-1)
          f2(3*i-0) = f2(3*i-0)
          potSum = potSum + p(i)
          potSum2 = potSum2 + p2(i)
          accDif = accDif  + (f(3*i-2) - f2(3*i-2)) * (f(3*i-2) - f2(3*i-2))&
               + (f(3*i-1) - f2(3*i-1)) * (f(3*i-1) - f2(3*i-1))&
               + (f(3*i-0) - f2(3*i-0)) * (f(3*i-0) - f2(3*i-0))
          accNrm = accNrm + f(3*i-2) * f(3*i-2) + f(3*i-1) * f(3*i-1) + f(3*i-0) * f(3*i-0)
          accNrm2 = accNrm2 + f2(3*i-2) * f2(3*i-2) + f2(3*i-1) * f2(3*i-1) + f2(3*i-0) * f2(3*i-0)
       endif
    enddo
    potSumGlob = 0
    potSumGlob2 = 0
    accDifGlob = 0
    accNrmGlob = 0
    accNrmGlob2 = 0
    call mpi_reduce(potSum,potSumGlob,1,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(potSum2,potSumGlob2,1,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(accDif,accDifGlob,1,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(accNrm,accNrmGlob,1,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(accNrm2,accNrmGlob2,1,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    potDifGlob = (potSumGlob - potSumGlob2) * (potSumGlob - potSumGlob2)
    potNrmGlob2 = potSumGlob2 * potSumGlob2
    pl2err = sqrt(potDifGlob/potNrmGlob2)
    fl2err = sqrt(accDifGlob/accNrmGlob2)
    enerf = potSumGlob*0.5
    enere = potSumGlob2*0.5
    grmsf = sqrt(accNrmGlob/3.0/nglobal)
    grmse = sqrt(accNrmGlob2/3.0/nglobal)
  end subroutine verify

  subroutine energy(nglobal,nat,nbonds,ntheta,ksize,&
       alpha,sigma,cutoff,cuton,pcycle,xold,&
       x,p,p2,f,f2,q,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
       ib,jb,it,jt,kt,atype,icpumap,numex,natex,etot,&
       pl2err,fl2err,ftotf,ftote,istep)
    implicit none
    integer nglobal,nat,nbonds,ntheta,ksize,i
    real(8) alpha,sigma,cutoff,cuton,etot,eb,et,efmm,evdw,pcycle
    real(8) pl2err,fl2err,enerf,enere,grmsf,grmse,ftotf,ftote
    integer,optional :: istep
    integer,allocatable,dimension(:) :: ib,jb,it,jt,kt,atype,icpumap,numex,natex
    real(8),allocatable,dimension(:) :: xold,x,p,p2,f,f2,q,gscale,fgscale,rscale
    real(8),allocatable,dimension(:,:) :: rbond,cbond
    real(8),allocatable,dimension(:,:,:) :: aangle,cangle

    call fmm_partition(nglobal,icpumap,x,q,xold,pcycle)
    p(1:nglobal)=0.0
    p2(1:nglobal)=0.0
    f(1:3*nglobal)=0.0
    f2(1:3*nglobal)=0.0
    call fmm_coulomb(nglobal,icpumap,x,q,p,f,pcycle)
    call ewald_coulomb(nglobal,icpumap,x,q,p2,f2,ksize,alpha,sigma,cutoff,pcycle)
    call verify(nglobal,icpumap,p,p2,f,f2,pl2err,fl2err,enerf,enere,grmsf,grmse)
    ftotf=0.0
    ftote=0.0
    do i = 1,nglobal
       if (icpumap(i) /= 1) cycle
       ftotf = ftotf+ f(3*i-2)+ f(3*i-1)+ f(3*i)
       ftote = ftote+f2(3*i-2)+f2(3*i-1)+f2(3*i)
    enddo
    call coulomb_exclusion(nglobal,icpumap,x,q,p,f,pcycle,numex,natex)
    efmm=0.0
    do i=1,nglobal
       if (icpumap(i) /= 1) cycle
       efmm=efmm+p(i)
    enddo
    efmm=efmm*0.5
    p(1:nglobal)=0.0
    call fmm_vanderwaals(nglobal,icpumap,atype,x,p,f,cuton,cutoff,&
         pcycle,nat,rscale,gscale,fgscale)
    call vanderwaals_exclusion(nglobal,icpumap,atype,x,p,f,cuton,cutoff,&
         pcycle,nat,rscale,gscale,fgscale,numex,natex)
    evdw=0.0
    do i = 1,nglobal
       if (icpumap(i) /= 1) cycle
       evdw=evdw+p(i)
    enddo
    evdw=evdw*0.5
    call bonded_terms(icpumap,atype,x,f,nbonds,ntheta,&
         ib,jb,it,jt,kt,rbond,cbond,aangle,cangle,eb,et)
    etot=eb+et+efmm+evdw
    return
  end subroutine energy

  subroutine force_testing(nglobal,nat,nbonds,ntheta,ksize,&
       alpha,sigma,cutoff,cuton,pcycle,xold,&
       x,p,p2,f,f2,q,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
       ib,jb,it,jt,kt,atype,icpumap,numex,natex)

    implicit none
    integer nglobal,nat,nbonds,ntheta,ksize,i,j,istart,iend
    real(8) alpha,sigma,cutoff,cuton,etot,pcycle
    real(8) xsave,e0,step,step2,eplus,eminus,nforce
    real(8) pl2err,fl2err,ftotf,ftote
    integer,allocatable,dimension(:) :: ib,jb,it,jt,kt,atype,icpumap,numex,natex
    real(8),allocatable,dimension(:) :: xold,x,p,p2,f,f2,q,gscale,fgscale,rscale
    real(8),allocatable,dimension(:,:) :: rbond,cbond
    real(8),allocatable,dimension(:,:,:) :: aangle,cangle
    real(8),allocatable,dimension(:) :: f_analytic
    call energy(nglobal,nat,nbonds,ntheta,ksize,&
         alpha,sigma,cutoff,cuton,pcycle,xold,&
         x,p,p2,f,f2,q,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
         ib,jb,it,jt,kt,atype,icpumap,numex,natex,etot,&
         pl2err,fl2err,ftotf,ftote)
    e0=etot
    allocate(f_analytic(3*nglobal))
    f_analytic(1:3*nglobal)=f(1:3*nglobal)
    istart=1
    iend=10
    step=0.0001
    step2=1.0/(2.0*step)
    write(*,*)'   I      J             Numeric           Analytic             Difference'
    do i = istart,iend
       do j = 1,3
          xsave=x(3*(i-1)+j) ! save it
          x(3*(i-1)+j)=x(3*(i-1)+j)-step
          call energy(nglobal,nat,nbonds,ntheta,ksize,&
               alpha,sigma,cutoff,cuton,pcycle,xold,&
               x,p,p2,f,f2,q,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
               ib,jb,it,jt,kt,atype,icpumap,numex,natex,eminus,&
               pl2err,fl2err,ftotf,ftote)
          x(3*(i-1)+j)=x(3*(i-1)+j)+2.0*step
          call energy(nglobal,nat,nbonds,ntheta,ksize,&
               alpha,sigma,cutoff,cuton,pcycle,xold,&
               x,p,p2,f,f2,q,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
               ib,jb,it,jt,kt,atype,icpumap,numex,natex,eplus,&
               pl2err,fl2err,ftotf,ftote)
          nforce=step2*(eplus-eminus)
          x(3*(i-1)+j)=xsave  ! restore
          write(*,'(2i7,3f20.9)')i,j,nforce,f_analytic(3*(i-1)+j),nforce-f_analytic(3*(i-1)+j)
       enddo
    enddo
    deallocate(f_analytic)
    return
  end subroutine force_testing

  subroutine print_energy(time,nglobal,f,v,mass,atype,icpumap,etot,pl2err,fl2err,ftotf,ftote)
    use mpi
    implicit none
    integer nglobal,i,ierr,mpirank
    real(8) time,temp,kboltz,ekinetic,ekineticGlob,grms,grmsGlob
    real(8) etot,etotGlob,ftotf,ftotfGlob,ftote,ftoteGlob
    real(8) pl2err,fl2err
    integer,allocatable,dimension(:) :: atype,icpumap
    real(8),allocatable,dimension(:) :: f,v,mass
    kboltz=1.987191d-03 !from CHARMM
    ekinetic=0.0
    grms=0.0
    do i=1,nglobal
       if(icpumap(i) /= 1) cycle
       ekinetic=ekinetic+mass(atype(i))*(v(3*i-2)**2+v(3*i-1)**2+v(3*i-0)**2)
       grms=grms+f(3*i-2)**2+f(3*i-1)**2+f(3*i-0)**2
    enddo
    call mpi_reduce(ekinetic,ekineticGlob,1,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(grms,grmsGlob,1,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(etot,etotGlob,1,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(ftotf,ftotfGlob,1,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(ftote,ftoteGlob,1,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    grmsGlob=sqrt(grmsGlob/3.0/real(nglobal))
    temp=ekineticGlob/3.0/nglobal/kboltz
    ekineticGlob=ekineticGlob/2.0
    ftotfGlob = ftotfGlob/nglobal
    ftoteGlob = ftoteGlob/nglobal
    call mpi_comm_rank(mpi_comm_world,mpirank,ierr)
    if(mpirank == 0) then
       write(*,'(''time:'',f9.3,'' Etotal:'',g14.5,'' Ekin:'',g14.5,'' Epot:'',g14.5,'' T:'',g12.3,'' Grms:'',g12.5)')&
            time,etotGlob+ekineticGlob,ekineticGlob,etotGlob,temp,grmsGlob
       write(2,'(f9.3,f14.5,f14.5,g14.5,g14.5,g14.5,g14.5)')&
            time,etotGlob,ekineticGlob,pl2err,fl2err,ftotfGlob,ftoteGlob
    endif
    return
  end subroutine print_energy

  subroutine run_dynamics(dynsteps,imcentfrq,printfrq,&
       nglobal,nat,nbonds,ntheta,ksize,&
       alpha,sigma,cutoff,cuton,pcycle,&
       x,p,p2,f,f2,q,v,mass,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
       ib,jb,it,jt,kt,atype,icpumap,numex,natex,nres,ires,time)
    use mpi
    implicit none
    integer dynsteps,nglobal,nat,nbonds,ntheta,ksize,imcentfrq,printfrq,nres
    integer i,istep,ierr,mpirank
    real(8) alpha,sigma,cutoff,cuton,etot,pcycle,tstep,tstep2,time
    real(8) pl2err,fl2err,ftotf,ftote
    real(8),parameter :: timfac=4.88882129D-02
    real(8),allocatable,dimension(:) :: x,v,mass,p,p2,f,f2,q,gscale,fgscale,rscale,xold
    real(8),allocatable,dimension(:,:) :: rbond,cbond
    real(8),allocatable,dimension(:,:,:) :: aangle,cangle
    integer,allocatable,dimension(:) :: ib,jb,it,jt,kt,atype,icpumap,numex,natex,ires
    real(8),allocatable,dimension(:) :: xnew,fac1,fac2

    call mpi_comm_rank(mpi_comm_world,mpirank,ierr)
    if(mpirank == 0)then
       open(unit=1,file='water.pdb',status='unknown')
       open(unit=2,file='verify.dat',status='unknown')
    endif

    tstep = 0.001/timfac !ps -> akma
    tstep2 = tstep**2

    allocate(xnew(3*nglobal),xold(3*nglobal),fac1(nglobal),fac2(nglobal))

    call energy(nglobal,nat,nbonds,ntheta,ksize,&
         alpha,sigma,cutoff,cuton,pcycle,xold,&
         x,p,p2,f,f2,q,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
         ib,jb,it,jt,kt,atype,icpumap,numex,natex,etot,&
         pl2err,fl2err,ftotf,ftote)

    call print_energy(time,nglobal,f,v,mass,atype,icpumap,etot,pl2err,fl2err,ftotf,ftote)

    ! precompute some constants and recalculate xold
    do i=1,nglobal
       fac1(i) = tstep2/mass(atype(i))
       fac2(i) = 0.5/tstep
       if(icpumap(i) /= 1)cycle
       xold(3*i-2) = v(3*i-2)*tstep-f(3*i-2)*fac1(i)*0.5
       xold(3*i-1) = v(3*i-1)*tstep-f(3*i-1)*fac1(i)*0.5
       xold(3*i-0) = v(3*i-0)*tstep-f(3*i-0)*fac1(i)*0.5
    enddo

    mainloop: do istep = 1,dynsteps

       do i = 1,nglobal
          if(icpumap(i) /= 1)cycle
          x(3*i-2) = x(3*i-2) + xold(3*i-2)
          x(3*i-1) = x(3*i-1) + xold(3*i-1)
          x(3*i-0) = x(3*i-0) + xold(3*i-0)
       enddo

       call energy(nglobal,nat,nbonds,ntheta,ksize,&
            alpha,sigma,cutoff,cuton,pcycle,xold,&
            x,p,p2,f,f2,q,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
            ib,jb,it,jt,kt,atype,icpumap,numex,natex,etot,&
            pl2err,fl2err,ftotf,ftote,istep)

       do i = 1,nglobal
          if(icpumap(i) /= 1)cycle
          xnew(3*i-2) = xold(3*i-2) - fac1(i)*f(3*i-2)
          xnew(3*i-1) = xold(3*i-1) - fac1(i)*f(3*i-1)
          xnew(3*i-0) = xold(3*i-0) - fac1(i)*f(3*i-0)
       enddo

       do i = 1,nglobal
          if(icpumap(i) /= 1)cycle
          v(3*i-2) = (xnew(3*i-2) + xold(3*i-2))*fac2(i)
          v(3*i-1) = (xnew(3*i-1) + xold(3*i-1))*fac2(i)
          v(3*i-0) = (xnew(3*i-0) + xold(3*i-0))*fac2(i)
       enddo

       do i = 1,nglobal
          if(icpumap(i) /= 1)cycle
          xold(3*i-2) = xnew(3*i-2)
          xold(3*i-1) = xnew(3*i-1)
          xold(3*i-0) = xnew(3*i-0)
       enddo

       if (mod(istep,imcentfrq) == 0) call image_center(nglobal,x,nres,ires,pcycle,icpumap)

       time=time+tstep*timfac ! for printing only
       if (mod(istep,printfrq) == 0) then
          call print_energy(time,nglobal,f,v,mass,atype,icpumap,etot,&
               pl2err,fl2err,ftotf,ftote)
          !call pdb_frame(1,time,nglobal,x,nres,icpumap)
       endif
    enddo mainloop

  end subroutine run_dynamics

  subroutine pdb_frame(unit,time,nglobal,x,nres,icpumap)
    use mpi
    implicit none
    integer nglobal,nres,ipt,i,unit,ierr,mpirank
    real(8) time
    integer,allocatable,dimension(:) :: icpumap
    real(8),allocatable,dimension(:) :: x
    call mpi_comm_rank(mpi_comm_world,mpirank,ierr)
    call bcast3(nglobal,icpumap,x)
    if (mpirank /= 0) return
    write(unit,'(''REMARK frame at time:'',f14.4,'' ps'')')time
    ipt=0
    do i = 1,nres
       ipt=ipt+1
       write(unit,'(a6,i5,2x,a3,1x,a3,1x,a1,i4,3x,3f8.3,2f6.2,11x,a1)')&
            'HETATM',ipt,'OH2','HOH','w',i,x(3*ipt-2),x(3*ipt-1),x(3*ipt),1.0,1.0,'O'
       ipt=ipt+1
       write(unit,'(a6,i5,2x,a3,1x,a3,1x,a1,i4,3x,3f8.3,2f6.2,11x,a1)')&
            'HETATM',ipt,'HO1','HOH','w',i,x(3*ipt-2),x(3*ipt-1),x(3*ipt),1.0,1.0,'H'
       ipt=ipt+1
       write(unit,'(a6,i5,2x,a3,1x,a3,1x,a1,i4,3x,3f8.3,2f6.2,11x,a1)')&
            'HETATM',ipt,'HO2','HOH','w',i,x(3*ipt-2),x(3*ipt-1),x(3*ipt),1.0,1.0,'H'
    enddo
    write(unit,'(''END'')')
    return

  end subroutine pdb_frame

  subroutine image_center(nglobal,x,nres,ires,pcycle,icpumap)
    implicit none
    integer nglobal,nres,flag,istart,iend,nsel,i,j
    real(8) pcycle,xmin,xmax,ymin,ymax,zmin,zmax,xcen,ycen,zcen
    integer,allocatable,dimension(:) :: ires,icpumap
    real(8),allocatable,dimension(:) :: x

    ! only simple translations, so velocities stay the same! "?"

    xmin=-pcycle/2.0
    ymin=-pcycle/2.0
    zmin=-pcycle/2.0
    xmax=pcycle/2.0
    ymax=pcycle/2.0
    zmax=pcycle/2.0

    residues: do i = 1,nres
       istart=ires(i)
       if(i == nres) then
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
       do j = istart,iend
          if(icpumap(j) == 1)flag=1
       enddo
       if (flag == 0) cycle residues
!!!
       xcen=0.0
       ycen=0.0
       zcen=0.0
       nsel=iend-istart+1
       do j = istart,iend
          xcen=xcen+x(3*j-2)
          ycen=ycen+x(3*j-1)
          zcen=zcen+x(3*j-0)
       enddo
       xcen=xcen/real(nsel)
       ycen=ycen/real(nsel)
       zcen=zcen/real(nsel)
       if(xcen < xmin)then
          do j = istart,iend
             x(3*j-2)=x(3*j-2)+pcycle
          enddo
       endif
       if(xcen > xmax)then
          do j = istart,iend
             x(3*j-2)=x(3*j-2)-pcycle
          enddo
       endif
       if(ycen < ymin)then
          do j = istart,iend
             x(3*j-1)=x(3*j-1)+pcycle
          enddo
       endif
       if(ycen > ymax)then
          do j = istart,iend
             x(3*j-1)=x(3*j-1)-pcycle
          enddo
       endif
       if(zcen < zmin)then
          do j = istart,iend
             x(3*j-0)=x(3*j-0)+pcycle
          enddo
       endif
       if(zcen > zmax)then
          do j = istart,iend
             x(3*j-0)=x(3*j-0)-pcycle
          enddo
       endif
    enddo residues

  end subroutine image_center

end module charmm_io

subroutine split_range(ista,iend,isplit,nsplit)
  implicit none
  integer ista,iend,isplit,nsplit
  integer increment,remainder,size
  size = iend - ista + 1
  increment = size / nsplit
  remainder = mod(size,nsplit)
  ista = ista + isplit * increment + min(isplit,remainder)
  iend = ista + increment - 1
  if(remainder > isplit) iend = iend + 1
  return
end subroutine split_range

program main
  use charmm_io
  implicit none
  include 'mpif.h'
  logical test_force
  character(len=128) infile,outfile,nstp
  integer dynsteps
  integer i,ierr,images,ista,iend,istat,ksize,lnam,mpirank,mpisize
  integer nat,nglobal,verbose,nbonds,ntheta,imcentfrq,printfrq,nres
  real(8) alpha,sigma,cuton,cutoff,average,pcycle,theta,time
  real(8) pl2err,fl2err,enerf,enere,grmsf,grmse
  integer,dimension (128) :: iseed
  integer,allocatable,dimension(:) :: icpumap,numex,natex,atype,ib,jb,it,jt,kt,ires
  real(8),parameter :: pi=3.14159265358979312d0
  real(8),allocatable,dimension(:) :: x,q,p,p2,f,f2,xc,v,mass
  real(8),allocatable,dimension(:) :: rscale,gscale,fgscale
  real(8),allocatable,dimension(:,:) :: rbond,cbond
  real(8),allocatable,dimension(:,:,:) :: aangle,cangle

  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world,mpisize,ierr)
  call mpi_comm_rank(mpi_comm_world,mpirank,ierr)
  nglobal = 1000
  images = 3
  theta = 0.3
  verbose = 0
  ksize = 11
  pcycle = 10 * pi
  sigma = .25 / pi
  cuton = 9.5
  cutoff = 10.0
  alpha = 10 / pcycle
  nat = 16
  time = 100. ! first 100ps was equilibration with standard CHARMM
  charmmio: if (command_argument_count() > 1) then
     call get_command_argument(1,infile,lnam,istat)
     call get_command_argument(2,outfile,lnam,istat)
     call charmm_cor_read(nglobal,x,q,pcycle,infile,numex,natex,nat,atype,&
          rscale,gscale,fgscale,nbonds,ntheta,ib,jb,it,jt,kt,rbond,cbond,&
          aangle,cangle,mass,xc,v,nres,ires,time)
     allocate( p(nglobal),f(3*nglobal),icpumap(nglobal) )
     allocate( p2(nglobal),f2(3*nglobal) )
     alpha = 10 / pcycle

  else
     allocate( x(3*nglobal),q(nglobal),v(3*nglobal) )
     allocate( p(nglobal),p2(nglobal),f(3*nglobal),f2(3*nglobal) )
     allocate( xc(3*nglobal),ires(nglobal),icpumap(nglobal) )
     allocate( numex(nglobal),natex(nglobal),atype(nglobal) )
     allocate( rscale(nat*nat),gscale(nat*nat),fgscale(nat*nat) )
     do i = 1,128
        iseed(i) = 0
     enddo
     call random_seed(put=iseed)
     call random_number(x)
     call random_number(q)
     average = 0
     do i = 1,nglobal
        x(3*i-2) = x(3*i-2) * pcycle - pcycle / 2
        x(3*i-1) = x(3*i-1) * pcycle - pcycle / 2
        x(3*i-0) = x(3*i-0) * pcycle - pcycle / 2
        p(i) = 0
        p2(i) = 0
        f(3*i-2) = 0
        f(3*i-1) = 0
        f(3*i-0) = 0
        f2(3*i-2) = 0
        f2(3*i-1) = 0
        f2(3*i-0) = 0
        icpumap(i) = 0
        average = average + q(i)
     enddo
     average = average / nglobal
     do i = 1,nglobal
        q(i) = q(i) - average
     enddo
     do i = 1,nglobal
        numex(i) = 1
        if(mod(i,2) == 1)then
           natex(i) = i+1
        else
           natex(i) = i-1
        endif
        atype(i) = 1
     enddo
     do i = 1,nat*nat
        rscale(i) = 1
        gscale(i) = 0.0001
        fgscale(i) = gscale(i)
     enddo
  endif charmmio
  if (mpirank == 0) print*,'I/O done'
  ista = 1
  iend = nglobal / 3
  call split_range(ista,iend,mpirank,mpisize)
  ista = (ista - 1) * 3 + 1
  iend = iend * 3
  do i = ista,iend
     icpumap(i) = 1
  enddo
  if (mpirank == 0) print*,'FMM init'
  call fmm_init(images,theta,verbose,nglobal)
  if (mpirank == 0) print*,'FMM partition'
  call fmm_partition(nglobal,icpumap,x,q,v,pcycle)
  if (mpirank == 0) print*,'FMM Coulomb'
  call fmm_coulomb(nglobal,icpumap,x,q,p,f,pcycle)
  do i = 1,nglobal
     p2(i) = 0
     f2(3*i-2) = 0
     f2(3*i-1) = 0
     f2(3*i-0) = 0
  enddo
  cutoff = 20
  if (mpirank == 0) print*,'Ewald Coulomb'
  call ewald_coulomb(nglobal,icpumap,x,q,p2,f2,ksize,alpha,sigma,cutoff,pcycle)
!  call direct_coulomb(nglobal,icpumap,x,q,p2,f2,pcycle)
  if(mpirank == 0) print*,'Coulomb exclusion'
  call coulomb_exclusion(nglobal,icpumap,x,q,p,f,pcycle,numex,natex)
  call coulomb_exclusion(nglobal,icpumap,x,q,p2,f2,pcycle,numex,natex)


  call verify(nglobal,icpumap,p,p2,f,f2,pl2err,fl2err,enerf,enere,grmsf,grmse)
  if (mpirank == 0) then
     print"(a)",'--- Coulomb FMM vs. Ewald -------'
     print"(a,f9.6)",'Rel. L2 Error (pot)  : ',pl2err
     print"(a,f9.6)",'Rel. L2 Error (acc)  : ',fl2err
     print"(a,f12.4)",'Energy (FMM)         : ',enerf
     print"(a,f12.4)",'Energy (Ewald)       : ',enere
     print"(a,f12.4)",'GRMS (FMM)           : ',grmsf
     print"(a,f12.4)",'GRMS (Ewald)         : ',grmse
  endif

  do i = 1,nglobal
     p(i) = 0
     p2(i) = 0
     f(3*i-2) = 0
     f(3*i-1) = 0
     f(3*i-0) = 0
     f2(3*i-2) = 0
     f2(3*i-1) = 0
     f2(3*i-0) = 0
  enddo
  call fmm_vanderwaals(nglobal,icpumap,atype,x,p,f,cuton,cutoff,&
       pcycle,nat,rscale,gscale,fgscale)
  call vanderwaals_exclusion(nglobal,icpumap,atype,x,p,f,cuton,cutoff,&
       pcycle,nat,rscale,gscale,fgscale,numex,natex)
  call direct_vanderwaals(nglobal,icpumap,atype,x,p2,f2,cuton,cutoff,&
       pcycle,nat,rscale,gscale,fgscale)
  call vanderwaals_exclusion(nglobal,icpumap,atype,x,p2,f2,cuton,cutoff,&
       pcycle,nat,rscale,gscale,fgscale,numex,natex)

  call verify(nglobal,icpumap,p,p2,f,f2,pl2err,fl2err,enerf,enere,grmsf,grmse)
  if (mpirank == 0) then
     print"(a)",'--- VdW FMM vs. Direct ----------'
     print"(a,f9.6)",'Rel. L2 Error (pot)  : ',pl2err
     print"(a,f9.6)",'Rel. L2 Error (acc)  : ',fl2err
     print"(a,f12.4)",'Energy (FMM)         : ',enerf
     print"(a,f12.4)",'Energy (Direct)      : ',enere
     print"(a,f12.4)",'GRMS (FMM)           : ',grmsf
     print"(a,f12.4)",'GRMS (Direct)        : ',grmse
  endif

  ! run dynamics if second command line argument specified
  if (command_argument_count() > 2) then
     call get_command_argument(3,nstp,lnam,istat)
     read(nstp,*)dynsteps
     if(mpirank == 0) write(*,*)'will run dynamics for ',dynsteps,' steps'
     ! for pure water systems there is no need for nbadd14() :-)
     test_force=.false.
     printfrq=1
     imcentfrq=10
     if (test_force) then
        call force_testing(nglobal,nat,nbonds,ntheta,ksize,&
             alpha,sigma,cutoff,cuton,pcycle,v,&
             xc,p,p2,f,f2,q,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
             ib,jb,it,jt,kt,atype,icpumap,numex,natex)
     else
        call run_dynamics(dynsteps,imcentfrq,printfrq,&
             nglobal,nat,nbonds,ntheta,ksize,&
             alpha,sigma,cutoff,cuton,pcycle,&
             xc,p,p2,f,f2,q,v,mass,gscale,fgscale,rscale,rbond,cbond,aangle,cangle,&
             ib,jb,it,jt,kt,atype,icpumap,numex,natex,nres,ires,time)
     endif
     call charmm_cor_write(nglobal,x,q,pcycle,trim(outfile),numex,natex,nat,atype,&
          rscale,gscale,fgscale,nbonds,ntheta,ib,jb,it,jt,kt,rbond,cbond,&
          aangle,cangle,mass,xc,v,time)
  endif

  deallocate( x,q,v,p,f,p2,f2,icpumap )
  deallocate( ires,numex,natex,rscale,gscale,fgscale,atype )

! from the end of run_dynamics():
!    deallocate(xnew,xold,fac1,fac2)
!    deallocate(ib,jb,it,jt,kt,rbond,cbond,mass,aangle,cangle,x)

  call fmm_finalize()
  call mpi_finalize(ierr)
end program main
