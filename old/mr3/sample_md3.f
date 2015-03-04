c********************************************************************
c C pre-processor macros are inserted by Testu Narumi (Jul.23, 1999)
c
c define the following line if you use Dr. yasuoka's program'
c for MDone
c#define USE_YASU_PROGRAM
c
c define the following line if you use WINE-2
c#define USE_WINE2
c
c********************************************************************
c
c*********************************************************
c* simulation.fort77 : molecular dynamics simulation *
c* Sample program for MD-GRAPE I *
c* coded by Kenji Yasuoka (C.S.L. RIKEN) *
c* (yasu@atlas.riken.go.jp) *
c* 2.0.2 ed. March 18, 1998 *
c (fixed at subroutine MR1calccoulomb_yasu)*
c (fixed at subroutine calc_force_eng_vdw) *
c (fixed at non-periodic condision) *
c*********************************************************
      program main
c %%%%%%%%%%%%%%%%%
c %%% main %%%%%%%%<< calls subroutines & controls data flow >>%%
c %%% 97/12/15 %%%
c %%%%%%%%%%%%%%%%%
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
c parameter (nmlcl=8)
      parameter (nmlcl=125)
c parameter (nmlcl=2197)
c      parameter (nmlcl=8000)
c parameter (nmlcl=15625)
c parameter (natot=1)
      parameter (natot=4)
      parameter (kmax=10)
c parameter (kmax=6)
c parameter (kmax=2)
c
      integer*4 natype(nmlcl)
      real*8 gseng(natot,natot),gsforc(natot,natot)
      real*8 rseng(natot,natot),rsforc(natot,natot)
      real*8 smass(natot)
      real*8 pos(3,nmlcl),vmom(3,nmlcl)
      real*8 engph(3,nmlcl),engpg(3,nmlcl)
      real*8 e_vdwg(3,nmlcl),e_vdwh(3,nmlcl)
      real*8 e_c_dg(3,nmlcl),e_c_dh(3,nmlcl)
      real*8 e_c_rg(3,nmlcl),e_c_rh(3,nmlcl)
      real*8 e_c_wg(3,nmlcl),e_c_wh(3,nmlcl),e_c_wg2(3,nmlcl)
      real*8 forceh(3,nmlcl),forceg(3,nmlcl)
      real*8 f_vdwg(3,nmlcl),f_vdwh(3,nmlcl)
      real*8 f_c_dg(3,nmlcl),f_c_dh(3,nmlcl)
      real*8 f_c_rg(3,nmlcl),f_c_rh(3,nmlcl)
      real*8 f_c_wg(3,nmlcl),f_c_wh(3,nmlcl)
      real*8 engt(nmlcl)
      real*8 charge(nmlcl)
      integer*4 kvect(3,4*kmax*kmax*kmax)
      real*8 rkvect(3,4*kmax*kmax*kmax)
      real*8 fac_k(4*kmax*kmax*kmax)
      real*8 cell(3,3),stress(3,3),cellsize(3),rkmax
      real*4 time0,time1,time2,time3,time4
      character*30 cuprnt
c
c data ncrpm, ncoum /100, 1/
        data ncrpm, ncoum /5, 1/
c        data ncrpm, ncoum /1, 1/
        data xvolum / 32.5d0/
c data xvolum / 10.0d0/
        data temp / 0.67d0 /
        data delt / 0.0050d0/
c
        data nuprnt / 11/
        data cuprnt /'sample_mdone.prnt'/
        data pi/3.1415926535897932385d0/
c
        eps=0.25d0/pi
c eps=0.1d0
        do i=1,3
           cell(i,i)=xvolum
           cellsize(i)=xvolum
        enddo
        rkmax=1.0d0*kmax
c
                open(nuprnt,file=cuprnt,status='unknown')
c
c*****
c system define
c
c Lennard-Jones Force style
c nsty_f = 0 : x**(-7)-x**(-4)
c nsty_f = 1 : 2*x**(-7)-x**(-4)
        nsty_f=1
        call sysdef(nmlcl,natot,natype,gseng,gsforc,rseng,rsforc,
     - smass,charge,alpha,nsty_f)
c call write_sys(nmlcl,natype,natot,gseng,gsforc,
c - rseng,rsforc,smass)
c*****
c set wave number and constant part of Ewald Sum in Wave number space
c
        call set_kvect(kmax,xvolum,eps,alpha,knum,kvect,fac_k)
c*
c change integer to real number for wave number vector
        call set_rkvect(knum,kvect,rkvect)
c
c*****
c make initial position and velocity for MD
c
        call start(pos,vmom,nmlcl,xvolum)
c*
c set temperature of system
        call tmpset(vmom,nmlcl,temp,smass,natype,natot)
c call write_start(pos,vmom,nmlcl)
c
c initialize MDone or MDGRAPE-2


        call MR3init
c        call MR1init
c call MR1setwaitsec(5)
c call MR1setwaitusec(500000)


        do 1100 imrpt=1,ncrpm,1
c*****
c initialize force and energy
c
        call values(nmlcl,e_vdwg,f_vdwg)
        call values(nmlcl,e_vdwh,f_vdwh)
        call values(nmlcl,e_c_dg,f_c_dg)
        call values(nmlcl,e_c_dh,f_c_dh)
        call values(nmlcl,e_c_rg,f_c_rg)
        call values(nmlcl,e_c_rh,f_c_rh)
        call values(nmlcl,e_c_wg,f_c_wg)
        call values(nmlcl,e_c_wg2,f_c_wg)
        call values(nmlcl,e_c_wh,f_c_wh)
                call cputime(time0)
c
c        goto 1098
c
c********
c calc Lennard-Jones force
c table number for calculating force
c notabl= 2 : Lennard Jones Force
c 3 : Lennard Jones Energy
c periodicflag
c nf_per= 1 : periodic boundary condition
c 0 : not periodic
c natchangeflag
c nf_chg= 1 : sorting the atom type
c 0 : not sorting without first step
c (at first time atom type is sorted)
        notabl=2
        nf_per=1
        nf_chg=0
c*
c in this code at host (force and energy)
        call calc_force_eng_vdw(pos,nmlcl,natype,natot,
     - gsforc,gseng,rseng,rsforc,xvolum,nf_per,
     - f_vdwh,e_vdwh)
c*
c
c*
c using RIKEN's library (by Dr.Tetsu Narumi) at MDone'
c        call MR1calcvdw(pos,nmlcl,natype,natot,gsforc,rsforc,notabl,
        call MR3calcvdw(pos,nmlcl,natype,natot,gsforc,rsforc,notabl,
     - xvolum,nf_per,nf_chg,f_vdwg)
c        call MR3calcvdw_ij(nmlcl,pos,natype,f_vdwg,
c     -       nmlcl,pos,natype,natot,gsforc,rsforc,notabl,
c     -       xvolum,nf_per)
c*
c using RIKEN's library (by Dr.Tetsu Narumi) at MDone'
c call MR1calcvdw_ij(nmlcl,pos,natype,f_vdwg,
c - nmlcl,pos,natype,natot,gsforc,rsforc,notabl,
c - xvolum,nf_per)
c*


        if(mod(imrpt,ncoum).eq.0) then
        notabl=3
c call MR1calcvdw(pos,nmlcl,natype,natot,gseng,rseng,notabl,
c - xvolum,nf_per,nf_chg,e_vdwg)
c call MR1calcvdw_ij(nmlcl,pos,natype,e_vdwg,
c - nmlcl,pos,natype,natot,gsforc,rsforc,notabl,
c - xvolum,nf_per)
c
c 010611 at PNL, MR1calcvdw_ij on potential mode does not work
c
c        call MR1calcvdw(pos,nmlcl,natype,natot,gseng,rseng,notabl,
        call MR3calcvdw(pos,nmlcl,natype,natot,gseng,rseng,notabl,
     - xvolum,nf_per,nf_chg,e_vdwg)
c        call MR3calcvdw_ij(nmlcl,pos,natype,e_vdwg,
c     -       nmlcl,pos,natype,natot,gseng,rseng,notabl,
c     -       xvolum,nf_per)

        endif
                call write_eng_forc(notabl,nmlcl,e_vdwh,e_vdwg,
     - f_vdwh,f_vdwg,natype)
c
                call cputime(time3)
c goto 1099
c
c********
c calc coulomb force
c table number for calculating force
c notabl= 0 : Coulomb Force
c 1 : Coulomb Energy
c (in this mode We don't need nf_chg)'
 1098   notabl=0
        nf_per=1
        nf_chg=0
        alphat=1.0d0
c***
c ver.2.0.2 (March 18, 1998)
c fixed for no-periodic conditon
       call calc_force_eng_coul(pos,nmlcl,charge,alphat,notabl,
     - xvolum,nf_per,f_c_dh,e_c_dh)
c       call MR1calccoulomb(pos,nmlcl,charge,alphat,notabl,xvolum,
       call MR3calccoulomb(pos,nmlcl,charge,alphat,notabl,xvolum,
     - nf_per,nf_chg,f_c_dg)
c       call MR3calccoulomb_ij(nmlcl,pos,charge,f_c_dg,
c     - nmlcl,pos,charge,
c     - alphat,notabl,xvolum,nf_per)


       do i=1,nmlcl
          do j=1,3
             f_c_dg(j,i)=f_c_dg(j,i)*charge(i)
          enddo
       enddo

       if(mod(imrpt,ncoum).eq.0) then
        notabl=notabl+1
       call MR3calccoulomb(pos,nmlcl,charge,alphat,notabl,xvolum,
     - nf_per,nf_chg,e_c_dg)
c       call MR1calccoulomb_ij(nmlcl,pos,charge,e_c_dg,
c       call MR3calccoulomb_ij(nmlcl,pos,charge,e_c_dg,
c     - nmlcl,pos,charge,
c     - alphat,notabl,xvolum,nf_per)


       do i=1,nmlcl
          e_c_dg(1,i)=e_c_dg(1,i)*charge(i)
       enddo
       endif
                call write_eng_forc(notabl,nmlcl,e_c_dh,e_c_dg,
     - f_c_dh,f_c_dg,natype)
c
                call cputime(time4)
c       goto 1099
c
c********
c calc Ewald real space force
c table number for calculating force
c notabl= 6 : Ewald Real space Force
c 7 : Ewald Real space Energy
c (in this mode We don't need nf_chg)'
        notabl=6
        nf_per=1
        nf_chg=0
c***
c ver.2.0.2 (March 18, 1998)
c fixed for no-periodic conditon
       call calc_force_eng_coul(pos,nmlcl,charge,alpha,notabl,
     - xvolum,nf_per,f_c_rh,e_c_rh)
c       call MR1calccoulomb(pos,nmlcl,charge,alpha,notabl,xvolum,
       call MR3calccoulomb(pos,nmlcl,charge,alpha,notabl,xvolum,
     - nf_per,nf_chg,f_c_rg)



       if(mod(imrpt,ncoum).eq.0) then
        notabl=notabl+1
c       call MR1calccoulomb(pos,nmlcl,charge,alpha,notabl,xvolum,
       call MR3calccoulomb(pos,nmlcl,charge,alpha,notabl,xvolum,
     - nf_per,nf_chg,e_c_rg)


       endif
                call write_eng_forc(notabl,nmlcl,e_c_rh,e_c_rg,
     - f_c_rh,f_c_rg,natype)
c
                call cputime(time1)
c*********
c calc coulomb force at Ewald sum of wave space part
c
  190 notabl=8
        call calc_force_coul_w(kmax,fac_k,kvect,knum,pos,nmlcl,
     - charge,alpha,eps,xvolum,f_c_wh,e_c_wh)
c        call MR3calcewald_host(kvect,knum,pos,nmlcl,charge,alpha,eps,
c     -       cell,f_c_wh,tpot,stress)
        tpot=0.0d0
        do i=1,nmlcl
           tpot=tpot+e_c_wh(1,i)
        enddo
        tpot=tpot*0.5d0
        write(*,*) 'tpot host ',tpot

        tpot=0.0d0
        do i=1,nmlcl
           do j=1,3
              f_c_wg(j,i)=0.0d0
           enddo
        enddo
        call MR3calcewald(kvect,knum,pos,nmlcl,charge,alpha,eps,
     -       cell,f_c_wg,tpot,stress)
        write(*,*) 'tpot grape',tpot

        if(mod(imrpt,ncoum).eq.0) then
c if(1.eq.0) then
        call MR3calcewald(kvect,-knum,pos,nmlcl,charge,alpha,eps,
     -          cell,e_c_wg,tpot,stress)
        endif
c
                call write_eng_forc(notabl,nmlcl,e_c_wh,e_c_wg,
     - f_c_wh,f_c_wg,natype)
c
 1099 call cputime(time2)
c
c********
c
        call sum_part_force_eng(nmlcl,forceg,f_c_rg,f_c_wg,f_vdwg,
     - engpg,e_c_rg,e_c_wg,e_vdwg)
        call sum_part_force_eng(nmlcl,forceh,f_c_rh,f_c_wh,f_vdwh,
     - engph,e_c_rh,e_c_wh,e_vdwh)

        call motion(pos,vmom,nmlcl,engt,forceg,smass,delt,xvolum,
     - natype,natot)
c call motion(pos,vmom,nmlcl,engt,forceh,smass,delt,xvolum,
c - natype,natot)

c if(mod(imrpt,ncoum).eq.0) then
c call output(nmlcl,imrpt,nuprnt,engt,engpg)
c call output(nmlcl,imrpt,nuprnt,engt,engph)
c endif
                write(*,*) time1,time2,time3,time4
c
c wait if other process wants to get MDone or MDGRAPE-2 board
c
c        call MR1wait
c
 1100 continue
c
c finalize MDone or MDGRAPE-2
c


c        call MR1free
        call MR3free


c
c
        end
c********************************************************************
      subroutine sum_part_force_eng(nmlcl,force,f_c_r,f_c_w,f_vdw,
     - engp,e_c_r,e_c_w,e_vdw)
c %%%%%%%%%%%%%%%%
c %%% sum forces %%%%%%%<< sum forces and potentials>>%%%%%%%%%%%%%%%
c %%% 97/12/25 %%%
c %%%%%%%%%%%%%%%%
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      real*8 force(3,nmlcl)
      real*8 f_c_r(3,nmlcl),f_c_w(3,nmlcl),f_vdw(3,nmlcl)
      real*8 engp(3,nmlcl)
      real*8 e_c_r(3,nmlcl),e_c_w(3,nmlcl),e_vdw(3,nmlcl)
c
      do 10 j=1,nmlcl,1
         engp(1,j)=e_c_r(1,j)+e_c_w(1,j)+e_vdw(1,j)
         do 10 i=1,3,1
            force(i,j)=f_c_r(i,j)+f_c_w(i,j)+f_vdw(i,j)
   10 continue
c

      return
      end
c


c********************************************************************
      subroutine set_rkvect(knum,kvect,rkvect)
c %%%%%%%%%%%%%%%
c % change kvect%%%%<< change kvect from int to real>>%%%%%%%%%%%%%%%
c % 98/01/05 %
c %%%%%%%%%%%%%%%
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      integer*4 kvect(3,knum)
      real*8 rkvect(3,knum)

        do 1000 i=1,knum,1
        do 1000 j=1,3,1
                rkvect(j,i)=kvect(j,i)
c write(*,*)'rkvect',rkvect(j,i)
 1000 continue

      return
      end
c********************************************************************
      subroutine set_kvect(kmax,xvolum,eps,alpha,knum,kvect,fac_k)
c %%%%%%%%%%%%%%%%
c %%% set_knum %%%%%%%%<< setting wave number of ewald sum >>%%%%%%%
c %%% 98/01/05 %%%
c %%%%%%%%%%%%%%%%
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      integer*4 kvect(3,4*kmax*kmax*kmax),knum
      real*8 fac_k(4*kmax*kmax*kmax)
c**********
c* summation over position
c* the values of l, m, & n are selected so that the symmetry of
c* reciprocal lattice is taken into account:
c* l : 0 --> lmax
c* m : 0 --> mmax ( l=0 )
c* -mmax --> mmax ( otherwise )
c* n : 1 --> nmax ( l=m=0 )
c* -nmax --> nmax ( otherwise )
c**********
        data pi/3.1415926535897932385d0/
c
      lmax=kmax
      mmax=kmax
      nmax=kmax
      kmaxsq=kmax*kmax
c
      ewfcgx=2d0*pi/xvolum
c
         knum=0
      do 1000 l=0,lmax,1
        if (l.eq.0) then
          mmin=0
        else
          mmin=-mmax
        endif
c
      do 1000 m=mmin,mmax,1
c
        if ((l.eq.0).and.(m.eq.0)) then
          nmin=1
        else
          nmin=-nmax
        endif
c
      do 1000 n=nmin,nmax,1
c
         ksq=l*l+m*m+n*n
         gsq=ksq*ewfcgx**2
c
        if (ksq.gt.kmaxsq) go to 1000
        knum=knum+1

        kvect(1,knum)=l
        kvect(2,knum)=m
        kvect(3,knum)=n
        fac_k(knum)=2d0/eps/xvolum**3*exp(-0.25d0*gsq/alpha**2)/gsq

c        kvect(1,knum)=abs(l)
c        kvect(2,knum)=abs(m)
c        kvect(3,knum)=abs(n)
c        kvect(1,knum)=1
c        kvect(2,knum)=1
c        kvect(3,knum)=1
c        ksq=kvect(1,knum)*kvect(1,knum)+
c     -       kvect(2,knum)*kvect(2,knum)+
c     -       kvect(3,knum)*kvect(3,knum)
c        gsq=ksq*ewfcgx**2
c        fac_k(knum)=2d0/eps/xvolum**3*exp(-0.25d0*gsq/alpha**2)/gsq

c write(*,*)'knum',knum
 1000 continue
c write(*,*) (kvect(1,i),kvect(2,i),kvect(3,i),i=2000,knum)
c write(*,*) 'fac_k',fac_k(1000)
        return
        end
c********************************************************************
      subroutine calc_force_coul_w(kmax,fac_k,kvect,knum,pos,nmlcl,
     - charge,alpha,eps,xvolum,f_c_wh,e_c_wh)
c %%%%%%%%%%%%%%%%
c %%% forces %%%%%%%%<< fourier component of ewalt sum >>%%%%%%%%%%
c %%% 97/12/25 %%%
c %%%%%%%%%%%%%%%%
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
c*** kkmax and nnn are temporary parameter
      parameter(kkmax=10,nnn=32000)
c
      integer*4 kvect(3,knum)
      real*8 pos(3,nmlcl),charge(nmlcl)
      real*8 f_c_wh(3,nmlcl),e_c_wh(3,nmlcl)
      real*8 fac_k(knum)
c
c*** local dimension **
      real*8 cosx(nnn*3,0:kkmax),sinx(nnn*3,0:kkmax),
     - cosy(nnn*3,0:kkmax),siny(nnn*3,0:kkmax),
     - cosz(nnn*3,0:kkmax),sinz(nnn*3,0:kkmax),
     - cosp(nnn*3),sinp(nnn*3)

      logical lsgnm,lsgnn
c
        data pi/3.1415926535897932385d0/
c
        if(nmlcl.gt.nnn) then
                write(*,'(a,i7,a,i7)')
     - 'nmlcl>nnn: nmlcl =',nmlcl,'  nnn',nnn
                stop '## error in calc_force_coul_w ##'
        elseif(kmax.gt.kkmax) then
                write(*,'(a,i7,a,i7)')
     - 'kmax>kkmax: kmax =',kmax,'  kkmax',kkmax
                stop '## error in calc_force_coul_w ##'
        endif
c
        ewfcgx=2d0*pi/xvolum
        ewfcgy=2d0*pi/xvolum
        ewfcgz=2d0*pi/xvolum
        lmax=kmax
        mmax=kmax
        nmax=kmax
c
c**********
c* sin & cos elements
c**********
c
      do 1000 ij=1,nmlcl,1
        cosx(ij,0)=1.0
        cosy(ij,0)=1.0
        cosz(ij,0)=1.0
        sinx(ij,0)=0.0
        siny(ij,0)=0.0
        sinz(ij,0)=0.0
        cosx(ij,1)=cos(ewfcgx*pos(1,ij))
        cosy(ij,1)=cos(ewfcgy*pos(2,ij))
        cosz(ij,1)=cos(ewfcgz*pos(3,ij))
        sinx(ij,1)=sin(ewfcgx*pos(1,ij))
        siny(ij,1)=sin(ewfcgy*pos(2,ij))
        sinz(ij,1)=sin(ewfcgz*pos(3,ij))
 1000 continue
c
      do 1001 k=2,lmax,1
      do 1001 ij=1,nmlcl,1
        cosx(ij,k)=cosx(ij,k-1)*cosx(ij,1)-sinx(ij,k-1)*sinx(ij,1)
        sinx(ij,k)=sinx(ij,k-1)*cosx(ij,1)+cosx(ij,k-1)*sinx(ij,1)
 1001 continue
      do 1002 k=2,mmax,1
      do 1002 ij=1,3*nmlcl,1
        cosy(ij,k)=cosy(ij,k-1)*cosy(ij,1)-siny(ij,k-1)*siny(ij,1)
        siny(ij,k)=siny(ij,k-1)*cosy(ij,1)+cosy(ij,k-1)*siny(ij,1)
 1002 continue
      do 1003 k=2,nmax,1
      do 1003 ij=1,3*nmlcl,1
        cosz(ij,k)=cosz(ij,k-1)*cosz(ij,1)-sinz(ij,k-1)*sinz(ij,1)
        sinz(ij,k)=sinz(ij,k-1)*cosz(ij,1)+cosz(ij,k-1)*sinz(ij,1)
 1003 continue

      do 1200 k=1,knum,1
           l=kvect(1,k)
           m=kvect(2,k)
           n=kvect(3,k)
        gx=ewfcgx*l
        gy=ewfcgy*m
        ma=abs(m)
        gz=ewfcgz*n
        na=abs(n)
c
        if (m.lt.0) then
          lsgnm=.false.
        else
          lsgnm=.true.
        endif
c
        if (n.lt.0) then
          lsgnn=.false.
        else
          lsgnn=.true.
        endif
c
        gsq=gx*gx+gy*gy+gz*gz
c
c fac_k=2.0*ewfcq*exp(gsq*ewasqi)/gsq
c g=2.0*(-1.0/gsq+ewasqi)
c vxfct=fac_k(k)*(1.0+g*gx*gx)
c vyfct=fac_k(k)*(1.0+g*gy*gy)
c vzfct=fac_k(k)*(1.0+g*gz*gz)
c
        if (lsgnm) then
          if (lsgnn) then
            do 1210 ij=1,nmlcl,1
              cmn= cosy(ij,ma)*cosz(ij,na)-siny(ij,ma)*sinz(ij,na)
              smn= siny(ij,ma)*cosz(ij,na)+cosy(ij,ma)*sinz(ij,na)
              cosp(ij)=charge(ij)*(cosx(ij,l)*cmn-sinx(ij,l)*smn)
              sinp(ij)=charge(ij)*(sinx(ij,l)*cmn+cosx(ij,l)*smn)
 1210 continue
          else
            do 1211 ij=1,nmlcl,1
              cmn= cosy(ij,ma)*cosz(ij,na)+siny(ij,ma)*sinz(ij,na)
              smn= siny(ij,ma)*cosz(ij,na)-cosy(ij,ma)*sinz(ij,na)
              cosp(ij)=charge(ij)*(cosx(ij,l)*cmn-sinx(ij,l)*smn)
              sinp(ij)=charge(ij)*(sinx(ij,l)*cmn+cosx(ij,l)*smn)
 1211 continue
          endif
        else
          if (lsgnn) then
            do 1212 ij=1,nmlcl,1
              cmn= cosy(ij,ma)*cosz(ij,na)+siny(ij,ma)*sinz(ij,na)
              smn=-siny(ij,ma)*cosz(ij,na)+cosy(ij,ma)*sinz(ij,na)
              cosp(ij)=charge(ij)*(cosx(ij,l)*cmn-sinx(ij,l)*smn)
              sinp(ij)=charge(ij)*(sinx(ij,l)*cmn+cosx(ij,l)*smn)
 1212 continue
          else
            do 1213 ij=1,nmlcl,1
              cmn= cosy(ij,ma)*cosz(ij,na)-siny(ij,ma)*sinz(ij,na)
              smn=-siny(ij,ma)*cosz(ij,na)-cosy(ij,ma)*sinz(ij,na)
              cosp(ij)=charge(ij)*(cosx(ij,l)*cmn-sinx(ij,l)*smn)
              sinp(ij)=charge(ij)*(sinx(ij,l)*cmn+cosx(ij,l)*smn)
 1213 continue
          endif
        endif
c
          csum=0.0
          ssum=0.0
        do 1219 ij=1,nmlcl,1
          csum=csum+cosp(ij)
          ssum=ssum+sinp(ij)
 1219 continue
c if(k.eq.1) then
c write(*,*)'c-ssum(bc-bs)',csum,ssum
c endif
c**********
c* potential & force
c**********

        do 1220 i=1,nmlcl,1
           e=csum*cosp(i)+ssum*sinp(i)
          e_c_wh(1,i)=e_c_wh(1,i)+fac_k(k)*e
          s=-fac_k(k)*(ssum*cosp(i)-csum*sinp(i))
           f_c_wh(1,i)=f_c_wh(1,i)+gx*s
           f_c_wh(2,i)=f_c_wh(2,i)+gy*s
           f_c_wh(3,i)=f_c_wh(3,i)+gz*s
c vrlx(i)=vrlx(i)+vxfct*e
c vrly(i)=vrly(i)+vyfct*e
c vrlz(i)=vrlz(i)+vzfct*e
 1220 continue
c
 1200 continue
c write(*,*)(frcxs(i,1),frcys(i,1),frczs(i,1),i=1,10)
c
c**********
c* convert site-site virial to molecular virial
c**********
c do 1300 i=1,nmlcl,1
c vrlx(i)=vrlx(i)-(stx(i,1)*frcxs(i,1)+stx(i,2)*frcxs(i,2)+
c - stx(i,3)*frcxs(i,3))
c vrly(i)=vrly(i)-(sty(i,1)*frcys(i,1)+sty(i,2)*frcys(i,2)+
c - sty(i,3)*frcys(i,3))
c vrlz(i)=vrlz(i)-(stz(i,1)*frczs(i,1)+stz(i,2)*frczs(i,2)+
c - stz(i,3)*frczs(i,3))
c1300 continue
        return
        end
c
c********************************************************************
      subroutine calc_force_eng_coul(pos,nmlcl,charge,alpha,notabl,
     - xvolum,nf_per,forceh,engph)
c %%%%%%%%%%%%%%%%
c %%% forces %%%%%%%%<< calulate coulom forces >>%%%%%%%%%%%%%%%%%
c %%% 97/12/25 %%%
c %%%%%%%%%%%%%%%%
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      real*8 pos(3,nmlcl)
      real*8 engph(3,nmlcl)
      real*8 forceh(3,nmlcl)
      real*8 charge(nmlcl)
        data rcut /8.0d0/
        data pi/3.1415926535897932385d0/

        zvolum=xvolum
        yvolum=xvolum
c***
c ver.2.0.2 (March 18, 1998)
c fixed for non-periodic conditon
        if(nf_per.eq.1) then
                xvhi=2.0d0/xvolum
                yvhi=2.0d0/yvolum
                zvhi=2.0d0/zvolum
        elseif(nf_per.eq.0) then
                xvhi=2.0d0/(xvolum*10)
                yvhi=2.0d0/(yvolum*10)
                zvhi=2.0d0/(zvolum*10)
        endif

c
        rcutsq=rcut*rcut
c
        if(notabl.eq.0) then
                goto 10
        elseif(notabl.eq.6) then
                goto 20
        else
                write(*,*) 'error in calc_force_coul: notabl =',notabl
                stop '## error in calc_force_eng_coul ##'
        endif
c
   10 do 1010 j=2,nmlcl
        do 1010 i=1,j-1
                z=pos(3,i)-pos(3,j)
                zd=z-int(z*zvhi)*zvolum
c if (abs(zd).gt.rcut) go to 1010
                y=pos(2,i)-pos(2,j)
                yd=y-int(y*yvhi)*yvolum
c if (abs(yd).gt.rcut) go to 1010
                x=pos(1,i)-pos(1,j)
                xd=x-int(x*xvhi)*xvolum
c
                rd=xd**2+yd**2+zd**2
c if (rd.gt.rcutsq) go to 1010
c
                rsqi=1.0d0/rd
                rsq1i=sqrt(rsqi)
                ep=charge(i)*charge(j)*rsq1i
                fc=ep*rsqi
                engph(1,i)=engph(1,i)+ep
                engph(1,j)=engph(1,j)+ep
                fx=xd*fc
                fy=yd*fc
                fz=zd*fc
                forceh(1,i)=forceh(1,i)+fx
                forceh(2,i)=forceh(2,i)+fy
                forceh(3,i)=forceh(3,i)+fz
                forceh(1,j)=forceh(1,j)-fx
                forceh(2,j)=forceh(2,j)-fy
                forceh(3,j)=forceh(3,j)-fz

 1010 continue
        return
c*****

   20 do 2010 j=2,nmlcl
        do 2010 i=1,j-1
                z=pos(3,i)-pos(3,j)
                zd=z-int(z*zvhi)*zvolum
c if (abs(zd).gt.rcut) go to 1010
                y=pos(2,i)-pos(2,j)
                yd=y-int(y*yvhi)*yvolum
c if (abs(yd).gt.rcut) go to 1010
                x=pos(1,i)-pos(1,j)
                xd=x-int(x*xvhi)*xvolum
c
                rd=xd**2+yd**2+zd**2
c if (rd.gt.rcutsq) go to 1010
c
                rsqi=1.0d0/rd
                rsq1i=sqrt(rsqi)
                rsq1=sqrt(rd)
c write(*,*) 'alpha',alpha,rsq1,derfc(rsq1*alpha)
                ep=charge(i)*charge(j)*rsq1i*derfc(rsq1*alpha)
                fc=(ep+2d0*charge(i)*charge(j)*alpha/sqrt(pi)
     - *exp(-rd*alpha**2))*rsqi
                engph(1,i)=engph(1,i)+ep
                engph(1,j)=engph(1,j)+ep
                fx=xd*fc
                fy=yd*fc
                fz=zd*fc
                forceh(1,i)=forceh(1,i)+fx
                forceh(2,i)=forceh(2,i)+fy
                forceh(3,i)=forceh(3,i)+fz
                forceh(1,j)=forceh(1,j)-fx
                forceh(2,j)=forceh(2,j)-fy
                forceh(3,j)=forceh(3,j)-fz

 2010 continue


c******
      return
      end


c********************************************************************
      subroutine calc_force_eng_vdw(pos,nmlcl,natype,natot,gsforc,
     - gseng,rseng,rsforc,xvolum,nf_per,forceh,engph)
c %%%%%%%%%%%%%%%%
c %%% forces %%%%%%%%<< calculate forces >>%%%%%%%%%%%%%%%%%%%%%%%%
c %%% 97/10/06 %%%
c %%%%%%%%%%%%%%%%
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      integer*4 natype(nmlcl)
      real*8 gseng(natot,natot),gsforc(natot,natot)
      real*8 rseng(natot,natot),rsforc(natot,natot)
      real*8 pos(3,nmlcl)
      real*8 engph(3,nmlcl)
      real*8 forceh(3,nmlcl)
c
        data rcut /8.0d0/
c
        zvolum=xvolum
        yvolum=xvolum
c***
c ver.2.0.2 (March 18, 1998)
c fixed for non-periodic conditon
        if(nf_per.eq.1) then
                xvhi=2.0d0/xvolum
                yvhi=2.0d0/yvolum
                zvhi=2.0d0/zvolum
        elseif(nf_per.eq.0) then
                xvhi=2.0d0/(xvolum*10)
                yvhi=2.0d0/(yvolum*10)
                zvhi=2.0d0/(zvolum*10)
        endif

c
        rcutsq=rcut*rcut

        do 1010 j=2,nmlcl
        do 1010 i=1,j-1
                z=pos(3,i)-pos(3,j)
                zd=z-int(z*zvhi)*zvolum
                y=pos(2,i)-pos(2,j)
                yd=y-int(y*yvhi)*yvolum
                x=pos(1,i)-pos(1,j)
                xd=x-int(x*xvhi)*xvolum
c
                rd=xd**2+yd**2+zd**2
c***
c ver.2.0.1 (March 12, 1998)
c fixed for cut off part
c
                rs1=rsforc(natype(i),natype(j))
        if (rd*rs1.gt.rcutsq) go to 1010
                rs=rseng(natype(i),natype(j))
c
                rsqi=1.0d0/rd
                r06=rsqi*rsqi*rsqi
                r12=r06*r06
c
                gs=gseng(natype(i),natype(j))
                        cett12=gs*rs**(-6)
                        cett06=-gs*rs**(-3)
c
                if (rd*rs.le.rcutsq) then
                                ep=cett12*r12+cett06*r06
                        engph(1,i)=engph(1,i)+ep
                        engph(1,j)=engph(1,j)+ep
                endif
c
                        cftt12=cett12*12
                        cftt06=cett06*6
                        fc=(cftt12*r12+cftt06*r06)*rsqi
                fx=xd*fc
                fy=yd*fc
                fz=zd*fc
                forceh(1,i)=forceh(1,i)+fx
                forceh(2,i)=forceh(2,i)+fy
                forceh(3,i)=forceh(3,i)+fz
                forceh(1,j)=forceh(1,j)-fx
                forceh(2,j)=forceh(2,j)-fy
                forceh(3,j)=forceh(3,j)-fz

 1010 continue


c******
      return
      end
c********************************************************************
      subroutine tmpset(vmom,nmlcl,temp,smass,natype,natot)
c %%%%%%%%%%%%%%%%%
c %%% tmpset %%%%%%%%%<< reset temperature >>%%%%%%%%%%%%%%%%%%%%%
c %%% 97/12/15 %%%
c %%%%%%%%%%%%%%%%%
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      integer*4 natype(nmlcl)
      real*8 smass(natot)
      real*8 vmom(3,nmlcl)
c*****<< initial >>*****
        ekt=0.0d0
        vx=0.0d0
        vy=0.0d0
        vz=0.0d0
      do 1010 i=1,nmlcl,1
        vx=vx+vmom(1,i)
        vy=vy+vmom(2,i)
        vz=vz+vmom(3,i)
 1010 continue
        vx=vx/nmlcl
        vy=vy/nmlcl
        vz=vz/nmlcl
      do 1030 i=1,nmlcl,1
        smassi=1d0/smass(natype(i))
        vmom(1,i)=vmom(1,i)-vx
        vmom(2,i)=vmom(2,i)-vy
        vmom(3,i)=vmom(3,i)-vz
        ekt=ekt+(vmom(1,i)**2+vmom(2,i)**2+vmom(3,i)**2)*smassi
 1030 continue
        ekt=0.5d0*ekt/nmlcl
c write(*,*) ekt
        vfc=sqrt(1.5d0*temp/ekt)
      do 1050 i=1,nmlcl,1
        vmom(1,i)=vmom(1,i)*vfc
        vmom(2,i)=vmom(2,i)*vfc
        vmom(3,i)=vmom(3,i)*vfc
 1050 continue
c write(*,*) vmom(1,1),vmom(2,1),vmom(3,1)

                                  ekt=0.0d0
                         do 3011 i=1,nmlcl,1
        smassi=1d0/smass(natype(i))
        ekt=ekt+(vmom(1,i)**2+vmom(2,i)**2+vmom(3,i)**2)*smassi
 3011 continue
                         ekt=0.5d0*ekt/nmlcl
c write(*,9500) ekt/1.5d0,temp
 9500 format(2f10.4)
c******
      return
      end

c

c********************************************************************
      subroutine output(nmlcl,imrpt,nuprnt,engt,engpg)
c %%%%%%%%%%%%%%%%
c %%% output %%%%%%%%<< output routine >>%%%%%%%%%%%%%%%%%%%%%%%%%%
c %%% 97/12/15 %%%
c %%%%%%%%%%%%%%%%
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      real*8 engt(nmlcl)
      real*8 engpg(3,nmlcl)
c
                tengt=0d0
                tengp=0d0
        do 1010 i=1,nmlcl,1
                tengt=tengt+engt(i)
                tengp=tengp+0.5d0*engpg(1,i)
 1010 continue
                teng=tengt+tengp
                write(nuprnt,9000) imrpt,tengt,tengp,teng
 9000 format(i6,3e15.7)

c******
      return
      end
c*********************************************************************
      subroutine motion(pos,vmom,nmlcl,engt,forceg,smass,delt,xvolum,
     - natype,natot)
c %%%%%%%%%%%%%%%%
c %%% motion %%%%%%%%<< motion under given forces >>%%%%%%%%%%%%%%%%
c %%% 97/12/15 %%%
c %%%%%%%%%%%%%%%%
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      integer*4 natype(nmlcl)
      real*8 smass(natot)
      real*8 pos(3,nmlcl),vmom(3,nmlcl)
      real*8 forceg(3,nmlcl)
      real*8 engt(nmlcl)
c
        do 1010 i=1,nmlcl,1
c write(*,*) 'smass',i,natype(i)
c stop
            smassi=1.0d0/smass(natype(i))
c write(*,*) i,natype(i),smassi
                vmomdx=forceg(1,i)*delt
                vmomdy=forceg(2,i)*delt
                vmomdz=forceg(3,i)*delt
          e=0.5d0*((vmom(1,i)+0.5d0*vmomdx)**2
     - +(vmom(2,i)+0.5d0*vmomdy)**2
     - +(vmom(3,i)+0.5d0*vmomdz)**2)*smassi
          engt(i)=e
c tengt=tengt+e
        vmom(1,i)=vmom(1,i)+vmomdx
        vmom(2,i)=vmom(2,i)+vmomdy
        vmom(3,i)=vmom(3,i)+vmomdz
          x=pos(1,i)+vmom(1,i)*delt*smassi
          y=pos(2,i)+vmom(2,i)*delt*smassi
          z=pos(3,i)+vmom(3,i)*delt*smassi
c
        if (x.lt.0) then
            pos(1,i)=x+xvolum
          else if (x.lt.xvolum) then
            pos(1,i)=x
          else
            pos(1,i)=x-xvolum
          endif
        if (y.lt.0) then
            pos(2,i)=y+xvolum
          else if (y.lt.xvolum) then
            pos(2,i)=y
          else
            pos(2,i)=y-xvolum
          endif
        if (z.lt.0) then
            pos(3,i)=z+xvolum
          else if (z.lt.xvolum) then
            pos(3,i)=z
          else
            pos(3,i)=z-xvolum
          endif
 1010 continue


c******
      return
      end



c*********************************************************************
      subroutine sysdef(nmlcl,natot,natype,gseng,gsforc,rseng,rsforc,
     - smass,charge,alpha,nsty_f)
c %%%%%%%%%%%%%%%%
c %%% sysdef %%%%%%%%<< defines system parameters >>%%%%%%%%%%%%%%%%
c %%% 97/12/15 %%%
c %%%%%%%%%%%%%%%%
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
c parameter (natot=1)
c parameter (natot=4)
c*** natot2 is temporary parameter
      parameter (natot2=20)
c
      integer*4 natype(nmlcl)
      real*8 gseng(natot,natot),gsforc(natot,natot)
      real*8 rseng(natot,natot),rsforc(natot,natot)
      real*8 smass(natot)
      real*8 charge(nmlcl)
c
c*** local dimension **
        real*8 eps(natot2),sgm(natot2)
c
        idum=100
c
        if(natot.gt.natot2) then
                write(*,'(a,i7,a,i7)')
     - 'natot>natot2: natot =',natot,'  natot2',natot2
                stop '## error in sysdef ##'
        endif
c
c
c*** set the parameter alpha in Ewald Sum
c
c        alpha=1d1
c        alpha=1d0
c        alpha=1d-1
        alpha=.8d-1
c        alpha=0.3d-1
c        alpha=2d-1
c****
c set atom type
c****
c
        do 1010 i=1,nmlcl
                natype(i)=int(ran2(idum)*natot+1)
                if(natype(i).gt.natot) natype(i)=natot
                if(natype(i).le.0) natype(i)=1
 1010 continue
c
c****
c set gscale and rscale
c****
        do 2000 i=1,natot,1
           sgm(i)=ran2(idum)+0.5
           eps(i)=ran2(idum)+0.5
 2000 continue
c
        do 2010 j=1,natot,1
        do 2010 i=1,natot,1
                rseng(i,j)=(0.5d0*(sgm(i)+sgm(j)))**(-2)
                gseng(i,j)=4.0d0*sqrt(eps(i)*eps(j))
        if(nsty_f.eq.0)then
                rsforc(i,j)=2d0**(-1d0/3d0)*rseng(i,j)
                gsforc(i,j)=3d0*gseng(i,j)*rsforc(i,j)
        elseif(nsty_f.eq.1)then
                rsforc(i,j)=rseng(i,j)
                gsforc(i,j)=6d0*gseng(i,j)*rsforc(i,j)
        endif
 2010 continue
c
c****
c set particle mass
c****
        do 3000 i=1,natot,1
           smass(i)=ran2(idum)+0.5
 3000 continue
c
c****
c set charge
c****
        sumch=0.0d0
        do 3010 i=1,nmlcl
           charge(i)=ran2(idum)+0.5
           sumch=sumch+charge(i)
 3010 continue
           sumch=sumch/nmlcl
        do 3020 i=1,nmlcl
           charge(i)=charge(i)-sumch
 3020 continue
c
c******

      return
      end

c*********************************************************************
      subroutine start(pos,vmom,nmlcl,xvolum)
c %%%%%%%%%%%%%%%%
c %%% start %%%%%%%%%<< setup initial configuration >>%%%%%%%%%%%%
c %%% 97/12/15 %%%
c %%%%%%%%%%%%%%%%
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      real*8 pos(3,nmlcl),vmom(3,nmlcl)
c
      idum=100
c
                num3=nmlcl**(1.0/3.0)
        do 1010 j=1,nmlcl,1
                do 1020 i=1,3,1
                   k=(j-1)/num3**(i-1)
                   pos(i,j)=xvolum/num3*(0.5*ran2(idum)+k-k/num3*num3)
c pos(i,j)=xvolum/num3*(k-k/num3*num3)
                   vmom(i,j)=2d0*ran2(idum)-1d0
 1020 continue
 1010 continue

c******
      return
      end
c*********************************************************************
      subroutine values(nmlcl,engp,force)
c %%%%%%%%%%%%%%%%
c %%% values %%%%%%%%<< physical quantities >>%%%%%%%%%%%%%%%%%%%%%
c %%% 97/12/15 %%%
c %%%%%%%%%%%%%%%%
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      real*8 engp(3,nmlcl)
      real*8 force(3,nmlcl)
c
        do 1010 j=1,nmlcl,1
                do 1020 i=1,3,1
                        engp(i,j)=0d0
                        force(i,j)=0d0
 1020 continue
 1010 continue
c******
      return
      end

c*********************************************************************
      subroutine write_sys(nmlcl,natype,natot,gseng,gsforc,rseng,
     - rsforc,smass)
c %%%%%%%%%%%%%%%%%%%%%
c %%% write_sysdef %%%%%%%%<< write system parameters >>%%%%%%%%%%
c %%% 97/12/15 %
c %%%%%%%%%%%%%%%%%%%%%
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
c
      integer*4 natype(nmlcl)
      real*8 gseng(natot,natot),gsforc(natot,natot)
      real*8 rseng(natot,natot),rsforc(natot,natot)
      real*8 smass(natot)

       write(*,9000) ('i=',i,'  natype=',natype(i),i=1,nmlcl)
 9000 format(a,i3,a,i3)
c
       do j=1,natot,1
       write(*,9010) (i,j,gseng(i,j),gsforc(i,j),
     - rseng(i,j),rsforc(i,j),i=1,natot)
       enddo
 9010 format(2i3,4f10.3)
c
       write(*,9020) (i,'smass= ',smass(i),i=1,natot)
 9020 format(i3,a,f10.3)

c******
      return
      end
c*********************************************************************
      subroutine write_start(pos,vmom,nmlcl)
c %%%%%%%%%%%%%%%%%%%%
c %%% write_start %%%%%%%%%<< write initial configuration >>%%%%%%%
c %%% 97/12/15 %%%
c %%%%%%%%%%%%%%%%%%%%
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      real*8 pos(3,nmlcl),vmom(3,nmlcl)
c
        write(*,9000)('pos',pos(1,j),pos(2,j),pos(3,j),j=1,nmlcl)
 9000 format(a,3f10.3)
c******
      return
      end
c*********************************************************************
      subroutine write_eng_forc(notabl,nmlcl,engph,engpg,
     - forceh,forceg,natype)
c %%%%%%%%%%%%%%%%%%%%%%%
c %%% write_eng_forc %%%%%%%%%<< write eng and force >>%%%%%%%%%%%%
c %%% 98/03/10 %%%
c %%%%%%%%%%%%%%%%%%%%%%%
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      real*8 engph(3,nmlcl),engpg(3,nmlcl)
      real*8 forceh(3,nmlcl),forceg(3,nmlcl)
      integer*4 natype(nmlcl)
        character*30 result
c
        if(notabl.eq.0.or.notabl.eq.1) then
                result=' at Coulomb part'
        elseif(notabl.eq.6.or.notabl.eq.7) then
                result=' at Ewald Real space part'
        elseif(notabl.eq.2.or.notabl.eq.3) then
                result=' at Lennard-Jones part'
        elseif(notabl.eq.8) then
                result=' at Ewald Wave space part'
        endif
c
c if(notabl.ne.8) then
                write(*,*)
                write(*,*)'******************************************'
                write(*,*)'RESULTS of ENERGY',result
                write(*,*)'******************************************'
                write(*,*)
                write(*,*)'   n ','   energy_host','   energy_MDone',
     - ' e_host/e_MDone',' atom type'
c
c do i=1,nmlcl,1
                do i=1,10,1
                if(engpg(1,i).ne.0d0) then
                        write(*,9000) i,engph(1,i),engpg(1,i),
     - engph(1,i)/engpg(1,i),natype(i)
                else
                        write(*,9005) i,engph(1,i),engpg(1,i),natype(i)
                endif
                enddo
c endif
c
                write(*,*)
                write(*,*)'******************************************'
                write(*,*)'RESULTS of FORCE',result
                write(*,*)'******************************************'
c        do i=1,nmlcl,1
        do i=1,10,1
                write(*,*)
        write(*,9010) i, '(' , natype(i) , ')' ,
     - '      F_x      ',
     - '      F_y      ',
     - '      F_z      '
        write(*,9020) 'host--> ',forceh(1,i),forceh(2,i),forceh(3,i)
        write(*,9020) 'grape-> ',forceg(1,i),forceg(2,i),forceg(3,i)
             if(abs(forceg(1,i)).gt.0d0.and.abs(forceg(3,i)).gt.0d0.
     - and.abs(forceg(3,i)).gt.0d0)then
             write(*,9025) 'host/grape->',
     - forceh(1,i)/forceg(1,i),
     - forceh(2,i)/forceg(2,i),
     - forceh(3,i)/forceg(3,i)
        r=forceh(1,i)*forceh(1,i)+forceh(2,i)*forceh(2,i)
     - +forceh(3,i)*forceh(3,i)
        r=sqrt(r)
c write(*,*) 'r         ->',r
        d=0.0d0
        do j=1,3
           d=d+(forceg(j,i)-forceh(j,i))**2
        enddo
        d=sqrt(d)
        write(*,9035) 'error -log10((host-grape)/host) ->',-log10(d/r)
c write(*,9030) 'error     ->',
c - -log10(1.0d0*abs((forceg(1,i)-forceh(1,i))/r)),
c - -log10(1.0d0*abs((forceg(2,i)-forceh(2,i))/r)),
c - -log10(1.0d0*abs((forceg(3,i)-forceh(3,i))/r))
c - (abs((forceg(1,i)-forceh(1,i))/r)),
c - (abs((forceg(2,i)-forceh(2,i))/r)),
c - (abs((forceg(3,i)-forceh(3,i))/r))
c - forceg(1,i)-forceh(1,i),
c - forceg(2,i)-forceh(2,i),
c - forceg(3,i)-forceh(3,i)
             endif
        enddo
        goto 8000
        s=0.0d0
        do i=1,nmlcl,1
           r=forceh(1,i)*forceh(1,i)+forceh(2,i)*forceh(2,i)
     - +forceh(3,i)*forceh(3,i)
           r=sqrt(r)
           d=0.0d0
           do j=1,3
              d=d+(forceg(j,i)-forceh(j,i))**2
           enddo
           d=sqrt(d)
           s=s-log10(d/r)
c do j=1,3
c s=s-log10(1.0d0*abs((forceg(j,i)-forceh(j,i))/r))
c enddo
        enddo
c write(*,*) 'average error = ',s/(3.0d0*nmlcl)
        write(*,*) 'average error = ',s/(1.0d0*nmlcl)

 8000 s=0.0d0
        d=0.0d0
        do i=1,nmlcl
           r=0.0d0
           do j=1,3
              r=r+forceh(j,i)**2
           enddo
           r=sqrt(r)
           s=s+r

           r=0.0d0
           do j=1,3
              r=r+(forceh(j,i)-forceg(j,i))**2
           enddo
           r=sqrt(r)
           d=d+r
        enddo
        s=s/(1.0d0*nmlcl)
        d=d/(1.0d0*nmlcl)
        write(*,*) 'average force by host =',s
        write(*,*) 'average error by grape=',d
        write(*,*) 'average error log     =',-log10(d/s)

 9000 format(i5,2e15.7,f15.10,i3)
 9005 format(i5,2e15.7,i3)
 9010 format(i5,a,i2,a,3a)
 9020 format(a,3e15.7)
 9025 format(a,3f15.10)
 9030 format(a,3f15.10)
 9035 format(a,f15.10)
c9020 format(a,3e15.7)
c9000 format(i5,3e15.7,i3)
c9010 format(i5,3a)
c******
      return
      end

c*********************************************************************
      subroutine cputime(tim)
c %%%%%%%%%%%%%%%%
c %%% cputime %%%%%%%%<< estimate cputime >>%%%%%%%%%%%%%%%%
c %%% 97/12/15 %%%
c %%%%%%%%%%%%%%%%
c
      real*4 dumm(2)
c

        tim=0



c
      return
      end
c********************************************************************
      real*8 FUNCTION ran2(idum)
c FROM NUMERICAL RECIPES in Fortran 2nd.ed.(Cambridge univ.pr.,1992)
c %%%%%%%%%%%%%%%%
c %%% ran2 %%%%%%%%<< make the random number >>%%%%%%%%%%%%%%%%
c %%%%%%%%%%%%%%%%
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
   11 continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END

c********************************************************************
      SUBROUTINE indexx(n,arr,indx)
c FROM NUMERICAL RECIPES in Fortran 2nd.ed.(Cambridge univ.pr.,1992)
c %%%%%%%%%%%%%%%%
c %%% indexx %%%%%%%%<< indexing >>%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c %%%%%%%%%%%%%%%%
      INTEGER n,indx(n),M,NSTACK
c REAL arr(n)
      INTEGER arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
c REAL a
      INTEGER a
      do 11 j=1,n
        indx(j)=j
   11 continue
      jstack=0
      l=1
      ir=n
    1 if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
   12 continue
          i=0
    2 indx(i+1)=indxt
   13 continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
    3 continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
    4 continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
    5 indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
