#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "mr3.h"

#if cpu
void MR3calccoulomb_ij(int ni, double xi[], double qi[], double force[],
                       int nj, double xj[], double qj[],
                       double rscale, 
                       int tblno, double xmax, int periodicflag)
{
  /* periodic flag bit 0 : 0 --- non periodic, 1 --- periodic
                   bit 1 : 0 --- no multiplication of qi
                           1 --- multiplication of qi
  */
  int i,j,k;
  double dr[3],dtmp,r2,x,factor,rsqrt;
  static int rcutflag=0,ini=0;
  static double rcut,rcut2;
  char *s;
  int multiplyq=0;

  if(ini==0){
    if((s=getenv("MR3_HOST_CUTOFF"))!=NULL){
      sscanf(s,"%le",&rcut);
      rcut2=rcut*rcut;
      printf("rcut=%e is read from MR3_HOST_CUTOFF\n",rcut);
      rcutflag=1;
    }
    ini=1;
  }

  if((periodicflag & 2)!=0) multiplyq=1;

  if((periodicflag & 1)==0){
    xmax*=2.0;
  }
  for(i=0;i<ni;i++){
    for(j=0;j<nj;j++){
      r2=0.0;
      for(k=0;k<3;k++){
        dr[k]=xi[i*3+k]-xj[j*3+k];
        if(dr[k]<-xmax/2.0){
          dr[k]+=xmax;
        }
        if(dr[k]>=xmax/2.0){
          dr[k]-=xmax;
        }
        r2+=dr[k]*dr[k];
      }
      x=r2*rscale*rscale;
      if(r2!=0.0 && (rcutflag==0 || (rcutflag==1 && x<rcut2))){
        if(tblno==0){
          rsqrt=1.0/sqrt(r2);
          dtmp=qj[j]*rsqrt*rsqrt*rsqrt;
          if(multiplyq) dtmp*=qi[i];
          for(k=0;k<3;k++){
            force[i*3+k]+=dtmp*dr[k];
          }
        }
        else if(tblno==1){
          rsqrt=1.0/sqrt(r2);
          dtmp=qj[j]*rsqrt;
          if(multiplyq) dtmp*=qi[i];
          for(k=0;k<3;k++){
            force[i*3+k]+=dtmp;
          }
        }
        else if(tblno==6){
          x=r2*rscale*rscale;
          dtmp=qj[j]*(M_2_SQRTPI*exp(-x)*pow(x,-1.0)
                      + erfc(sqrt(x))*pow(x,-1.5));
          factor=rscale*rscale*rscale*qi[i];
          for(k=0;k<3;k++){
            force[i*3+k]+=dtmp*dr[k]*factor;
          }
        }
        else if(tblno==7){
          x=r2*rscale*rscale;
          dtmp=qj[j]*erfc(sqrt(x))*pow(x,-0.5);
          factor=rscale*qi[i];
          for(k=0;k<3;k++){
            force[i*3+k]+=dtmp*factor;
          }
        }
      }
    }
  }
}


void MR3calcvdw_ij(int ni, double xi[], int atypei[], double force[],
                   int nj, double xj[], int atypej[],
                   int nat, double gscale[], double rscale[],
                   int tblno, double xmax, int periodicflag)
{
  /* periodic flag 0 --- non periodic
                   1 --- periodic
  */
  int i,j,k;
  double dr[3],r,dtmp,rs,gs,rrs,r2,r1,r6;
  //  double r2min=0.25,r2max=64.0;
  double r2min=MD_LJ_R2MIN,r2max=MD_LJ_R2MAX;

  if((periodicflag & 1)==0){
    xmax*=2.0;
  }
  if(tblno==5){ r2min=0.01;r2max=1e6;}
  for(i=0;i<ni;i++){
    for(j=0;j<nj;j++){
      r2=0.0;
      for(k=0;k<3;k++){
        dr[k]=xi[i*3+k]-xj[j*3+k];
        if(dr[k]<-xmax/2.0){
          dr[k]+=xmax;
        }
        if(dr[k]>=xmax/2.0){
          dr[k]-=xmax;
        }
        r2+=dr[k]*dr[k];
      }
      if(r2!=0.0){
        r=sqrt(r2);
        rs=rscale[atypei[i]*nat+atypej[j]];
        gs=gscale[atypei[i]*nat+atypej[j]];
        rrs=r2*rs;
        if(rrs>=r2min && rrs<r2max){
          r1=1.0/rrs;
          r6=r1*r1*r1;
          if(tblno==2){
            dtmp=gs*r6*r1*(2.0*r6-1.0);
            for(k=0;k<3;k++){
              force[i*3+k]+=dtmp*dr[k];
            }
          }
          else if(tblno==3){
            dtmp=gs*r6*(r6-1.0);
            for(k=0;k<3;k++){
              force[i*3+k]+=dtmp;
            }
          }
          else if(tblno==5){
            dtmp=-gs*r6;
            for(k=0;k<3;k++){
              force[i*3+k]+=dtmp*dr[k];
            }
          }
        }
      }
    }
  }
}


void MR3calcewald_dft(int k[], int knum, double x[], int n,
                      double chg[], double cellsize[3],
                      double bs[], double bc[])
{
    int i,i3,j,j3,c;
    double th,cellsize_1[3];

    for(i=0;i<3;i++) cellsize_1[i]=1.0/cellsize[i];
    for(i=i3=0;i<knum;i++,i3+=3){
        bs[i]=bc[i]=0.0;
        for(j=j3=0;j<n;j++,j3+=3){
            th=0.0;
            for(c=0;c<3;c++) th+=k[i3+c]*x[j3+c]*cellsize_1[c];
            th*=2.0*M_PI;
            bs[i]+=chg[j]*sin(th);
            bc[i]+=chg[j]*cos(th);
        }
    }
}

void MR3calcewald_idft_eng(int k[], double bs[], double bc[],
                           int knum, double x[], int n,
                           double cellsize[3],
                           double force[])
{
  int i,i3,j,j3,c;
  double th,sth,cth;
  double cellsize_1[3];
  double fstmp,fctmp;
  
  for(i=0;i<3;i++) cellsize_1[i]=1.0/cellsize[i];
  for(i=i3=0;i<n;i++,i3+=3){
    fstmp=fctmp=0.0;
    for(j=j3=0;j<knum;j++,j3+=3){
      th=0.0;
      for(c=0;c<3;c++) th+=k[j3+c]*x[i3+c]*cellsize_1[c];
      th*=2.0*M_PI;
      sth=sin(th);
      cth=cos(th);
      fstmp+=bs[j]*sth;
      fctmp+=bc[j]*cth;
    }
    force[i3]+=fstmp+fctmp;
  }
}


void MR3calcewald_idft_force(int k[], double bs[], double bc[],
                             int knum, double x[], int n,
                             double cellsize[3],
                             double force[])
{
  int i,i3,j,j3,c;
  double th,sth,cth,cellsize_1[3];
  double fst[3],fct[3];
  double fstmp[3],fctmp[3];

  for(i=0;i<3;i++) cellsize_1[i]=1.0/cellsize[i];
  for(i=i3=0;i<n;i++,i3+=3){
    for(c=0;c<3;c++) fstmp[c]=fctmp[c]=0.0;
    for(j=j3=0;j<knum;j++,j3+=3){
      th=0.0;
      for(c=0;c<3;c++) th+=k[j3+c]*x[i3+c]*cellsize_1[c];
      th*=2.0*M_PI;
      sth=sin(th);
      cth=cos(th);
      for(c=0;c<3;c++){
        fst[c]=bc[j]*sth*k[j3+c];
        fct[c]=bs[j]*cth*k[j3+c];
      }
      for(c=0;c<3;c++){
        fstmp[c]+=fst[c];
        fctmp[c]+=fct[c];
      }


    }
    for(c=0;c<3;c++){
      force[i3+c]+=fstmp[c]-fctmp[c];
    }
  }
}


void MR3calcewald(int *k, int knum_org, double *x, int n, double *chg,
                  double alpha, double epsilon, double cell[3][3],
                  double *force, double *tpot, double(*)[3]) {
  double *bs,*bc,cellsize[3],cellsize_1[3];
  int knum,i,j,c,i3,j3;
  double factor1_tmp,vol1,eps1,alpha4,r2,kvtmp;

  knum=knum_org<0 ? -knum_org:knum_org;
  if((bs=(double *)malloc(sizeof(double)*knum))==NULL){
    fprintf(stderr,"** error : can't malloc bs **\n");
    exit(1);
  }
  if((bc=(double *)malloc(sizeof(double)*knum))==NULL){
    fprintf(stderr,"** error : can't malloc bc **\n");
    exit(1);
  }
  for(i=0;i<3;i++) cellsize[i]=cell[i][i];
  for(i=0;i<3;i++) cellsize_1[i]=1.0/cellsize[i];
  for(i=0,vol1=1.0;i<3;i++) vol1*=cellsize_1[i];
  eps1=1.0/epsilon;
  alpha4=1.0/(4.0*alpha*alpha);
  for(i=0;i<n*3;i++) force[i]=0.0;

  /* DFT */
  MR3calcewald_dft(k,knum,x,n,chg,cellsize,bs,bc);

  /* multiply factor etc */  
  *tpot=0.0;
  for(j=j3=0;j<knum;j++,j3+=3){
    for(c=0,r2=0.0;c<3;c++){
      kvtmp=2.0*M_PI*k[j3+c]*cellsize_1[c];
      r2+=kvtmp*kvtmp;
    }
    factor1_tmp=2.0*eps1*vol1*exp(-r2*alpha4)/r2;
    *tpot+=0.5*factor1_tmp*(bs[j]*bs[j]+bc[j]*bc[j]);
    bs[j]*=factor1_tmp;
    bc[j]*=factor1_tmp;
  }

  /* IDFT */
  if(knum_org<0){ /* potential energy */
    MR3calcewald_idft_eng(k,bs,bc,knum,x,n,cellsize,force);
    for(i=i3=0;i<n;i++,i3+=3) force[i3]*=chg[i];
  }
  else{           /* force */
    MR3calcewald_idft_force(k,bs,bc,knum,x,n,cellsize,force);
    for(i=i3=0;i<n;i++,i3+=3){
      for(c=0;c<3;c++) force[i3+c]*=2.0*M_PI*chg[i]*cellsize_1[c];
    }
  }

  free(bs);
  free(bc);
}

int get_knum(double ksize) {
  double kmaxsq = ksize * ksize;
  int kmax = ksize;
  int knum = 0;    
  for( int l=0; l<=kmax; l++ ) {
    int mmin = -kmax;
    if( l==0 ) mmin = 0;
    for( int m=mmin; m<=kmax; m++ ) {
      int nmin = -kmax;
      if( l==0 && m==0 ) nmin=1;
      for( int n=nmin; n<=kmax; n++ ) {
        double ksq = l * l + m * m + n * n;
        if( ksq <= kmaxsq ) {
          knum++;
        }
      }
    }
  }
  return knum;
}

void init_kvec(double ksize, int *kvec) {
  double kmaxsq = ksize * ksize;
  int kmax = ksize;
  int knum = 0;
  for( int l=0; l<=kmax; l++ ) {
    int mmin = -kmax;
    if( l==0 ) mmin = 0;
    for( int m=mmin; m<=kmax; m++ ) {
      int nmin = -kmax;
      if( l==0 && m==0 ) nmin=1;
      for( int n=nmin; n<=kmax; n++ ) {
        double ksq = l * l + m * m + n * n;
        if( ksq <= kmaxsq ) {
          kvec[3*knum+0] = l;
          kvec[3*knum+1] = m;
          kvec[3*knum+2] = n;
          knum++;
        }
      }
    }
  }
}

#endif
