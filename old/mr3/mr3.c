#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "matrixMul.h"

#define MD_LJ_R2MIN 0.25f
#define MD_LJ_R2MAX 64.0f

#define COULOMB_VDW_FACTOR
static double Coulomb_vdw_factor=1.0;

#define MR3_malloc_pointer(size,string) malloc(size)
#define MR3_free_pointer(p,string)      free(p)

static double Eps2=0.0;
int Deviceid=0;


void vg_exit(int ret)
{
  exit(ret);
}


void *m3_get_unit(void)
{
  return NULL;
}


void m3_set_softening(void *unit, double eps)
{
  Eps2=eps*eps;
}


void MR3_exit(int ret)
{
  exit(ret);
}


void MR3init(void)
{
  char *s;

  s=getenv("VG_DEVICEID");
  if(s!=NULL){
    sscanf(s,"%d",&Deviceid);
    printf("VG_DEVICEID is set %d\n",Deviceid);
  }
}


void mr3init_(void)
{
  MR3init();
}


void mr3init__(void)
{
  mr3init_();
}


void MR3free(void)
{
}


void mr3free_(void)
{
  MR3free();
}


void mr3free__(void)
{
  mr3free_();
}


#if 0
void MR3SetTable(char *filename, int tblno, int flag)
{
}
#endif


void MR3calccoulomb_ij(int ni, double xi[], double qi[], double force[],
		       int nj, double xj[], double qj[],
		       double rscale, 
		       int tblno, double xmax, int periodicflag)
{
  int i,j,natchangeflag=0;

  switch(tblno){
  case 0:
    gpucoulombforce_ij_(xi,ni,qj,xj,nj,rscale,tblno,xmax,
			periodicflag,natchangeflag,force,Eps2);
    if((periodicflag & 2)!=0){
      for(i=0;i<ni;i++) for(j=0;j<3;j++) force[i*3+j]*=qi[i];
    }
    break;
  default:
    fprintf(stderr,"** error : not supported tblno = %d **\n",tblno);
    MR3_exit(1);
    break;
  }
}


void mr3calccoulomb_ij_(int *ni, double xi[], double qi[], double force[],
			int *nj, double xj[], double qj[],
			double *rscale, 
			int *tblno, double *xmax, int *periodicflag)
{
  MR3calccoulomb_ij(*ni,xi,qi,force,
		    *nj,xj,qj,
		    *rscale, 
		    *tblno,*xmax,*periodicflag);
}


void mr3calccoulomb_ij__(int *ni, double xi[], double qi[], double force[],
			 int *nj, double xj[], double qj[],
			 double *rscale, 
			 int *tblno, double *xmax, int *periodicflag)
{
  mr3calccoulomb_ij_(ni,xi,qi,force,
		     nj,xj,qj,
		     rscale, 
		     tblno,xmax,periodicflag);
}


void MR3calccoulomb(double x[], int n, double q[], double rscale,
		    int tblno, double xmax, int periodicflag,
		    int natchangeflag, double force[])
{
  int i,j;

  switch(tblno){
  case 0:
    gpucoulombforce_(x,n,q,rscale,tblno,xmax,
		     periodicflag,natchangeflag,force,Eps2);
    if((periodicflag & 2)!=0){
      for(i=0;i<n;i++) for(j=0;j<3;j++) force[i*3+j]*=q[i];
    }
    break;
  case 1:
    gpucoulombpot_(x,n,q,rscale,tblno,xmax,
		   periodicflag,natchangeflag,force);
    if((periodicflag & 2)!=0){
      for(i=0;i<n;i++) force[i*3]*=q[i];
    }
    break;
  case 6:
    gpuewreforce_(x,n,q,rscale,tblno,xmax,
		  periodicflag,natchangeflag,force);
    break;
  case 7:
    gpuewrepot_(x,n,q,rscale,tblno,xmax,
		periodicflag,natchangeflag,force);
    break;
  default:
    fprintf(stderr,"** error : not supported tblno = %d **\n",tblno);
    MR3_exit(1);
    break;
  }
}


void mr3calccoulomb_(double x[], int *n, double q[], double *rscale,
		     int *tblno, double *xmax, int *periodicflag,
		     int *natchangeflag, double force[])
{
  MR3calccoulomb(x,*n,q,*rscale,*tblno,*xmax,*periodicflag,
		 *natchangeflag,force);
}


void mr3calccoulomb__(double x[], int *n, double q[], double *rscale,
		      int *tblno, double *xmax, int *periodicflag,
		      int *natchangeflag, double force[])
{
  mr3calccoulomb_(x,n,q,rscale,
		  tblno,xmax,periodicflag,
		  natchangeflag,force);
}


void MR3calcvdw_ij(int ni, double xi[], int atypei[], double force[],
		   int nj, double xj[], int atypej[],
		   int nat, double gscale[], double rscale[],
		   int tblno, double xmax, int periodicflag)
{
#if 0
  MR3calcvdw_ij_host(ni,xi,atypei,force,
		     nj,xj,atypej,
		     nat,gscale,rscale,
		     tblno,xmax,periodicflag);
#endif
}


void MR3calcvdw(double x[], int n, int atype[], int nat,
		double gscale[], double rscale[],
		int tblno, double xmax, int periodicflag,
		int natchangeflag, double force[])
{
  //  double r2mind=0.25,r2maxd=64.0;
  double r2mind=MD_LJ_R2MIN,r2maxd=MD_LJ_R2MAX;
  int i;

#if 0
  for(i=0;i<n;i++) x[i*3]=x[i*3+1]=x[i*3+2]=0.0;
  x[0]=x[5*3]=1.0;
  printf("tblno=%d nat=%d xmax=%e\n",tblno,nat,xmax);
  for(i=0;i<10;i++){
    printf("pos[%d]=%e %e %e\n",i,x[i*3],x[i*3+1],x[i*3+2]);
  }
#endif
  if((nat+1)*(nat+1)*2>THD){
    fprintf(stderr,"** too many atom types = %d **\n",nat);
    MR3_exit(1);
  }
  switch(tblno){
  case 2:
    gpuvdwforce_(x,n,atype,nat,gscale,rscale,tblno,xmax,
		 periodicflag,natchangeflag,force,r2mind,r2maxd);
    break;
  case 3:
#if 1
    gpuvdwpot_(x,n,atype,nat,gscale,rscale,tblno,xmax,
	       periodicflag,natchangeflag,force,r2mind,r2maxd);
#else
    printf("gpuvdwpot is canged\n");
    gpuvdwpot_(x,n,atype,1.0,3,xmax,periodicflag,natchangeflag,force);
#if 0
    {
      double *q;
      if((q=(double *)MR3_malloc_pointer(sizeof(double)*n,"MR3calcvdw"))==NULL){
	fprintf(stderr,"** error : can't malloc q **\n");
	MR3_exit(1);
      }
      for(i=0;i<n;i++) q[i]=0.0;
      gpuvdwpot_(x,n,q,1.0,3,xmax,periodicflag,natchangeflag,force);
      MR3_free_pointer(q,"MR3calcvdw");
    }
#endif
#endif
    break;
  default:
    fprintf(stderr,"** error : not supported tblno = %d **\n",tblno);
    MR3_exit(1);
    break;
  }
}


void mr3calcvdw_(double x[], int *n, int atype[], int *nat,
		 double gscale[], double rscale[],
		 int *tblno, double *xmax, int *periodicflag,
		 int *natchangeflag, double force[])
{
  int *atype2;
  int i;

  if((atype2=(int *)MR3_malloc_pointer(sizeof(int)*(*n),"atypei2 in mr3calcvdw_"))==NULL){
    fprintf(stderr,
	    "** error at malloc atype2 in mr3calcvdw_ **\n");
    MR3_exit(1);
  }
  for(i=0;i<*n;i++){
    if(atype[i]>0){
      atype2[i]=atype[i]-1;
    }
    else{
      printf("  warning : atype[%d]=%d should be positive\n",
	     i,atype[i]);
      atype2[i]=0;
    }
  }
  MR3calcvdw(x,*n,atype2,*nat,
	     gscale,rscale,
	     *tblno,*xmax,*periodicflag,
	     *natchangeflag,force);
  MR3_free_pointer(atype2,"atype2 in mr3calcvdw_");
}


void mr3calcvdw__(double x[], int *n, int atype[], int *nat,
		  double gscale[], double rscale[],
		  int *tblno, double *xmax, int *periodicflag,
		  int *natchangeflag, double force[])
{
  mr3calcvdw_(x,n,atype,nat,
	      gscale,rscale,
	      tblno,xmax,periodicflag,
	      natchangeflag,force);
}


void MR3calcewald(int *k, int knum_org, double *x, int n, double *chg,
		  double alpha, double epsilon, double cell[3][3],
		  double *force, double *tpot, double stress[3][3])
{
  double *gvec,cell_1[3],*xq,vol1,eps1,alpha4,kvtmp,r2;
  int knum,i,j;
  
  knum=(knum_org>0 ? knum_org:-knum_org);
  if((gvec=(double *)MR3_malloc_pointer(sizeof(double)*knum*4,"gvec in MR3calcewald"))==NULL){
    fprintf(stderr,"** can't malloc gvec in MR3calcewald **\n");
    MR3_exit(1);
  }
  if((xq=(double *)MR3_malloc_pointer(sizeof(double)*n*4,"xq in MR3calcewald"))==NULL){
    fprintf(stderr,"** can't malloc xq in MR3calcewald **\n");
    MR3_exit(1);
  }
  for(i=0;i<3;i++) cell_1[i]=1.0/cell[i][i];
  for(i=0;i<n;i++){
    for(j=0;j<3;j++) xq[i*4+j]=x[i*3+j];
    xq[i*4+3]=chg[i];
  }
  for(i=0,vol1=1.0;i<3;i++) vol1*=cell_1[i];
  eps1=1.0/epsilon;
  alpha4=1.0/(4.0*alpha*alpha);
  for(i=0;i<knum;i++){
    r2=0.0;
    for(j=0;j<3;j++){
      kvtmp=gvec[i*4+j]=k[i*3+j]*2.0*M_PI*cell_1[j];
      r2+=kvtmp*kvtmp;
    }
    gvec[i*4+3]=2.0*eps1*vol1*exp(-r2*alpha4)/r2;
  }
  if(knum_org>0){
    gpuewwvforce_(xq,n,gvec,knum_org,force,tpot);
  }
  else if(knum_org<0){
    gpuewwvpot_(xq,n,gvec,-knum_org,force,tpot);
    *tpot=0.0;
    for(i=0;i<n;i++) *tpot+=force[i*3];
    *tpot*=0.5;
  }
  MR3_free_pointer(gvec,"gvec in MR3calcewald");
  MR3_free_pointer(xq,"xq in MR3calcewald");
}


void mr3calcewald_(int *k, int *knum, double *x, int *n,
		   double *chg, double *alpha, double *epsilon,
		   double cell[3][3], double *force,
		   double *tpot, double stress[3][3])
{
    MR3calcewald(k,*knum,x,*n,chg,*alpha,*epsilon,cell,
		 force,tpot,stress);
}


void mr3calcewald__(int *k, int *knum, double *x, int *n,
		    double *chg, double *alpha, double *epsilon,
		    double cell[3][3], double *force,
		    double *tpot, double stress[3][3])
{
    mr3calcewald_(k,knum,x,n,chg,alpha,epsilon,cell,
		  force,tpot,stress);
}


void MR3calcnacl(void)
{
  fprintf(stderr,"** error : MR3calcnacl is not supported **\n");
}


void vtgrape_force(double xj[][3], // position of j-th particles
		   double mj[],    // mass of j-th particles
		   double xi[][3], // position of i-th particles
		   double eps2,    // softening parameter
		   double a[][3], // force of i-th particles
		   int ni,         // number of i-th particles
		   int nj)
{
  double size=10.0;
  static int ini=0;

  if(ini==0){
    MR3init();
    ini=1;
  }
  gpucoulombforce_ij_(xi,ni,mj,xj,nj,1.0,0,size,
		      0,0,a,eps2);
}
#include "mr3_host.c"


void vg_get_r1result(int n, float r[], float r1[])
{
  printf("** warning : this routine is not supported **\n");
}
