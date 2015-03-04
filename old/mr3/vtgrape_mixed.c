static void sort_natex2(int ni, int numex2[], int natex2[])
{
  int i,j,offset,k,flag;

  // sort natex2
  for(i=offset=flag=0;i<ni;i++){
    if(numex2[i]>=2){
      for(j=0;j<numex2[i]-1;j++){
	for(k=j+1;k<numex2[i];k++){
	  if(natex2[offset+j]>natex2[offset+k]){ // swap
	    int d;
	    d=natex2[offset+k];
	    natex2[offset+k]=natex2[offset+j];
	    natex2[offset+j]=d;
	  }
	}
      }
    }
    offset+=numex2[i];
  }

  // check order of natex2
  for(i=offset=flag=0;i<ni;i++){
    if(numex2[i]>0){
      //      printf("numex2[%d]=%d natex2=%d",i,numex2[i],natex2[offset]);
      k=natex2[offset];
      for(j=1;j<numex2[i];j++){
	//	printf(" %d",natex2[offset+j]);
	if(natex2[offset+j]<k){
	  //	  printf("order is not correct: i=%d j=%d natex2[%d+%d]=%d\n",i,j,offset,j,natex2[offset+j]);
	  flag=1;
	}
	else{
	  k=natex2[offset+j];
	}
      }
      //      printf("\n");
    }
    offset+=numex2[i];
  }
  if(flag!=0){
    fprintf(stderr,"** order of natex2 is not good **\n");
    exit(1);
  }
}


static void calc_ij_host_mixedaccuracy_setup(float *largeval, int method)
{
  if(method=='N' || method=='P') *largeval=LARGE_VAL;
}


static void calc_ij_host_mixedaccuracy_clearforce(float forcef[3], float forcef2[3],
						  int forcei[3], float largeval)
{
  forcef[0]=forcef[1]=forcef[2]=largeval;
  forcef2[0]=forcef2[1]=forcef2[2]=0.0;
  forcei[0]=forcei[1]=forcei[2]=0;
}


static void calc_ij_host_mixedaccuracy_updateforce(float forcef[3], float forcef2[3],
						   int forcei[3], double forced[3],
						   int method)
{
  int k;

  if(method=='P' && 1){
    for(k=0;k<3;k++){
#if 1 // only single precision
#define LOWER_VAL_SHIFT2 (LOWER_VAL_SHIFT)
      int itmp;
      unsigned int ui;
      itmp=forcei[k] & (MASK(LOWER_VAL_SHIFT2)<<(32-LOWER_VAL_SHIFT2));
      //      if(itmp>0) printf("    k=%d forcei=%08x itmp=%08x ",k,forcei[k],itmp);
      forcei[k]&=MASK(32-LOWER_VAL_SHIFT2);
#if 1
      //      if(itmp>0) printf("forcei_new=%08x forcef=%e(%08x) ",forcei[k],forcef[k],ui=((unsigned int *)forcef)[k]);
      forcef[k]+=((float)itmp)*((float)LOWER_VAL_FACTOR_1)*1.0f;
      //      if(itmp>0) printf("forcef_new=%e(%08x) itmp*LOWER_VAL_FACTOR_1=%e\n",forcef[k],ui=((unsigned int *)forcef)[k],((float)itmp)*((float)LOWER_VAL_FACTOR_1)*1.0f);
#else
      {
	float fnew;
	ui=((unsigned int *)forcef)[k];
	if(itmp>0) printf("forcei_new=%08x forcef=%e(%08x) ",forcei[k],(double)(forcef[k]),ui);
	fnew=forcef[k]+((float)itmp)*((float)LOWER_VAL_FACTOR_1)*1.0f;
	ui=*((unsigned int *)&fnew);
	if(itmp>0) printf("forcef_new=%e(%08x) itmp*LOWER_VAL_FACTOR_1=%e\n",(double)fnew,ui,((float)itmp)*((float)LOWER_VAL_FACTOR_1)*1.0f);
	forcef[k]=fnew;
      }	
#endif
#else // use double precision
      forced[k]+=forcei[k]*LOWER_VAL_FACTOR_1;
      forcei[k]=0;
#endif
    }
  }
}


static void calc_ij_host_mixedaccuracy_addforce(float forcef[3], float forcef2[3],
						int forcei[3], double forced[3],
						float largeval, char method)
{
  int k;

  for(k=0;k<3;k++){
    if(method=='f'){
      forced[k]+=forcef[k];
    }
    else if(method=='i' || method=='n' || method=='N'){
      forced[k]+=(double)forcef[k]-(double)largeval+(double)forcef2[k];
      //        printf("forcef2[%d]=%e\n",k,forcef2[k]);
    }
    else if(method=='P'){
      forced[k]+=(double)forcef[k]-(double)largeval+(forcei[k]*LOWER_VAL_FACTOR_1);
    }
  }
}
							   

static void calc_ij_host_mixedaccuracy_calcforce(double xid[3], double xjd[3],
						 double qid, double qjd,
						 int atypei, int atypej, 
						 float forcef[3], float forcef2[3],
						 int forcei[3], double forced[3],
						 int nat,
						 double *gscaled, double *rscaled,
						 int tblno, double xmaxd,
						 double factord, 
						 float largeval, char method)
{
  float xi[3],xj[3],xmax=xmaxd,factor=factord;
  float dr[3],r,dtmp,rs,gs,rrs,qi,qj;
  int k;

  if(method=='D'){
    double ddr[3],dr,ddtmp,drrs;
    dr=0.0;
    for(k=0;k<3;k++){
      ddr[k]=xid[k]-xjd[k];
      if(ddr[k]<-xmaxd/2.0){
	ddr[k]+=xmaxd;
      }
      if(ddr[k]>=xmaxd/2.0){
	ddr[k]-=xmaxd;
      }
      dr+=ddr[k]*ddr[k];
    }
    dr=sqrt(dr);
    if(dr!=0.0){
      if(tblno==2 || tblno==3){
	drrs=dr*dr*rscaled[atypei*nat+atypej];
	if(tblno==2){
	  double drrs6;
	  drrs=1.0f/drrs;
	  drrs6=drrs*drrs*drrs;
	  ddtmp=gscaled[atypei*nat+atypej]*drrs6*drrs*(2.0*drrs6-1.0);
	  ddtmp*=factord;
	  for(k=0;k<3;k++) ddr[k]*=ddtmp;
	}
	else if(tblno==3){
	  ddtmp=gscaled[atypei*nat+atypej]*(pow(drrs,-6.0)-pow(drrs,-3.0));
	  ddtmp*=factord;
	  for(k=0;k<3;k++) ddr[k]=ddtmp;
	}
      }
      else if(tblno==0 || tblno==1 || tblno==6 || tblno==7){
	qi=qid;
	qj=qjd;
	if(tblno==0){
	  drrs=1.0/dr;
	  ddtmp=qjd*drrs*drrs*drrs;
	  ddtmp*=factord;
	  for(k=0;k<3;k++) ddr[k]*=ddtmp;
	}
	else if(tblno==1){
	  drrs=1.0/dr;
	  ddtmp=qjd*drrs;
	  ddtmp*=factord;
	  for(k=0;k<3;k++) ddr[k]=ddtmp;
	}
      }
      for(k=0;k<3;k++) forced[k]+=ddr[k];
    }
    return;
  }
  r=0.0;
  for(k=0;k<3;k++){
    xi[k]=xid[k];
    xj[k]=xjd[k];
    dr[k]=xi[k]-xj[k];
    if(dr[k]<-xmax/2.0){
      dr[k]+=xmax;
    }
    if(dr[k]>=xmax/2.0){
      dr[k]-=xmax;
    }
    r+=dr[k]*dr[k];
  }
  r=sqrtf(r);
  if(r!=0.0){
    if(tblno==2 || tblno==3){
      rs=rscaled[atypei*nat+atypej];
      gs=gscaled[atypei*nat+atypej];
      rrs=r*r*rs;
      if(tblno==2){
	float rrs6;
	rrs=1.0f/rrs;
	rrs6=rrs*rrs*rrs;
	dtmp=gs*rrs6*rrs*(2.0f*rrs6-1.0f);
	dtmp*=factor;
	for(k=0;k<3;k++) dr[k]*=dtmp;
      }
      else if(tblno==3){
	dtmp=gs*(powf(rrs,-6.0f)-powf(rrs,-3.0f));
	dtmp*=factor;
	for(k=0;k<3;k++) dr[k]=dtmp;
      }
    }
    else if(tblno==0 || tblno==1 || tblno==6 || tblno==7){
      qi=qid;
      qj=qjd;
      if(tblno==0){
	rrs=1.0f/r;
	dtmp=qj*rrs*rrs*rrs;
	dtmp*=factor;
	for(k=0;k<3;k++) dr[k]*=dtmp;
      }
      else if(tblno==1){
	rrs=1.0f/r;
	dtmp=qj*rrs;
	dtmp*=factor;
	for(k=0;k<3;k++) dr[k]=dtmp;
      }
    }
    if(method=='F'){
      for(k=0;k<3;k++) forced[k]+=dr[k];
    }
    else if(method=='f'){
      for(k=0;k<3;k++){
	forcef[k]+=dr[k];
      }
    }
    else if(method=='i'){
      float z,w,z1,v,z2,zz;
      for(k=0;k<3;k++){         // ah->forcef[k], al->forcef2[k]
	z=forcef[k]+dr[k];      // x->forcef[k], y->dr[k], z->z, zz->zz
	w=z-forcef[k];
	z1=dr[k]-w;
	v=z-w;
	z2=v-forcef[k];
	zz=z1-z2;               // add higher part with 6 ops.
	forcef2[k]+=zz;         // add lower part with 1 op.
	forcef[k]=z+forcef2[k]; // x->z, y->forcef2[k], z->forcef[k] zz->forcef2[k]
	w=forcef[k]-z;
	forcef2[k]=forcef2[k]-w;// balance higher and lower part with 3 ops.
      }
    }
    else if(method=='n' || method=='N'){
      float z,w,zz;
      for(k=0;k<3;k++){
	z=forcef[k]+dr[k];      // x->forcef[k], y->dr[k], z->z, zz->zz
	w=forcef[k]-z;
	zz=w+dr[k];
	forcef2[k]+=zz;         // 4 ops.
	forcef[k]=z;
      }
    }
    else if(method=='P'){
      float z,w,zz;
      for(k=0;k<3;k++){
	z=forcef[k]+dr[k];      // x->forcef[k], y->dr[k], z->z, zz->zz
	w=forcef[k]-z;
	zz=w+dr[k];
	forcef[k]=z;
	forcei[k]+=(int)(zz*((float)LOWER_VAL_FACTOR));
#if 0
	if((j % (1<<LOWER_VAL_SHIFT))==(1<<LOWER_VAL_SHIFT)-1){
	  forced[k]+=forcei[k]*LOWER_VAL_FACTOR_1;
	  forcei[k]=0;
	}
#endif
      }
    }
  }
}


/* 
     method        D -- all double precision (not included in this routine)
                   f -- all single precision
                   F -- all but accumulation is single precision
                   i -- Iitaka's method: most accurate (10 operation/add)
                   n -- Nitadori's method (4 operation/add)
                   N -- Nitadori's method and always accumulated value is larger (4 op/add)
                   P -- proposed method (5 op/add)
 */


void MR3calc_ij_host_mixedaccuracy(int ni, double xid[], double *qid, int *atypei, 
				   double forced[],
				   int nj, double xjd[], double *qjd, int *atypej,
				   int nat, double *gscaled, double *rscaled,
				   int tblno, double xmaxd, int periodicflag,
				   char method)
{
  /* periodic flag bit 0 : 0 --- non periodic
                           1 --- periodic
     method ------------ see upper comments
  */
  int i,j,k,ati,atj;
  int forcei[3];
  float forcef[3],forcef2[3],largeval=0.0f;
  double qi,qj;

  if((periodicflag & 6)!=0){
    fprintf(stderr,"** error : bit 1 or 2 of periodicflag is not supported in MR3calc_ij_host_mixedaccuracy **\n");
    vg_exit(1);
  }
  //  printf("in calc_ij_host_mixedaccracy\n");
  //  printf("ni=%d nat=%d xid[0]=%e %e %e tblno=%d xmaxd=%e periodicflag=%d method=%c\n",ni,nat,xid[0],xid[1],xid[2],tblno,xmaxd,periodicflag,method);
  if((periodicflag & 1)==0){
    xmaxd*=2.0;
  }
  calc_ij_host_mixedaccuracy_setup(&largeval,method);
  for(i=0;i<ni;i++){
    calc_ij_host_mixedaccuracy_clearforce(forcef,forcef2,forcei,largeval);
    /*    forcef[0]=forcef[1]=forcef[2]=largeval;
    forcef2[0]=forcef2[1]=forcef2[2]=0.0;
    forcei[0]=forcei[1]=forcei[2]=0;*/
    for(j=0;j<nj;j++){
      if(tblno==2 || tblno==3){	ati=atypei[i];atj=atypej[j];qi=qj=0.0;}
      else if(tblno==0||tblno==1||tblno==6||tblno==7){ati=atj=0;qi=qid[i];qj=qjd[j];}
      calc_ij_host_mixedaccuracy_calcforce(xid+i*3,xjd+j*3,qi,qj,ati,atj,
					   forcef,forcef2,forcei,forced+i*3,
					   nat,gscaled,rscaled,tblno,xmaxd,
					   1.0,largeval,method);
      if((j % (1<<LOWER_VAL_SHIFT))==(1<<LOWER_VAL_SHIFT)-1){
	calc_ij_host_mixedaccuracy_updateforce(forcef,forcef2,forcei,forced+i*3,
					       method);
      }
    }
    calc_ij_host_mixedaccuracy_addforce(forcef,forcef2,forcei,forced+i*3,
					largeval,method);
  }
}


void MR3calc_ij_exlist_host_mixedaccuracy(int ni, double xid[], double *qid, 
					  int *atypei, double forced[],
					  int nj, double xjd[], double *qjd, int *atypej,
					  int nat, double *gscaled, double *rscaled,
					  int tblno, double xmaxd, int periodicflag,
					  int numex[], int natex[],
					  char method)
{
  /* periodicflag bit 0: 0 --- non periodic
                         1 --- periodic
                  bit 1: 0 --- natex does not include duplicate list
                         1 --- natex includes duplicate list
                 (bit 2: 0 --- qi is not multiplied to force)
                 (       1 --- qi is multiplied to force    )
  */
  int *numex2,*natex2,*natex_offset,*natex2_offset;
  int forcei[3];
  float forcef[3],forcef2[3],largeval=0.0f;
  int i,j,k,iexcl,jj,jskip,ati,atj;
  double qi,qj,*ftmp,*fmalloc;

  //  printf("in MR3calc_ij_exlist_host_mixedaccuracy: tblno=%d periodicflag=%d\n",tblno,periodicflag);
  if((periodicflag & 4)!=0){
    fprintf(stderr,"** error : bit 2 of periodicflag is not supported in MR3calc_ij_exlist_host_mixedaccuracy **\n");
    vg_exit(1);
  }
  if((periodicflag & 2)==0){  // convert natex to natex2
    MR3_make_natex2_from_natex(ni,natex,numex,
			       &natex_offset,&natex2,&numex2,&natex2_offset,1);
  }
  else{
    numex2=numex;
    natex2=natex;
  }
  sort_natex2(ni,numex2,natex2);

  if((periodicflag & 4)!=0){
    if((fmalloc=(double *)MR3_malloc_pointer(sizeof(double)*ni*3,"MR3calc_ij_exlist_host_mixedaccuracy"))==NULL){
      fprintf(stderr,"** error : can't malloc ftmp in MR3calc_ij_exlist_host_mixedaccuracy **\n");
      vg_exit(1);
    }
    ftmp=fmalloc;
  }
  else{
    ftmp=forced;
  }
  if((periodicflag & 1)==0){
    xmaxd*=2.0;
  }
  calc_ij_host_mixedaccuracy_setup(&largeval,method);
#if 0 // calc N^2 and subtract exlist
  for(i=iexcl=0;i<ni;i++){
    calc_ij_host_mixedaccuracy_clearforce(forcef,forcef2,forcei,largeval);
    /*    forcef[0]=forcef[1]=forcef[2]=largeval;
    forcef2[0]=forcef2[1]=forcef2[2]=0.0;
    forcei[0]=forcei[1]=forcei[2]=0;*/
    for(j=0;j<nj;j++){
      if(tblno==2 || tblno==3){ati=atypei[i];atj=atypej[j];qi=qj=0.0;}
      else if(tblno==0||tblno==1||tblno==6||tblno==7){ati=atj=0;qi=qid[i];qj=qjd[j];}
      calc_ij_host_mixedaccuracy_calcforce(xid+i*3,xjd+j*3,qi,qj,ati,atj,
					   forcef,forcef2,forcei,ftmp+i*3,
					   nat,gscaled,rscaled,tblno,xmaxd,
					   1.0,largeval,method);
      if((j % (1<<LOWER_VAL_SHIFT))==(1<<LOWER_VAL_SHIFT)-1){
	calc_ij_host_mixedaccuracy_updateforce(forcef,forcef2,forcei,ftmp+i*3,
					       method);
      }
    }

    for(jj=0;jj<numex2[i];jj++){
      j=natex2[iexcl+jj];
      if(j>=0){
	if(tblno==2 || tblno==3){ati=atypei[i];atj=atypej[j];qi=qj=0.0;}
	else if(tblno==0||tblno==1||tblno==6||tblno==7){ati=atj=0;qi=qid[i];qj=qjd[j];}
	calc_ij_host_mixedaccuracy_calcforce(xid+i*3,xjd+j*3,qi,qj,ati,atj,
					     forcef,forcef2,forcei,ftmp+i*3,
					     nat,gscaled,rscaled,tblno,xmaxd,
					     -1.0,largeval,method);
      }
    }
    calc_ij_host_mixedaccuracy_addforce(forcef,forcef2,forcei,ftmp+i*3,
					largeval,method);
    iexcl+=numex2[i];
  }
#else
  for(i=iexcl=0;i<ni;i++){
    calc_ij_host_mixedaccuracy_clearforce(forcef,forcef2,forcei,largeval);
    /*    forcef[0]=forcef[1]=forcef[2]=largeval;
    forcef2[0]=forcef2[1]=forcef2[2]=0.0;
    forcei[0]=forcei[1]=forcei[2]=0;*/
    if(numex2[i]>0) jskip=natex2[iexcl];
    else            jskip=-1;
    for(j=jj=0;j<nj;j++){
      if(j==jskip){
	if(jj++<numex2[i]-1) jskip=natex2[iexcl+jj];
	else                 jskip=-1;
      }
      else{
	if(tblno==2 || tblno==3){ati=atypei[i];atj=atypej[j];qi=qj=0.0;}
	else if(tblno==0||tblno==1||tblno==6||tblno==7){ati=atj=0;qi=qid[i];qj=qjd[j];}
	calc_ij_host_mixedaccuracy_calcforce(xid+i*3,xjd+j*3,qi,qj,ati,atj,
					     forcef,forcef2,forcei,ftmp+i*3,
					     nat,gscaled,rscaled,tblno,xmaxd,
					     1.0,largeval,method);
      }
      if((j % (1<<LOWER_VAL_SHIFT))==(1<<LOWER_VAL_SHIFT)-1){
	calc_ij_host_mixedaccuracy_updateforce(forcef,forcef2,forcei,ftmp+i*3,
					       method);
      }
    }
    calc_ij_host_mixedaccuracy_addforce(forcef,forcef2,forcei,ftmp+i*3,
					largeval,method);
    iexcl+=numex2[i];
  }
#endif

  if((periodicflag & 2)==0){  // fr,""ee
    MR3_free_pointer(natex_offset,"MR3calc_ij_exlist_host_mixedaccuracy");
    MR3_free_pointer(natex2,"MR3calc_ij_exlist_host_mixedaccuracy");
    MR3_free_pointer(numex2,"MR3calc_ij_exlist_host_mixedaccuracy");
    MR3_free_pointer(natex2_offset,"MR3calc_ij_exlist_host_mixedaccuracy");
  }
  if((periodicflag & 4)!=0){
    for(i=0;i<ni;i++){
      for(j=0;j<3;j++){
	forced[i*3+j]+=ftmp[i*3+j]*qid[i];
      }
    }
    //    for(i=0;i<3;i++) printf("ftmp[%d]=%e %e %e forced=%e %e %e\n",i,ftmp[i*3],ftmp[i*3+1],ftmp[i*3+2],forced[i*3],forced[i*3+1],forced[i*3+2]);
    MR3_free_pointer(fmalloc,"MR3calc_ij_exlist_host_mixedaccuracy");
  }
}


void MR3calc_ij_nlist_host_mixedaccuracy(int ni, double xid[], double *qid, int *atypei, 
					 double forced[],
					 int nj, double xjd[], double *qjd, int *atypej,
					 int nat, double *gscaled, double *rscaled,
					 int tblno, double xmaxd, int periodicflag,
					 int numex[], int natex[],
					 char method)
{
  /* periodicflag bit 0: 0 --- non periodic
                         1 --- periodic
                  bit 1: 0 --- natex does not include duplicate list
                         1 --- natex includes duplicate list
  */
  int *numex2,*natex2,*natex_offset,*natex2_offset;
  int forcei[3];
  float forcef[3],forcef2[3],largeval=0.0f;
  int i,j,k,iexcl,jj,jskip,ati,atj;
  double qi,qj;

  if((periodicflag & 4)!=0){
    fprintf(stderr,"** error : bit 2 of periodicflag is not supported in MR3calc_ij_nlist_host_mixedaccuracy **\n");
    vg_exit(1);
  }
  if((periodicflag & 2)==0){  // convert natex to natex2
    MR3_make_natex2_from_natex(ni,natex,numex,
			       &natex_offset,&natex2,&numex2,&natex2_offset,1);
  }
  else{
    numex2=numex;
    natex2=natex;
  }
  sort_natex2(ni,numex2,natex2);

  if((periodicflag & 1)==0){
    xmaxd*=2.0;
  }
  calc_ij_host_mixedaccuracy_setup(&largeval,method);
  for(i=iexcl=0;i<ni;i++){
    calc_ij_host_mixedaccuracy_clearforce(forcef,forcef2,forcei,largeval);
    /*    forcef[0]=forcef[1]=forcef[2]=largeval;
    forcef2[0]=forcef2[1]=forcef2[2]=0.0;
    forcei[0]=forcei[1]=forcei[2]=0;*/
    for(jj=0;jj<numex2[i];jj++){
      j=natex2[iexcl+jj];
      if(j>=0){
	if(tblno==2 || tblno==3){ati=atypei[i];atj=atypej[j];qi=qj=0.0;}
	else if(tblno==0||tblno==1||tblno==6||tblno==7){ati=atj=0;qi=qid[i];qj=qjd[j];}
	calc_ij_host_mixedaccuracy_calcforce(xid+i*3,xjd+j*3,qi,qj,ati,atj,
					     forcef,forcef2,forcei,forced+i*3,
					     nat,gscaled,rscaled,tblno,xmaxd,
					     -1.0,largeval,method);
      }
      if((jj % (1<<LOWER_VAL_SHIFT))==(1<<LOWER_VAL_SHIFT)-1){
	calc_ij_host_mixedaccuracy_updateforce(forcef,forcef2,forcei,forced+i*3,
					       method);
      }
    }
    calc_ij_host_mixedaccuracy_addforce(forcef,forcef2,forcei,forced+i*3,
					largeval,method);
    iexcl+=numex2[i];
  }
}


// Flowwing program is for HPC Asia 2009 paper
/*
#include <stdlib.h>
#include <stdio.h>
#include "rettest.h"


//--------------------------------------------------------------------
typedef struct {
    float hs;
    float ls;
} SS;


SS add_knuth_and_dekker1(float xs, float ys)
{
    float ws;
    SS z;

    z.hs = xs   + ys;
    ws   = z.hs - xs;
    z.ls = ys   - ws;

    return z;
}
//--------------------------------------------------------------------


SS add_knuth_and_dekker2(float xs, float ys)
{
    float ws, vs, z1s, z2s;
    SS z;

    z.hs = xs   + ys;
    ws   = z.hs - xs;
    z1s  = ys   - ws;
    vs   = z.hs - ws;
    z2s  = vs   - xs;
    z.ls = z1s  - z2s;

    return z;
}
//--------------------------------------------------------------------


SS add_nitadori(SS a, float ys)
{
    SS b;

    b    = add_knuth_and_dekker1(a.hs, ys);
    b.ls = a.ls + b.ls;

    return b;
}
//--------------------------------------------------------------------


#if 1 // change to 0 for modified Nitadori's method
#define LARGE_VAL        0.0f
#else
#define LARGE_VAL_SHIFT  21
#define LARGE_VAL  \
        (float)(3 << (LARGE_VAL_SHIFT-1))
#endif

void accumulate_nitadori(float ri[3], float rj[][3], 
                         int nj, double fd[3])
{
  int j, k;
  SS a[3];
  float fs[3];

  for(k=0; k<3; k++){
    a[k].hs = LARGE_VAL;
    a[k].ls = 0.0f;
  }
  for(j=0; j<nj; j++){
    pairwiseforce(ri, rj[j], fs);
    for(k=0; k<3; k++){
      a[k] = add_nitadori(a[k], fs[k]);
    }
  }
  for(k=0;k<3;k++){
    fd[k] = a[k].hs - LARGE_VAL;
    fd[k] = fd[k] + a[k].ls;
  }
}
//--------------------------------------------------------------------


void pairwiseforce(float ris[3], float rjs[3], 
                   float fs[3])
{
  int k;
  float ds[3], rs, qs, ps;

  for(k=0; k<3; k++) ds[k] = ris[k] - rjs[k];
  rs  = ds[0] * ds[0] + ds[1] * ds[1] + ds[2] * ds[2];
  if(rs != 0.0f){
    rs  = 1.0f / rs;
    qs  = rs * rs * rs;
    ps  = rs * qs * (2.0f * qs - 1.0f);
    for(k=0; k<3; k++) fs[k] = ps * ds[k];
  }
  else{
    for(k=0; k<3; k++) fs[k] = 0.0f;
  }
}

void calcvdw(int ni, int nj, float ri[][3], 
	     float rj[][3], double fd[][3])
{
  int i;

  for(i=0; i<ni; i++){
    accumulate_nitadori(ri[i], rj, nj, fd[i]);
    // accumulate_narumi(ri[i], rj, nj, fd[i]);
  }
}
//--------------------------------------------------------------------


typedef struct {
    float hs;
    int   li;
} SI;

SI add_narumi(SI a, float ys)
{
  SS b;
  SI c;
  
  b    = add_knuth_and_dekker1(a.hs, ys);
  c.hs = b.hs;
  c.li = a.li+(int)(b.ls * LOWER_VAL_FACTOR);
  
  return c;
}
 
#define LARGE_VAL_SHIFT  21
#define LARGE_VAL  \
        (float)(3 << (LARGE_VAL_SHIFT-1))
#define LOWER_VAL_SHIFT  7
#define LOWER_VAL_LOOP   ( 1 << LOWER_VAL_SHIFT )
#define LOWER_VAL_FACTOR   \
        (float)( 1LL << ( 23 - LARGE_VAL_SHIFT \
                        + 32 - LOWER_VAL_SHIFT ) )
#define LOWER_VAL_FACTOR_1 \
        ( 1.0f / LOWER_VAL_FACTOR ) 
#define MASK(n)          ((0x1<<(n)) -1)

void accumulate_narumi(float ri[3], float rj[][3],
                       int nj, double fd[3])
{
  int j, jj, k;
  SI a[3];
  float fs[3];

  for(k=0; k<3; k++){
    a[k].hs = LARGE_VAL;
    a[k].li = 0;
  }
  for(jj=0; jj<nj; jj+=LOWER_VAL_LOOP){
    for(j=jj;j-jj<LOWER_VAL_LOOP && j<nj;j++){
      pairwiseforce(ri, rj[j], fs);
      for(k=0; k<3; k++){
	a[k] = add_narumi(a[k], fs[k]);
      }
    }
    for(k=0; k<3; k++){
      a[k].hs = a[k].hs + (float)(a[k].li & 
        (MASK(LOWER_VAL_SHIFT)<<(32-LOWER_VAL_SHIFT)))
        *LOWER_VAL_FACTOR_1;
      a[k].li = a[k].li & MASK(32-LOWER_VAL_SHIFT);
    }
  }
  for(k=0; k<3; k++){
    fd[k] = a[k].hs - LARGE_VAL;
    fd[k] = fd[k] + a[k].li * (double)LOWER_VAL_FACTOR_1;
  }
}
//--------------------------------------------------------------------


int main(int argc, char **argv)
{
  float ri[10][3],rj[10][3];
  int i,k,ni,nj;
  double fd[10][3];

  ni=nj=10;
  for(i=0;i<ni;i++) for(k=0;k<3;k++) ri[i][k]=drand48();
  for(i=0;i<nj;i++) for(k=0;k<3;k++) rj[i][k]=drand48();
  calcvdw(ni,nj,ri,rj,fd);
  for(i=0;i<ni;i++){
    printf("fd[%d]=%e %e %e\n",i,fd[i][0],fd[i][1],fd[i][2]);
  }
}

*/
