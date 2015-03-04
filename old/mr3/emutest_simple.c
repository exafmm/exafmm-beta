#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>
#include <string.h>
#include "vtgrape.h"
#include "md.h"
#include "vtgrapeproto.h"

#define N (1<<23)
//#define N (1<<22)

#if defined(CUDA_SDK_2) && 0
int Argc=0;
char **Argv=NULL;
#endif


void get_cputime(double *laptime, double *sprittime)
{
  struct timeval tv;
  struct timezone tz;
  double sec,microsec;

  gettimeofday(&tv, &tz);
  sec=tv.tv_sec;
  microsec=tv.tv_usec;

  *sprittime = sec + microsec * 1e-6 - *laptime;
  *laptime = sec + microsec * 1e-6;
  //  printf("    vg_get_cputime is called: ltime=%e stime=%e\n",*laptime,*sprittime);
}


int main(int argc, char **argv)
{
  if(argc<2){
    printf("usage: %s option1 option2 ...\n",argv[0]);
    printf("               option1: ni(.nj)\n");
    printf("               option2: number of loops\n");
    printf("               option3: size of a cell. default is 100.0\n");
    printf("               option4: function_type\n");
    printf("                        0 -- Coulomb force\n");
    printf("                        1 -- Coulomb potential\n");
    printf("                        2 -- van der Waals force\n");
    printf("                             number of atom types can be specified by 2.nat\n");
    printf("                             default number of atoms is 16\n");
    printf("                        3 -- van der Waals potential\n");
    printf("                             number of atom types can be specified by 3.nat\n");
    printf("               option5: periodic(1) or non-periodic(0). default is non-periodic\n");
    printf("     ex. %s 1000.5000 100 200.0 2 1\n",argv[0]);
  }
  else{
    double *xi,*xj,*qi,*qj,*fi,*fh;
    int *atypei,*atypej;
    double gscalesf[64*64],gscalesp[64*64],rscales[64*64],rscale=0.3;
    int nat=16,ni,nj=0;
    int i,j,tblno,periodicflag=0,loop=1;
    double size=100.0,s,emax,eavr,favr,e,ltime,stime;
    double nfi[40]={38,38,32,32,0,0,55,55,0,0,
		    0,0,0,0,0,0,0,0,0,0,
		    0,0,0,0,0,0,0,0,0,0,
		    0,0,0,0,0,0,32+55,32+55,0,0};
    double *gscalesfp=NULL;
    double sum=0.0,cellindexfactor=1.0;
    ni=nj=0;
    
    if(argc>=2){
      sscanf(argv[1],"%d.%d",&ni,&nj);
      if(nj==0) nj=ni;
      printf("ni=%d nj=%d\n",ni,nj);
    }
    if(argc>=3){
      sscanf(argv[2],"%d",&loop);
    }
    printf("loop=%d\n",loop);
    if(argc>=4){
      sscanf(argv[3],"%lf",&size);
    }
    printf("size is %f\n",size);
    if(argc>=5){
      sscanf(argv[4],"%d.%d",&tblno,&nat);
      if(nat>64){
	fprintf(stderr,"** error : too many atom types **\n");
	exit(1);
      }
    }
    printf("tblno=%d, nat=%d\n",tblno,nat);
    if(argc>=6){
      sscanf(argv[5],"%d",&periodicflag);
    }
    printf("periodicflag=%d\n",periodicflag);
    if((xi=(double *)MR3_malloc_pointer(sizeof(double)*ni*3,"main"))==NULL){
      fprintf(stderr,"** error : can't malloc xi **\n");
      exit(1);
    }
    if((fi=(double *)MR3_malloc_pointer(sizeof(double)*ni*3,"main"))==NULL){
      fprintf(stderr,"** error : can't malloc fi **\n");
      exit(1);
    }
    if((fh=(double *)MR3_malloc_pointer(sizeof(double)*ni*3,"main"))==NULL){
      fprintf(stderr,"** error : can't malloc fh **\n");
      exit(1);
    }
    if((qi=(double *)MR3_malloc_pointer(sizeof(double)*ni,"main"))==NULL){
      fprintf(stderr,"** error : can't malloc qi **\n");
      exit(1);
    }
    if((atypei=(int *)MR3_malloc_pointer(sizeof(int)*ni,"main"))==NULL){
      fprintf(stderr,"** error : can't malloc atypei **\n");
      exit(1);
    }
    if((xj=(double *)MR3_malloc_pointer(sizeof(double)*nj*3,"main"))==NULL){
      fprintf(stderr,"** error : can't malloc xj **\n");
      exit(1);
    }
    if((qj=(double *)MR3_malloc_pointer(sizeof(double)*nj,"main"))==NULL){
      fprintf(stderr,"** error : can't malloc qj **\n");
      exit(1);
    }
    if((atypej=(int *)MR3_malloc_pointer(sizeof(int)*nj,"main"))==NULL){
      fprintf(stderr,"** error : can't malloc atypej **\n");
      exit(1);
    }
    for(i=0;i<ni;i++){
      for(j=0;j<3;j++){
	xi[i*3+j]=drand48()*size;
	fi[i*3+j]=fh[i*3+j]=0.0;
      }
      qi[i]=drand48()*2.0-1.0;
      atypei[i]=drand48()*nat;
    }
    for(i=0;i<nj;i++){
      for(j=0;j<3;j++){
	xj[i*3+j]=drand48()*size;
      }
      qj[i]=drand48()*2.0-1.0;
      atypej[i]=drand48()*nat;
    }
    for(i=0;i<nat;i++){ // only diagonal part is activated
      gscalesf[i*nat+i]=drand48();
      gscalesp[i*nat+i]=drand48();
      rscales[i*nat+i]=drand48();
    }
    for(i=0;i<nat;i++) for(j=0;j<nat;j++) if(i!=j){
      gscalesf[i*nat+j]=sqrt(gscalesf[i*nat+i]*gscalesf[j*nat+j]);
      gscalesp[i*nat+j]=sqrt(gscalesp[i*nat+i]*gscalesp[j*nat+j]);
      rscales[i*nat+j]=(sqrt(rscales[i*nat+i])+sqrt(rscales[j*nat+j]))*0.5;
      rscales[i*nat+j]*=rscales[i*nat+j];
    }
    
    MR3init();

    // GPU calculation
    get_cputime(&ltime,&stime);
    for(i=0;i<loop;i++){
      bzero(fi,sizeof(double)*ni*3);
      switch(tblno){
      case 0:
      case 1:
      case 6:
      case 7:
	MR3calccoulomb_ij(ni,xi,qi,fi,nj,xj,qj,rscale,tblno,size,periodicflag);
	break;
      case 2:
	gscalesfp=gscalesf;
      case 3:
	if(gscalesfp==NULL) gscalesfp=gscalesp;
	MR3calcvdw_ij(ni,xi,atypei,fi,nj,xj,atypej,nat,gscalesfp,rscales,
		      tblno,size,periodicflag);
	break;
      default:
	fprintf(stderr,"** error : not supported tblno=%d **\n",tblno);
	exit(1);
	break;
      }
    }
    get_cputime(&ltime,&stime);
      
    // host calculation
    switch(tblno){
      case 0:
      case 1:
      case 6:
      case 7:
	MR3calccoulomb_ij_host(ni,xi,qi,fh,nj,xj,qj,rscale,tblno,size,periodicflag);
	break;
      case 2:
      case 3:
	MR3calcvdw_ij_host(ni,xi,atypei,fh,nj,xj,atypej,nat,gscalesfp,rscales,tblno,size,periodicflag);
	break;
    }
    
    favr=emax=0.0;
    for(i=0;i<ni;i++){
      for(j=0,s=0.0;j<((tblno % 2)==0 ? 3:1);j++) s+=fh[i*3+j]*fh[i*3+j];
      favr+=sqrt(s);
    }
    favr/=ni;
    for(i=0;i<ni;i++){
      for(j=0,s=0.0;j<((tblno % 2)==0 ? 3:1);j++){
	s+=(fi[i*3+j]-fh[i*3+j])*(fi[i*3+j]-fh[i*3+j]);
      }
      s=sqrt(s);
      e=s/favr;
      if(e>1e-4){
	printf("i=%d e=%e fi=%20.12e %20.12e %20.12e fh=%20.12e %20.12e %20.12e\n",i,e,fi[i*3],fi[i*3+1],fi[i*3+2],fh[i*3],fh[i*3+1],fh[i*3+2]);
      }
      eavr+=e;
      if(e>emax) emax=e;
    }
    eavr/=ni;
    printf("average_force=%e, max_relative_err=%e avr_relative_err=%e\n",
	   favr,emax,eavr);
    for(i=0,sum=0.0;i<ni;i++){
      if((tblno % 2)==0){ // force
	for(j=0;j<3;j++) sum+=fi[i*3+j]*fi[i*3+j];
      }
      else{ // potential
	sum+=fi[i*3];
      }
    }
    printf("sum of force or potential=%20.17e\n",sum);
    printf("ni %d nj=%d loop=%d time/step=%e nfi=%d effective_speed= %f Gflops\n",
	   ni,nj,loop,stime/loop,(int)nfi[tblno],
	   cellindexfactor*(double)ni*((double)nj)*nfi[tblno]*((double)loop)/stime/1e9);
    MR3free();
  }
   
  return 0;
}
