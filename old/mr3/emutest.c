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
#if defined(CUDA_SDK_2) && 0
  Argc=argc;
  Argv=argv;
#endif
  if(argc<2){
    printf("usage: %s option1 option2 ...\n",argv[0]);
    printf("option1: e --- emulator test\n");
    printf("         E --- emulator test for multiplication\n");
    printf("         c,C,v,V\n");
    printf("           --- coulomb or vdw speed test\n");
    printf("               when 'C', calculation speed is measured with changing N\n");
    printf("               when 'v', emulator is used for host calculation\n");
    printf("               when 'V', emulator(nlist) is used for host calculation\n");
    printf("               option2: ni.nj\n");
    printf("               option3: number of loops\n");
    printf("                        if minus, host calculation is skipped\n");
    printf("                        loops.size\n");
    printf("                        size can be specified. default is 100.0\n");
    printf("               option4: function_type\n");
    printf("                        0 -- Coulomb force\n");
    printf("                        1 -- Coulomb potential\n");
    printf("                        2 -- van der Waals force\n");
    printf("                        3 -- van der Waals potential\n");
    printf("                        number of cells per dimension can be specified\n");
    printf("                        by like, function_type.number_of_cells_per_dimension.\n");
    printf("                        then cellindex routine is used.\n");
    printf("               ex. %s C 65536 2 0 |grep Gflops|awk '{print $2,$8}'\n",argv[0]);
    printf("                   %s c 16384 1.4 2\n",argv[0]);
    //    printf("-- max_err:1.3e-4, avr_err:1.1e-7->1.2e-8(pair float)\n",argv[0]);
    printf("                      MD_SORT_ATYPEI, no ACCUMULATE_PAIR_FLOAT\n");
    printf("                      MD_PERIODIC_FIXED=1, MD_LJ_05SIGMA_8SIGMA, periodicflag=0\n");
    printf("                        max_relative_err=1.126121e-06 avr_relative_err=1.150693e-07\n");
    printf("                      ACCUMULATE_PAIR_FLOAT=2,\n");
    printf("                        max_relative_err=9.955265e-09 avr_relative_err=1.171249e-09\n");
    printf("                      ACCUMULATE_PAIR_FLOAT=1,\n");
    printf("                        max_relative_err=9.957317e-09 avr_relative_err=1.171236e-09\n");
    printf("                   %s c 16384 1.4 0\n",argv[0]);
    printf("                      similar to the above, and no ACCUMULATE_PAIR_FLOAT\n");
    printf("                        max_relative_err=9.317259e-07 avr_relative_err=2.115614e-09\n");
    printf("                      ACCUMULATE_PAIR_FLOAT=1,\n");
    printf("                        max_relative_err=6.081641e-07 avr_relative_err=1.959465e-10\n");
    printf("                      ACCUMULATE_PAIR_FLOAT=2,\n");
    printf("                        max_relative_err=9.317259e-07 avr_relative_err=2.115614e-09\n");
    printf("                   %s c 1000 1.100 36.3\n",argv[0]);
    printf("                      cellindex virial test for charmm\n");
    printf("          A,a -- accumulation test\n");
    printf("                 if 'A', emulator is used.\n");
    printf("               option1: function table no\n");
    printf("          x ---- exlist test\n");
    printf("               option1: number of particles. ex.1000.200\n");
    printf("                 ex.%s x 100\n",argv[0]);
    printf("          w ---- write emu file\n");
    printf("          r ---- generate emu table internaly\n");
    printf("          n ---- nlist emulator performance test\n");
    printf("               option1: number of particles\n");
    exit(1);
  }
  else{
    switch(argv[1][0]){
    case 'w':
      MR3init();
      vg_write_emufile(0);
      MR3free();
      break;
    case 'r':
      MR3init();
      vg_write_emufile(1);
      MR3free();
      break;
    case 'x':
    case 'n':
      {
	double *xi,*xj,*qi,*qj,*fi,*fh;
	int *atypei,*atypej,*numex,*natex,npa=3;
	double gscalesf[32*32],gscalesp[32*32],rscales[32*32],rscale=1.0;
	int i,j,k,l,tblno,periodicflag=2,loop=1,nmin,hostflag=1;
	int nat,ni,nj=0,potflag=1,changeflag=0,njtmp=0;
	double size=10.0,potc=0.0,potv=0.0,favr,emax,eavr,s,e;
	double potch=0.0,potvh=0.0,volume[3];

	if(argv[1][0]=='n'){
	  npa=10;
	  printf("npa=%d\n",npa);
	}
	nat=16;
	printf("periodicflag=%d\n",periodicflag);
	if(argc>=3){
	  sscanf(argv[2],"%d.%d.%d",&ni,&nj,&njtmp);
	  if(nj==0) nj=ni;
	  printf("ni=%d nj=%d nat=%d njtmp=%d\n",ni,nj,nat,njtmp);
	}
	else{
	  printf("number of atoms must be specified\n");
	  exit(1);
	}
	size=pow((double)nj,1.0/3.0)*2.0;
	printf("size=%e\n",size);
	for(i=0;i<3;i++) volume[i]=size;
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
	if((numex=(int *)MR3_malloc_pointer(sizeof(int)*ni,"main"))==NULL){
	  fprintf(stderr,"** error : can't malloc numex **\n");
	  exit(1);
	}
	if((natex=(int *)MR3_malloc_pointer(sizeof(int)*ni*npa,"main"))==NULL){
	  fprintf(stderr,"** error : can't malloc natex **\n");
	  exit(1);
	}
	for(i=0;i<ni;i++){
	  for(j=0;j<3;j++){
	    xi[i*3+j]=drand48()*size;
	    fi[i*3+j]=fh[i*3+j]=0.0;
	  }
	  qi[i]=drand48()*2.0-1.0;
	  //	  qi[i]=1.0;printf("qi is set 1\n");
	  atypei[i]=drand48()*nat;
	  //	  atypei[i]=0;/* ************************************** */
	}
	for(i=0;i<nj;i++){
	  for(j=0;j<3;j++){
	    xj[i*3+j]=drand48()*size;
	  }
	  qj[i]=drand48()*2.0-1.0;
	  //	  qj[i]=1.0;printf("qj is set 1\n");
	  atypej[i]=drand48()*nat;
	  //	  atypej[i]=1;/* ************************************** */
	}
	for(i=0;i<nat;i++){ // only diagonal part is activated
	  gscalesf[i*nat+i]=drand48()+1.0;
	  gscalesp[i*nat+i]=drand48()+1.0;
	  rscales[i*nat+i]=drand48()+1.0;
	  //	  gscalesf[i*nat+i]=gscalesp[i*nat+i]=rscales[i*nat+i]=1.0;printf("** gscales and rscales are set 1\n");
	  //	  rscales[i*nat+i]=1.0;printf("** rscales are set 1\n");
	}
	for(i=0;i<nat;i++) for(j=0;j<nat;j++) if(i!=j){
	  gscalesf[i*nat+j]=sqrt(gscalesf[i*nat+i]*gscalesf[j*nat+j]);
	  gscalesp[i*nat+j]=sqrt(gscalesp[i*nat+i]*gscalesp[j*nat+j]);
	  //	  rscales[i*nat+j]=(rscales[i*nat+i]+rscales[j*nat+j])*0.5;
	  rscales[i*nat+j]=(sqrt(rscales[i*nat+i])+sqrt(rscales[j*nat+j]))*0.5;
	  rscales[i*nat+j]*=rscales[i*nat+j];
	  //	  printf("g %e %e r %e\n",gscalesf[i*nat+j],gscalesp[i*nat+j],rscales[i*nat+j]);
	}
	// atom type 0 do not interact with other atoms
	//	for(i=0;i<nat;i++) gscalesf[i*nat]=gscalesp[i*nat]=gscalesf[i]=gscalesp[i]=0.0;

	// all atom type is set to zero
	//	for(i=0;i<ni;i++) atypei[i]=0;

	// charges are set to zero
	//	for(i=0;i<ni;i++) qi[i]=0.0;
	//	for(i=0;i<nj;i++) qj[i]=0.0;

	// charges are set to one
	for(i=0;i<ni;i++) qi[i]=1.0;
	for(i=0;i<nj;i++) qj[i]=1.0;

	for(i=k=0;i<ni;i++){
	  for(j=0;j<npa;j++){
	    int flag;
	    do{
	      flag=0;
	      natex[k+j]=drand48()*nj;
	      for(l=0;l<j-1;l++){
		if(natex[k+j]==natex[k+l]) flag=1;
	      }
	    }
	    while(flag==1);
	    //	    printf("i=%d exclude_index=%d\n",i,natex[k+j]);
	  }
	  numex[i]=npa;
	  //	  numex[i]=0;printf("** numex is set zero\n");
	  k+=npa;
	}
	//	for(i=0;i<ni;i++) printf("numex[%d]=%d\n",i,numex[i]);
	//	for(i=0;i<ni*npa;i++) printf("natex[%d]=%d\n",i,numex[i]);

	MR3init();

	if(argv[1][0]=='n'){
	  double ltime,stime;
	  int l,loop=100;
	  get_cputime(&ltime,&stime);
	  for(l=0;l<loop;l++){
	    MR3calccoulomb_vdw_nlist_ij_emu(ni,xi,qi,atypei,fi,
					    nj,xj,qj,atypej,
					    nat,gscalesf,rscales,
					    rscale,tblno,volume,periodicflag,
					    numex,natex,-1.0);
	  }
	  get_cputime(&ltime,&stime);
	  printf("stime=%f\n",stime);
	  exit(0);
	}

	for(l=0;l<=5;l++){
	  for(i=0;i<ni;i++) for(j=0;j<3;j++) fi[i*3+j]=fh[i*3+j]=0.0;
	  potc=potv=potch=potvh=0.0;
	  if(l==0 || l==2 || l==4) tblno=0;
	  else if(l==5)            tblno=6;
	  else                     tblno=1;
	  if(l==0){
	    MR3calccoulomb_ij(ni,xi,qi,fi,nj,xj,qj,rscale,tblno,size,(periodicflag & 1)+2);
	    MR3calccoulomb_ij_host(ni,xi,qi,fh,nj,xj,qj,rscale,tblno,size,(periodicflag & 1)+2);
	  }
	  else if(l==1){
	    MR3calccoulomb_ij(ni,xi,qi,fi,nj,xj,qj,rscale,tblno,size,(periodicflag & 1)+2);
	    MR3calccoulomb_ij_host(ni,xi,qi,fh,nj,xj,qj,rscale,tblno,size,(periodicflag & 1)+2);
	  }
	  else if(l==2){
	    MR3calcvdw_ij(ni,xi,atypei,fi,nj,xj,atypej,nat,
			  gscalesf,rscales,(tblno & 1)+2,size,(periodicflag & 1));
	    MR3calcvdw_ij_host(ni,xi,atypei,fh,nj,xj,atypej,nat,
			       gscalesf,rscales,(tblno & 1)+2,size,(periodicflag & 1));
	  }
	  else if(l==3){
	    int njj=nj;
	    if(njtmp!=0) njj=njtmp;
	    MR3calcvdw_ij(ni,xi,atypei,fi,njj,xj,atypej,nat,
			  gscalesf,rscales,(tblno & 1)+2,size,(periodicflag & 1));
	    MR3calcvdw_ij_host(ni,xi,atypei,fh,njj,xj,atypej,nat,
			       gscalesf,rscales,(tblno & 1)+2,size,(periodicflag & 1));
	  }
	  else if(l==4){// non-cutoff coulomb and cutoff vdw 
	    MR3calccoulomb_vdw_ij_exlist(ni,xi,qi,atypei,fi,nj,xj,qj,atypej,
					 nat,gscalesf,gscalesp,rscales,rscale,
					 tblno,size,&potc,&potv,(periodicflag & 3),
					 potflag,changeflag,numex,natex);
	    MR3calccoulomb_vdw_ij_exlist_host2(ni,xi,qi,atypei,fh,nj,xj,qj,atypej,nat,
					       gscalesf,gscalesp,rscales,rscale,tblno,
					       size,&potch,&potvh,periodicflag,potflag,
					       changeflag,numex,natex);
	    printf("potc=%e potch=%e rerr=%e potv=%e potvh=%e rerr=%e\n",
		   potc,potch,fabs((potc-potch)/potch),
		   potv,potvh,fabs((potv-potvh)/potvh));
	  }
	  else if(l==5){// ewald real-space and cutoff vdw (cutoff is used for host ewald real-space)
	    MR3calccoulomb_vdw_ij_exlist(ni,xi,qi,atypei,fi,nj,xj,qj,atypej,
					 nat,gscalesf,gscalesp,rscales,rscale,
					 tblno,size,&potc,&potv,(periodicflag & 3),
					 potflag,changeflag,numex,natex);
	    MR3calccoulomb_vdw_ij_exlist_host(ni,xi,qi,atypei,fh,nj,xj,qj,atypej,nat,
					      gscalesf,gscalesp,rscales,rscale,tblno,
					      size,&potch,&potvh,periodicflag,potflag,
					      changeflag,numex,natex);
	    printf("potc=%e potch=%e rerr=%e potv=%e potvh=%e rerr=%e\n",
		   potc,potch,fabs((potc-potch)/potch),
		   potv,potvh,fabs((potv-potvh)/potvh));
	  }
#if 1
	  {
	    double rerr;
	    //	    for(i=0;i<(ni<3 ? ni:3);i++){
	    for(i=0;i<ni;i++){
	      rerr=((tblno % 2)==0 ? 
		    sqrt(fabs(((fi[i*3]-fh[i*3])*(fi[i*3]-fh[i*3])+(fi[i*3+1]-fh[i*3+1])*(fi[i*3+1]-fh[i*3+1])+(fi[i*3+2]-fh[i*3+2])*(fi[i*3+2]-fh[i*3+2]))/
			      (fh[i*3]*fh[i*3]+fh[i*3+1]*fh[i*3+1]+fh[i*3+2]*fh[i*3+2])))
		    :sqrt(fabs((fi[i*3]-fh[i*3])*(fi[i*3]-fh[i*3])/
			       (fh[i*3]*fh[i*3])))
		    );
	      if(rerr>1e-4){
		//		if(1){
		printf("l=%d fi[%d]=%20.12e %20.12e %20.12e fh=%20.12e %20.12e %20.12e rerr=%e\n",l,i,fi[i*3],fi[i*3+1],fi[i*3+2],fh[i*3],fh[i*3+1],fh[i*3+2],rerr);
	      }
	    }
	  }
#endif

	  // compare with host
	    favr=emax=eavr=0.0;
	    for(i=0;i<ni;i++){
	      for(j=0,s=0.0;j<((tblno % 2)==0 ? 3:1);j++) s+=fh[i*3+j]*fh[i*3+j];
	      favr+=sqrt(s);
	    }
	    favr/=ni;
	    for(i=0;i<ni;i++){
	      for(j=0,s=0.0;j<((tblno % 2)==0 ? 3:1);j++) s+=(fi[i*3+j]-fh[i*3+j])*(fi[i*3+j]-fh[i*3+j]);
	      s=sqrt(s);
	      e=s/favr;
	      eavr+=e;
	      if(e>emax) emax=e;
	    }
	    eavr/=ni;
	    printf("l=%d average_force=%e, max_relative_err=%e avr_relative_err=%e\n",
		   l,favr,emax,eavr);
	}

	MR3free();
      }
      break;
    case 'e':
      {
	int i,n,count,errcount,l;
	float *r,*r1,*r1test,r1host,r1est;
	VG_UNION_FI fihost,figpu;

	//#if VG_MINIMUM_PARTICLE_BLOCK_I==128
	//	fprintf(stderr,"** VG_MINIMUM_PARTICLE_BLOCK_I=128 does't work, which is very strange\n");
	//#endif
	if((r=(float *)MR3_malloc_pointer(sizeof(float)*N*2,"main"))==NULL){
	  fprintf(stderr,"** error : can't malloc r **\n");
	  exit(1);
	}
	if((r1=(float *)MR3_malloc_pointer(sizeof(float)*N*2,"main"))==NULL){
	  fprintf(stderr,"** error : can't malloc r1 **\n");
	  exit(1);
	}
	if((r1test=(float *)MR3_malloc_pointer(sizeof(float)*N*2,"main"))==NULL){
	  fprintf(stderr,"** error : can't malloc r1test **\n");
	  exit(1);
	}
	//	for(i=0;i<N;i++) r[i]=1.0+i/(double)N;
	MR3init();
	//	vg_initialize_emu();
	for(l=0;l<2;l++){
	//	for(l=1;l<2;l++){
	  if(l==0) n=N;
	  else     n=N*2;
	  //	  else     n=N;
	  if(l==0) vg_get_r1result(n,r,r1);
	  else     vg_get_rsqrtresult(n,r,r1);
	  for(i=0;i<n;i++){
	    if(l==0) r1host=1.0f/r[i];
	    else     r1host=1.0/sqrt((double)r[i]);
	    //    printf("i=%d r=%f r1host=%f r1gpu=%f\n",i,r[i],r1host,r1[i]);
#if 0
	    if((i % 10000)==0){
	      if(r1[i]!=r1host){
		fihost.f=r1host;
		figpu.f=r1[i];
		printf("i=%d r=%f r1host=%f(%08x) r1gpu=%f(%08x)\n",
		       i,r[i],r1host,fihost.i,r1[i],figpu.i);
	      }
	    }
#endif
	  }
	  n=100000;
	  if(n>N) n=N;
	  for(i=0;i<n;i++) r[i]=drand48()*100.0;
	  if(l==0) vg_get_r1result(-n,r,r1test);
	  else     vg_get_rsqrtresult(-n,r,r1test);
	  for(i=count=errcount=0;i<n;i++){
	    float exp,man,one=1.0f;
	    unsigned int ui;
#if 1
	    if(l==0) r1est=vg_r1emu(r[i]);
	    else     r1est=vg_rsqrtemu(r[i]);
#else
	    if(i==0) printf("** vg_r1emu rsqrtemu is not used **\n");
	    if(l==0){
	      ui=(((unsigned int *)r)[i] & 0x7fffff) | *((unsigned int *)&one);
	      man=*((float *)&ui);
	      ui=(((unsigned int *)r)[i] & 0xff800000);
	      exp=*((float *)&ui);
	      ui=(((unsigned int *)r)[i] & 0x7fffff);
	      r1est=r1[ui]/exp;
	    }
	    else{
	      float r2;
	      ui=(((unsigned int *)r)[i] & 0x7fffff) | *((unsigned int *)&one);
	      man=*((float *)&ui);
	      ui=(((unsigned int *)r)[i] & 0xff800000);
	      exp=*((float *)&ui);
	      ui=(((unsigned int *)r)[i] & 0x7fffff);
	      if((*((unsigned int *)&exp)>>23) % 2==0){
		man*=2.0f;
		exp/=2.0f;
		ui+=(1<<23);
	      }
	      r1est=r1[ui]/sqrtf(exp);
	    }
#endif
	    if(r1est!=r1test[i]){
	      printf("r[%d]=%e man=%e(%08x) exp=%e adr=%x r1est=%e r1test=%e diff=%e\n",
		     i,r[i],man,*((unsigned int *)&man),exp,ui,r1est,r1test[i],r1est-r1test[i]);
	      errcount++;
	    }
	    if((i % 10000)==0) printf("i=%d r1est=%e r1test=%e\n",i,r1est,r1test[i]);
	    count++;
	  }
	  if(l==0) printf("r1: ");
	  else     printf("rsqrt: ");
	  printf("count=%d errcount=%d\n",count,errcount);
	}
	MR3free();
      }
      break;
    case 'E':
      {
	int i,n,count,errcount,l;
	float *r,*rq,*rqtest,rqhost,rqest,*q;
	VG_UNION_FI fihost,figpu;

	if((r=(float *)MR3_malloc_pointer(sizeof(float)*N*2,"main"))==NULL){
	  fprintf(stderr,"** error : can't malloc r **\n");
	  exit(1);
	}
	if((q=(float *)MR3_malloc_pointer(sizeof(float)*N*2,"main"))==NULL){
	  fprintf(stderr,"** error : can't malloc q **\n");
	  exit(1);
	}
	if((rq=(float *)MR3_malloc_pointer(sizeof(float)*N*2,"main"))==NULL){
	  fprintf(stderr,"** error : can't malloc rq **\n");
	  exit(1);
	}
	if((rqtest=(float *)MR3_malloc_pointer(sizeof(float)*N*2,"main"))==NULL){
	  fprintf(stderr,"** error : can't malloc rqtest **\n");
	  exit(1);
	}
	MR3init();
	n=100000;
	if(n>N) n=N;
	for(i=0;i<n;i++){
	  r[i]=drand48()*100.0;
	  q[i]=0.0f;
	  //	  q[i]=drand48();
	  //	  r[i]=2.0;
	  //	  q[i]=2.0;
	}
	for(l=0;l<10000;l++){
	  q[0]=drand48();
	  vg_get_mulresult(-n,r,q,rqtest);
	  for(i=count=errcount=0;i<n;i++){
	    double dq,dr,drqest;
	    unsigned long long ull;
#if 1
	    dr=r[i];
	    dr*=dr;
	    ull=*((unsigned long long *)&dr);
	    ull&=0xffffffffe0000000LL;
	    dr=*((double *)&ull);
	    rqest=q[0]+(float)dr;
#else
	    rqest=q[0]+r[i]*r[i];
#endif
	    if(rqest!=rqtest[i]){
	      unsigned int man,manq;
	      man=*((unsigned int *)&r[i]) & 0x7fffff;
	      manq=*((unsigned int *)&q[0]) & 0x7fffff;
	      printf("r[%d]=%e(%06x) q=%e(%06x) rqest=%e(%08x) rqtest=%e(%08x) diff=%e\n",
		     i,r[i],man,q[0],manq,rqest,*((unsigned int *)&rqest),
		     rqtest[i],*((unsigned int *)&rqtest[i]),rqest-rqtest[i]);
	      printf("  rqest=%08x rqtest=%08x\n",
		     (*((unsigned int *)&rqest) & 0x7fffff)<<1,
		     (*((unsigned int *)&rqtest[i]) & 0x7fffff)<<1);
	      dq=q[0];dr=r[i];dr*=dr;drqest=dq+dr;
	      printf("  dr*dr=%20.14e(%016llx) dq=%20.14e(%016llx) drqest=%20.14e(%016llx)\n",
		     dr,*((unsigned long long *)&dr),dq,*((unsigned long long *)&dq),
		     drqest,*((unsigned long long *)&drqest));
	      errcount++;
	      exit(1);
	    }
	    //	    if((i % 10000)==0) printf("i=%d rqest=%e rqtest=%e\n",i,rqest,rqtest[i]);
	    count++;
	  }
	  printf("rq: ");
	  printf("count=%d errcount=%d\n",count,errcount);
	}
	MR3free();
      }
      break;
    case 'A':
    case 'a':
      {
	double xi[512*3],xj[512*3],qi[512],qj[512],fi[512*3],size=99.0;
	double rscale=1.0,fh[512*3],gscales[1]={1.0},rscales[1]={1.0};
	int i,j,l,n=512,tblno=0,ni,nj,periodicflag=0;
	int atypei[512],atypej[512],nat=1,hostflag=1;
	if(argc>=3){
	  sscanf(argv[2],"%d",&tblno);
	  printf("tblno=%d\n",tblno);
	}
	if(argv[1][0]=='A'){
	  hostflag=2;
	  printf("emulator is used\n");
	}
	ni=nj=n;
	ni=3;
	for(i=0;i<n;i++){
	  for(j=0;j<3;j++) xi[i*3+j]=xj[i*3+j]=fi[i*3+j]=fh[i*3+j]=0.0;
	  qi[i]=qj[i]=0.0;
	  atypei[i]=atypej[i]=0;
	}
	if(tblno==2 || tblno==3){
	  xj[0]=0.5;xj[1*3]=2.0;xj[2*3]=2.5;xj[3*3]=3;
	  xi[1*3]=2.5;xi[2*3]=3.0;
	}
	else{
	  xj[0]=0.01;xj[1*3]=1.0;xj[2*3]=1.01;xj[3*3]=1.02;
	  xi[1*3]=1.01;xi[2*3]=2.0;
	}
	qi[0]=qi[1]=qi[2]=qj[0]=qj[1]=qj[2]=qj[3]=1.0;
	MR3init();
	for(l=0;l<5;l++){
	  printf("l=%d nj=%d\n",l,nj);
	  for(i=0;i<ni;i++) for(j=0;j<3;j++) fi[i*3+j]=fh[i*3+j]=0.0;
	  if(tblno==0 || tblno==1 || tblno==6 || tblno==7){
	    MR3calccoulomb_ij(ni,xi,qi,fi,nj,xj,qj,rscale,
			      tblno,size,periodicflag);
	    if(hostflag==1) MR3calccoulomb_ij_host(ni,xi,qi,fh,nj,xj,qj,rscale,
						   tblno,size,periodicflag);
	    else            MR3calccoulomb_ij_emu(ni,xi,qi,fh,nj,xj,qj,rscale,
						  tblno,size,periodicflag);
	  }
	  else if(tblno==2 || tblno==3){
	    MR3calcvdw_ij(ni,xi,atypei,fi,nj,xj,atypej,nat,gscales,rscales,
			  tblno,size,periodicflag);
	    if(hostflag==1) MR3calcvdw_ij_host(ni,xi,atypei,fh,nj,xj,atypej,nat,
					       gscales,rscales,tblno,size,periodicflag);
	    else            MR3calcvdw_ij_emu(ni,xi,atypei,fh,nj,xj,atypej,nat,
					      gscales,rscales,tblno,size,periodicflag);
	  }
	  for(i=0;i<3;i++){
	    printf("fi[%d]=%20.12f fh=%20.12f fi/fh=%20.12f\n",
		   i,fi[i*3],fh[i*3],fi[i*3]/fh[i*3]);
	  }
	  nj=4-l;
	}
	MR3free();
      }
      break;
    case 'c':
    case 'C':
    case 'v':
    case 'V':
      {
	double *xi,*xj,*qi,*qj,*fi,*fh;
	int *atypei,*atypej;
	double gscalesf[32*32],gscalesp[32*32],rscales[32*32],rscale=0.3;
	int nat,ni,nj=0;
	int i,j,k,tblno,periodicflag=1,loop=1,nmin,hostflag=1;
	double size=100.0,s,emax,eavr,favr,e,ltime,stime;
	double nfi[40]={38,38,32,32,0,0,55,55,0,0,
	                0,0,0,0,0,0,0,0,0,0,
	                0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,32+55,32+55,0,0};
	double *gscalesfp=NULL;
	double rcut=10.0,skinnb=0.0,volume[3];
	int ldim[3]={1,1,1};
	double sum=0.0,cellindexfactor=1.0;
	ni=nj=0;
	//	for(i=0;i<40;i++)printf("nif[%d]=%e\n",i,nfi[i]);
#if 1
	//	nat=4;
	nat=16;
#else
	nat=1;
	printf("**** nat is set %d ****\n",nat);
#endif
	if(argc>=3){
	  sscanf(argv[2],"%d.%d",&ni,&nj);
	  if(nj==0) nj=ni;
	  printf("ni=%d nj=%d nat=%d\n",ni,nj,nat);
	}
	else{
	  exit(1);
	}
	if(argv[1][0]=='v'){
	  hostflag=2;
	  printf("emulator is used\n");
	}
	else if(argv[1][0]=='V'){
	  hostflag=3;
	  printf("emulator (nlist) is used\n");
	}
	printf("periodicflag=%d\n",periodicflag);
	if(argc>=4){
	  sscanf(argv[3],"%d.%lf",&loop,&size);
	  if(loop<0){
	    loop=-loop;
	    hostflag=0;
	    printf("host calculation is skipped\n");
	  }
	  printf("loop=%d\n",loop);
	}
	printf("size is %f\n",size);
	if(argc>=5){
	  sscanf(argv[4],"%d.%d",&tblno,&ldim[0]);
	  printf("tblno=%d\n",tblno);
	  if(ldim[0]>1){
	    ldim[1]=ldim[2]=ldim[0];
	    cellindexfactor=ldim[0]*ldim[1]*ldim[2];
	    cellindexfactor=27.0/cellindexfactor;
	  }
	}
	if(argv[1][0]=='C'){
	  hostflag=0;
	  nj=ni;
	  nmin=512;
	}
	else{
	  nmin=ni;
	}
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
	  //	  atypei[i]=0;/* ************************************** */
	}
	for(i=0;i<nj;i++){
	  for(j=0;j<3;j++){
	    xj[i*3+j]=drand48()*size;
	  }
	  qj[i]=drand48()*2.0-1.0;
	  atypej[i]=drand48()*nat;
	  //	  atypej[i]=1;/* ************************************** */
	}
	for(i=0;i<nat;i++){ // only diagonal part is activated
	  gscalesf[i*nat+i]=drand48();
	  gscalesp[i*nat+i]=drand48();
	  rscales[i*nat+i]=drand48();
	  //	  gscalesf[i*nat+i]=gscalesp[i*nat+i]=rscales[i*nat+i]=1.0;printf("** gscales and rscales are set 1\n");
	  //	  rscales[i*nat+i]=1.0;printf("** rscales are set 1\n");
	}
	for(i=0;i<nat;i++) for(j=0;j<nat;j++) if(i!=j){
	  gscalesf[i*nat+j]=sqrt(gscalesf[i*nat+i]*gscalesf[j*nat+j]);
	  gscalesp[i*nat+j]=sqrt(gscalesp[i*nat+i]*gscalesp[j*nat+j]);
	  //	  rscales[i*nat+j]=(rscales[i*nat+i]+rscales[j*nat+j])*0.5;
	  rscales[i*nat+j]=(sqrt(rscales[i*nat+i])+sqrt(rscales[j*nat+j]))*0.5;
	  rscales[i*nat+j]*=rscales[i*nat+j];
	  //	  printf("g %e %e r %e\n",gscalesf[i*nat+j],gscalesp[i*nat+j],rscales[i*nat+j]);
	}
	MR3init();
	//	vg_initialize_emu();
	// dummy to eliminate initial latency
	if(ldim[0]>1){
	  printf("** warning: dummy routine to remove initial latency is skipped for cellindex **\n");
	}
	else{
	  MR3calccoulomb_ij(ni,xi,qi,fi,nj,xj,qj,rscale,0,size,periodicflag);
	}
	//	tblno=0;
	for(;ni>=nmin;ni/=2,nj/=2,loop*=3){
	  get_cputime(&ltime,&stime);
	  for(i=0;i<loop;i++){
	    bzero(fi,sizeof(double)*ni*3);
	    switch(tblno){
	    case 0:
	    case 1:
	    case 6:
	    case 7:
#if 1
	      if(ldim[0]>1){
		printf("MR3calccoulomb_ij_ci is used: periodicflag=%d, rcut=%f, skinnb=%f, ldim=(%d,%d,%d)\n",
		       periodicflag,rcut,skinnb,ldim[0],ldim[1],ldim[2]);
		volume[0]=volume[1]=volume[2]=size;
		MR3calccoulomb_ij_ci(ni,xi,qi,fi,nj,xj,qj,rscale,tblno,rcut,skinnb,volume,ldim);
	      }
	      else{
		MR3calccoulomb_ij(ni,xi,qi,fi,nj,xj,qj,rscale,tblno,size,periodicflag);
	      }
	      /*
#if 1 // use cell-index routine for debug
	      {
		double volume[3]={size,size,size},rcut=10.0,skinnb=0.1;
		int ldim[3]={3,3,3};
	      }
#else
	      MR3calccoulomb_ij(ni,xi,qi,fi,nj,xj,qj,rscale,
				tblno,size,periodicflag);
#endif
	      */
#else
	      //	  for(j=0;j<ni;j++) for(k=0;k<3;k++) fi[j*3+k]=0.0;
	      //	  bzero(fi,sizeof(double)*ni*3);
	      if(i!=0){
		MR3calccoulomb_ij(ni,xi,qi,fh,nj,xj,qj,rscale,tblno,size,periodicflag);
	      }
	      else{
		MR3calccoulomb_ij(ni,xi,qi,fi,nj,xj,qj,rscale,tblno,size,periodicflag);
	      }
	      for(j=0;j<ni;j++) for(k=0;k<3;k++) if(fi[j*3+k]!=fh[j*3+k]) printf("fh[%d][%d]=%e fi=%e\n",j,k,fh[j*3+k],fi[j*3+k]);
#endif
	      break;
	    case 2:
	      gscalesfp=gscalesf;
	    case 3:
	      if(gscalesfp==NULL) gscalesfp=gscalesp;
	      if(ldim[0]>1){
		printf("MR3calcvdw_ij_ci is used: periodicflag=%d, rcut=%f, skinnb=%f, ldim=(%d,%d,%d)\n",
		       periodicflag,rcut,skinnb,ldim[0],ldim[1],ldim[2]);
		volume[0]=volume[1]=volume[2]=size;
		MR3calcvdw_ij_ci(ni,xi,atypei,fi,nj,xj,atypej,nat,gscalesfp,rscales,tblno,rcut,skinnb,volume,ldim);
	      }
	      else{
		MR3calcvdw_ij(ni,xi,atypei,fi,nj,xj,atypej,nat,gscalesfp,rscales,
			      tblno,size,periodicflag);
	      }
	      break;
	    case 36:
	      gscalesfp=gscalesf;
	    case 37:
	    case 47:
	      if(gscalesfp==NULL) gscalesfp=gscalesp;
	      if(ldim[0]>1){
		printf("MR3calccoulomb_vdw_ij_ci is used: periodicflag=%d, rcut=%f, skinnb=%f, ldim=(%d,%d,%d)\n",
		       periodicflag,rcut,skinnb,ldim[0],ldim[1],ldim[2]);
		volume[0]=volume[1]=volume[2]=size;
#if 0 // 
		printf("** particles positions is explicitely specified.\n");
		xi[0]=3.0;xi[1]=15.0;xi[2]=15.0;
		xi[3]=28.0;xi[4]=15.0;xi[5]=15.0;
#endif
#if 1 // copy i particles to j-particles
		printf("** i-particles are copyed to j-particles.\n");
		if(ni!=nj){ fprintf(stderr,"** error : ni != nj **\n");exit(1);}
		for(i=0;i<ni;i++){
		  for(j=0;j<3;j++) xj[i*3+j]=xi[i*3+j];
		  qj[i]=qi[i];
		  atypej[i]=atypei[i];
		}
#endif
		MR3calccoulomb_vdw_ij_ci_old_09(ni,xi,qi,atypei,fi,nj,xj,qj,atypej,nat,gscalesfp,rscales,rscale,(tblno % 10),rcut,skinnb,volume,ldim);
#if 1
		{
		  int ni2max;
		  double *fi2,*xi2,*qi2,*fi3,stress[3][3];
		  int *atypei2,*index;
		  int cid[3],ni2=0,cxyz[3];
		  ni2max=ni*(ldim[0]+2.0)*(ldim[1]+2.0)*(ldim[2]+2.0)/(ldim[0]*ldim[1]*ldim[2])*2.0;
		  xi2=(double *)MR3_malloc_pointer(sizeof(double)*ni2max*3,"main");
		  fi2=(double *)MR3_malloc_pointer(sizeof(double)*ni2max*3,"main");
		  fi3=(double *)MR3_malloc_pointer(sizeof(double)*ni2max*3,"main");
		  qi2=(double *)MR3_malloc_pointer(sizeof(double)*ni2max,"main");
		  atypei2=(int *)MR3_malloc_pointer(sizeof(int)*ni2max,"main");
		  index=(int *)MR3_malloc_pointer(sizeof(int)*ni2max,"main");
		  for(i=0;i<ni;i++){
		    for(j=0;j<3;j++){
		      cid[j]=xi[i*3+j]/size*3.0;
		    }		    
		    for(cxyz[0]=(cid[0]==ldim[0]-1 ?-1:0);cxyz[0]<=(cid[0]==0 ? 1:0);cxyz[0]++){
		      for(cxyz[1]=(cid[1]==ldim[1]-1 ?-1:0);cxyz[1]<=(cid[1]==0 ? 1:0);cxyz[1]++){
			for(cxyz[2]=(cid[2]==ldim[2]-1 ?-1:0);cxyz[2]<=(cid[2]==0 ? 1:0);cxyz[2]++){
			  for(j=0;j<3;j++) xi2[ni2*3+j]=xi[i*3+j]+cxyz[j]*size;
			  qi2[ni2]=qi[i];
			  atypei2[ni2]=atypei[i];
			  index[ni2]=i;
			  //			  printf("i=%d ni2=%d cxyz=(%d,%d,%d) xi=(%e,%e,%e)(%d,%d,%d)\n",i,ni2,cxyz[0],cxyz[1],cxyz[2],xi2[ni2*3],xi2[ni2*3+1],xi2[ni2*3+2],cid[0],cid[1],cid[2]);
			  ni2++;
			}
		      }
		    }
		  }
		  if(ni2>ni2max){
		    fprintf(stderr,"** error : ni2 overflow **\n",ni2);
		    exit(1);
		  }
		  bzero(fi2,sizeof(double)*ni2*3);
		  bzero(fi3,sizeof(double)*ni2*3);
		  printf("ni2=%d\n",ni2);
#if 1
		  MR3calccoulomb_vdw_ij_ci(ni2,xi2,qi2,atypei2,fi3,nj,xj,qj,atypej,nat,gscalesfp,rscales,rscale,(tblno % 10),rcut,skinnb,volume,ldim);
		  //		  MR3calccoulomb_vdw_ij_ci(ni2,xi2,qi2,atypei2,fi3,nj,xj,qj,atypej,nat,gscalesfp,rscales,rscale,(tblno % 10),-rcut,skinnb,volume,ldim);
		  for(i=0;i<ni2;i++){
		    for(j=0;j<3;j++){
		      fi2[index[i]*3+j]+=fi3[i*3+j];
		      //		      if(index[i]==0 && j==2) printf("index[%d]=%d fi3=%e %e %e fi2=%e %e %e\n",i,index[i],fi3[i*3],fi3[i*3+1],fi3[i*3+2],fi2[index[i]*3],fi2[index[i]*3+1],fi2[index[i]*3+2]);
		    }
		  }
#else
		  MR3calccoulomb_vdw_ij_ci(ni,xi,qi,atypei,fi2,nj,xj,qj,atypej,nat,gscalesfp,rscales,rscale,(tblno % 10),-rcut,skinnb,volume,ldim);
#endif
		  for(i=0;i<ni;i++){
		    for(j=0;j<3;j++) cid[j]=xi[i*3+j]/size*3.0;
		    for(j=0;j<3;j++){
		      if(fi[i*3+j]!=fi2[i*3+j]){
			printf("i=%d j=%d fi=%24.17e fi2=%24.17e\n",i,j,fi[i*3+j],fi2[i*3+j]);
			printf("  cid=(%d,%d,%d)\n",cid[0],cid[1],cid[2]);
		      }
		      else{
			printf("i=%d j=%d agree (fi=%e)\n",i,j,fi[i*3+j]);
		      }
		    }
		  }
		}
#endif
	      }
	      else{
		fprintf(stderr,"** error : non-cellindex routine is not supported **\n");
		exit(1);
	      }
	      break;
	    default:
	      fprintf(stderr,"** error : not supported tblno=%d **\n",tblno);
	      exit(1);
	      break;
	    }
	  }
	  get_cputime(&ltime,&stime);
	  if(hostflag){
	    switch(tblno){
	    case 0:
	    case 1:
	    case 6:
	    case 7:
	      if(hostflag==1)      MR3calccoulomb_ij_host(ni,xi,qi,fh,nj,xj,qj,rscale,tblno,size,periodicflag);
	      else if(hostflag==2) MR3calccoulomb_ij_emu(ni,xi,qi,fh,nj,xj,qj,rscale,tblno,size,periodicflag);
	      else{
		int *natex,*numex;
		if((natex=(int *)MR3_malloc_pointer(sizeof(int)*ni*nj,"main"))==NULL){
		  fprintf(stderr,"** error : can't malloc natex **\n");
		  vg_exit(1);
		}
		if((numex=(int *)MR3_malloc_pointer(sizeof(int)*ni,"main"))==NULL){
		  fprintf(stderr,"** error : can't malloc numex **\n");
		  vg_exit(1);
		}
		for(i=0;i<ni;i++) for(j=0;j<nj;j++) natex[i*nj+j]=j;
		for(i=0;i<ni;i++) numex[i]=nj;
		MR3calccoulomb_nlist_ij_emu(ni,xi,qi,fh,nj,xj,qj,rscale,tblno,size,periodicflag+2+4,natex,numex,1.0);
		MR3_free_pointer(natex,"main");
		MR3_free_pointer(numex,"main");
	      }
	      break;
	    case 2:
	    case 3:
	      if(hostflag==1)      MR3calcvdw_ij_host(ni,xi,atypei,fh,nj,xj,atypej,nat,gscalesfp,rscales,tblno,size,periodicflag);
	      else if(hostflag==2) MR3calcvdw_ij_emu(ni,xi,atypei,fh,nj,xj,atypej,nat,gscalesfp,rscales,tblno,size,periodicflag);
	      else{
		int *natex,*numex;
		if((natex=(int *)MR3_malloc_pointer(sizeof(int)*ni*nj,"main"))==NULL){
		  fprintf(stderr,"** error : can't malloc natex **\n");
		  vg_exit(1);
		}
		if((numex=(int *)MR3_malloc_pointer(sizeof(int)*ni,"main"))==NULL){
		  fprintf(stderr,"** error : can't malloc numex **\n");
		  vg_exit(1);
		}
		for(i=0;i<ni;i++) for(j=0;j<nj;j++) natex[i*nj+j]=j;
		for(i=0;i<ni;i++) numex[i]=nj;
		MR3calcvdw_nlist_ij_emu(ni,xi,atypei,fh,nj,xj,atypej,nat,gscalesfp,rscales,tblno,size,periodicflag+2,numex,natex,1.0);
		MR3_free_pointer(natex,"main");
		MR3_free_pointer(numex,"main");
	      }
	      break;
	    }
#if 0
	    for(i=0;i<10;i++){
	      printf("fi=%20.12e %20.12e %20.12e fh=%20.12e %20.12e %20.12e\n",fi[i*3],fi[i*3+1],fi[i*3+2],fh[i*3],fh[i*3+1],fh[i*3+2]);
	    }
#endif
	    
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
	  }
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
	}
	MR3free();
      }
      break;
    default:
      fprintf(stderr,"** error : not supported option **\n");
      break;
    }
  }
  return 0;
}
