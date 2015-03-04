//static double Coulomb_vdw_factor=1.0;
//#define COULOMB_VDW_FACTOR

__inline__ double vg_mul_double(double x, double y){ return x*y;}
__inline__ float  vg_mul_float (float  x, float  y){ return x*y;}
__inline__ double vg_add_double(double x, double y){ return x+y;}
__inline__ float  vg_add_float (float  x, float  y){ return x+y;}


void MR3calccoulomb_ij_host(int ni, double xi[], double qi[], double force[],
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
#if defined(COULOMB_SHIFT) || MD_LJ_05SIGMA_8SIGMA==2
  VG_UNIT *unit;
  unit=m3_get_unit();
#endif
#ifdef COULOMB_SHIFT
  double rcut21;
  rcut21=1.0/unit->rcut2;
#endif

//  printf("in MR3calccoulomb_ij_host: rscale=%e\n",rscale);
  if(ini==0){
#if MD_LJ_05SIGMA_8SIGMA==2
    //    rcut2=MD_LJ_R2MAX;
    rcut2=unit->rcut2;
    rcutflag=1;
#endif
    if((s=getenv("MR3_HOST_CUTOFF"))!=NULL){
      sscanf(s,"%le",&rcut);
      rcut2=rcut*rcut;
      printf("rcut=%e is read from MR3_HOST_CUTOFF\n",rcut);
      rcutflag=1;
    }
    ini=1;
  }

#if 0
  rcut2=10.0;
  rcutflag=1;
  printf("**** rcut2 is set %e **\n",rcut2);
#endif

  if((periodicflag & 2)!=0) multiplyq=1;

#if 0
  if(rscale!=1.0 || tblno==6 || tblno==7){
    fprintf(stderr,"** error : not supported parameter in MR3calccoulomb_ij_emu **\n");
    MR3_exit(1);
  }
#endif

  if((periodicflag & 1)==0){
    xmax*=2.0;
  }
  for(i=0;i<ni;i++){
    //    if((i % 1000)==0) printf("MR3calccoulomb_ij_host i=%d\n",i);
#if 0
    for(k=0;k<3;k++){
      force[i*3+k]=0.0;
    }
#endif
    //    printf("i=%d xi=%e %e %e qi=%e\n",i,xi[i*3],xi[i*3+1],xi[i*3+2],qi[i]);
    for(j=0;j<nj;j++){
#if 0
      if(i==3 && j<10) printf("MR3calccoulomb_ij_host : force[%d]=%e %e %e\n",i,force[i*3],force[i*3+1],force[i*3+2]);
#endif
      r2=0.0;
      for(k=0;k<3;k++){
	dr[k]=xi[i*3+k]-xj[j*3+k];
#if 1
	if(dr[k]<-xmax/2.0){
	  dr[k]+=xmax;
	  /*	  printf("  i=%d j=%d k=%d shifted xmax\n",i,j,k);*/
	}
	if(dr[k]>=xmax/2.0){
	  dr[k]-=xmax;
	  /*printf("  i=%d j=%d k=%d shifted -xmax\n",i,j,k);*/
	}
#endif
	r2+=dr[k]*dr[k];
      }
      //      r=sqrt(r2);
      x=r2*rscale*rscale;
#if MD_LJ_05SIGMA_8SIGMA==2
      if(r2!=0.0 && (rcutflag==0 || (rcutflag==1 && r2<rcut2))){
#else
      if(r2!=0.0 && (rcutflag==0 || (rcutflag==1 && x<rcut2))){
#endif
	if(tblno==0){
	  rsqrt=1.0/sqrt(r2);
#ifdef COULOMB_SHIFT
	  dtmp=qj[j]*rsqrt;
	  x=1.0-r2*rcut21;
	  dtmp=dtmp*x*x*rsqrt*rsqrt+4.0*rcut21*dtmp*x;
#else
	  dtmp=qj[j]*rsqrt*rsqrt*rsqrt;
#endif
	  if(multiplyq) dtmp*=qi[i]*Coulomb_vdw_factor;
	  for(k=0;k<3;k++){
	    force[i*3+k]+=dtmp*dr[k];
	  }
	  //	  printf(" i=%d j=%d dtmp=%e force=%e %e %e\n",i,j,dtmp,force[i*3],force[i*3+1],force[i*3+2]);
#if 0
	  if(i==0) printf("  force[%d]=%e %e %e\n",
			  i,force[i*3]*qi[i],force[i*3+1]*qi[i],force[i*3+2]*qi[i]);
#endif
	}
	else if(tblno==1){
	  rsqrt=1.0/sqrt(r2);
	  dtmp=qj[j]*rsqrt;
#ifdef COULOMB_SHIFT
	  x=1.0-r2*rcut21;
	  dtmp*=x*x;
#endif
	  if(multiplyq) dtmp*=qi[i]*Coulomb_vdw_factor;
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
#if 0
	  if(i<20 && j<20){
	    printf("    i=%d j=%d rscale=%e r2=%e qj=%e dtmp=%e\n",
		   i,j,rscale,r2,qj[j],dtmp);
	  }
#endif
#if 0
	  if(i==3 && x<10.0){
	    printf("host i=%d j=%d r2=%e rrs=%e xi=%e %e %e xj=%e %e %e\n",i,j,r2,x,xi[i*3],xi[i*3+1],xi[i*3+2],xj[j*3],xj[j*3+1],xj[j*3+2]);
	    printf("     xmax=%e dr=%e %e %e dtmp=%e\n",xmax,dr[0],dr[1],dr[2],dtmp);
	    printf("     factor=%e dtmp*factor=%e\n",factor,dtmp*factor);
	    printf("     force=%e %e %e\n",force[i*3],force[i*3+1],force[i*3+2]);
	  }
#endif
	  for(k=0;k<3;k++){
	    force[i*3+k]+=dtmp*factor;
	  }
#if 0
	  if(i==3 && x<10.0) printf("     force=%e %e %e\n",force[i*3],force[i*3+1],force[i*3+2]);
#endif
	}
      }
#if 0
      if(i==0 && j<4){
	printf("i=%d j=%d xi=%e %e %e\n",i,j,xi[i*3],xi[i*3+1],xi[i*3+2]);
	printf("  xj=%e %e %e r2=%e\n",xj[j*3],xj[j*3+1],xj[j*3+2],r2);
	printf("  qi=%e qj=%e fcore=%e\n",qi[i],qj[j],dtmp*qi[i]);
	printf("  force=%e %e %e\n",force[i*3]*qi[i],force[i*3+1]*qi[i],force[i*3+2]*qi[i]);
      }
#endif
    }
  }
}


void mr3calccoulomb_ij_host_(int *ni, double xi[], double qi[], double force[],
			     int *nj, double xj[], double qj[],
			     double *rscale, 
			     int *tblno, double *xmax, int *periodicflag)
{
  MR3calccoulomb_ij_host(*ni,xi,qi,force,*nj,xj,qj,
			 *rscale,*tblno,*xmax,*periodicflag);
}


void mr3calccoulomb_ij_host__(int *ni, double xi[], double qi[], double force[],
			      int *nj, double xj[], double qj[],
			      double *rscale, 
			      int *tblno, double *xmax, int *periodicflag)
{
  mr3calccoulomb_ij_host_(ni,xi,qi,force,
			  nj,xj,qj,
			  rscale, 
			  tblno,xmax,periodicflag);
}


static double vdw_shift(double r2, 
			double rs, double gs, double rcut21_org, int potflag)
{
  /* potflag 0 --- force
             1 --- potential
  */
  double SIGFAC=pow(2.0,1.0/3.0);
  double A,B;
  static double CACX,CBCX,CADX,CBDX,CCDX;
  double CCA,CCB,CCC,CCD,CA,CC;
  double ENEVDW,ENN,TFVDW;
  static double rcut21=0.0;

  if(rcut21!=rcut21_org){
    rcut21=rcut21_org;
    CACX=-2.0*pow(rcut21,9);
    CBCX=pow(rcut21,6);
    CADX=-pow(rcut21,6);
    CBDX=pow(rcut21,3);
    CCDX=pow(rcut21,-3);
  }
  A=SIGFAC/rs;
  if(potflag) B=gs/4.0;
  else        B=gs*A/(24.0*SIGFAC);
  CCA=B*pow(A,6);
  CCB=2.0*B*pow(A,3);
  CCC=CACX*CCA+CBCX*CCB;
  CCD=CADX*CCA+CBDX*CCB+CCDX*CCC;
  CA=CCA*pow(r2,-6);
  CC=pow(r2,3)*CCC;
  ENEVDW=CA-CCB*pow(r2,-3);
  ENN=ENEVDW-CC+CCD;
  if(potflag){
    return ENN;
  }
  else{
    TFVDW=-6.0*(ENEVDW+CA+CC)/r2;
    return TFVDW;
  }
}


void MR3calcvdw_ij_host(int ni, double xi[], int atypei[], double force[],
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
  double r2min=MD_LJ_R2MIN,r2max=MD_LJ_R2MAX,cuton=CUTON_DEFAULT;
  double tmp1,tmp2;
  VG_UNIT *unit;
  unit=m3_get_unit();
  if(unit!=NULL){
    r2max=unit->rcut2;
    cuton=unit->cuton;
  }
#if VDW_SHIFT>=1
  double rcut21,shift;
  rcut21=1.0/unit->rcut2;
#endif

#if 0
  printf("  ni=%d nj=%d nat=%d tblno=%d xmax=%e periodicflag=%d\n",
	 ni,nj,nat,tblno,xmax,periodicflag);
  for(i=0;i<nat;i++){
    for(j=0;j<nat;j++){
      printf("  i=%d j=%d gscale=%e rscale=%e\n",
	     i,j,gscale[i*nat+j],rscale[i*nat+j]);
    }
  }
  for(i=0;i<ni;i++){
    printf("  xi[%d]=%e %e %e atypei=%d\n",
	   i,xi[i*3],xi[i*3+1],xi[i*3+2],atypei[i]);
  }
  for(i=0;i<nj;i++){
    printf("  xj[%d]=%e %e %e atypej=%d\n",
	   i,xj[i*3],xj[i*3+1],xj[i*3+2],atypej[i]);
  }
#endif
  if((periodicflag & 1)==0){
    xmax*=2.0;
  }
  for(i=0;i<ni;i++){
#if 0
    for(k=0;k<3;k++){
      force[i*3+k]=0.0;
    }
#endif
    for(j=0;j<nj;j++){
      r2=0.0;
      for(k=0;k<3;k++){
#if 0 
	dr[k]=xj[j*3+k]-xi[i*3+k];
#else /* old one until 050512 , work OK */
	dr[k]=xi[i*3+k]-xj[j*3+k];
#endif
#if 1
	if(dr[k]<-xmax/2.0){
	  dr[k]+=xmax;
	  /*	  printf("  i=%d j=%d k=%d shifted xmax\n",i,j,k);*/
	}
	if(dr[k]>=xmax/2.0){
	  dr[k]-=xmax;
	  /*printf("  i=%d j=%d k=%d shifted -xmax\n",i,j,k);*/
	}
#endif
	r2+=dr[k]*dr[k];
      }
      if(r2!=0.0){
	r=sqrt(r2);
#if 1
	rs=rscale[atypei[i]*nat+atypej[j]];
	gs=gscale[atypei[i]*nat+atypej[j]];
#else
	rs=(rscale[atypei[i]*nat+atypei[i]]+rscale[atypej[j]*nat+atypej[j]])*0.5;
	gs=sqrt(gscale[atypei[i]*nat+atypei[i]]*gscale[atypej[j]*nat+atypej[j]]);
#endif
	rrs=r2*rs;
#if MD_LJ_05SIGMA_8SIGMA==2
	if(r2>=r2min && r2<r2max){
#else
	if(rrs>=r2min && rrs<r2max){
#endif
	  r1=1.0/rrs;
	  r6=r1*r1*r1;
	  if(tblno==2){
#if 1 // use separate routine
 #if VDW_SHIFT==0
	    dtmp=VDW_FORCE(r1,r6,1.0,2.0,vg_mul_double,vg_add_double);
 #elif VDW_SHIFT==1
	    shift=(float)(pow(rscale[atypei[i]*nat+atypej[j]]*r2max,-3.0));
	    dtmp=VDW_FORCE_SHIFT(r1,r6,r2,1.0/r2max,shift,1.0,2.0,3.0,vg_mul_double,vg_add_double,tmp1);
 #elif VDW_SHIFT==2
	    shift=cuton*cuton;
	    if(r2>=shift){
	      dtmp=VDW_FORCE_SWITCH(r1,r6,r2,1.0/r2max,shift,1.0,2.0,3.0,vg_mul_double,vg_add_double,tmp1,tmp2);
	    }
	    else{
	      dtmp=VDW_FORCE(r1,r6,1.0,2.0,vg_mul_double,vg_add_double);
	    }
 #endif
#else // original
 #if VDW_SHIFT==1
  #if 0
	    dtmp=r6*r1*(2.0f*r6-1.0f);
	    r1=shift;
	    r2*=rs;
	    r2*=r2;
	    dtmp+=r1*r1*(1.0f-r1)*r2;
  #else
	    dtmp=r1;
	    r1=shift;
   #ifdef VDW_SHIFT_CORRECT
	    dtmp*=r6*(2.0f*r6-1.0f)-r1*(3.0f*r1-2.0f);
   #else
	    r2*=rcut21;
	    r2*=r2*r2;
	    dtmp*=r6*(2.0f*r6-1.0f)-r2*r1*(2.0f*r1-1.0f);
   #endif
  #endif
 #else
	    dtmp=r6*r1*(2.0*r6-1.0);
 #endif
#endif // end of use separate routine
	    dtmp*=gs;
            
	    for(k=0;k<3;k++){
	      force[i*3+k]+=dtmp*dr[k];
	    }
	  }
	  else if(tblno==3){
#if VDW_SHIFT==1
#if 0
	    dtmp=r6*(r6-1.0f);
	    r1=shift;
	    r2*=rs;
	    r2*=r2;
	    //	    r2*=rrs;
	    dtmp+=r1*r1*(r1-1.0f)*r2+2.0f*r1*(1.0f-r1);
	    //	    dtmp+=r1*r1*(r1-1.0f)*rrs+2.0f*r1*(1.0f-r1);
	    dtmp*=gs;
#else // charmm compatible routine
	    {
	      VG_SCALER *scal=unit->scalers;
	      double rcut21=scal->rcut21;
	      dtmp=vdw_shift(r2,rs,gs,rcut21,1);
	    }
#endif

#else
	    dtmp=gs*r6*(r6-1.0);
	    //	    dtmp=gs*(pow(rrs,-6.0)-pow(rrs,-3.0));
#endif
	    for(k=0;k<3;k++){
	      force[i*3+k]+=dtmp;
	    }
	  }
	}
      }
    }
  }
}


void mr3calcvdw_ij_host_(int *ni, double xi[], int atypei[], double force[],
			 int *nj, double xj[], int atypej[],
			 int *nat, double gscale[], double rscale[],
			 int *tblno, double *xmax, int *periodicflag)
{
    int *atypei2,*atypej2;
    int i;

    if((atypei2=(int *)MR3_malloc_pointer(sizeof(int)*(*ni),"atypei2 in mr3calcvdw_ij_host_"))==NULL){
        fprintf(stderr,
                "** error at malloc atypei2 in mr3calcvdw_ij_ **\n");
        MR3_exit(1);
    }
    if((atypej2=(int *)MR3_malloc_pointer(sizeof(int)*(*nj),"atypej2 in mr3calcvdw_ij_host_"))==NULL){
        fprintf(stderr,
                "** error at malloc atypej2 in mr3calcvdw_ij_ **\n");
        MR3_exit(1);
    }
    for(i=0;i<*ni;i++){
      if(atypei[i]>0){
	atypei2[i]=atypei[i]-1;
      }
      else{
	printf("  warning : atypei[%d]=%d should be positive\n",
	       i,atypei[i]);
	atypei2[i]=0;
      }
    }
    for(i=0;i<*nj;i++){
      if(atypej[i]>0){
	atypej2[i]=atypej[i]-1;
      }
      else{
	printf("  warning : atypej[%d]=%d should be positive\n",
	       i,atypej[i]);
      }
    }
    
    MR3calcvdw_ij_host(*ni,xi,atypei2,force,
		      *nj,xj,atypej2,
		      *nat,gscale,rscale,
		      *tblno,*xmax,*periodicflag);
    MR3_free_pointer(atypei2,"atypei2 in mr3calcvdw_ij_host_");
    MR3_free_pointer(atypej2,"atypej2 in mr3calcvdw_ij_host_");
}


void mr3calcvdw_ij_host__(int *ni, double xi[], int atypei[], double force[],
			 int *nj, double xj[], int atypej[],
			 int *nat, double gscale[], double rscale[],
			 int *tblno, double *xmax, int *periodicflag)
{
  mr3calcvdw_ij_host_(ni,xi,atypei,force,
		      nj,xj,atypej,
		      nat,gscale,rscale,
		      tblno,xmax,periodicflag);
}


void MR3calccoulomb_nlist_ij_host_notwork(int ni, double xi[], double qi[], double force[],
				  int nj, double xj[], double qj[], 
				  double rscale, int tblno, double xmax, int periodicflag,
				  int numex[], int natex[], double factorall)
{
  /*
      charge of i-particle is always multiplied 

      periodicflag bit 0 ---- 0 : non periodic
                              1 : periodic
                   bit 1 ---- 0 : use Newton's third law
                                  do not include duplicate list
                              1 : do not use Newton's third law
                   bit 2 ---- 0 : multiply charge and rscale
                              1 : do not multiply charge nor rscale
  */
  int i,j,k;
  double dr[3],r,dtmp,r2,x,dtmp1,factor;
  static int rcutflag=0,ini=0;
  static double rcut,rcut2;
  char *s;
  int multiplyq=0,thirdlaw=1;
  int iexcl=0,jj;

  if(ini==0){
    if((s=getenv("MR3_HOST_CUTOFF"))!=NULL){
      sscanf(s,"%le",&rcut);
      rcut2=rcut*rcut;
      printf("rcut=%e is read from MR3_HOST_CUTOFF\n",rcut);
      rcutflag=1;
    }
    ini=1;
  }

  if((periodicflag & 4)==0) multiplyq=1;
  if((periodicflag & 2)!=0) thirdlaw=0;

  if(periodicflag==0){
    xmax*=2.0;
  }
  for(i=0;i<ni;i++){
    if((i % 10000)==0) printf("MR3calccoulomb_nlist_host i=%d\n",i);
    for(jj=iexcl;jj<iexcl+numex[i];jj++){
      j=natex[jj];
      if(j<0) continue;
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
      r=sqrt(r2);
      if(r!=0.0 && (rcutflag==0 || (rcutflag==1 && r2<rcut2))){
	if(tblno==0){
	  dtmp=qj[j]/(r*r*r);
	  dtmp*=factorall;
	  if(multiplyq) dtmp*=qi[i]*Coulomb_vdw_factor;
	  for(k=0;k<3;k++){
	    force[i*3+k]+=dtmp*dr[k];
	    if(thirdlaw) force[j*3+k]-=dtmp*dr[k];
	  }
	}
	else if(tblno==1){
	  dtmp=qj[j]/(r);
	  dtmp*=factorall;
	  if(multiplyq) dtmp*=qi[i]*Coulomb_vdw_factor;
	  for(k=0;k<3;k++){
	    force[i*3+k]+=dtmp;
	    if(thirdlaw) force[j*3+k]+=dtmp;
	  }
	}
	else if(tblno==6){
	  x=r2*rscale*rscale;
	  dtmp1=(M_2_SQRTPI*exp(-x)*pow(x,-1.0)
	              + erfc(sqrt(x))*pow(x,-1.5));
	  dtmp=qj[j]*dtmp1;
	  factor=rscale*rscale*rscale*qi[i];
	  for(k=0;k<3;k++){
	    force[i*3+k]+=dtmp*dr[k]*factor;
	    if(thirdlaw) force[j*3+k]-=dtmp*dr[k]*factor;
	  }
	}
	else if(tblno==7){
	  x=r2*rscale*rscale;
	  dtmp=qj[j]*erfc(sqrt(x))*pow(x,-0.5);
	  factor=rscale*qi[i];
	  for(k=0;k<3;k++){
	    force[i*3+k]+=dtmp*factor;
	    if(thirdlaw) force[j*3+k]+=dtmp*factor;
	  }
	}
      }
    }
    iexcl+=numex[i];
  }
}


void MR3calccoulomb_nlist_ij_host(int ni, double xid[], double qid[], double force[],
				  int nj, double xjd[], double qjd[], 
				  double rscale, int tblno, double xmaxd, int periodicflag,
				  int numex[], int natex[], double factord)
{
  /*
      charge of i-particle is always multiplied 

      periodicflag bit 0 ---- 0 : non periodic
                              1 : periodic
                   bit 1 ---- 0 : use Newton's third law
                                  do not include duplicate list
                              1 : do not use Newton's third law
                   bit 2 ---- 0 : multiply charge and rscale
                              1 : do not multiply charge nor rscale
  */
  int i,j,k,jj;
  static int rcutflag=0,ini=0;
  static double rcut;
  static float rcut2=MD_LJ_R2MAX;
  char *s;
  int multiplyq=0,thirdlaw=1,iexcl=0;
  double xi[3],xj[3],qi,qj,dr[3],r2,dtmp,x,rsqrt,xmax=xmaxd,factor=factord;
  double rscale2=(rscale*rscale),dtmp2;
  VG_UNIT *unit;
  unit=m3_get_unit();
  if(unit!=NULL) rcut2=unit->rcut2;
#ifdef COULOMB_SHIFT
  double rcut21;
  rcut21=1.0/unit->rcut2;
#endif

  if(ini==0){
#if MD_LJ_05SIGMA_8SIGMA==2
    rcutflag=1;
#endif
    if((s=getenv("MR3_HOST_CUTOFF"))!=NULL){
      sscanf(s,"%le",&rcut);
      rcut2=rcut*rcut;
      printf("rcut=%e is read from MR3_HOST_CUTOFF\n",rcut);
      rcutflag=1;
    }
    ini=1;
  }

  if((periodicflag & 4)==0) multiplyq=1;
  if((periodicflag & 2)!=0) thirdlaw=0;
  //  printf("multiplyq=%d thirdlaw=%d\n",multiplyq,thirdlaw);
  if((periodicflag & 1)==0){
    xmax*=2.0;
  }
  for(i=0;i<ni;i++){
    for(k=0;k<3;k++) xi[k]=xid[i*3+k];
    qi=qid[i];
    //    printf("i=%d iexcl=%d numex[%d]=%d\n",i,iexcl,i,numex[i]);
    for(jj=iexcl;jj<iexcl+numex[i];jj++){
      j=natex[jj];
      //      printf("i=%d jj=%d iexcl=%d numex[%d]=%d j=%d\n",i,jj,iexcl,i,numex[i],j);
      if(j<0) continue;
      for(k=0;k<3;k++) xj[k]=xjd[j*3+k];
      qj=qjd[j];
      r2=0.0;
      for(k=0;k<3;k++){
	dr[k]=xi[k]-xj[k];
	if(dr[k]<-xmax/2.0){
	  dr[k]+=xmax;
	}
	if(dr[k]>=xmax/2.0){
	  dr[k]-=xmax;
	}
	r2+=dr[k]*dr[k];
      }
      x=r2*rscale2;
      //      printf("  r2=%e\n",r2);
      if(r2!=0.0 && (rcutflag==0 || (rcutflag==1 && x<rcut2))){
	if(tblno==0){
	  //	  rsqrt=vg_rsqrtemu(r2);
	  //	  rsqrt=rsqrtf(r2);
	  rsqrt=1.0/sqrt(r2);
#ifdef COULOMB_SHIFT
	  dtmp=qj*rsqrt;
	  x=1.0-r2*rcut21;
	  dtmp=dtmp*x*x*rsqrt*rsqrt+4.0*rcut21*dtmp*x;
#else
	  dtmp=qj*rsqrt*rsqrt*rsqrt;
#endif
	  dtmp*=factor;
	  if(multiplyq) dtmp*=qi*Coulomb_vdw_factor;
	  for(k=0;k<3;k++){
	    force[i*3+k]+=dtmp*dr[k];
	    if(thirdlaw) force[j*3+k]-=dtmp*dr[k];
	  }
	  //	  printf("  i=%d j=%d f=%e %e %e\n",i,j,dtmp*dr[0],dtmp*dr[1],dtmp*dr[2]);
	}
	else if(tblno==1){
	  //rsqrt=vg_rsqrtemu(r2);
	  //	  rsqrt=rsqrtf(r2);
	  rsqrt=1.0/sqrt(r2);
	  dtmp=qj*rsqrt;
#ifdef COULOMB_SHIFT
	  x=1.0-r2*rcut21;
	  dtmp*=x*x;
#endif
	  dtmp*=factor;
	  if(multiplyq) dtmp*=qi*Coulomb_vdw_factor;
	  for(k=0;k<3;k++){
	    force[i*3+k]+=dtmp;
	    if(thirdlaw) force[j*3+k]+=dtmp;
	  }
	}
	else if(tblno==6){
	  x=r2*rscale2;
	  dtmp=qj*(M_2_SQRTPI*exp(-x)*pow(x,-1.0)
	              + erfc(sqrt(x))*pow(x,-1.5));
	  dtmp*=factor*rscale2*rscale;
	  if(multiplyq) dtmp*=qi*Coulomb_vdw_factor;
	  for(k=0;k<3;k++){
	    force[i*3+k]+=dtmp*dr[k];
	    if(thirdlaw) force[j*3+k]-=dtmp*dr[k];
	  }
	}
	else if(tblno==7){
	  x=r2*rscale2;
	  dtmp=qj*erfc(sqrt(x))*pow(x,-0.5);
	  dtmp*=factor*rscale;
	  if(multiplyq) dtmp*=qi*Coulomb_vdw_factor;
	  for(k=0;k<3;k++){
	    force[i*3+k]+=dtmp;
	    if(thirdlaw) force[j*3+k]+=dtmp;
	  }
	}
	else{
	  fprintf(stderr,"** error : tblno=%d is not supported in MR3calccoulomb_nlist_ij_host\n",tblno);
	  vg_exit(1);
	}
      }
    }
    iexcl = iexcl + numex[i];
  }
}


void MR3calccoulomb_nlist_host(double xd[], int n, double qd[], double rscale,
			       int tblno, double xmaxd, int periodicflag,
			       int numex[], int natex[], double factord,
			       double force[])
{
  MR3calccoulomb_nlist_ij_host(n,xd,qd,force,n,xd,qd,rscale,tblno,xmaxd,periodicflag,
			       numex,natex,factord);
}


void MR3calcvdw_nlist_ij_host_notwork(int ni, double xi[], int atypei[], double force[],
			      int nj, double xj[], int atypej[], 
			      int nat, double gscale[], double rscale[], int tblno,
			      double xmax, int periodicflag,
			      int numex[], int natex[], double factor)
{
  /*
      periodicflag bit 0 ---- 0 : non periodic
                              1 : periodic
                   bit 1 ---- 0 : use Newton's third law
                                  do not include duplicate list
                              1 : do not use Newton's third law
  */
  int i,j,k;
  double dr[3],r,dtmp,rs,gs,rrs;
  int iexcl=0,jj,thirdlaw=1;
  double r2min=MD_LJ_R2MIN,r2max=MD_LJ_R2MAX;
  VG_UNIT *unit;
  unit=m3_get_unit();
  if(unit!=NULL) r2max=unit->rcut2;

  if(periodicflag==0){
    xmax*=2.0;
  }
  if((periodicflag & 2)!=0) thirdlaw=0;
  for(i=0;i<ni;i++){
    for(jj=iexcl;jj<iexcl+numex[i];jj++){
      j=natex[jj];
      if(j<0) continue;
      r=0.0;
      for(k=0;k<3;k++){
	dr[k]=xi[i*3+k]-xj[j*3+k];
#if 1
	if(dr[k]<-xmax/2.0){
	  dr[k]+=xmax;
	  /*	  printf("  i=%d j=%d k=%d shifted xmax\n",i,j,k);*/
	}
	if(dr[k]>=xmax/2.0){
	  dr[k]-=xmax;
	  /*printf("  i=%d j=%d k=%d shifted -xmax\n",i,j,k);*/
	}
#endif
	r+=dr[k]*dr[k];
      }
      r=sqrt(r);
      //      if(r!=0.0){
      if(rrs>=r2min && rrs<r2max){
	rs=rscale[atypei[i]*nat+atypej[j]];
	gs=gscale[atypei[i]*nat+atypej[j]];
	rrs=r*r*rs;
#if 0
	if(rrs<0.25){
	  printf("  i=%d j=%d rrs=%e\n",i,j,rrs);
	}
#endif
	if(tblno==2){
	  dtmp=gs*(2.0*pow(rrs,-7.0)-pow(rrs,-4.0));
	  for(k=0;k<3;k++){
	    force[i*3+k]+=dtmp*dr[k];
	    if(thirdlaw) force[j*3+k]-=dtmp*dr[k];
	  }
#if 0
	  if(i<10 || j<10){
	    printf("  i=%d j=%d MR3nlisthost pairforce=%20.17e %20.17e %20.17e\n",
		   i,j,dtmp*dr[0],dtmp*dr[1],dtmp*dr[2]);
	  }
#endif
	}
	else if(tblno==3){
	  dtmp=gs*(pow(rrs,-6.0)-pow(rrs,-3.0));
	  force[i*3]+=dtmp;
	  if(thirdlaw) force[j*3]+=dtmp;
	}
      }
    }
    iexcl+=numex[i];
  }
}


void MR3calcvdw_nlist_ij_host(int ni, double xid[], int atypei[], double force[],
			      int nj, double xjd[], int atypej[], 
			      int nat, double gscale[], double rscale[], int tblno,
			      double xmaxd, int periodicflag,
			      int numex[], int natex[], double factord)
{
  /*
      periodicflag bit 0 ---- 0 : non periodic
                              1 : periodic
                   bit 1 ---- 0 : use Newton's third law
                                  do not include duplicate list
                              1 : do not use Newton's third law
   */
  int i,j,k,iexcl=0,jj,thirdlaw=1;
  double r2min=MD_LJ_R2MIN,r2max=MD_LJ_R2MAX;
  double factor=factord;
  double dtmp,r1,r6,r2,rs,gs,rrs,dr[3],xi[3],xj[3],xmax=xmaxd;
  VG_UNIT *unit;
  unit=m3_get_unit();
  if(unit!=NULL) r2max=unit->rcut2;
#if VDW_SHIFT>=1
  double rcut21,shift;
  rcut21=1.0/unit->rcut2;
#endif

  if((periodicflag & 1)==0){
    xmax*=2.0;
  }
  if((periodicflag & 2)!=0) thirdlaw=0;
  for(i=0;i<ni;i++){
    for(k=0;k<3;k++) xi[k]=xid[i*3+k];
    for(jj=iexcl;jj<iexcl+numex[i];jj++){
      j=natex[jj];
      if(j<0) continue;
      r2=0.0;
      for(k=0;k<3;k++) xj[k]=xjd[j*3+k];
      for(k=0;k<3;k++){
	dr[k]=xi[k]-xj[k];
	if(dr[k]<-xmax/2.0){
	  dr[k]+=xmax;
	}
	if(dr[k]>=xmax/2.0){
	  dr[k]-=xmax;
	}
	r2+=dr[k]*dr[k];
      }
      if(r2!=0.0){
#if VDW_SHIFT==1
	shift=(float)(pow(rscale[atypei[i]*nat+atypej[j]]*r2max,-3.0));
#endif
	rs=rscale[atypei[i]*nat+atypej[j]];
	gs=gscale[atypei[i]*nat+atypej[j]];
	rrs=r2*rs;
#if MD_LJ_05SIGMA_8SIGMA==2
	if(gs!=0.0f && r2>=r2min && r2<r2max){
#else
	if(gs!=0.0f && rrs>=r2min && rrs<r2max){
#endif
	  //r1=vg_r1emu(rrs);
	  r1=1.0/rrs;
	  r6=r1*r1*r1;
	  if(tblno==2){
#if VDW_SHIFT==1
#if 0
	    dtmp=r6*r1*(2.0f*r6-1.0f);
	    r1=shift;
	    r2*=rs;
	    r2*=r2;
	    dtmp+=r1*r1*(1.0f-r1)*r2;
#else
	    dtmp=r1;
	    r1=shift;
#ifdef VDW_SHIFT_CORRECT
	    dtmp*=r6*(2.0f*r6-1.0f)-r1*(3.0f*r1-2.0f);
#else
	    r2*=rcut21;
	    r2*=r2*r2;
	    dtmp*=r6*(2.0f*r6-1.0f)-r2*r1*(2.0f*r1-1.0f);
#endif
#endif
	    dtmp*=gs;
#else
	    dtmp=gs*r6*r1*(2.0*r6-1.0);
#endif
	    dtmp*=factor;
	    for(k=0;k<3;k++){
	      force[i*3+k]+=dtmp*dr[k];
	      if(thirdlaw) force[j*3+k]-=dtmp*dr[k];
	    }
	  }
	  else if(tblno==3){
#if VDW_SHIFT==1
#if 0
	    dtmp=r6*(r6-1.0f);
	    r1=shift;
	    r2*=rs;
	    r2*=r2;
	    //	    r2*=rrs;
	    dtmp+=r1*r1*(r1-1.0f)*r2+2.0f*r1*(1.0f-r1);
	    //	    dtmp+=r1*r1*(r1-1.0f)*rrs+2.0f*r1*(1.0f-r1);
	    dtmp*=gs;
#else // charmm compatible routine
	    {
	      VG_SCALER *scal=unit->scalers;
	      double rcut21=scal->rcut21;
	      dtmp=vdw_shift(r2,rs,gs,rcut21,1);
	    }
#endif
#else
	    dtmp=gs*r6*(r6-1.0);
	    //	    dtmp=gs*(pow(rrs,-6.0)-pow(rrs,-3.0));
#endif
	    dtmp*=factor;
	    for(k=0;k<3;k++){
	      force[i*3+k]+=dtmp;
	      if(thirdlaw) force[j*3+k]+=dtmp;
	    }
	  }
	}
      }
    }
    iexcl = iexcl + numex[i];
  }
}


void MR3calcvdw_nlist_host2(double x[], int n, int atype[], int nat,
			    double gscale[], double rscale[], int tblno,
			    double xmax, int periodicflag,
			    int numex[], int natex[], double factor,
			    double force[])
{
  MR3calcvdw_nlist_ij_host(n,x,atype,force,n,x,atype,
			   nat,gscale,rscale,tblno,xmax,periodicflag,
			   numex,natex,factor);
}


void MR3calccoulomb_ij_exlist_host(int ni, double xi[], double qi[], double force[],
				   int nj, double xj[], double qj[],
				   double rscale, 
				   int tblno, double xmax, int periodicflag,
				   int numex[], int natex[])
{
  /* periodicflag bit 0: 0 --- non periodic
                         1 --- periodic
                  bit 1: 0 --- natex does not include duplicate list
                         1 --- natex includes duplicate list
		  bit 2: 0 --- no multiplication of qi
                         1 --- multiplication of qi
  */
  int i;

  if((periodicflag & 4)==0){
    fprintf(stderr,"** error : bit2 of periodicflag must be activated **\n");
    vg_exit(1);
  }
  MR3calccoulomb_nlist_ij_host(ni,xi,qi,force,nj,xj,qj,rscale,tblno,xmax,(periodicflag & 3),
			       numex,natex,-1.0);
  MR3calccoulomb_ij_host(ni,xi,qi,force,nj,xj,qj,rscale,tblno,xmax,(periodicflag & 1)+2);
}


void MR3calcvdw_ij_exlist_host(int ni, double xi[], int atypei[], double force[],
			       int nj, double xj[], int atypej[],
			       int nat, double gscale[], double rscale[],
			       int tblno, double xmax, int periodicflag,
			       int numex[], int natex[])
{
  /* periodicflag bit 0: 0 --- non periodic
                         1 --- periodic
                  bit 1: 0 --- natex does not include duplicate list
                         1 --- natex includes duplicate list
  */
  int i,j;
  MR3calcvdw_nlist_ij_host(ni,xi,atypei,force,nj,xj,atypej,nat,gscale,rscale,tblno,
			   xmax,(periodicflag & 3),numex,natex,-1.0);
  MR3calcvdw_ij_host(ni,xi,atypei,force,nj,xj,atypej,nat,gscale,rscale,tblno,
		     xmax,(periodicflag & 1));
}


void MR3calcewald_dft_host(int k[], int knum, double x[], int n,
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

void MR3calcewald_idft_eng_host(int k[], double bs[], double bc[],
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

    //    printf("HOST s*cos c*sin %d %e %e\n",i,fctmp,fstmp);

  }
}


void MR3calcewald_idft_force_host(int k[], double bs[], double bc[],
                                         int knum, double x[], int n,
                                         double cellsize[3],
                                         double force[])
{
  int i,i3,j,j3,c;
  double th,sth,cth,cellsize_1[3];
  double fst[3],fct[3];
  double fstmp[3],fctmp[3];
  /*
  double cmin,cmax;
  double smin,smax;
  */

  for(i=0;i<3;i++) cellsize_1[i]=1.0/cellsize[i];
  for(i=i3=0;i<n;i++,i3+=3){
    /*
    {
      cmin=1e40;
      cmax=0.0;
      smin=1e40;
      smax=0.0;
    }
    */
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
      /*
      {
	if(fabs(fst[0])>smax)smax=fabs(fst[0]);
	if(fabs(fst[0])<smin)smin=fabs(fst[0]);
	if(fabs(fct[0])>cmax)cmax=fabs(fst[0]);
	if(fabs(fct[0])<cmin)cmin=fabs(fst[0]);
      }
      */
      for(c=0;c<3;c++){
        fstmp[c]+=fst[c];
        fctmp[c]+=fct[c];
      }


    }
    for(c=0;c<3;c++){
      force[i3+c]+=fstmp[c]-fctmp[c];
    }
    /*
    {
#define WARN_RANGE 1.0e-4
      if(smax/fabs(fstmp[0])>WARN_RANGE)printf("%d fsc max/sum %e / %e = %e\n",i,smax,fstmp[0],smax/fstmp[0]);
      if(cmax/fabs(fctmp[0])>WARN_RANGE)printf("%d fcs max/sum %e / %e = %e\n",i,cmax,fctmp[0],cmax/fctmp[0]);
      if(fabs(fstmp[c]/force[i3])>WARN_RANGE)printf("%d fsc/force %e / %e = %e\n",i,fstmp[c],force[i3],fstmp[c]/force[i3]);
      if(fabs(fctmp[c]/force[i3])>WARN_RANGE)printf("%d fcs/force %e / %e = %e\n",i,fctmp[c],force[i3],fctmp[c]/force[i3]);
    }
    */
  }
}


void MR3calcewald_host(int *k, int knum_org, double *x, int n, double *chg,
                       double alpha, double epsilon, double cell[3][3],
                       double *force, double *tpot, double stress[3][3])
{
  double *bs,*bc,cellsize[3],cellsize_1[3];
  int knum,i,j,c,i3,j3;
  double factor1_tmp,vol1,eps1,alpha4,r2,kvtmp;

  knum=knum_org<0 ? -knum_org:knum_org;
  if((bs=(double *)MR3_malloc_pointer(sizeof(double)*knum,"bs in MR3calcewald_host"))==NULL){
    fprintf(stderr,"** error : can't malloc bs **\n");
    MR3_exit(1);
  }
  if((bc=(double *)MR3_malloc_pointer(sizeof(double)*knum,"bc in MR3calcewald_host"))==NULL){
    fprintf(stderr,"** error : can't malloc bc **\n");
    MR3_exit(1);
  }
  for(i=0;i<3;i++) cellsize[i]=cell[i][i];
  for(i=0;i<3;i++) cellsize_1[i]=1.0/cellsize[i];
  for(i=0,vol1=1.0;i<3;i++) vol1*=cellsize_1[i];
  eps1=1.0/epsilon;
  alpha4=1.0/(4.0*alpha*alpha);
  for(i=0;i<n*3;i++) force[i]=0.0;

  /* DFT */
  MR3calcewald_dft_host(k,knum,x,n,chg,cellsize,bs,bc);

#if 0
  for(i=0;i<10;i++){
    printf("bs bc %d %e %e\n",i,bs[i],bc[i]);
  }
#endif

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
    MR3calcewald_idft_eng_host(k,bs,bc,knum,x,n,cellsize,force);
    for(i=i3=0;i<n;i++,i3+=3) force[i3]*=chg[i];
  }
  else{           /* force */
    MR3calcewald_idft_force_host(k,bs,bc,knum,x,n,cellsize,force);
    for(i=i3=0;i<n;i++,i3+=3){
      for(c=0;c<3;c++) force[i3+c]*=2.0*M_PI*chg[i]*cellsize_1[c];
    }
  }

  MR3_free_pointer(bs,"bs in MR3calcewald_host");
  MR3_free_pointer(bc,"bc in MR3calcewald_host");
}


void mr3calcewald_host__(k,knum,x,n,chg,alpha,epsilon,cell,
                         force,tpot,stress)
int *n,*k,*knum;
double *x,*chg,cell[3][3],*force,*tpot,*alpha,*epsilon,stress[3][3];
{
    MR3calcewald_host(k,*knum,x,*n,chg,*alpha,*epsilon,cell,
                      force,tpot,stress);
}

void mr3calcewald_host_(k,knum,x,n,chg,alpha,epsilon,cell,
			force,tpot,stress)
int *n,*k,*knum;
double *x,*chg,cell[3][3],*force,*tpot,*alpha,*epsilon,stress[3][3];
{
    MR3calcewald_host(k,*knum,x,*n,chg,*alpha,*epsilon,cell,
                      force,tpot,stress);
}


void MR3SetTable(char *filename, int tblno, int flag)
{
  fprintf(stderr,"*** warning : this function is not implemented ***\n");
}


void mr3settable_nf_(char *fname, int *tblno, int *flag, int *n)
{
}


void mr3settable_nf__(char *fname, int *tblno, int *flag, int *n)
{
  mr3settable_nf_(fname,tblno,flag,n);
}


void MR3_setup_exlist(int n, int *bckptr)
{
}


void mr3_setup_exlist_(int *n, int *bckptr)
{
#if 0
  int i;
  int *p=NULL;

  if((p=(int *)MR3_malloc_pointer(sizeof(int)*(*n),"p in mr3_setup_exlist_"))==NULL){
    fprintf(stderr,"** error : can't malloc p **\n");
    MR3_exit(1);
  }
  for(i=0;i<*n;i++){
    p[i]=bckptr[i]-1;
  }
  MR3_setup_exlist(*n,p);
  MR3_free_pointer(p,"p in mr3_setup_exlist_");
#endif
}


void MR3_make_natex2_from_natex(int n, int natex[], int numex[],
				int **natex_offset,
				int **natex2, int **numex2,
				int **natex2_offset, int flag)
{
  /* flag 0 -- do not malloc
          1 -- malloc
  */
  int i,s,j,e;
  int dispflag=0;

  for(i=s=0;i<n;i++) s+=numex[i];

  if(flag==1){
    if((*natex_offset=(int *)MR3_malloc_pointer(sizeof(int)*n,"*natex_offset in MR3_make_natex2_from_natex"))==NULL){
      fprintf(stderr,"** error : can't malloc natex_offset **\n");
      MR3_exit(1);
    }
    if((*natex2=(int *)MR3_malloc_pointer(sizeof(int)*s*2,"*natex2 in MR3_make_natex2_from_natex"))==NULL){
      fprintf(stderr,"** error : can't malloc natex2 **\n");
      MR3_exit(1);
    }
    if((*numex2=(int *)MR3_malloc_pointer(sizeof(int)*n,"*numex2 in MR3_make_natex2_from_natex"))==NULL){
      fprintf(stderr,"** error : can't malloc numex2 **\n");
      MR3_exit(1);
    }
    if((*natex2_offset=(int *)MR3_malloc_pointer(sizeof(int)*n,"*natex2_offset in MR3_make_natex2_from_natex"))==NULL){
      fprintf(stderr,"** error : can't malloc natex2_offset **\n");
      MR3_exit(1);
    }
  }

  /* make natex_offset */
  (*natex_offset)[0]=0;
  for(i=1;i<n;i++){
    (*natex_offset)[i]=(*natex_offset)[i-1]+numex[i-1];
  }
#if 1
  if(dispflag){
    printf("original natex and numex\n");
    for(i=0;i<n;i++){
      for(e=0;e<numex[i];e++){
	j=natex[(*natex_offset)[i]+e];
	if(dispflag){
	  if(e==0) printf("natex[%d]=",i);
	  printf("%d ",j);
	  if(e==numex[i]-1) printf("\n");
	}
      }
    }
  }
#endif

  /* make numex2 */
  for(i=0;i<n;i++) (*numex2)[i]=0;
  for(i=0;i<n;i++){
    for(e=0;e<numex[i];e++){
      j=natex[(*natex_offset)[i]+e];
      if(dispflag){
	if(e==0) printf("natex[%d]=",i);
	printf("%d ",j);
	if(e==numex[i]-1) printf("\n");
      }
      if(j>=0){
	(*numex2)[i]++;
	(*numex2)[j]++;
      }
#if 0
      else(j<0){
	fprintf(stderr,"** error : j=%d < 0 **\n",j);
	MR3_exit(1);
      }
#endif
    }
  }
  if(dispflag){
    for(i=0;i<n;i++){
      printf("numex2[%d]=%d\n",i,(*numex2)[i]);
    }
  }

  /* make natex2_offset */
  (*natex2_offset)[0]=0;
  for(i=1;i<n;i++){
    (*natex2_offset)[i]=(*natex2_offset)[i-1]+(*numex2)[i-1];
  }

  /* make natex2 */
  for(i=0;i<n;i++) (*numex2)[i]=0;
  for(i=0;i<n;i++){
    for(e=0;e<numex[i];e++){
      j=natex[(*natex_offset)[i]+e];
      if(j>=0){
	(*natex2)[(*natex2_offset)[i]+(*numex2)[i]]=j;
	(*natex2)[(*natex2_offset)[j]+(*numex2)[j]]=i;
	(*numex2)[i]++;
	(*numex2)[j]++;
      }
    }
  }
  if(dispflag){
    for(i=0;i<n;i++){
      for(e=0;e<(*numex2)[i];e++){
	j=(*natex2)[(*natex2_offset)[i]+e];
	if(j<0){
	  fprintf(stderr,"** error : j=%d < 0 **\n",j);
	  MR3_exit(1);
	}
	if(dispflag){
	  if(e==0) printf("natex2[%d]=",i);
	  printf("%d ",j);
	  if(e==(*numex2)[i]-1) printf("\n");
	}
      }
    }
  }
}


void calccoulomb_vdw_ij_exlist_host_sub(int i, int j,
					       double xi[], double xj[],
					       double qj,
					       double gsf, double gsp,
					       double rs, double rscale,
					       double cutoff,
					       double *fc, double *fv,
					       double *pc, double *pv,
					       int tblno, double size[3],
					       int potflag)
{
  double r,dr[3],rrs,dtmp,r2,rrs_1,rrs_3,shift,tmp1,tmp2;
  int k;
  double r2min=MD_LJ_R2MIN,r2max=MD_LJ_R2MAX,cuton=CUTON_DEFAULT;
  VG_UNIT *unit;
  unit=m3_get_unit();
  if(unit!=NULL){
    r2max=unit->rcut2;
    cuton=unit->cuton;
  }

  r2=0.0;
  for(k=0;k<3;k++){
    dr[k]=xi[k]-xj[k];
#if 1
    if(dr[k]<-size[k]/2.0){
      dr[k]+=size[k];
    }
    if(dr[k]>=size[k]/2.0){
      dr[k]-=size[k];
    }
#endif
    r2+=dr[k]*dr[k];
  }
  r=sqrt(r2);
  if(r!=0.0 && r<cutoff){
    rrs=r2*rs;
    rrs_1=1.0/rrs;
    rrs_3=rrs_1*rrs_1*rrs_1;
    if(tblno==0 || tblno==6){ /* vdw force */
#if MD_LJ_05SIGMA_8SIGMA==2
      if(r2>=r2min && r2<r2max){
#else
      if(rrs>=r2min && rrs<r2max){
#endif
#if VDW_SHIFT==0
	//	dtmp=gsf*(2.0*pow(rrs,-7.0)-pow(rrs,-4.0));
	dtmp=gsf*VDW_FORCE(rrs_1,rrs_3,1.0,2.0,vg_mul_double,vg_add_double);
#elif VDW_SHIFT==1
	shift=pow(rs*r2max,-3.0);
	dtmp=gsf*VDW_FORCE_SHIFT(rrs_1,rrs_3,r2,1.0/r2max,shift,1.0,2.0,3.0,vg_mul_double,vg_add_double,tmp1);
	//	printf("rr=%20.10f tfvdw=%20.10f\n",r2,-dtmp);
#elif VDW_SHIFT==2
	shift=cuton*cuton;
	if(r2>=shift){
	  dtmp=gsf*VDW_FORCE_SWITCH(rrs_1,rrs_3,r2,1.0/r2max,shift,1.0,2.0,3.0,vg_mul_double,vg_add_double,tmp1,tmp2);
	}
	else{
	  dtmp=gsf*VDW_FORCE(rrs_1,rrs_3,1.0,2.0,vg_mul_double,vg_add_double);
	}
#endif
	for(k=0;k<3;k++){
	  fv[i*3+k]+=dtmp*dr[k];
	}
      }
    }

    if((tblno==1 || tblno==7) ||
       ((tblno==0 || tblno==6) && potflag)){ /* vdw potential */
#if MD_LJ_05SIGMA_8SIGMA==2
      if(r2>=r2min && r2<r2max){
#else
      if(rrs>=r2min && rrs<r2max){
#endif
#if VDW_SHIFT==0
	//	dtmp=gsp*(pow(rrs,-6.0)-pow(rrs,-3.0));
	dtmp=gsp*VDW_ENERGY(rrs_3,1.0,vg_mul_double,vg_add_double);
#elif VDW_SHIFT==1
	shift=pow(rs*r2max,-3.0);
	dtmp=gsp*VDW_ENERGY_SHIFT(rrs_1,rrs_3,r2,1.0/r2max,shift,1.0,2.0,3.0,vg_mul_double,vg_add_double,tmp1);
	//	printf("rr=%20.10f vdw=%20.10f\n",r2,dtmp);
#elif VDW_SHIFT==2
	shift=cuton*cuton;
	if(r2>=shift){
	  dtmp=gsp*VDW_ENERGY_SWITCH(rrs_1,rrs_3,r2,shift,1.0,2.0,3.0,vg_mul_double,vg_add_double,tmp1);
	}
	else{
	  dtmp=gsp*VDW_ENERGY(rrs_3,1.0,vg_mul_double,vg_add_double);
	}
#endif
	pv[i]+=dtmp;
      }
    }
    
    if(tblno==0){ /* Coulomb force */
      dtmp=qj*pow(r,-3.0);
      for(k=0;k<3;k++){
	fc[i*3+k]+=dtmp*dr[k];
      }
    }
    if(tblno==1 ||
       (tblno==0 && potflag)){ /* Coulomb potential */
      pc[i]+=qj*pow(r,-1.0);
    }
    if(tblno==6){ /* real-space of Coulomb force */
      rrs=r*r*rscale*rscale;
      dtmp=qj*(M_2_SQRTPI*exp(-rrs)*pow(rrs,-1.0)
		  + erfc(sqrt(rrs))*pow(rrs,-1.5));
      for(k=0;k<3;k++){
	fc[i*3+k]+=dtmp*dr[k];
      }
    }
    if(tblno==7 ||
       (tblno==6 && potflag)){ /* real-space of Coulomb potential */
      //      rrs=r*r*rscale*rscale;
      rrs=r2*rscale*rscale;
      dtmp=qj*erfc(sqrt(rrs))*pow(rrs,-0.5);
      //      {double cele=332.0716;printf("r=%20.10f rk=%20.10f q=%20.10f e=%20.10f\n",r,r*rscale,qj,rscale*cele*dtmp);}
      pc[i]+=dtmp;
    }
  }
}


void MR3calccoulomb_vdw_ij_exlist_cutoff_host(int ni, double xi[], double qi[], 
					      int atypei[], double force[], 
					      int nj, double xj[], double qj[], 
					      int atypej[],
					      int nat,
					      double gscalesf[], 
					      double gscalesp[], double rscales[],
					      double rscale, int tblno,
					      double cutoff,
					      double size[3], double *potc, double *potv,
					      int periodicflag, 
					      int potflag, int changeflag,
					      int numex[], int natex[])
{
  int i,j,k;
  double dr[3],r,dtmp,rs,gsf,gsp,rrs;
  int iexcl=0,jj;
  double *fc,*fv,*pc,*pv;
  int *natex_offset,*natex2,*numex2,*natex2_offset;
  int nistart=0;
  int nidisp=20,dispflag=0;

  if((periodicflag & 2)==0){ /* duplicate list must be generated */
    MR3_make_natex2_from_natex(ni,natex,numex,&natex_offset,
			       &natex2,&numex2,&natex2_offset,1);
  }
  else{
    natex2=natex;
    numex2=numex;
    natex_offset=natex2_offset=NULL;
    /*
    if((natex2_offset=(int *)MR3_malloc_pointer(sizeof(int)*ni,"natex2_offset in MR3calccoulomb_vdw_ij_exlist_cutoff_host"))==NULL){
      fprintf(stderr,"** error : can't malloc natex2_offset **\n");
      MR3_exit(1);
    }
    natex2_offset[0]=0;
    for(i=1;i<ni;i++) natex2_offset[i]=natex2_offset[i-1]+numex[i];
    */
  }

  if((fc=(double *)MR3_malloc_pointer(sizeof(double)*ni*3,"fc in MR3calccoulomb_vdw_ij_exlist_cutoff_host"))==NULL){
    fprintf(stderr,"** error : can't malloc fc **\n");
    MR3_exit(1);
  }
  if((fv=(double *)MR3_malloc_pointer(sizeof(double)*ni*3,"fv in MR3calccoulomb_vdw_ij_exlist_cutoff_host"))==NULL){
    fprintf(stderr,"** error : can't malloc fv **\n");
    MR3_exit(1);
  }
  if((pc=(double *)MR3_malloc_pointer(sizeof(double)*ni,"pc in MR3calccoulomb_vdw_ij_exlist_cutoff_host"))==NULL){
    fprintf(stderr,"** error : can't malloc pc **\n");
    MR3_exit(1);
  }
  if((pv=(double *)MR3_malloc_pointer(sizeof(double)*ni,"pv in MR3calccoulomb_vdw_ij_exlist_cutoff_host"))==NULL){
    fprintf(stderr,"** error : can't malloc pv **\n");
    MR3_exit(1);
  }
  for(i=0;i<ni*3;i++) fc[i]=fv[i]=0.0;
  for(i=0;i<ni;i++)   pc[i]=pv[i]=0.0;

#if 0
  nistart=80670;
  ni=80673;
  printf("ni or/and nistart are changed to %d, %d\n",ni,nistart);
#endif

  if((periodicflag & 1)==0){
    for(j=0;j<3;j++) size[j]*=2.0;
  }
  for(i=0;i<nistart;i++) iexcl+=numex2[i];
  for(i=nistart;i<ni;i++){
    //    printf("xi[%d]=%e %e %e\n",i,xi[i*3],xi[i*3+1],xi[i*3+2]);
    for(jj=iexcl;jj<iexcl+numex2[i];jj++){
      j=natex2[jj];
      if(j<0) continue;
      rs=rscales[atypei[i]*nat+atypej[j]];
      gsf=gscalesf[atypei[i]*nat+atypej[j]];
      gsp=gscalesp[atypei[i]*nat+atypej[j]];
      calccoulomb_vdw_ij_exlist_host_sub(i,j,xi+i*3,xj+j*3,qj[j],
					 gsf,gsp,rs,rscale,
					 cutoff,
					 fc,fv,pc,pv,tblno,size,
					 potflag & 1);
    }
    iexcl+=numex2[i];
    //    if(i==0) printf("host pc=%e pv=%e fc=%e %e %e fv=%e %e %e\n",pc[0],pv[0],fc[0],fc[1],fc[2],fv[0],fv[1],fv[2]);
  }

  if(dispflag){
    printf("bonded part\n");
    //    for(i=0;i<nidisp;i++){
    for(i=nistart;i<ni;i++){
      //      if((i>=121 && i<=123) || (i>=247 && i<=249)){
      //      if(i>=123 && i<=126){
      if(1){
	printf("fc[%d]=%e %e %e pc=%e pcq=%e\n",
	       i,fc[i*3],fc[i*3+1],fc[i*3+2],pc[i],pc[i]*qi[i]);
	//	printf("fv[%d]=%e %e %e pv=%e\n",
	//	       i,fv[i*3],fv[i*3+1],fv[i*3+2],pv[i]);
      }
    }
  }
#if 0
  for(i=nistart;i<ni;i++){
    if(i>2530){
	printf("fc[%d]=%e %e %e pc=%e pcq=%e\n",
	       i,fc[i*3],fc[i*3+1],fc[i*3+2],pc[i],pc[i]*qi[i]);
	printf("fv[%d]=%e %e %e pv=%e\n",
	       i,fv[i*3],fv[i*3+1],fv[i*3+2],pv[i]);
    }
  }
#endif

  /* reverse bonded part */
  for(i=nistart;i<ni;i++){
    for(j=0;j<3;j++){
      fc[i*3+j]*=-1.0;
      fv[i*3+j]*=-1.0;
    }
    pc[i]*=-1.0;
    pv[i]*=-1.0;
  }

  /* i-j interaction */
  for(i=nistart;i<ni;i++){
    if((i % 1000)==0) printf("host calculation: i=%d ni=%d\n",i,ni);
    for(j=0;j<nj;j++){
    //    {static int ini=0;if(ini==0){ fprintf(stderr,"*** j is reduced ***\n");ini=1;}} for(j=i+1;j<nj;j++){
      rs=rscales[atypei[i]*nat+atypej[j]];
      gsf=gscalesf[atypei[i]*nat+atypej[j]];
      gsp=gscalesp[atypei[i]*nat+atypej[j]];
      calccoulomb_vdw_ij_exlist_host_sub(i,j,xi+i*3,xj+j*3,qj[j],
					 gsf,gsp,rs,rscale,
					 cutoff,
					 fc,fv,pc,pv,tblno,size,
					 potflag & 1);
      //      {static double dtmpbak=0.0;printf("qi=%20.10f qj=%20.10f qij=%20.10f pc[%d]=%20.10f\n",qi[i],qj[j],qi[i]*qj[j],i,pc[i]);printf("e=%20.10f\n",qi[i]*(pc[i]-dtmpbak)*rscale*332.0716);printf("aqi=%20.10f e=%20.10f cvf=%20.10f\n",qi[i],qi[i]*(pc[i]-dtmpbak)*rscale*332.0716,Coulomb_vdw_factor); if(j==nj-1) dtmpbak=pc[i+1];else dtmpbak=pc[i];}
    }
    //    if(i==0) printf("host pc=%e pv=%e fc=%e %e %e fv=%e %e %e\n",pc[0],pv[0],fc[0],fc[1],fc[2],fv[0],fv[1],fv[2]);
  }

  if(dispflag){
    printf("bonded part + non-bonded part\n");
    //    for(i=0;i<nidisp;i++){
    for(i=nistart;i<ni;i++){
      //      if(fabs(pc[i]*qi[i])>0.1){
      //      if((i>=121 && i<=123) || (i>=247 && i<=249)){
      //      if(i>=123 && i<=126){
      //      if(1){
      if(fabs(pc[i]*qi[i])>0.01){
	printf("fc[%d]=%e %e %e pc=%e pcq=%e\n",
	       i,fc[i*3],fc[i*3+1],fc[i*3+2],pc[i],pc[i]*qi[i]);
	//	printf("fv[%d]=%e %e %e pv=%e\n",
	//	       i,fv[i*3],fv[i*3+1],fv[i*3+2],pv[i]);
      }
    }
  }

  /* multiply qi */
#if !defined(COULOMB_VDW_FACTOR) && 0
  fprintf(stderr,"** erorr : this routine need to define COULOMB_VDW_FACTOR **\n");
  MR3_exit(1);
#endif
  for(i=nistart;i<ni;i++){
    if(tblno==0){
      for(j=0;j<3;j++) fc[i*3+j]*=qi[i]*Coulomb_vdw_factor;
      if(potflag & 1) pc[i]*=qi[i]*Coulomb_vdw_factor;
    }
    else if(tblno==6){
      for(j=0;j<3;j++) fc[i*3+j]*=rscale*rscale*rscale*qi[i]*Coulomb_vdw_factor;
      if(potflag & 1) pc[i]*=rscale*qi[i]*Coulomb_vdw_factor;
    }
  }

  /* add to return array */
  for(i=nistart;i<ni;i++){
    for(j=0;j<3;j++) force[i*3+j]+=fc[i*3+j]+fv[i*3+j];
  }
  for(i=nistart;i<ni;i++){
    *potc+=pc[i];
    *potv+=pv[i];
#if 0
    if(i>2530) printf("  pc[%d]=%e pv=%e\n",i,pc[i],pv[i]);
#endif
  }
  
  MR3_free_pointer(fc,"fc in MR3calccoulomb_vdw_ij_exlist_cutoff_host");
  MR3_free_pointer(fv,"fv in MR3calccoulomb_vdw_ij_exlist_cutoff_host");
  MR3_free_pointer(pc,"pc in MR3calccoulomb_vdw_ij_exlist_cutoff_host");
  MR3_free_pointer(pv,"pv in MR3calccoulomb_vdw_ij_exlist_cutoff_host");
  MR3_free_pointer(natex_offset,"natex_offset in MR3calccoulomb_vdw_ij_exlist_cutoff_host");
  if((periodicflag & 2)==0){
    MR3_free_pointer(natex2,"natex2 in MR3calccoulomb_vdw_ij_exlist_cutoff_host");
    MR3_free_pointer(numex2,"numex2 in MR3calccoulomb_vdw_ij_exlist_cutoff_host");
  }
  MR3_free_pointer(natex2_offset,"natex2_offset in MR3calccoulomb_vdw_ij_exlist_cutoff_host");
}


void MR3calccoulomb_vdw_ij_exlist_host(int ni, double xi[], double qi[], 
				       int atypei[], double force[], 
				       int nj, double xj[], double qj[], 
				       int atypej[],
				       int nat,
				       double gscalesf[], 
				       double gscalesp[], double rscales[],
				       double rscale, int tblno,
				       double xmax, double *potc, double *potv,
				       int periodicflag, 
				       int potflag, int changeflag,
				       int numex[], int natex[])
{
  /*
    periodicflag bit 0 : 0 --- non periodic
                         1 --- periodic
                 bit 1 : 0 --- natex and numex do not have duplicate list
                         1 ---                 have duplicate list
   */
  int i,j,k;
  double dr[3],r,dtmp,rs,gsf,gsp,rrs;
  int iexcl=0,jj;
  double *fc,*fv,*pc,*pv;
  int *natex_offset,*natex2,*numex2,*natex2_offset;
  int nidisp=20,dispflag=0;
  double cutoff=9999.9,size[3];
  double *ftmp;
  VG_UNIT *unit=m3_get_unit();

  size[0]=size[1]=size[2]=xmax;
#if MD_LJ_05SIGMA_8SIGMA==2
  cutoff=sqrt(MD_LJ_R2MAX);
  if(unit!=NULL){
    if(unit->rcut2!=0.0) cutoff=sqrt(unit->rcut2);
  }
#endif
  MR3calccoulomb_vdw_ij_exlist_cutoff_host(ni,xi,qi,atypei,force, 
					   nj,xj,qj,atypej,
					   nat,gscalesf,gscalesp,rscales,
					   rscale,tblno,cutoff,
					   size,potc,potv,
					   periodicflag,potflag,changeflag,
					   numex,natex);
}


void MR3calccoulomb_vdw_ij_exlist_host2(int ni, double xi[], double qi[], 
					int atypei[], double force[], 
					int nj, double xj[], double qj[], 
					int atypej[],
					int nat,
					double gscalesf[], 
					double gscalesp[], double rscales[],
					double rscale, int tblno,
					double xmax, double *potc, double *potv,
					int periodicflag, 
					int potflag, int changeflag,
					int numex[], int natex[])
{
  /*
    periodicflag bit 0 : 0 --- non periodic
                         1 --- periodic
                 bit 1 : 0 --- natex and numex do not have duplicate list
                         1 ---                 have duplicate list
   */
  double cutoff=9999.9,size[3];
  double *ftmp;
  int i;

  size[0]=size[1]=size[2]=xmax;
  if(potflag!=0){
    if((ftmp=(double *)MR3_malloc_pointer(sizeof(double)*ni*3,"MR3calccoulomb_vdw_ij_exlist_host2"))==NULL){
      fprintf(stderr,"** error : can't malloc ftmp **\n");
      vg_exit(1);
    }
    bzero(ftmp,sizeof(double)*ni*3);
    MR3calccoulomb_ij_exlist_host(ni,xi,qi,ftmp,nj,xj,qj,
				  rscale,tblno+1,xmax,(periodicflag & 3)+4,
				  numex,natex);
    for(i=0;i<ni;i++) (*potc)+=ftmp[i*3];
    bzero(ftmp,sizeof(double)*ni*3);
    MR3calcvdw_ij_exlist_host(ni,xi,atypei,ftmp,nj,xj,atypej,
			      nat,gscalesp,rscales,3,xmax,
			      (periodicflag & 3),
			      numex,natex);
    for(i=0;i<ni;i++) (*potv)+=ftmp[i*3];
    MR3_free_pointer(ftmp,"MR3calccoulomb_vdw_ij_exlist_host2");
  }
  MR3calccoulomb_ij_exlist_host(ni,xi,qi,force,nj,xj,qj,
				rscale,tblno,xmax,(periodicflag & 3)+4,
				numex,natex);
  MR3calcvdw_ij_exlist_host(ni,xi,atypei,force,nj,xj,atypej,
			    nat,gscalesf,rscales,2,xmax,(periodicflag & 3),
			    numex,natex);
}


void mr3calccoulomb_vdw_ij_exlist_host_(int *ni, double xi[], double qi[], 
					int atypei[], double force[], 
					int *nj, double xj[], double qj[], 
					int atypej[],
					int *nat,
					double gscalesf[], 
					double gscalesp[], double rscales[],
					double *rscale, int *tblno,
					double *xmax, double *potc, double *potv,
					int *periodicflag, 
					int *potflag, int *changeflag,
					int numex[], int natex[])
{
  int *atypei2,*atypej2;
  int *natex2,s;
  int i;

  for(i=s=0;i<*ni;i++){
    s+=numex[i];
  }
  if((atypei2=(int *)MR3_malloc_pointer(sizeof(int)*(*ni),"atypei2 in mr3calccoulomb_vdw_ij_exlist_host_"))==NULL){
    fprintf(stderr,
	    "** error at malloc atypei2 in mr3calccoulomb_vdw_ij_exlist_host_ **\n");
    MR3_exit(1);
  }
  if((atypej2=(int *)MR3_malloc_pointer(sizeof(int)*(*nj),"atypej2 in mr3calccoulomb_vdw_ij_exlist_host_"))==NULL){
    fprintf(stderr,
	    "** error at malloc atypej2 in mr3calccoulomb_vdw_ij_exlist_host_ **\n");
    MR3_exit(1);
  }
  if((natex2=(int *)MR3_malloc_pointer(sizeof(int)*s,"natex2 in mr3calccoulomb_vdw_ij_exlist_host_"))==NULL){
    fprintf(stderr,
	    "** error at malloc natex2 in mr3calccoulomb_vdw_ij_exlist_host_ **\n");
    MR3_exit(1);
  }

  for(i=0;i<*ni;i++){
    if(atypei[i]>0){
      atypei2[i]=atypei[i]-1;
    }
    else{
      printf("  warning : atypei[%d]=%d should be positive\n",
	     i,atypei[i]);
      atypei2[i]=0;
    }
  }
  for(i=0;i<*nj;i++){
    if(atypej[i]>0){
      atypej2[i]=atypej[i]-1;
    }
    else{
      printf("  warning : atypej[%d]=%d should be positive\n",
	     i,atypej[i]);
    }
  }
  for(i=0;i<s;i++) natex2[i]=natex[i]-1;

  MR3calccoulomb_vdw_ij_exlist_host(*ni,xi,qi,atypei2,force, 
				    *nj,xj,qj,atypej2,
				    *nat,gscalesf,gscalesp,rscales,*rscale,
				    *tblno,*xmax,potc,potv,
				    *periodicflag,*potflag,*changeflag,
				    numex,natex2);
  
  MR3_free_pointer(atypei2,"atypei2 in mr3calccoulomb_vdw_ij_exlist_host_");
  MR3_free_pointer(atypej2,"atypej2 in mr3calccoulomb_vdw_ij_exlist_host_");
  MR3_free_pointer(natex2,"natex2 in mr3calccoulomb_vdw_ij_exlist_host_");
}


void mr3calccoulomb_vdw_ij_exlist_host__(int *ni, double xi[], double qi[], 
					 int atypei[], double force[], 
					 int *nj, double xj[], double qj[], 
					 int atypej[],
					 int *nat,
					 double gscalesf[], 
					 double gscalesp[], double rscales[],
					 double *rscale, int *tblno,
					 double *xmax, double *potc, double *potv,
					 int *periodicflag, 
					 int *potflag, int *changeflag,
					 int numex[], int natex[])
{
  mr3calccoulomb_vdw_ij_exlist_host_(ni,xi,qi, 
				     atypei,force, 
				     nj,xj,qj, 
				     atypej,
				     nat,
				     gscalesf, 
				     gscalesp,rscales,
				     rscale,tblno,
				     xmax,potc,potv,
				     periodicflag, 
				     potflag,changeflag,
				     numex,natex);
}


/*
void MR3_get_forces_and_virial_overlap(double force[], double virial[3][3])
{
}


void mr3_get_forces_and_virial_overlap_(double force[], double virial[3][3])
{
  MR3_get_forces_and_virial_overlap(force,virial);
}
*/

