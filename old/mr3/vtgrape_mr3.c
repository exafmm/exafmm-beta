static int Overlapflag=0;
static double *F2compare=NULL,Potccompare=0.0,Potvcompare=0.0;
static int    Ncompare=0;


void MR3_get_cputime(double *laptime, double *sprittime)
{
  vg_get_cputime(laptime,sprittime);
}

static void setup_volume(double xmax, double volume[3])
{
  static int volumeenvflag[3]={0,0,0};
  static int ini=0;
  char *s;
  int i;

  if(ini==0){
    if((s=getenv("VG_CELLSIZEX"))!=NULL){
      sscanf(s,"%lf",&volume[0]);
      printf("VG_CELLSIZEX is set %f\n",volume[0]);
      volumeenvflag[0]=1;
    }
    if((s=getenv("VG_CELLSIZEY"))!=NULL){
      sscanf(s,"%lf",&volume[1]);
      printf("VG_CELLSIZEY is set %f\n",volume[1]);
      volumeenvflag[1]=1;
    }
    if((s=getenv("VG_CELLSIZEZ"))!=NULL){
      sscanf(s,"%lf",&volume[2]);
      printf("VG_CELLSIZEZ is set %f\n",volume[2]);
      volumeenvflag[2]=1;
    }
    ini=1;
  }
  for(i=0;i<3;i++){
    if(volumeenvflag[i]==0) volume[i]=xmax;
  }
}


void MR3_set_ldim(int ldim[3])
{
  int i;
  for(i=0;i<3;i++) Ldim[i]=ldim[i];
}


void mr3_set_ldim_(int ldim[3])
{
  MR3_set_ldim(ldim);
}


void mr3_set_ldim__(int ldim[3])
{
  MR3_set_ldim(ldim);
}


void MR3_set_coulomb_vdw_factor(double cgf)
{
  Coulomb_vdw_factor=cgf;
}


void mr3_set_coulomb_vdw_factor_(double *cgf)
{
  MR3_set_coulomb_vdw_factor(*cgf);
}


void mr3_set_coulomb_vdw_factor__(double *cgf)
{
  mr3_set_coulomb_vdw_factor_(cgf);
}


void MR3calccoulomb_ij_emu(int ni, double xid[], double qid[], double force[],
			   int nj, double xjd[], double qjd[],
			   double rscale, 
			   int tblno, double xmaxd, int periodicflag)
{
  /* periodic flag bit 0 : 0 --- non periodic, 1 --- periodic
                   bit 1 : 0 --- no multiplication of qi
                           1 --- multiplication of qi
  */
  int i,j,k;
  static int rcutflag=0,ini=0;
  static double rcut;
  static float rcut2;
  float r2min=MD_LJ_R2MIN,r2max;
  char *s;
  int multiplyq=0;
  float xi[3],xj[3],qi,qj,dr[3],r2,dtmp,x,factor,rsqrt,xmax=xmaxd;
  float rscale2=(float)(rscale*rscale);
  VG_UNIT *unit;
#ifdef COULOMB_SHIFT
  float rcut21;
#endif

  //  printf("xmaxd=%e\n",xmaxd);
  unit=m3_get_unit();
  //  r2max=MD_LJ_R2MAX;
  r2max=unit->rcut2;
#ifdef MD_PRINT_WARN
  {
    static int ini=0;
    if(ini==0){
      if(r2max!=MD_LJ_R2MAX) printf("r2max is modified in MR3calccoulomb_ij_emu\n",r2max);
      ini=1;
    }
  }
#endif
#ifdef COULOMB_SHIFT
  rcut21=1.0/unit->rcut2;
#endif
  if(unit->r1==NULL || unit->rsqrt==NULL) vg_initialize_emu();
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
    //    if(i<3) printf("xid[%d]=%e %e %e\n",i,xid[i*3],xid[i*3+1],xid[i*3+2]);
    for(k=0;k<3;k++) xi[k]=xid[i*3+k];
    qi=qid[i];
    for(j=0;j<nj;j++){
      //      if(j<3) printf("xjd[%d]=%e %e %e\n",j,xjd[j*3],xjd[j*3+1],xjd[j*3+2]);
      for(k=0;k<3;k++) xj[k]=xjd[j*3+k];
      qj=qjd[j];
      r2=0.0;
      for(k=0;k<3;k++){
	dr[k]=xi[k]-xj[k];
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
      x=r2*rscale2;
#if MD_LJ_05SIGMA_8SIGMA==2
      if(r2>=r2min && r2<r2max){
#else
      if(r2!=0.0 && (rcutflag==0 || (rcutflag==1 && x<rcut2))){
#endif
	if(tblno==0){
	  rsqrt=vg_rsqrtemu(r2);
	  //	  rsqrt=rsqrtf(r2);
	  //rsqrt=1.0f/sqrtf(r2);
#ifdef COULOMB_SHIFT
	  dtmp=qj*rsqrt;
	  x=1.0f-r2*rcut21;
	  dtmp=dtmp*x*x*rsqrt*rsqrt+4.0f*rcut21*dtmp*x;
#else
	  dtmp=qj*rsqrt*rsqrt*rsqrt;
#endif
	  if(multiplyq) dtmp*=qi*Coulomb_vdw_factor;
	  for(k=0;k<3;k++){
	    force[i*3+k]+=dtmp*dr[k];
	  }
#if 0
	  if(i<3){
	    printf("host i=%d j=%d r2=%e rrs=%e xi=%e %e %e xj=%e %e %e\n",i,j,r2,x,xi[0],xi[1],xi[2],xj[0],xj[1],xj[2]);
	    printf("     xmax=%e dr=%e %e %e dtmp=%e rsqrt=%e\n",xmax,dr[0],dr[1],dr[2],dtmp,rsqrt);
	    printf("     factor=%e dtmp=%e dtmp*factor=%e\n",factor,dtmp,dtmp*factor);
	    printf("     force=%e %e %e\n",force[i*3],force[i*3+1],force[i*3+2]);
	  }
#endif
	}
	else if(tblno==1){
	  rsqrt=vg_rsqrtemu(r2);
	  //	  rsqrt=rsqrtf(r2);
	  //rsqrt=1.0f/sqrtf(r2);
	  dtmp=qj*rsqrt;
#ifdef COULOMB_SHIFT
	  x=1.0f-r2*rcut21;
	  dtmp*=x*x;
#endif
	  if(multiplyq) dtmp*=qi*Coulomb_vdw_factor;
	  for(k=0;k<3;k++){
	    force[i*3+k]+=dtmp;
	  }
	}
#if 0
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
#endif
	else{
	  fprintf(stderr,"** error : tblno=%d is not supported in MR3calccoulomb_ij_emu\n",tblno);
	  vg_exit(1);
	}
      }
    }
  }
}


void MR3calccoulomb_nlist_ij_emu(int ni, double xid[], double qid[], double force[],
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
  int i,j,k,jj,n;
  static int rcutflag=0,ini=0;
  static double rcut;
  static float rcut2;
  float r2min=MD_LJ_R2MIN,r2max;
  char *s;
  int multiplyq=0,thirdlaw=1,iexcl=0;
  float xi[3],xj[3],qi,qj,dr[3],r2,dtmp,x,rsqrt,xmax=xmaxd,factor=factord,dtmp0;
  float rscale2=(float)(rscale*rscale),m_2_sqrtpi=(float)M_2_SQRTPI;
  VG_UNIT *unit;
  double *ftmp,*fmixed;
#ifdef COULOMB_SHIFT
  float rcut21;
#endif

  //  printf("in MR3calccoulomb_nlist_ij_emu : ni=%d nj=%d rscale=%e tblno=%d xmax=%e periodicflag=%d Coulomb_vdw_factor=%e\n",ni,nj,rscale,tblno,xmaxd,periodicflag,Coulomb_vdw_factor);
  if(Coulomb_vdw_factor!=1.0){
    n=(ni>nj ? ni:nj);
    if((ftmp=(double *)MR3_malloc_pointer(sizeof(double)*n*3,"MR3calccoulomb_nlist_ij_emu"))==NULL){
      fprintf(stderr,"** error : can't malloc ftmp **\n");
      vg_exit(1);
    }
    bzero(ftmp,sizeof(double)*n*3);
    fmixed=ftmp;
  }
  else{
    fmixed=force;
  }
  unit=m3_get_unit();
  //  r2max=MD_LJ_R2MAX;
  r2max=unit->rcut2;
#ifdef MD_PRINT_WARN
  {
    static int ini=0;
    if(ini==0){
      if(r2max!=MD_LJ_R2MAX) printf("r2max is modified in MR3calccoulomb_nlist_ij_emu\n",r2max);
      ini=1;
    }
  }
#endif
#ifdef COULOMB_SHIFT
  rcut21=1.0/unit->rcut2;
#endif
  if(unit->r1==NULL || unit->rsqrt==NULL) vg_initialize_emu();
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
	//	if(i==53) printf("    dr[%d]=%20.14e xi=%20.14e xj=%20.14e xi-xj=%20.14e\n",k,dr[k],xi[k],xj[k],xi[k]-xj[k]);
	if(dr[k]<-xmax/2.0f){
	  dr[k]+=xmax;
	}
	if(dr[k]>=xmax/2.0f){
	  dr[k]-=xmax;
	}
	//	if(i==53) printf("    after dr[%d]=%e\n",k,dr[k]);
	r2+=dr[k]*dr[k];
      }
      x=r2*rscale2;
      //      printf("  r2=%e\n",r2);
#if MD_LJ_05SIGMA_8SIGMA==2
      if(r2>=r2min && r2<r2max){
#else
      if(r2!=0.0 && (rcutflag==0 || (rcutflag==1 && x<rcut2))){
#endif
	if(tblno==0){
	  rsqrt=vg_rsqrtemu(r2);
	  //	  rsqrt=rsqrtf(r2);
	  //rsqrt=1.0f/sqrtf(r2);
#ifdef COULOMB_SHIFT
	  dtmp=qj*rsqrt;
	  x=1.0f-r2*rcut21;
	  dtmp=dtmp*x*x*rsqrt*rsqrt+4.0f*rcut21*dtmp*x;
#else
	  dtmp=qj*rsqrt*rsqrt*rsqrt;
#endif
	  dtmp*=factor;
	  if(multiplyq) dtmp*=qi;
	  //	  if(1) dtmp*=qi*Coulomb_vdw_factor;
	  for(k=0;k<3;k++){
	    fmixed[i*3+k]+=dtmp*dr[k];
	    if(thirdlaw){
	      if(multiplyq) fmixed[j*3+k]-=dtmp*dr[k];
	      else          fmixed[j*3+k]-=qi*dtmp0*dr[k];
	    }
	  }
	  //	  printf("  i=%d j=%d f=%e %e %e\n",i,j,dtmp*dr[0],dtmp*dr[1],dtmp*dr[2]);
	}
	else if(tblno==1){
	  rsqrt=vg_rsqrtemu(r2);
	  //	  rsqrt=rsqrtf(r2);
	  //rsqrt=1.0f/sqrtf(r2);
	  dtmp0=rsqrt;
	  //	  dtmp0*=factor;
	  dtmp=qj*dtmp0;
#ifdef COULOMB_SHIFT
	  x=1.0f-r2*rcut21;
	  dtmp*=x*x;
#endif
	  dtmp*=factor;
	  if(multiplyq) dtmp*=qi;
	  //	  if(1) dtmp*=qi*Coulomb_vdw_factor;
	  //	  for(k=0;k<3;k++){
	  for(k=0;k<1;k++){
	    fmixed[i*3+k]+=dtmp;
	    if(thirdlaw){
	      if(multiplyq) fmixed[j*3+k]+=dtmp;
	      else          fmixed[j*3+k]+=qi*dtmp0;
	    }
	  }
	}
	else if(tblno==6){
	  x=r2*rscale2;
	  dtmp=qj*(m_2_sqrtpi*expf(-x)*powf(x,-1.0)
	              + erfcf(sqrtf(x))*powf(x,-1.5));
	  dtmp*=factor*rscale2*rscale;
	  if(1) dtmp*=qi;
	  for(k=0;k<3;k++){
	    fmixed[i*3+k]+=dtmp*dr[k];
	    if(thirdlaw) fmixed[j*3+k]-=dtmp*dr[k];
	  }
	  /*
	  if(i==53){
	    printf("i=%d j=%d dtmp=%e fmixed[i]=%e %e %e j=%e %e %e\n",i,j,dtmp,fmixed[i*3],fmixed[i*3+1],fmixed[i*3+2],fmixed[j*3],fmixed[j*3+1],fmixed[j*3+2]);
	    printf("  r2=%e rscale2=%e qj=%e x=%e factor=%e qi=%e\n",
		   r2,rscale2,qj,x,factor,qi);
	    printf("  xi=%e %e %e xj=%e %e %e\n",xi[0],xi[1],xi[2],xj[0],xj[1],xj[2]);
	    printf("  dr=%e %e %e\n",dr[0],dr[1],dr[2]);
	  }
	  */
	}
	else if(tblno==7){
	  x=r2*rscale2;
	  dtmp=qj*erfcf(sqrtf(x))*powf(x,-0.5);
	  dtmp*=factor*rscale;
	  if(1) dtmp*=qi;
	  //	  for(k=0;k<3;k++){
	  for(k=0;k<1;k++){
	    fmixed[i*3+k]+=dtmp;
	    if(thirdlaw) fmixed[j*3+k]+=dtmp;
	  }
	}
	else{
	  fprintf(stderr,"** error : tblno=%d is not supported in MR3calccoulomb_nlist_ij_emu\n",tblno);
	  vg_exit(1);
	}
      }
    }
    iexcl = iexcl + numex[i];
  }
  if(Coulomb_vdw_factor!=1.0){
    for(i=0;i<ni;i++){
      for(j=0;j<3;j++) force[i*3+j]+=ftmp[i*3+j]*Coulomb_vdw_factor;
    }
    MR3_free_pointer(ftmp,"MR3calccoulomb_nlist_ij_emu");
  }
}


void MR3calccoulomb_nlist_emu(double xd[], int n, double qd[], double rscale,
			      int tblno, double xmaxd, int periodicflag,
			      int numex[], int natex[], double factord,
			      double force[])
{
  MR3calccoulomb_nlist_ij_emu(n,xd,qd,force,n,xd,qd,rscale,tblno,xmaxd,periodicflag,
			      numex,natex,factord);
}


static float vdw_shift_emu(float r2, 
			   float rs, float gs, float rcut21_org, int potflag)
{
  /* potflag 0 --- force
             1 --- potential
  */
  float SIGFAC=pow(2.0,1.0/3.0);
  float A,B;
  static float CACX,CBCX,CADX,CBDX,CCDX;
  float CCA,CCB,CCC,CCD,CA,CC;
  float ENEVDW,ENN,TFVDW;
  static float rcut21=0.0;

  if(rcut21!=rcut21_org){
    rcut21=rcut21_org;
    CACX=-2.0*pow((double)rcut21,9);
    CBCX=pow((double)rcut21,6);
    CADX=-pow((double)rcut21,6);
    CBDX=pow((double)rcut21,3);
    CCDX=pow((double)rcut21,-3);
  }
  A=SIGFAC/rs;
  if(potflag) B=gs/4.0f;
  else        B=gs*A/(24.0f*SIGFAC);
  CCA=B*powf(A,6.0f);
  CCB=2.0f*B*powf(A,3.0f);
  CCC=CACX*CCA+CBCX*CCB;
  CCD=CADX*CCA+CBDX*CCB+CCDX*CCC;
  CA=CCA*powf(r2,-6.0f);
  CC=powf(r2,3.0f)*CCC;
  ENEVDW=CA-CCB*powf(r2,-3.0f);
  ENN=ENEVDW-CC+CCD;
  if(potflag){
    return ENN;
  }
  else{
    TFVDW=-6.0f*(ENEVDW+CA+CC)/r2;
    return TFVDW;
  }
}


void MR3calcvdw_ij_emu(int ni, double xid[], int atypei[], double force[],
			int nj, double xjd[], int atypej[],
			int nat, double gscale[], double rscale[],
			int tblno, double xmax, int periodicflag) __attribute__((optimize(0)));

#pragma optimize("",off)
void MR3calcvdw_ij_emu(int ni, double xid[], int atypei[], double force[],
			int nj, double xjd[], int atypej[],
			int nat, double gscale[], double rscale[],
			int tblno, double xmax, int periodicflag)
{
  /* periodic flag 0 --- non periodic
                   1 --- periodic
  */
  int i,j,k;
  float r2min=MD_LJ_R2MIN,r2max;
  float dtmp,r1,r6,r2,rs,gs,rrs,dr[3],xi[3],xj[3];
#if VDW_SHIFT>=1
  float rcut21,shift;
  VG_SCALER *scal;
#endif
  VG_UNIT *unit;

  unit=m3_get_unit();
  //  r2max=MD_LJ_R2MAX;
  r2max=unit->rcut2;
#ifdef MD_PRINT_WARN
  {
    static int ini=0;
    if(ini==0){
      if(r2max!=MD_LJ_R2MAX) printf("r2max is modified in MR3calcvdw_ij_emu\n",r2max);
      ini=1;
    }
  }
#endif
#if VDW_SHIFT>=1
  scal=unit->scalers;
  rcut21=scal->rcut21;
#endif
  if(unit->r1==NULL || unit->rsqrt==NULL) vg_initialize_emu();
  if((periodicflag & 1)==0){
    xmax*=2.0;
  }
  for(i=0;i<ni;i++){
    for(k=0;k<3;k++) xi[k]=xid[i*3+k];
    //    printf("i=%d xi=%e %e %e\n",i,xi[0],xi[1],xi[2]);
    for(j=0;j<nj;j++){
      r2=0.0;
      for(k=0;k<3;k++) xj[k]=xjd[j*3+k];
      //      printf(" j=%d xj=%e %e %e\n",j,xj[0],xj[1],xj[2]);
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
	rs=rscale[atypei[i]*nat+atypej[j]];
	gs=gscale[atypei[i]*nat+atypej[j]];
#if VDW_SHIFT==1
	shift=(float)(pow(rscale[atypei[i]*nat+atypej[j]]*(unit->rcut2),-3.0));
#endif
	rrs=r2*rs;
	//	printf("  i=%d j=%d xi=%e %e %e r2=%e rrs=%e\n",i,j,xi[0],xi[1],xi[2],r2,rrs);
#if MD_LJ_05SIGMA_8SIGMA==2
	if(gs!=0.0f && r2>=r2min && r2<r2max){
#else
	if(rrs>=r2min && rrs<r2max){
#endif
	  r1=vg_r1emu(rrs);
	  //	  printf("   i=%d j=%d r1=%e\n",i,j,r1);
	  //r1=1.0f/rrs;
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
	    dtmp=gs*r6*r1*(2.0f*r6-1.0f);
#endif
	    for(k=0;k<3;k++){
	      force[i*3+k]+=dtmp*dr[k];
	    }
	  }
	  else if(tblno==3){
#if VDW_SHIFT==1
#if 1
	    dtmp=r6*(r6-1.0f);
#if 0
	    //	    r1=powf(rs*rcut2,-3.0f);
	    r1=shift;
	    r2*=r2*rs;
	    //	    r2*=r2;
	    dtmp+=r1*r1*(r1-1.0f)*r2+2.0f*r1*(1.0f-r1);
#endif
	    r2*=rcut21;
	    r2*=r2*r2;
	    dtmp+=r2*shift*(2.0f*shift-1.0f)+shift*(2.0f-3.0f*shift);
	    dtmp*=gs;
	    //	    dtmp=r6*10000.0f;
#else
	    {
	      VG_SCALER *scal=unit->scalers;
	      float rcut21=scal->rcut21;
	      dtmp=vdw_shift_emu(r2,rs,gs,rcut21,1);
	    }
#endif
#else
	    dtmp=gs*r6*(r6-1.0f);
	    //	    dtmp=r2;
#endif
	    for(k=0;k<3;k++){
	      force[i*3+k]+=dtmp;
	    }
	  }
	  else{
	    fprintf(stderr,"** error : tblno=%d is not supported in MR3calcvdw_ij_emu\n",tblno);
	    vg_exit(1);
	  }
	}
      }
    }
  }
}
#pragma optimize("",on)


static __inline__ float madd_emu(float a, float x, float y)
  //static float madd_emu(float a, float x, float y)
{
  /*
     c = a * x + y
   */
#if 0
  unsigned long long ull;
  double ddr;
  ddr=a;
  //  printf("   madd_emu : a=%016llx x=%016llx",*((unsigned long long *)&ddr),*((unsigned long long *)&x));
  ddr*=(double)x;
  //  printf(" ax=%016llx",*((unsigned long long *)&ddr));
  ull=*((unsigned long long *)&ddr);
  ull&=0xffffffffe0000000LL;
  ddr=*((double *)&ull);
  return (y+(float)ddr);
#else
#ifdef MD_PRINT_WARN
  static int ini=0;
  if(ini==0){
    fprintf(stderr,"** madd_emu does not use FMAD fusion emulation **\n");
    ini=1;
  }
#endif
  return (a*x+y);
#endif
}


float calc_dr_r2(
#if MD_PERIODIC_FIXED==1
                 int xi[3], int xj[3], float al2[3],
#elif MD_PERIODIC_FIXED==2
                 float xi[3], float xj[3], float volumef[3],
#else
                 float xi[3], float xj[3],
#endif
                 float dr[]) __attribute__((optimize(0)));

#pragma optimize("",off)
__inline__ float calc_dr_r2(
#if MD_PERIODIC_FIXED==1
                 int xi[3], int xj[3], float al2[3],
#elif MD_PERIODIC_FIXED==2
                 float xi[3], float xj[3], float volumef[3],
#else
                 float xi[3], float xj[3],
#endif
                 float dr[])
{
  float r2=0.0;
  int k;

  for(k=0;k<3;k++){
    dr[k]=(float)(xi[k]-xj[k]);
#if MD_PERIODIC_FIXED==1
    dr[k]*=al2[k];
#elif MD_PERIODIC_FIXED==2
    if(dr[k]<-volumef[k]/2.0f){
      dr[k]+=volumef[k];
    }
    if(dr[k]>=volumef[k]/2.0f){
      dr[k]-=volumef[k];
    }
#endif
#if 0 // round off error occurrs because of MADD accuracy is different from GPU
    r2+=dr[k]*dr[k];
#endif
  }
#if 1      
  r2=dr[0]*dr[0];
  r2=madd_emu(dr[1],dr[1],r2);
  r2=madd_emu(dr[2],dr[2],r2);
#endif
  return r2;
}

__inline__ 
float calc_dtmp_vdw_force(float r1, float r2, float r6, 
#if VDW_SHIFT>=1
			  float shift, float rcut21, 
#endif
			  float gs, float factorv)
{
  float dtmp;

#if VDW_SHIFT==1
  dtmp=r1;
  r1=shift;
#ifdef VDW_SHIFT_CORRECT
  dtmp*=r6*(2.0f*r6-1.0f)-r1*(3.0f*r1-2.0f);
#else
  r2*=rcut21;
  r2*=r2*r2;
  dtmp*=r6*(2.0f*r6-1.0f)-r2*r1*(2.0f*r1-1.0f);
#endif
  dtmp*=gs;
#else
  dtmp=gs*r6*r1*(2.0f*r6-1.0f);
#endif
  dtmp*=factorv;
  
  return dtmp;
}


#pragma optimize("",on)


static void MR3calccoulomb_vdw_nlist_ij_emu_work(int ni, double xid[], double qid[], int atypei[], double force[],
						 int nj, double xjd[], double qjd[], int atypej[],
						 int nat, double gscales[], double rscales[],
						 double rscale, int tblno, double volume[], int periodicflag,
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
  int i,j,k,jj,n;
  static int rcutflag=0,ini=0;
  static double rcut;
  static float rcut2;
  float r2min=MD_LJ_R2MIN,r2max;
  char *s;
  int multiplyq=0,thirdlaw=1,iexcl=0;
  float qi,qj,dr[3],r2,dtmp,x,rsqrt,volumef[3]={volume[0],volume[1],volume[2]},dtmp0;
  float rscale2=(float)(rscale*rscale),m_2_sqrtpi=(float)M_2_SQRTPI;
  VG_UNIT *unit;
  double *ftmp,*fmixed;
  float factorc,factorv,rs,gs,rrs,r1,r6;
  int tblnoc=tblno,tblnov=2+(tblno % 2);
#if MD_PERIODIC_FIXED==1
  int xi[3],xj[3];
  float al2[3];
  double volume_1[3];
  DI2 di2;
#else
  float xi[3],xj[3];
#endif
#if VDW_SHIFT>=1
  float rcut21,shift;
  VG_SCALER *scal;
#endif
#ifdef COULOMB_SHIFT
  float rcut21c;
#endif

#ifdef MD_MEASURE_TIME
  if(tblno<0 && tblno>=10){
    fprintf(stderr,"** error : tblno is outof range for MD_MEASURE_TIME **\n");
    vg_exit(1);
  }
  vg_start_timer(20+tblno);
#endif
#if 0
  {
    double sum=0.0;
    for(i=0;i<ni;i++) sum+=force[i*3];
    printf("before MR3calccoulomb_vdw_nlist_ij_emu_work is called: sum=%24.17e\n",sum);
  }
#endif
  //  printf("periocicflag=%d\n",periodicflag);
  fmixed=force;
  factorc=(float)(factord*Coulomb_vdw_factor);
  factorv=(float)factord;
  unit=m3_get_unit();
  //  r2max=MD_LJ_R2MAX;
  r2max=unit->rcut2;
#ifdef MD_PRINT_WARN
  {
    static int ini=0;
    if(ini==0){
      if(r2max!=MD_LJ_R2MAX) printf("r2max is modified in MR3calccoulomb_vdw_nlist_ij_emu_work\n",r2max);
      ini=1;
    }
  }
#endif
#ifdef COULOMB_SHIFT
  rcut21c=1.0/unit->rcut2;
#endif
#if VDW_SHIFT>=1
  scal=unit->scalers;
  rcut21=scal->rcut21;
#endif
  if(unit->r1==NULL || unit->rsqrt==NULL) vg_initialize_emu();
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
  if((periodicflag & 1)==0){
    volumef[0]*=volumef[1]*=volumef[2]*=2.0;
  }
#if MD_PERIODIC_FIXED==1
  for(i=0;i<3;i++) volume_1[i]=1.0/volume[i];
  for(i=0;i<3;i++) al2[i]=scalbnf(volumef[i],-32);
#endif
  //#pragma omp parallel for private(k,xi,qi,jj,j,qj,r2,dr,x,rsqrt,dtmp,dtmp0,rs,gs,shift,rrs,r1,r6)
  for(i=0;i<ni;i++){
    for(k=0;k<3;k++){
      COPY_POS_INT(xi[k],xid[i*3+k],volume_1[k]);
    }
    qi=qid[i];
    for(jj=iexcl;jj<iexcl+numex[i];jj++){
      j=natex[jj];
      if(j<0) continue;
      for(k=0;k<3;k++){
	COPY_POS_INT(xj[k],xjd[j*3+k],volume_1[k]);
      }
      qj=qjd[j];
#if 1 // separate funciton
#if MD_PERIODIC_FIXED==1
      r2=calc_dr_r2(xi,xj,al2,dr);
#elif MD_PERIODIC_FIXED==2
      r2=calc_dr_r2(xi,xj,volumef,dr);
#else
      r2=calc_dr_r2(xi,xj,dr);
#endif
#else // original
      r2=0.0;
      for(k=0;k<3;k++){
	dr[k]=(float)(xi[k]-xj[k]);
#if MD_PERIODIC_FIXED==1
	dr[k]*=al2[k];
#elif MD_PERIODIC_FIXED==2
	if(dr[k]<-volumef[k]/2.0f){
	  dr[k]+=volumef[k];
	}
	if(dr[k]>=volumef[k]/2.0f){
	  dr[k]-=volumef[k];
	}
#endif
#if 0 // round off error occurrs because of MADD accuracy is different from GPU
	r2+=dr[k]*dr[k];
#endif
      }
#if 1      
      r2=dr[0]*dr[0];
      r2=madd_emu(dr[1],dr[1],r2);
      r2=madd_emu(dr[2],dr[2],r2);
#endif
#endif

      // for coulomb and real space
      x=r2*rscale2;
#if MD_LJ_05SIGMA_8SIGMA==2
      if(r2>=r2min && r2<r2max){
#else
      if(r2!=0.0 && (rcutflag==0 || (rcutflag==1 && x<rcut2))){
#endif
	if(tblnoc==0){
	  rsqrt=vg_rsqrtemu(r2);
#ifdef COULOMB_SHIFT
	  dtmp=qj*rsqrt;
	  x=1.0f-r2*rcut21c;
	  dtmp=dtmp*x*x*rsqrt*rsqrt+4.0f*rcut21c*dtmp*x;
#else
	  dtmp=qj*rsqrt*rsqrt*rsqrt;
#endif
	  dtmp*=factorc;
	  if(multiplyq) dtmp*=qi;
	  for(k=0;k<3;k++){
	    fmixed[i*3+k]+=dtmp*dr[k];
	    if(thirdlaw){
	      if(multiplyq) fmixed[j*3+k]-=dtmp*dr[k];
	      else          fmixed[j*3+k]-=qi*dtmp0*dr[k];
	      // dtmp0 seems not initialized
	    }
	  }
	}
	else if(tblnoc==1){
	  rsqrt=vg_rsqrtemu(r2);
	  dtmp0=rsqrt;
	  dtmp=qj*dtmp0;
#ifdef COULOMB_SHIFT
	  x=1.0f-r2*rcut21c;
	  dtmp*=x*x;
#endif
	  dtmp*=factorc;
	  if(multiplyq) dtmp*=qi;
          for(k=0;k<1;k++){
	    fmixed[i*3+k]+=dtmp;
	    if(thirdlaw){
	      if(multiplyq) fmixed[j*3+k]+=dtmp;
	      else          fmixed[j*3+k]+=qi*dtmp0;
	    }
	  }
	}
	else if(tblnoc==6){
	  x=r2*rscale2;
	  dtmp=qj*(m_2_sqrtpi*expf(-x)*powf(x,-1.0)
	              + erfcf(sqrtf(x))*powf(x,-1.5));
	  dtmp*=factorc*rscale2*rscale;
	  if(1) dtmp*=qi;
	  for(k=0;k<3;k++){
	    fmixed[i*3+k]+=dtmp*dr[k];
	    if(thirdlaw) fmixed[j*3+k]-=dtmp*dr[k];
	  }
	}
	else if(tblnoc==7){
	  x=r2*rscale2;
	  dtmp=qj*erfcf(sqrtf(x))*powf(x,-0.5);
	  dtmp*=factorc*rscale;
	  if(1) dtmp*=qi;
          for(k=0;k<1;k++){
	    fmixed[i*3+k]+=dtmp;
	    if(thirdlaw) fmixed[j*3+k]+=dtmp;
	  }
	}
      }

      // for vdw
      if(r2!=0.0){
	rs=rscales[atypei[i]*nat+atypej[j]];
	gs=gscales[atypei[i]*nat+atypej[j]];
#if VDW_SHIFT==1
	shift=(float)(pow(rscales[atypei[i]*nat+atypej[j]]*(unit->rcut2),-3.0));
#endif
	rrs=r2*rs;
#if MD_LJ_05SIGMA_8SIGMA==2
	if(gs!=0.0f && r2>=r2min && r2<r2max){
#else
	if(gs!=0.0f && rrs>=r2min && rrs<r2max){
#endif
	  r1=vg_r1emu(rrs);
	  //r1=1.0f/rrs;
	  r6=r1*r1*r1;
	  if(tblnov==2){
#if 1 // use separate routine
            dtmp=calc_dtmp_vdw_force(r1,r2,r6,
#if VDW_SHIFT>=1
				     shift,rcut21,
#endif
				     gs,factorv);
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
	    dtmp*=gs;
#else
	    dtmp=gs*r6*r1*(2.0f*r6-1.0f);
#endif
	    dtmp*=factorv;
#endif
#if 1 // work
	    for(k=0;k<3;k++){
	      force[i*3+k]+=dtmp*dr[k];
	      if(thirdlaw) force[j*3+k]-=dtmp*dr[k];
	    }
#else // debug
	    //	    force[i*3]+=rrs;
	    {
	      float ftmp;
	      //	      ftmp=dr[0]*dr[0];
	      force[i*3]=madd_emu(dr[0],dr[0],force[i*3]);
	      force[i*3]=madd_emu(dr[1],dr[1],force[i*3]);
	      //	      ftmp=dr[0]*dr[0]+dr[1]*dr[1];
	      //	      force[i*3]+=ftmp;
	    }
#endif
	  }
	  else if(tblnov==3){
#if VDW_SHIFT==1
#if 1
	    dtmp=r6*(r6-1.0f);
#if 0
	    //	    r1=powf(rs*rcut2,-3.0f);
	    r1=shift;
	    r2*=r2*rs;
	    //	    r2*=r2;
	    dtmp+=r1*r1*(r1-1.0f)*r2+2.0f*r1*(1.0f-r1);
#endif
	    r2*=rcut21;
	    r2*=r2*r2;
	    dtmp+=r2*shift*(2.0f*shift-1.0f)+shift*(2.0f-3.0f*shift);
	    dtmp*=gs;
#else
	    {
	      VG_SCALER *scal=unit->scalers;
	      float rcut21=scal->rcut21;
	      dtmp=vdw_shift_emu(r2,rs,gs,rcut21,1);
	    }
#endif
#else
	    dtmp=gs*r6*(r6-1.0f);
	    //	    dtmp=gs*(pow(rrs,-6.0)-pow(rrs,-3.0));
#endif
	    dtmp*=factorv;
	    //	    for(k=0;k<3;k++){
	    for(k=1;k<2;k++){
              //	    for(k=0;k<2;k++){
              //	    for(k=0;k<1;k++){
	      force[i*3+k]+=dtmp;
	      if(thirdlaw) force[j*3+k]+=dtmp;
	    }
	  }
	}
      }
    }
    iexcl = iexcl + numex[i];
  }
#if 0
  {
    double sum=0.0;
    for(i=0;i<ni;i++) sum+=force[i*3];
    printf("after MR3calccoulomb_vdw_nlist_ij_emu_work is called: sum=%24.17e\n",sum);
  }
#endif
#ifdef MD_MEASURE_TIME
  vg_stop_and_accumulate_timer(20+tblno);
#endif
}


#ifdef VG_GCC
static v4sf erfcfv(v4sf xv)
{
  v4sf retv;
  int i;
  
  for(i=0;i<4;i++){
    ((float *)&retv)[i]=erfcf(((float *)&xv)[i]);
  }

  return retv;
}


static v4sf expfv(v4sf xv)
{
  v4sf retv;
  int i;
  
  for(i=0;i<4;i++){
    ((float *)&retv)[i]=expf(((float *)&xv)[i]);
  }

  return retv;
}


static v4sf rintfv(v4sf xv)
{
  v4sf retv;
  int i;
  
  for(i=0;i<4;i++){
    ((float *)&retv)[i]=rintf(((float *)&xv)[i]);
  }

  return retv;
}


static v4sf powfv(v4sf xv, float p)
{
  v4sf retv;
  int i;
  
  for(i=0;i<4;i++){
    ((float *)&retv)[i]=powf(((float *)&xv)[i],p);
  }

  return retv;
}


static v4sf sqrtfv(v4sf xv)
{
  v4sf retv;
  int i;
  
  for(i=0;i<4;i++){
    ((float *)&retv)[i]=sqrtf(((float *)&xv)[i]);
  }

  return retv;
}


static void MR3calccoulomb_vdw_nlist_ij_emu_vec_sub(v4sf fiv[3], v4sf fjv[3],
#if MD_PERIODIC_FIXED==1
						    v4si xiv[3], v4si xjv[3], float al2[],
#else
						    v4sf xiv[3], v4sf xjv[3],
						    //float xia[3][4], v4sf xja[3][4],
#endif
						    v4sf qiv, v4sf qjv, v4sf gsv, v4sf rsv, v4sf shv,
						    //float qia[], float qja[], float gsa[], float rsa[], float sha[],
						    float al[], float l[],
						    float r2max, float rcut21, float rcut21c, float factorc, float factorv,
						    int nat, float rscale, float rscale2, int tblnoc, int tblnov, float volume[],
						    int multiplyq, int thirdlaw)
{
  float r2,dr[3],x,dtmp,rsqrt,qi,qj,gs,rs,shift,rrs,dtmp0;
  float r1,r6;
  float r2min=MD_LJ_R2MIN;
  float m_2_sqrtpi=(float)M_2_SQRTPI;
  int k,v;
  float zv[VG_EMU_VECSIZE];
#if MD_PERIODIC_FIXED==1
  int xi[3],xj[3];
#else
  float xi[3],xj[3];
#endif

  for(v=0;v<VG_EMU_VECSIZE;v++) zv[v]=0.0f;
  for(k=0;k<3;k++) fiv[k]=fjv[k]=*(v4sf *)zv;

  if(tblnoc==6 && tblnov==2 && thirdlaw==1){
    v4sf drv[3],r2v,al2v[3],rscale2v,xv,rrsv,r2zerov,r1v,r6v;
    v4sf rcut21v,dtmpv,m_2_sqrtpiv,rscalev,factorcv,factorvv;
    v4sf alv[3],lv[3],twov={2.0f,2.0f,2.0f,2.0f},onev={1.0f,1.0f,1.0f,1.0f};
    float dum[4];
    dum[0]=rscale2;rscale2v=__builtin_ia32_shufps(*(v4sf *)dum,*(v4sf *)dum, 0x00);
    dum[0]=rscale;rscalev=__builtin_ia32_shufps(*(v4sf *)dum,*(v4sf *)dum, 0x00);
    dum[0]=factorc;factorcv=__builtin_ia32_shufps(*(v4sf *)dum,*(v4sf *)dum, 0x00);
    dum[0]=rcut21;rcut21v=__builtin_ia32_shufps(*(v4sf *)dum,*(v4sf *)dum, 0x00);
    dum[0]=factorv;factorvv=__builtin_ia32_shufps(*(v4sf *)dum,*(v4sf *)dum, 0x00);
    dum[0]=m_2_sqrtpi;m_2_sqrtpiv=__builtin_ia32_shufps(*(v4sf *)dum,*(v4sf *)dum, 0x00);
#if MD_PERIODIC_FIXED==1
    al2v[0]=__builtin_ia32_shufps(*(v4sf *)al2,*(v4sf *)al2, 0x00);
    al2v[1]=__builtin_ia32_shufps(*(v4sf *)al2,*(v4sf *)al2, 0x55);
    al2v[2]=__builtin_ia32_shufps(*(v4sf *)al2,*(v4sf *)al2, 0xaa);
#else
    alv[0]=__builtin_ia32_shufps(*(v4sf *)al,*(v4sf *)al, 0x00);
    alv[1]=__builtin_ia32_shufps(*(v4sf *)al,*(v4sf *)al, 0x55);
    alv[2]=__builtin_ia32_shufps(*(v4sf *)al,*(v4sf *)al, 0xaa);
    lv[0]=__builtin_ia32_shufps(*(v4sf *)l,*(v4sf *)l, 0x00);
    lv[1]=__builtin_ia32_shufps(*(v4sf *)l,*(v4sf *)l, 0x55);
    lv[2]=__builtin_ia32_shufps(*(v4sf *)l,*(v4sf *)l, 0xaa);
#endif
    r2v=*(v4sf *)zv;
    for(k=0;k<3;k++){
      drv[k]=xiv[k]-xjv[k];
#if MD_PERIODIC_FIXED==1
      drv[k]*=al2v[k];
#elif MD_PERIODIC_FIXED==2
      drv[k]-=rintfv(drv[k]*alv[k])*lv[k];
#endif
    }
    r2v=drv[0]*drv[0]+drv[1]*drv[1]+drv[2]*drv[2];
    /* is r2 zero? */
    for(k=0;k<4;k++){
      r2=((float *)&r2v)[k];
      if(r2<r2min || r2>=r2max) dum[k]=0.0f;
      else                      dum[k]=1.0f;
    }
    r2zerov=*(v4sf *)dum;
    /* real space */
    xv=r2v*rscale2v;
    dtmpv=qjv*(m_2_sqrtpiv*expfv(-xv)*powfv(xv,-1.0)
		 + erfcfv(sqrtfv(xv))*powfv(xv,-1.5));
    dtmpv*=factorcv*rscale2v*rscalev;
    dtmpv*=qiv;
    dtmpv*=r2zerov;
    for(k=0;k<3;k++) fiv[k]+=dtmpv*drv[k];
    for(k=0;k<3;k++) fjv[k]-=dtmpv*drv[k];
    /* vdw */    
    rrsv=r2v*rsv;
    r1v=vg_r1emu_vec(rrsv);
    r6v=r1v*r1v*r1v;
    dtmpv=r1v;
    r1v=shv;
    r2v*=rcut21v;
    r2v*=r2v*r2v;
    dtmpv*=r6v*(twov*r6v-onev)-r2v*r1v*(twov*r1v-onev);
    dtmpv*=gsv;
    dtmpv*=factorvv;
    dtmpv*=r2zerov;
    for(k=0;k<3;k++) fiv[k]+=dtmpv*drv[k];
    for(k=0;k<3;k++) fjv[k]-=dtmpv*drv[k];

    return;
  }
  for(v=0;v<VG_EMU_VECSIZE;v++){
    for(k=0;k<3;k++) xi[k]=((float *)(xiv+k))[v];
    for(k=0;k<3;k++) xj[k]=((float *)(xjv+k))[v];
    qi=((float *)&qiv)[v];
    qj=((float *)&qjv)[v];
    r2=0.0;
    for(k=0;k<3;k++){
      dr[k]=(float)(xi[k]-xj[k]);
#if MD_PERIODIC_FIXED==1
      dr[k]*=al2[k];
#elif MD_PERIODIC_FIXED==2
      dr[k]-=rintf(dr[k]*al[k])*l[k];
#endif
    }
    r2=dr[0]*dr[0];
    r2=madd_emu(dr[1],dr[1],r2);
    r2=madd_emu(dr[2],dr[2],r2);
    
    // for coulomb and real space
    x=r2*rscale2;
    if(r2>=r2min && r2<r2max){
      if(tblnoc==0){
	rsqrt=vg_rsqrtemu(r2);
	dtmp=rsqrt;
	x=1.0f-r2*rcut21c;
	dtmp=dtmp*x*x*rsqrt*rsqrt+4.0f*rcut21c*dtmp*x;
	dtmp*=factorc;
	dtmp0=dtmp;
	dtmp*=qj;
	if(multiplyq) dtmp*=qi;
	for(k=0;k<3;k++){
	  ((float *)(fiv+k))[v]+=dtmp*dr[k];
	  if(thirdlaw){
	    if(multiplyq) ((float *)(fjv+k))[v]-=dtmp*dr[k];
	    else          ((float *)(fjv+k))[v]-=qi*dtmp0*dr[k];
	  }
	}
      }
      else if(tblnoc==1){
	rsqrt=vg_rsqrtemu(r2);
	dtmp=rsqrt;
	x=1.0f-r2*rcut21c;
	dtmp*=x*x;
	dtmp*=factorc;
	dtmp0=dtmp;
	dtmp*=qj;
	if(multiplyq) dtmp*=qi;
	for(k=0;k<1;k++){
	  ((float *)(fiv+k))[v]+=dtmp;
	  if(thirdlaw){
	    if(multiplyq) ((float *)(fjv+k))[v]+=dtmp;
	    else          ((float *)(fjv+k))[v]+=qi*dtmp0;
	  }
	}
      }
      else if(tblnoc==6){
	x=r2*rscale2;
	dtmp=qj*(m_2_sqrtpi*expf(-x)*powf(x,-1.0)
		 + erfcf(sqrtf(x))*powf(x,-1.5));
	dtmp*=factorc*rscale2*rscale;
	if(1) dtmp*=qi;
	for(k=0;k<3;k++){
	  ((float *)(fiv+k))[v]+=dtmp*dr[k];
	  if(thirdlaw) ((float *)(fjv+k))[v]-=dtmp*dr[k];
	}
      }
      else if(tblnoc==7){
	x=r2*rscale2;
	dtmp=qj*erfcf(sqrtf(x))*powf(x,-0.5);
	dtmp*=factorc*rscale;
	if(1) dtmp*=qi;
	for(k=0;k<1;k++){
	  ((float *)(fiv+k))[v]+=dtmp;
	  if(thirdlaw) ((float *)(fjv+k))[v]+=dtmp;
	}
      }
    }

    // for vdw
    if(r2!=0.0){
      rs=((float *)&rsv)[v];
      gs=((float *)&gsv)[v];
      shift=((float *)&shv)[v];
      rrs=r2*rs;
      if(gs!=0.0f && r2>=r2min && r2<r2max){
	r1=vg_r1emu(rrs);
	r6=r1*r1*r1;
	if(tblnov==2){
	  dtmp=r1;
	  r1=shift;
	  r2*=rcut21;
	  r2*=r2*r2;
	  dtmp*=r6*(2.0f*r6-1.0f)-r2*r1*(2.0f*r1-1.0f);
	  dtmp*=gs;
	  dtmp*=factorv;
	  for(k=0;k<3;k++){
	    ((float *)(fiv+k))[v]+=dtmp*dr[k];
	    if(thirdlaw) ((float *)(fjv+k))[v]-=dtmp*dr[k];
	  }
	}
	else if(tblnov==3){
	  dtmp=r6*(r6-1.0f);
	  r2*=rcut21;
	  r2*=r2*r2;
	  dtmp+=r2*shift*(2.0f*shift-1.0f)+shift*(2.0f-3.0f*shift);
	  dtmp*=gs;
	  dtmp*=factorv;
	  for(k=0;k<1;k++){
	    ((float *)(fiv+k))[v]+=dtmp;
	    if(thirdlaw) ((float *)(fjv+k))[v]+=dtmp;
	  }
	}
      }
    }
  }
}


static void MR3calccoulomb_vdw_nlist_ij_emu_vec_sub2(int iindex[],int jindex[], 
						     v4sf fiv[3], v4sf fjv[3], double force[])
{
  int v,k;

  for(v=0;v<VG_EMU_VECSIZE;v++){
    if(iindex[v]>=0 && jindex[v]>=0){
      for(k=0;k<3;k++){
	force[iindex[v]*3+k]+=((float *)(fiv+k))[v];
      }
      for(k=0;k<3;k++){
	force[jindex[v]*3+k]+=((float *)(fjv+k))[v];
      }
    }
  }
}


static void MR3calccoulomb_vdw_nlist_ij_emu_vec(int ni, double xid[], double qid[], int atypei[], double force[],
						int nj, double xjd[], double qjd[], int atypej[],
						int nat, double gscales[], double rscales[],
						double rscale, int tblno, double volume[], int periodicflag,
						int numex[], int natex[], double factord)
{
  /*
    MR3_HOST_CUTOFF is not supported 
   */
  int i,j,k,jj,n;
  char *s;
  int multiplyq=0,thirdlaw=1,iexcl=0;
  float dr[3],r2,dtmp,x,rsqrt,volumef[3]={volume[0],volume[1],volume[2]},dtmp0;
  float rscale2=(float)(rscale*rscale),m_2_sqrtpi=(float)M_2_SQRTPI;
  VG_UNIT *unit;
  double *ftmp,*fmixed;
  float factorc,factorv,rs,gs,rrs,r1,r6;
  int tblnoc=tblno,tblnov=2+(tblno % 2);
  float l[4]={(float)volume[0],(float)volume[1],(float)volume[2],0.0f};
  float al[4]={(float)(1.0/volume[0]),(float)(1.0/volume[1]),(float)(1.0/volume[2]),0.0f};
#if MD_PERIODIC_FIXED==1
  float al2[4]={0.0f,0.0f,0.0f,0.0f};
  double volume_1[3];
  DI2 di2;
  v4si xiv[3],xjv[3];
  int xia[3][VG_EMU_VECSIZE],xja[3][VG_EMU_VECSIZE];
#else
  v4sf xiv[3],xjv[3];
  float xia[3][VG_EMU_VECSIZE],xja[3][VG_EMU_VECSIZE];
#endif
  float r2max,rcut21,shift;
  VG_SCALER *scal;
  float rcut21c;
  v4sf fiv[3],fjv[3];
  v4sf qiv,qjv;
  float qia[VG_EMU_VECSIZE],qja[VG_EMU_VECSIZE];
  v4sf gsv,rsv,shv;
  float gsa[VG_EMU_VECSIZE],rsa[VG_EMU_VECSIZE],sha[VG_EMU_VECSIZE];
  int v,vv,iindex[VG_EMU_VECSIZE],jindex[VG_EMU_VECSIZE];

  fmixed=force;
  factorc=(float)(factord*Coulomb_vdw_factor);
  factorv=(float)factord;
  unit=m3_get_unit();
  //  r2max=MD_LJ_R2MAX;
  r2max=unit->r2max;
#ifdef MD_PRINT_WARN
  if(r2max!=MD_LJ_R2MAX) printf("r2max is modified in MR3calccoulomb_vdw_nlist_ij_emu_vec\n",r2max);
#endif
  rcut21c=1.0/unit->rcut2;
  scal=unit->scalers;
  rcut21=scal->rcut21;
  if(unit->r1==NULL || unit->rsqrt==NULL) vg_initialize_emu();

  if((periodicflag & 4)==0) multiplyq=1;
  if((periodicflag & 2)!=0) thirdlaw=0;
  if((periodicflag & 1)==0){
    volumef[0]*=volumef[1]*=volumef[2]*=2.0;
  }
#if MD_PERIODIC_FIXED==1
  for(i=0;i<3;i++) volume_1[i]=1.0/volume[i];
  for(i=0;i<3;i++) al2[i]=scalbnf(volumef[i],-32);
#endif
  for(i=v=0;i<ni;i++){
    for(jj=iexcl;jj<iexcl+numex[i];jj++){
      j=natex[jj];
      if(j<0) continue;
      iindex[v]=i;
      jindex[v]=j;
      for(k=0;k<3;k++){ 
	COPY_POS_INT(xia[k][v],xid[i*3+k],volume_1[k]);
      }
      qia[v]=qid[i];
      for(k=0;k<3;k++){
	COPY_POS_INT(xja[k][v],xjd[j*3+k],volume_1[k]);
      }
      qja[v]=qjd[j];
      rsa[v]=rscales[atypei[i]*nat+atypej[j]];
      gsa[v]=gscales[atypei[i]*nat+atypej[j]];
      sha[v]=(float)(pow(rscales[atypei[i]*nat+atypej[j]]*(unit->rcut2),-3.0));
      v++;
      if(v==VG_EMU_VECSIZE){
	for(k=0;k<3;k++){
#if MD_PERIODIC_FIXED==1
	  xiv[k]=*(v4si *)(xia[k]);
	  xjv[k]=*(v4si *)(xja[k]);
#else
	  xiv[k]=*(v4sf *)(xia[k]);
	  xjv[k]=*(v4sf *)(xja[k]);
#endif
	}
	qiv=*(v4sf *)qia;
	qjv=*(v4sf *)qja;
	gsv=*(v4sf *)gsa;
	rsv=*(v4sf *)rsa;
	shv=*(v4sf *)sha;
	MR3calccoulomb_vdw_nlist_ij_emu_vec_sub(fiv,fjv,xiv,xjv,
#if MD_PERIODIC_FIXED==1
						al2,
#endif
						qiv,qjv,gsv,rsv,shv,
						al,l,r2max,rcut21,rcut21c,factorc,factorv,
						nat,rscale,rscale2,tblnoc,tblnov,volumef,
						multiplyq,thirdlaw);
	MR3calccoulomb_vdw_nlist_ij_emu_vec_sub2(iindex,jindex,fiv,fjv,force);
	v=0;
      }
    }
    iexcl = iexcl + numex[i];
  }
  for(vv=v;vv<VG_EMU_VECSIZE;vv++){
    iindex[vv]=jindex[vv]=-1;
    for(k=0;k<3;k++) xia[k][vv]=xja[k][vv]=0;
    qia[vv]=qja[vv]=gsa[vv]=sha[vv]=0.0f;
    rsa[vv]=1.0f;
  }
  for(k=0;k<3;k++){
#if MD_PERIODIC_FIXED==1
    xiv[k]=*(v4si *)(xia[k]);
    xjv[k]=*(v4si *)(xja[k]);
#else
    xiv[k]=*(v4sf *)(xia[k]);
    xjv[k]=*(v4sf *)(xja[k]);
#endif
  }
  qiv=*(v4sf *)qia;
  qjv=*(v4sf *)qja;
  gsv=*(v4sf *)gsa;
  rsv=*(v4sf *)rsa;
  shv=*(v4sf *)sha;
  MR3calccoulomb_vdw_nlist_ij_emu_vec_sub(fiv,fjv,xiv,xjv,
#if MD_PERIODIC_FIXED==1
					  al2,
#endif
					  qiv,qjv,gsv,rsv,shv,
					  al,l,rcut21,rcut21c,factorc,factorv,
					  nat,rscale,rscale2,tblnoc,tblnov,volumef,
					  multiplyq,thirdlaw);
  MR3calccoulomb_vdw_nlist_ij_emu_vec_sub2(iindex,jindex,fiv,fjv,force);
}
#endif


void MR3calccoulomb_vdw_nlist_ij_emu(int ni, double xid[], double qid[], int atypei[], double force[],
				     int nj, double xjd[], double qjd[], int atypej[],
				     int nat, double gscales[], double rscales[],
				     double rscale, int tblno, double volume[], int periodicflag,
				     int numex[], int natex[], double factord)
/*
  Only MR3calccoulomb_vdw_nlist_ij_emu_work supports separate potential result
  in force. force[i*3] contains Coulomb. force[i*3+1] contains only vdW (from 101104).
  in force. force[i*3] contains Coulomb + vdW. force[i*3+1] contains only vdW (until 1000826).
 */
{
#if MD_LJ_05SIGMA_8SIGMA==2 && defined(COULOMB_SHIFT) && VDW_SHIFT>=1 && defined(VG_GCC)
  // not work well
  MR3calccoulomb_vdw_nlist_ij_emu_vec(ni,xid,qid,atypei,force,
				      nj,xjd,qjd,atypej,
				      nat,gscales,rscales,
				      rscale,tblno,volume,periodicflag,
				      numex,natex,factord);
#else
  MR3calccoulomb_vdw_nlist_ij_emu_work(ni,xid,qid,atypei,force,
				       nj,xjd,qjd,atypej,
				       nat,gscales,rscales,
				       rscale,tblno,volume,periodicflag,
				       numex,natex,factord);
#endif
}


void MR3calcvdw_nlist_ij_emu(int ni, double xid[], int atypei[], double force[],
			     int nj, double xjd[], int atypej[], 
			     int nat, double gscale[], double rscale[], int tblno,
			     double xmaxd, int periodicflag,
			     int numex[], int natex[], double factord) __attribute__((optimize(0)));

#pragma optimize("",off)
void MR3calcvdw_nlist_ij_emu(int ni, double xid[], int atypei[], double force[],
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
  int i,j,k,iexcl=0,jj,thirdlaw=1,n;
  float r2min=MD_LJ_R2MIN,r2max;
  float factor=factord;
  float dtmp,r1,r6,r2,rs,gs,rrs,dr[3],xmax;
#if VDW_SHIFT>=1
  float rcut21,shift;
  VG_SCALER *scal;
#endif
  VG_UNIT *unit;
#if MD_PERIODIC_FIXED==1
  int xi[3],xj[3];
  float al2[3];
  //  double xmax_1;
  double volume_1[3];
  DI2 di2;
#else
  float xi[3],xj[3];
#endif

#if defined(MD_PRINT_WARN) 
  {
    static int ini=0;
    if(ini==0){
      printf("** warning : MR3calcvdw_nlist_ij_emu might be wrong with compiler optimization **\n");
      ini=1;
    }
  }
#endif
#if 0
  {
    double sum=0.0;
    for(i=0;i<ni;i++) sum+=force[i*3];
    printf("before MR3calc_vdw_nlist_ij_emu is called: sum=%24.17e\n",sum);
  }
#endif
  unit=m3_get_unit();
  //  r2max=MD_LJ_R2MAX;
  r2max=unit->rcut2;
#ifdef MD_PRINT_WARN
  {
    static int ini=0;
    if(ini==0){
      if(r2max!=MD_LJ_R2MAX) printf("r2max is modified in MR3calcvdw_nlist_ij_emu\n",r2max);
      ini=1;
    }
  }
#endif
#if VDW_SHIFT>=1
  scal=unit->scalers;
  rcut21=scal->rcut21;
#endif
  if(unit->r1==NULL || unit->rsqrt==NULL) vg_initialize_emu();
  if((periodicflag & 1)==0){
    xmaxd*=2.0;
  }
  xmax=xmaxd;
#if MD_PERIODIC_FIXED==1
  //xmax_1=1.0/xmaxd;
  for(i=0;i<3;i++) volume_1[i]=1.0/xmaxd;
  for(i=0;i<3;i++) al2[i]=scalbnf(xmax,-32);
#endif
  if((periodicflag & 2)!=0) thirdlaw=0;
  for(i=0;i<ni;i++){
    for(k=0;k<3;k++){ 
      COPY_POS_INT(xi[k],xid[i*3+k],volume_1[k]);
    }
    for(jj=iexcl;jj<iexcl+numex[i];jj++){
      j=natex[jj];
      if(j<0) continue;
      r2=0.0;
      for(k=0;k<3;k++){
	COPY_POS_INT(xj[k],xjd[j*3+k],volume_1[k]);
      }
      for(k=0;k<3;k++){
	dr[k]=(float)(xi[k]-xj[k]);
#if MD_PERIODIC_FIXED==1
	//	printf("  xi[%d]=%e(%08x) xj=%e(%08x) dr=%e\n",k,xid[i*3+k],xi[k],xjd[i*3+k],xj[k],dr[k]);
	dr[k]*=al2[k];
#elif MD_PERIODIC_FIXED==2
	if(dr[k]<-xmax/2.0){
	  dr[k]+=xmax;
	}
	if(dr[k]>=xmax/2.0){
	  dr[k]-=xmax;
	}
#endif
	//	printf("i=%d j=%d dr[%d]=%e\n",i,j,k,dr[k]);
#if 0 // round off error occurrs because of MADD accuracy is different from GPU
	r2+=dr[k]*dr[k];
#elif 1
	r2=dr[0]*dr[0];
	r2=madd_emu(dr[1],dr[1],r2);
	r2=madd_emu(dr[2],dr[2],r2);
#else // work 
	if(k==0){
	  r2+=dr[k]*dr[k];
	}
	else{
	  unsigned long long ull;
	  double ddr;
	  ddr=dr[k];
	  ddr*=ddr;
	  ull=*((unsigned long long *)&ddr);
	  ull&=0xffffffffe0000000LL;
	  ddr=*((double *)&ull);
	  r2+=(float)ddr;
	}
#endif
      }
      if(r2!=0.0){
	rs=rscale[atypei[i]*nat+atypej[j]];
	gs=gscale[atypei[i]*nat+atypej[j]];
#if VDW_SHIFT==1
	shift=(float)(pow(rscale[atypei[i]*nat+atypej[j]]*(unit->rcut2),-3.0));
#endif
	rrs=r2*rs;
#if MD_LJ_05SIGMA_8SIGMA==2
	if(gs!=0.0f && r2>=r2min && r2<r2max){
#else
	if(gs!=0.0f && rrs>=r2min && rrs<r2max){
#endif
	  r1=vg_r1emu(rrs);
	  //r1=1.0f/rrs;
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
	    dtmp=gs*r6*r1*(2.0f*r6-1.0f);
#endif
	    dtmp*=factor;
#if 1 // work
	    for(k=0;k<3;k++){
	      force[i*3+k]+=dtmp*dr[k];
	      if(thirdlaw) force[j*3+k]-=dtmp*dr[k];
	    }
#else // debug
	    //	    force[i*3]+=rrs;
	    {
	      float ftmp;
	      //	      ftmp=dr[0]*dr[0];
	      force[i*3]=madd_emu(dr[0],dr[0],force[i*3]);
	      force[i*3]=madd_emu(dr[1],dr[1],force[i*3]);
	      //	      ftmp=dr[0]*dr[0]+dr[1]*dr[1];
	      //	      force[i*3]+=ftmp;
	    }
#endif
	  }
	  else if(tblno==3){
#if VDW_SHIFT==1
#if 1
	    dtmp=r6*(r6-1.0f);
#if 0
	    //	    r1=powf(rs*rcut2,-3.0f);
	    r1=shift;
	    r2*=r2*rs;
	    //	    r2*=r2;
	    dtmp+=r1*r1*(r1-1.0f)*r2+2.0f*r1*(1.0f-r1);
#endif
	    r2*=rcut21;
	    r2*=r2*r2;
	    dtmp+=r2*shift*(2.0f*shift-1.0f)+shift*(2.0f-3.0f*shift);
	    dtmp*=gs;
#else
	    {
	      VG_SCALER *scal=unit->scalers;
	      float rcut21=scal->rcut21;
	      dtmp=vdw_shift_emu(r2,rs,gs,rcut21,1);
	    }
#endif
#else
	    dtmp=gs*r6*(r6-1.0f);
	    //	    dtmp=gs*(pow(rrs,-6.0)-pow(rrs,-3.0));
#endif
	    dtmp*=factor;
	    //	    for(k=0;k<3;k++){
	    for(k=0;k<1;k++){
	      force[i*3+k]+=dtmp;
	      if(thirdlaw) force[j*3+k]+=dtmp;
	    }
	  }
	}
      }
    }
    iexcl = iexcl + numex[i];
  }
#if 0
  {
    double sum=0.0;
    for(i=0;i<ni;i++) sum+=force[i*3];
    printf("after MR3calc_vdw_nlist_ij_emu is called: sum=%24.17e\n",sum);
  }
#endif
}
#pragma optimize("",on)


void MR3calcvdw_nlist_emu2(double xd[], int n, int atype[], int nat,
			   double gscale[], double rscale[], int tblno,
			   double xmaxd, int periodicflag,
			   int numex[], int natex[], double factord,
			   double force[])
{
  MR3calcvdw_nlist_ij_emu(n,xd,atype,force,n,xd,atype,
			  nat,gscale,rscale,tblno,xmaxd,periodicflag,
			  numex,natex,factord);
  /*
  if(periodicflag==0){
    xmax*=2.0;
  }
  if((periodicflag & 2)!=0) thirdlaw=0;
  for(i=0;i<n;i++){
    for(k=0;k<3;k++) xi[k]=xd[i*3+k];
    for(jj=iexcl;jj<iexcl+numex[i];jj++){
      j=natex[jj];
      if(j<0) continue;
      r2=0.0;
      for(k=0;k<3;k++) xj[k]=xd[j*3+k];
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
	rs=rscale[atype[i]*nat+atype[j]];
	gs=gscale[atype[i]*nat+atype[j]];
	rrs=r2*rs;
	if(gs>0.0f && rrs>=r2min && rrs<r2max){
	  r1=vg_r1emu(rrs);
	  //r1=1.0f/rrs;
	  r6=r1*r1*r1;
	  if(tblno==2){
	    dtmp=gs*r6*r1*(2.0f*r6-1.0f);
	    dtmp*=factor;
	    for(k=0;k<3;k++){
	      force[i*3+k]+=dtmp*dr[k];
	      if(thirdlaw) force[j*3+k]-=dtmp*dr[k];
	    }
	  }
	  else if(tblno==3){
	    dtmp=gs*r6*(r6-1.0f);
	    //	    dtmp=gs*(pow(rrs,-6.0)-pow(rrs,-3.0));
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
  */
}


void MR3calccoulomb_ij_exlist(int ni, double xi[], double qi[], double force[],
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
  int i,multiplyq=0;
  VG_UNIT *unit=m3_get_unit();

#if 1
  {
    char *s;
    static char method=0;
    static int ini=0;
    int j;
    if(ini==0){
      if((s=getenv("MR3_MIXEDACCURACY_METHOD"))!=NULL){
	sscanf(s,"%c",&method);
	printf("MR3_MIXEDACCURACY_METHOD is defined and method=%c is used for MR3calccoulomb_ij_exlist\n",method);
      }
      ini=1;
    }
    if(method!=0){
#if 0 // not work currently
      MR3calc_ij_exlist_host_mixedaccuracy(ni,xi,qi,NULL,force,
					   nj,xj,qj,NULL,
					   0,NULL,NULL,
					   tblno,xmax,periodicflag,
					   numex,natex,method);
#else
      MR3calc_ij_exlist_host_mixedaccuracy(ni,xi,qi,NULL,force,
					   nj,xj,qj,NULL,
					   0,NULL,NULL,
					   tblno,xmax,periodicflag & 3,
					   numex,natex,method);
      for(i=0;i<ni;i++) for(j=0;j<3;j++) force[i*3+j]*=qi[i];
#endif
      //      for(i=0;i<3;i++) printf("after MR3calccoulomb_ij_exlist : force[%d]=%e %e %e\n",i,force[i*3],force[i*3+1],force[i*3+2]);
      return;
    }
  }
#endif

  //  printf("in MR3calccoulomb_ij_exlist : tblno=%d periodicflag=%d\n",tblno,periodicflag);
#if 0
  if(ni!=nj){
    fprintf(stderr,"** warning : ni should be equal to nj **\n");
    //    vg_exit(1);
  }
#endif
  //  printf("in MR3calccoulomb_ij_exlist : ni=%d nj=%d rscale=%e tblno=%d xmax=%e periodicflag=%d\n",ni,nj,rscale,tblno,xmax,periodicflag);
  if((periodicflag & 4)!=0) multiplyq=1;
  if(tblno==6 || tblno==7) multiplyq=0;
  /*
  if((periodicflag & 4)==0){
    fprintf(stderr,"** error : bit2 of periodicflag must be activated **\n");
    vg_exit(1);
    }*/

#if 1
  //  printf("tblno is modified\n");if(tblno>=6) tblno-=6;
  MR3calccoulomb_nlist_ij_emu(ni,xi,qi,force,nj,xj,qj,rscale,
			      tblno,xmax,(periodicflag & 3)+(1-multiplyq)*4,
			      numex,natex,-1.0);
  //  for(i=0;i<3;i++) printf("after coulomb_nlist_ij_emu: force[%d]=%e %e %e\n",i,force[i*3],force[i*3+1],force[i*3+2]);
#else
  MR3calccoulomb_nlist_ij_host(ni,xi,qi,force,nj,xj,qj,rscale,
			       tblno,xmax,(periodicflag & 3)+(1-multiplyq)*4,
			       numex,natex,-1.0);
  printf("** MR3calccoulomb_nlist_ij_host is called in MR3calccoulomb_ij_exlist **\n");
#endif
  //  for(i=0;i<3;i++) printf("after MR3calccoulomb_nlist_ij_emu force[%d]=%e %e %e\n",i,force[i*3],force[i*3+1],force[i*3+2]);
#if 1 // use GPU
#ifdef MD_PRINT_WARN
  printf("before calling MR3calccoulomb_ij_ci: Ldim=(%d,%d,%d)\n",Ldim[0],Ldim[1],Ldim[2]);
#endif
#if defined(MD_CELLINDEX) // use cell-index routine
  if(Ldim[0]>1 || Ldim[1]>1 || Ldim[2]>1){
    double volume[3]={xmax,xmax,xmax},rcut=sqrt(unit->rcut2),skinnb=0.0;
    MR3calccoulomb_ij_ci(ni,xi,qi,force,nj,xj,qj,rscale,tblno,rcut,skinnb,volume,Ldim);
    //    MR3calccoulomb_ij_ci_old(ni,xi,qi,force,nj,xj,qj,rscale,tblno,rcut,skinnb,volume,Ldim);printf("** MR3calccoulomb_ij_ci_old is called\n");
  }
  else{
    MR3calccoulomb_ij(ni,xi,qi,force,nj,xj,qj,rscale,tblno,xmax,(periodicflag & 1)+multiplyq*2);
  }
#else
  MR3calccoulomb_ij(ni,xi,qi,force,nj,xj,qj,rscale,tblno,xmax,(periodicflag & 1)+multiplyq*2);
#endif
#elif 1 // use host
  MR3calccoulomb_ij_host(ni,xi,qi,force,nj,xj,qj,rscale,tblno,xmax,(periodicflag & 1)+multiplyq*2);
  printf("** MR3calccoulomb_ij_host is called in MR3calccoulomb_ij_exlist **\n");
#else // use emu
  MR3calccoulomb_ij_emu(ni,xi,qi,force,nj,xj,qj,rscale,tblno,xmax,(periodicflag & 1)+multiplyq*2);
  printf("** MR3calccoulomb_ij_emu is called in MR3calccoulomb_ij_exlist **\n");
#endif
  //  MR3calccoulomb_ij_emu(ni,xi,qi,force,nj,xj,qj,rscale,tblno,xmax,(periodicflag & 1)+multiplyq*2);printf("** MR3calccoulomb_ij_emu is called **\n");
  //  MR3calccoulomb_ij_host(ni,xi,qi,force,nj,xj,qj,rscale,tblno,xmax,(periodicflag & 1)+multiplyq*2);printf("** MR3calccoulomb_ij_host is called **\n");
  //  for(i=0;i<3;i++) printf("after MR3calccoulomb_ij force[%d]=%e %e %e\n",i,force[i*3],force[i*3+1],force[i*3+2]);

#if 0
  for(i=0;i<2;i++) printf("after MR3calccoulomb_ij force[%d]=%e %e %e\n",i,force[i*3],force[i*3+1],force[i*3+2]);
  MR3calccoulomb_nlist_emu(xi,ni,qi,rscale,tblno,xmax,periodicflag & 3,
			   numex,natex,-1.0,force);
#endif
  //  for(i=0;i<3;i++) printf("after MR3calccoulomb_ij_exlist : force[%d]=%e %e %e\n",i,force[i*3],force[i*3+1],force[i*3+2]);
}



void mr3calccoulomb_ij_exlist_(int *ni, double xi[], double qi[], 
			       double force[],
			       int *nj, double xj[], double qj[],
			       double *rscale, 
			       int *tblno, double *xmax, int *periodicflag,
			       int numexorg[], int natexorg[])
{
  /*
    periodicflag bit 0 : 0 --- non periodic
                         1 --- periodic
                 bit 4 : 0 --- no overlap
                         1 --- overlap
  */
  int i,*natex2,s;
  int *numex,*natex;
  static int *numextmp=NULL,*natextmp=NULL;

  if((*periodicflag & (1<<4))!=0){
    Overlapflag=1;
  }

#ifdef MD_PRINT_WARN
  printf("in mr3calccoulomb_ij_exlist_: Overlapflag=%d\n",Overlapflag);fflush(stdout);
#endif
  
  for(i=s=0;i<*ni;i++){
    s+=numexorg[i];
  }
  if(s==0){
    mr3calccoulomb_ij_(ni,xi,qi,force,nj,xj,qj,
		       rscale,tblno,xmax,periodicflag);
    return;
  }

#if 0
  {
    int j,jpoint,jpointorg;
    printf("  numextmp and natextmp are used\n");
    if(numextmp==NULL){
      if((numextmp=(int *)MR3_malloc_pointer(sizeof(int)*(*ni),"numextmp in mr3calccoulomb_ij_exlist_"))==NULL){
	fprintf(stderr,"** error : can't malloc numextmp **\n");
	MR3_exit(1);
      }
      if((natextmp=(int *)MR3_malloc_pointer(sizeof(int)*s,"natextmp in mr3calccoulomb_ij_exlist_"))==NULL){
	fprintf(stderr,"** error : can't malloc natextmp **\n");
	MR3_exit(1);
      }
    }
    for(i=0;i<*ni;i++) numextmp[i]=numexorg[i];
    jpointorg=jpoint=0;
#if 1
    numextmp[0]=2;
    for(j=0;j<numextmp[0];j++){
      natextmp[jpoint++]=natexorg[j];
    }
    jpointorg+=numexorg[0];
#else
    for(j=0;j<numexorg[0];j++){
      natextmp[jpoint++]=natexorg[jpointorg++];
    }
#endif
    for(i=1;i<*ni;i++){
      for(j=0;j<numexorg[i];j++){
	natextmp[jpoint++]=natexorg[jpointorg++];
      }
    }
    for(i=0;i<*ni;i++) numexorg[i]=numextmp[i];
    jpointorg=jpoint=0;
    for(i=0;i<*ni;i++){
      for(j=0;j<numextmp[i];j++){
	natexorg[jpoint++]=natextmp[jpointorg++];
      }
    }
    
    numex=numextmp;
    natex=natextmp;
  }
#else
  numex=numexorg;
  natex=natexorg;
#endif

  if((natex2=(int *)MR3_malloc_pointer(sizeof(int)*s,"natex2 in mr3calccoulomb_ij_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc natex2 in mr3calccoulomb_ij_exlist_ **\n");
    MR3_exit(1);
  }
  for(i=0;i<s;i++)  natex2[i]=natex[i]-1;
  MR3calccoulomb_ij_exlist(*ni,xi,qi,force,
			   *nj,xj,qj,
			   *rscale,*tblno,*xmax,*periodicflag,
			   numex,natex2);
  MR3_free_pointer(natex2,"natex2 in mr3calccoulomb_ij_exlist_");
  Overlapflag=0;
}


void mr3calccoulomb_ij_exlist__(int *ni, double xi[], double qi[], 
			       double force[],
			       int *nj, double xj[], double qj[],
			       double *rscale, 
			       int *tblno, double *xmax, int *periodicflag,
			       int numexorg[], int natexorg[])
{
  mr3calccoulomb_ij_exlist_(ni,xi,qi, 
			    force,
			    nj,xj,qj,
			    rscale, 
			    tblno,xmax,periodicflag,
			    numexorg,natexorg);
}



void MR3calcvdw_ij_exlist(int ni, double xi[], int atypei[], double force[],
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
  VG_UNIT *unit=m3_get_unit();

#if 1
  {
    char *s;
    static char method=0;
    static int ini=0;
    if(ini==0){
      if((s=getenv("MR3_MIXEDACCURACY_METHOD"))!=NULL){
	sscanf(s,"%c",&method);
	printf("MR3_MIXEDACCURACY_METHOD is defined and method=%c is used for MR3calcvdw_ij_exlist\n",method);
      }
      ini=1;
    }
    if(method!=0){
      MR3calc_ij_exlist_host_mixedaccuracy(ni,xi,NULL,atypei,force,
					   nj,xj,NULL,atypej,
					   nat,gscale,rscale,
					   tblno,xmax,periodicflag,
					   numex,natex,method);
      return;
    }
  }
#endif

  //  printf("in MR3calcvdw_ij_exlist : ni=%d nj=%d nat=%d tblno=%d xmax=%e periodicflag=%d\n",ni,nj,nat,tblno,xmax,periodicflag);
#if 1
  MR3calcvdw_nlist_ij_emu(ni,xi,atypei,force,nj,xj,atypej,nat,gscale,rscale,tblno,
			  xmax,(periodicflag & 3),numex,natex,-1.0);
#else
  MR3calcvdw_nlist_ij_host(ni,xi,atypei,force,nj,xj,atypej,nat,gscale,rscale,tblno,
			   xmax,(periodicflag & 3),numex,natex,-1.0);
  printf("MR3calcvdw_nlist_ij_host is called in MR3calcvdw_ij_exlist\n");
#endif
  //  for(i=0;i<ni;i++) for(j=0;j<3;j++) if(force[i*3+j]!=0) printf("emu force[%d*3+%d]=%e\n",i,j,force[i*3+j]);
  //  for(i=0;i<3;i++) printf(" after MR3calcvdw_nlist_ij_emu force[%d]=%e(%016llx) %e(%016llx) %e(%016llx)\n",i,force[i*3],((unsigned long long *)(force+i*3))[0],force[i*3+1],((unsigned long long *)(force+i*3))[1],force[i*3+2],((unsigned long long *)(force+i*3))[2]);
#if 1 // use GPU
#if defined(MD_CELLINDEX) // use cell-index routine
  if(Ldim[0]>1 || Ldim[1]>1 || Ldim[2]>1){
    double volume[3]={xmax,xmax,xmax},rcut=sqrt(unit->rcut2),skinnb=0.0;
    MR3calcvdw_ij_ci(ni,xi,atypei,force,nj,xj,atypej,nat,gscale,rscale,tblno,rcut,skinnb,volume,Ldim);
  }
  else{
    MR3calcvdw_ij(ni,xi,atypei,force,nj,xj,atypej,nat,gscale,rscale,tblno,xmax,(periodicflag & 1));
  }
#else
  MR3calcvdw_ij(ni,xi,atypei,force,nj,xj,atypej,nat,gscale,rscale,tblno,xmax,(periodicflag & 1));
#endif
#elif 0
  MR3calcvdw_ij_emu(ni,xi,atypei,force,nj,xj,atypej,nat,gscale,rscale,tblno,xmax,(periodicflag & 1));printf("** MR3calcvdw_ij_emu is called **\n");
  printf("*** MR3calcvdw_ij_emu is called ***\n");
#else
  printf("*** MR3calcvdw_ij_host is called ***\n");
  MR3calcvdw_ij_host(ni,xi,atypei,force,nj,xj,atypej,nat,gscale,rscale,tblno,
		     xmax,(periodicflag & 1));
#endif
  //  for(i=0;i<ni;i++) for(j=0;j<3;j++) if(force[i*3+j]!=0) printf("force[%d*3+%d]=%e\n",i,j,force[i*3+j]);
  //  for(i=0;i<3;i++) printf(" after MR3calcvdw_ij_exlist force[%d]=%e %e %e\n",i,force[i*3],force[i*3+1],force[i*3+2]);
}


void mr3calcvdw_ij_exlist_(int *ni, double xi[], int atypei[], double force[],
			   int *nj, double xj[], int atypej[],
			   int *nat, double gscale[], double rscale[],
			   int *tblno, double *xmax, int *periodicflag,
			   int numex[], int natex[])
{
  /*
    periodicflag bit 0 : 0 --- non periodic
                         1 --- periodic
                 bit 4 : 0 --- no overlap
                         1 --- overlap
  */
  int *atypei2,*atypej2;
  int i,*natex2,s;

  for(i=s=0;i<*ni;i++){
    s+=numex[i];
  }
  if(s==0){
    mr3calcvdw_ij_(ni,xi,atypei,force,nj,xj,atypej,nat,gscale,rscale,
		  tblno,xmax,periodicflag);
    return;
  }
  if((*periodicflag & (1<<4))!=0){
    /*    printf("Overlapflag is set in md3calcvdw_ij_exlist_\n");*/
    Overlapflag=1;
  }

  if((atypei2=(int *)MR3_malloc_pointer(sizeof(int)*(*ni),"atypei2 in mr3calcvdw_ij_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc atypei2 in mr3calcvdw_ij_ **\n");
    MR3_exit(1);
  }
  if((atypej2=(int *)MR3_malloc_pointer(sizeof(int)*(*nj),"atypej2 in mr3calcvdw_ij_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc atypej2 in mr3calcvdw_ij_ **\n");
    MR3_exit(1);
  }
  for(i=0;i<*ni;i++) atypei2[i]=atypei[i]-1;
  for(i=0;i<*nj;i++) atypej2[i]=atypej[i]-1;
  if((natex2=(int *)MR3_malloc_pointer(sizeof(int)*s,"natex2 in mr3calcvdw_ij_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc natex2 in mr3calcvdw_ij_exlist_ **\n");
    MR3_exit(1);
  }
  for(i=0;i<s;i++)  natex2[i]=natex[i]-1;
  
  MR3calcvdw_ij_exlist(*ni,xi,atypei2,force,
		       *nj,xj,atypej2,
		       *nat,gscale,rscale,
		       *tblno,*xmax,*periodicflag,
		       numex,natex2);
  MR3_free_pointer(atypei2,"atypei2 in mr3calcvdw_ij_exlist_");
  MR3_free_pointer(atypej2,"atypej2 in mr3calcvdw_ij_exlist_");
  MR3_free_pointer(natex2,"natex2 in mr3calcvdw_ij_exlist_");
  Overlapflag=0;
}


void mr3calcvdw_ij_exlist__(int *ni, double xi[], int atypei[], double force[],
			    int *nj, double xj[], int atypej[],
			    int *nat, double gscale[], double rscale[],
			    int *tblno, double *xmax, int *periodicflag,
			    int numex[], int natex[])
{
  mr3calcvdw_ij_exlist_(ni,xi,atypei,force,
			nj,xj,atypej,
			nat,gscale,rscale,
			tblno,xmax,periodicflag,
			numex,natex);
}


static void make_natex2(int n, int numex[], int natex[],
			int **numex2, int **natex2, int flag)
{
  /* flag 0 -- do not malloc
          1 -- malloc
  */
  int i,j,k,s,s2,iexcl,*natex2_offset;

  for(i=s=0;i<n;i++){
    s+=numex[i];
  }
  if(flag==1){
    if(((*numex2)=(int *)MR3_malloc_pointer(sizeof(int)*n,"numex2 in make_natex2"))==NULL){
      fprintf(stderr,"** error : can't malloc numex2 in make_natex2 **\n");
      MR3_exit(1);
    }
  }
  if((natex2_offset=(int *)MR3_malloc_pointer(sizeof(int)*n,"natex2_offset in make_natex2"))==NULL){
    fprintf(stderr,"** error : can't malloc natex2_offset in make_natex2 **\n");
    MR3_exit(1);
  }
  for(i=0;i<n;i++) (*numex2)[i]=0;
  iexcl=s2=0;
  for(i=0;i<n;i++){
    for(k=iexcl;k<iexcl+numex[i];k++){
      j=natex[k];
      if (j>=0){
	(*numex2)[i]++;
	(*numex2)[j]++;
	s2++;
      }
      else{
	//	printf("j=%d is minus. i=%d\n",j,i);
      }
    }
    iexcl += numex[i];
  }
  //  printf("total sum of natex=%d, for j>=0 sum is %d\n",s,s2);

  for(i=s=0;i<n;i++) s+=(*numex2)[i];
  //  printf("s2=%d s2/2=%d\n",s,s/2);

  if(s!=s2*2){
    fprintf(stderr,"** error : s=%d and s2*2=%d should be same **\n",s,s2*2);
    MR3_exit(1);
  }

  if(flag==1){
    if(((*natex2)=(int *)MR3_malloc_pointer(sizeof(int)*s,"natex2 in make_natex2"))==NULL){
      fprintf(stderr,"** error : can't malloc natex2 in make_natex2 **\n");
      MR3_exit(1);
    }
  }

  natex2_offset[0]=0;
  for(i=1;i<n;i++){
    natex2_offset[i]=natex2_offset[i-1]+(*numex2)[i-1];
  }
  for(i=0;i<n;i++) (*numex2)[i]=0;

  iexcl=0;
  for(i=0;i<n;i++){
    for(k=iexcl;k<iexcl+numex[i];k++){
      j=natex[k];
      if (j>=0){
	(*natex2)[natex2_offset[i]+(*numex2)[i]]=j;
	(*numex2)[i]++;
	(*natex2)[natex2_offset[j]+(*numex2)[j]]=i;
	(*numex2)[j]++;
      }
    }
    iexcl += numex[i];
  }

#if 0
  for(i=s=0;i<n;i++) s+=(*numex2)[i];
  printf("s=%d\n",s);
#endif

  MR3_free_pointer(natex2_offset,"natex2_offset in make_natex2");
}


void MR3calccoulomb_vdw_ij_exlist_core_old(int ni, double xi[], double qi[], 
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
  //#define DEBUG_CVIEC  // to debug this routine
  double *ftmp;
  int i,j,multiplyq=1;

  printf("** warning : differenet rcut is not used in this routine **\n");
#if 0
  printf("in MR3calccoulomb_vdw_ij_exlist, host2 routine is called.\n");
  /*
  printf("periodicflag=%d potflag=%d\n",periodicflag,potflag);
  MR3calccoulomb_vdw_ij_exlist_host(ni,xi,qi,atypei,force, 
				    nj,xj,qj,atypej,
				    nat,gscalesf,gscalesp,rscales,rscale,tblno,
				    xmax,potc,potv,0,1,0,
				    numex,natex);
  */
  MR3calccoulomb_vdw_ij_exlist_host2(ni,xi,qi,atypei,force,nj,xj,qj,atypej,
				     nat,gscalesf,gscalesp,rscales,rscale,
				     tblno,xmax,potc,potv,
				     periodicflag,potflag,changeflag,
				     numex,natex);
  return;
#endif  


#if 1 // work at 0812
#ifdef DEBUG_CVIEC
  {
    static int ini=0;
    if(ini==0){
      printf("** DEBUG_CVIEC is defined **\n");
      ini=1;
    }
  }
  if(1){
#else
  if(potflag!=0){
#endif
    if((ftmp=(double *)MR3_malloc_pointer(sizeof(double)*ni*3,"MR3calccoulomb_vdw_ij_exlist_core"))==NULL){
      fprintf(stderr,"** error : can't malloc ftmp **\n");
      vg_exit(1);
    }
  }
  // Coulomb
  if(tblno==6) multiplyq=0;
#if defined(DEBUG_CVIEC) && 0
  if(1){
#else
  if(potflag!=0){
#endif
    bzero(ftmp,sizeof(double)*ni*3);
    MR3calccoulomb_ij_exlist(ni,xi,qi,ftmp,nj,xj,qj,
			     rscale,tblno+1,xmax,(periodicflag & 3)+4*multiplyq,
			     numex,natex);
    for(i=0;i<ni;i++) (*potc)+=ftmp[i*3];
  }
#if defined(DEBUG_CVIEC) && 1
  {
#if 0
    double *ftmp2=(double *)MR3_malloc_pointer(sizeof(double)*ni*3);
    bzero(ftmp2,sizeof(double)*ni*3);// force
    MR3calccoulomb_ij_exlist(ni,xi,qi,ftmp2,nj,xj,qj,
			     rscale,tblno,xmax,(periodicflag & 3)+4*multiplyq,
			     numex,natex);
#endif
#if 1
    bzero(ftmp,sizeof(double)*ni*3);// energy
    MR3calccoulomb_ij_exlist(ni,xi,qi,ftmp,nj,xj,qj,
			     rscale,tblno+1,xmax,(periodicflag & 3)+4*multiplyq,
			     numex,natex);
    for(i=0;i<ni;i++) (*potc)+=ftmp[i*3];
#endif
#if 0
    bzero(ftmp,sizeof(double)*ni*3);// force
    MR3calccoulomb_ij_exlist(ni,xi,qi,ftmp,nj,xj,qj,
			     rscale,tblno,xmax,(periodicflag & 3)+4*multiplyq,
			     numex,natex);
    for(i=0;i<ni*3;i++) force[i]+=ftmp[i];
#endif
#if 1 // force
    MR3calccoulomb_ij_exlist(ni,xi,qi,force,nj,xj,qj,
			     rscale,tblno,xmax,(periodicflag & 3)+4*multiplyq,
			     numex,natex);
#endif
#if 0
    for(i=0;i<ni*3;i++){
      if(ftmp[i]!=ftmp2[i]){
	printf("ftmp[%d][%d]=%20.14f ftmp2=%20.14f\n",i/3,i% 3,ftmp[i],ftmp2[i]);
      }
    }
    MR3_free_pointer(ftmp2,"MR3calccoulomb_vdw_ij_exlist_core");
#endif
  }
#else // else of (defined(DEBUG_CVIEC) && 1)
  MR3calccoulomb_ij_exlist(ni,xi,qi,force,nj,xj,qj,
			   rscale,tblno,xmax,(periodicflag & 3)+4*multiplyq,
			   numex,natex);
#endif // end of (defined(DEBUG_CVIEC) && 1)
  // vdw
#if defined(DEBUG_CVIEC) && 0
  if(1){
#else
  if(potflag!=0){
#endif
    bzero(ftmp,sizeof(double)*ni*3);
    MR3calcvdw_ij_exlist(ni,xi,atypei,ftmp,nj,xj,atypej,
			 nat,gscalesp,rscales,3,xmax,(periodicflag & 3),
			 numex,natex);
    for(i=0;i<ni;i++) (*potv)+=ftmp[i*3];
  }
  MR3calcvdw_ij_exlist(ni,xi,atypei,force,nj,xj,atypej,
		       nat,gscalesf,rscales,2,xmax,(periodicflag & 3),
		       numex,natex);

#ifdef DEBUG_CVIEC
  if(1){
#else
  if(potflag!=0){
#endif
    MR3_free_pointer(ftmp,"MR3calccoulomb_vdw_ij_exlist_core");
  }
  return;
#endif // end of if 1 work at 0812

  // Coulomb
  if(tblno==6) multiplyq=0;
  if(potflag!=0){
    if((ftmp=(double *)MR3_malloc_pointer(sizeof(double)*ni*3,"MR3calccoulomb_vdw_ij_exlist_core"))==NULL){
      fprintf(stderr,"** error : can't malloc ftmp **\n");
      vg_exit(1);
    }
    bzero(ftmp,sizeof(double)*ni*3);
    MR3calccoulomb_ij_exlist(ni,xi,qi,ftmp,nj,xj,qj,
			     rscale,tblno+1,xmax,(periodicflag & 3)+4*multiplyq,
			     numex,natex);
    for(i=0;i<ni;i++) (*potc)+=ftmp[i*3];
    //    printf("  potc=%e in MR3coulomb_vdw_ij_exlist_core\n",*potc);
    bzero(ftmp,sizeof(double)*ni*3);
#if 1
    MR3calcvdw_ij_exlist(ni,xi,atypei,ftmp,nj,xj,atypej,
			 nat,gscalesp,rscales,3,xmax,(periodicflag & 3),
			 numex,natex);
#else
    printf("MR3calcvdw_ij_exlist_host is called for potential\n");
    MR3calcvdw_ij_exlist_host(ni,xi,atypei,ftmp,nj,xj,atypej,
			 nat,gscalesp,rscales,3,xmax,(periodicflag & 3),
			 numex,natex);
#endif
    for(i=0;i<ni;i++) (*potv)+=ftmp[i*3];
    //    printf("  potv=%e in MR3coulomb_vdw_ij_exlist_core\n",*potv);
    MR3_free_pointer(ftmp,"MR3calccoulomb_vdw_ij_exlist_core");
  }

#if 0 // calc separately (note that tblno=6 and 7 may not work. 
      // please see multiplyq=0 in MR3calccoulomb_ij_exlist)
  printf("** calc separately in MR3calccoulomb_vdw_ij_exlist_core\n");
  MR3calccoulomb_nlist_ij_emu(ni,xi,qi,force,nj,xj,qj,rscale,tblno,
			      xmax,(periodicflag & 3),numex,natex,-1.0);
  MR3calccoulomb_ij(ni,xi,qi,force,nj,xj,qj,rscale,tblno,
		    xmax,(periodicflag & 1)+2*multiplyq);
  MR3calcvdw_nlist_ij_emu(ni,xi,atypei,force,nj,xj,atypej,
			  nat,gscalesf,rscales,2,
			  xmax,(periodicflag & 3),numex,natex,-1.0);
  MR3calcvdw_ij(ni,xi,atypei,force,nj,xj,atypej,nat,gscalesf,rscales,2,
		xmax,(periodicflag & 1));
#else // else of calc separately

#if 1 // use Coulomb_vdw_factor in MR3calccoulomb_ij_exlist
  MR3calccoulomb_ij_exlist(ni,xi,qi,force,nj,xj,qj,
			   rscale,tblno,xmax,(periodicflag & 3)+4*multiplyq,
			   numex,natex);
#else // not use Coulomb_vdw_factor in MR3calccoulomb_ij_exlist
  {
    double cvf=Coulomb_vdw_factor;
    printf("debugging in MR3calccoulomb_vdw_ij_exlist\n");
    if((ftmp=(double *)MR3_malloc_pointer(sizeof(double)*ni*3,"MR3calccoulomb_vdw_ij_exlist_core"))==NULL){
      fprintf(stderr,"** error : can't malloc ftmp **\n");
      vg_exit(1);
    }
    bzero(ftmp,sizeof(double)*ni*3);
    Coulomb_vdw_factor=1.0;
    MR3calccoulomb_ij_exlist(ni,xi,qi,ftmp,nj,xj,qj,
			     rscale,tblno,xmax,(periodicflag & 3)+4*multiplyq,
			     numex,natex);
    for(i=0;i<ni*3;i++) force[i]+=ftmp[i]*cvf;
    Coulomb_vdw_factor=cvf;
    MR3_free_pointer(ftmp,"MR3calccoulomb_vdw_ij_exlist_core");
  }
#endif
#if 1
#if 1
  MR3calcvdw_ij_exlist(ni,xi,atypei,force,nj,xj,atypej,
		       nat,gscalesf,rscales,2,xmax,(periodicflag & 3),
		       numex,natex);
#else
  printf("debugging in MR3calccoulomb_vdw_ij_exlist\n");
  if((ftmp=(double *)MR3_malloc_pointer(sizeof(double)*ni*3,"MR3calccoulomb_vdw_ij_exlist_core"))==NULL){
    fprintf(stderr,"** error : can't malloc ftmp **\n");
    vg_exit(1);
  }
  bzero(ftmp,sizeof(double)*ni*3);
  MR3calcvdw_ij_exlist(ni,xi,atypei,ftmp,nj,xj,atypej,
		       nat,gscalesf,rscales,2,xmax,(periodicflag & 3),
		       numex,natex);
  for(i=0;i<ni;i++){
    if(i>=219 && i<224){
      printf("before_and_after_vdw %d %e %e\n",i+1,force[i*3],ftmp[i*3]);
    }
    for(j=0;j<3;j++) force[i*3+j]+=ftmp[i*3+j];
  }
  MR3_free_pointer(ftmp,"MR3calccoulomb_vdw_ij_exlist_core");
#endif
#else
  printf("MR3calcvdw_ij_exlist_host is called for force\n");
  MR3calcvdw_ij_exlist_host(ni,xi,atypei,force,nj,xj,atypej,
		       nat,gscalesf,rscales,2,xmax,(periodicflag & 3),
		       numex,natex);
#endif
#endif // end of calc separately

#if 0
  {
    static int ini=0;
    if(ini==0){
      printf("** warning : force is scaled by 0.1 **\n");
      ini=1;
    }
    for(i=0;i<ni;i++) for(j=0;j<3;j++) if(isnan(force[i*3+j])) printf("f[%d][%d] is NaN\n",i,j);
    //    for(i=0;i<ni;i++) for(j=0;j<3;j++) if(isnan(force[i*3+j]) force[i*3+j]*=0.1;
  }
#endif
}


void MR3calccoulomb_vdw_ij_exlist_core(int ni, double xi[], double qi[], 
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
#if !defined(MD_CELLINDEX) || !defined(MD_QAUNION_ATYPEBIT)
  MR3calccoulomb_vdw_ij_exlist_core_old(ni,xi,qi,atypei,force, 
					nj,xj,qj,atypej,
					nat,gscalesf,gscalesp,
					rscales,rscale,tblno,
					xmax,potc,potv,periodicflag, 
					potflag,changeflag,
					numex,natex);
#else
  double *ftmp,*forcep;
  int i,j,multiplyq=1,tblnop=tblno+1,overlapflag=0;
  double rcut,skinnb=0.0;
  VG_UNIT *unit=m3_get_unit();
  double volume[3]={0.0,0.0,0.0};
#define __USE_COULOMBVDW // define this to use MR3calccoulomb_vdw_ij_ci routine

  //  rcut=sqrt(MD_LJ_R2MAX);
  rcut=sqrt(unit->rcut2);
  setup_volume(xmax,volume);
  if((periodicflag & (1<<6))!=0) overlapflag=1;

  if((ftmp=(double *)MR3_malloc_pointer(sizeof(double)*ni*3,"MR3calccoulomb_vdw_ij_exlist_core"))==NULL){
    fprintf(stderr,"** error : can't malloc ftmp **\n");
    vg_exit(1);
  }

  if(tblno==6 || tblno==7) multiplyq=0;
  if(potflag!=0){
#ifdef MD_MEASURE_TIME
    vg_start_timer(53);
#endif
    //    *potc=*potv=0.0;
#if 1 // separate emu and 2-body
    bzero(ftmp,sizeof(double)*ni*3);
#if 1 // use coulomb and vdw emulator instead of separate emulator
    //    printf("using non-separate emulator in MR3calccoulomb_vdw_ij_exlist_core\n");
    MR3calccoulomb_vdw_nlist_ij_emu(ni,xi,qi,atypei,ftmp,nj,xj,qj,atypej,
                                    nat,gscalesp,rscales,rscale,
                                    tblnop,volume,(periodicflag & 3)+(1-multiplyq)*4,
                                    numex,natex,-1.0);
    //    for(i=0;i<ni;i++) (*potc)+=ftmp[i*3]-ftmp[i*3+1];
    for(i=0;i<ni;i++) (*potc)+=ftmp[i*3];
    for(i=0;i<ni;i++) (*potv)+=ftmp[i*3+1];
#else // use separate emulator
    MR3calccoulomb_nlist_ij_emu(ni,xi,qi,ftmp,nj,xj,qj,rscale,
    				tblnop,xmax,(periodicflag & 3)+(1-multiplyq)*4,
    				numex,natex,-1.0);
    //    for(i=0;i<ni;i++) printf("grape1: coulomb_emu[%d]=%e %e %e\n",i,ftmp[i*3],ftmp[i*3+1],ftmp[i*3+2]);
    for(i=0;i<ni;i++) (*potc)+=ftmp[i*3];
    bzero(ftmp,sizeof(double)*ni*3);
    MR3calcvdw_nlist_ij_emu(ni,xi,atypei,ftmp,nj,xj,atypej,nat,gscalesp,rscales,3,
			    xmax,(periodicflag & 3),numex,natex,-1.0);
    //    for(i=0;i<ni;i++) printf("grape1: vdw_emu[%d]=%e %e %e\n",i,ftmp[i*3],ftmp[i*3+1],ftmp[i*3+2]);
    //    for(i=0;i<3;i++) printf("emu potv[%d]=%e\n",i,ftmp[i*3]);
    for(i=0;i<ni;i++) (*potv)+=ftmp[i*3];
#endif // end of separate emulator
    //    printf("potc=%e potv=%e\n",*potc,*potv);
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(53);
#endif
    if(Ldim[0]>1 || Ldim[1]>1 || Ldim[2]>1){
#ifdef MD_MEASURE_TIME
      vg_start_timer(54);
#endif
#ifndef __USE_COULOMBVDW  // non coulomb_vdw
      bzero(ftmp,sizeof(double)*ni*3);
      MR3calccoulomb_ij_ci(ni,xi,qi,ftmp,nj,xj,qj,rscale,tblnop,rcut,skinnb,volume,Ldim);
      //      for(i=0;i<3;i++) printf("OK cpot[%d]=%e\n",i,ftmp[i*3]);
      for(i=0;i<ni;i++) (*potc)+=ftmp[i*3];
      bzero(ftmp,sizeof(double)*ni*3);
      MR3calcvdw_ij_ci(ni,xi,atypei,ftmp,nj,xj,atypej,nat,gscalesp,rscales,3,rcut,skinnb,volume,Ldim);
      //      for(i=0;i<3;i++) printf("OK vpot[%d]=%e\n",i,ftmp[i*3]);
      for(i=0;i<ni;i++) (*potv)+=ftmp[i*3];
#else
      bzero(ftmp,sizeof(double)*ni*3);
      MR3calccoulomb_vdw_ij_ci(ni,xi,qi,atypei,ftmp,nj,xj,qj,atypej,nat,gscalesp,rscales,rscale,tblnop,rcut,skinnb,volume,Ldim);
      //      for(i=0;i<ni;i++) printf("grape1: coulomb_vdw_ci[%d]=%e %e %e\n",i,ftmp[i*3],ftmp[i*3+1],ftmp[i*3+2]);
      //      for(i=0;i<3;i++) printf("NG cpot[%d]=%e vpot=%e\n",i,ftmp[i*3],ftmp[i*3+1]);
      for(i=0;i<ni;i++){
      	(*potc)+=ftmp[i*3];
	(*potv)+=ftmp[i*3+1];
      }
      //      printf("potc=%e potv=%e\n",*potc,*potv);
      //      bzero(ftmp,sizeof(double)*ni*3);
      //      MR3calcvdw_ij_ci(ni,xi,atypei,ftmp,nj,xj,atypej,nat,gscalesp,rscales,3,rcut,skinnb,volume,Ldim);
      //      for(i=0;i<ni;i++) (*potv)+=ftmp[i*3];
      //      printf("MR3calcvdw_ij_ci potc=%e potv=%e\n",*potc,*potv);
#endif
#ifdef MD_MEASURE_TIME
      vg_stop_and_accumulate_timer(54);
#endif
    }
    else{
      bzero(ftmp,sizeof(double)*ni*3);
      MR3calccoulomb_ij(ni,xi,qi,ftmp,nj,xj,qj,rscale,tblnop,xmax,(periodicflag & 1)+multiplyq*2);
      for(i=0;i<ni;i++) (*potc)+=ftmp[i*3];
      bzero(ftmp,sizeof(double)*ni*3);
      MR3calcvdw_ij(ni,xi,atypei,ftmp,nj,xj,atypej,nat,gscalesp,rscales,3,xmax,(periodicflag & 1));
      for(i=0;i<ni;i++) (*potv)+=ftmp[i*3];
    }
#else // use _exlist 
    bzero(ftmp,sizeof(double)*ni*3);
    MR3calccoulomb_ij_exlist(ni,xi,qi,ftmp,nj,xj,qj,
			     rscale,tblno+1,xmax,(periodicflag & 3)+4*multiplyq,
			     numex,natex);
    for(i=0;i<ni;i++) (*potc)+=ftmp[i*3];
    bzero(ftmp,sizeof(double)*ni*3);
    MR3calcvdw_ij_exlist(ni,xi,atypei,ftmp,nj,xj,atypej,
			 nat,gscalesp,rscales,3,xmax,(periodicflag & 3),
			 numex,natex);
    for(i=0;i<ni;i++) (*potv)+=ftmp[i*3];
#endif
  }

#if 0 // dummy routine
  bzero(ftmp,sizeof(double)*ni*3);
  MR3calccoulomb_nlist_ij_emu(ni,xi,qi,ftmp,nj,xj,qj,rscale,
			      tblno,xmax,(periodicflag & 3)+(1-multiplyq)*4,
			      numex,natex,-1.0);
  MR3calcvdw_nlist_ij_emu(ni,xi,atypei,ftmp,nj,xj,atypej,nat,gscalesf,rscales,2,
			  xmax,(periodicflag & 3),numex,natex,-1.0);
#endif

#ifdef MD_MEASURE_TIME
  vg_start_timer(55);
#endif
  bzero(ftmp,sizeof(double)*ni*3);
  if(Ldim[0]>1 || Ldim[1]>1 || Ldim[2]>1){
#ifndef __USE_COULOMBVDW // non coulomb_vdw
    MR3calccoulomb_ij_ci(ni,xi,qi,ftmp,nj,xj,qj,rscale,tblno,rcut,skinnb,volume,Ldim);
    MR3calcvdw_ij_ci(ni,xi,atypei,ftmp,nj,xj,atypej,nat,gscalesf,rscales,2,rcut,skinnb,volume,Ldim);
#else
    if(overlapflag!=0){
      unit->gpuoverlapflag=1;
      unit->ni_overlap=ni;
    }
    else{
      unit->gpuoverlapflag=0;
    }
    //    printf("unit->gpuoverlap is set %d in MR3calccoulomb_vdw_ij_exlist_core\n",unit->gpuoverlapflag);
    MR3calccoulomb_vdw_ij_ci(ni,xi,qi,atypei,ftmp,nj,xj,qj,atypej,nat,gscalesf,rscales,rscale,tblno,rcut,skinnb,volume,Ldim);
#endif
  }
  else{
    MR3calccoulomb_ij(ni,xi,qi,ftmp,nj,xj,qj,rscale,tblno,xmax,(periodicflag & 1)+multiplyq*2);
    MR3calcvdw_ij(ni,xi,atypei,ftmp,nj,xj,atypej,nat,gscalesf,rscales,2,xmax,(periodicflag & 1));
  }
#ifdef MD_MEASURE_TIME
  vg_stop_and_accumulate_timer(55);
  vg_start_timer(56);
#endif

#if 0 // dummy emulator
  {
    double *ftmp2=MR3_malloc_pointer(sizeof(double)*ni*3,"MR3calccoulomb_vdw_ij_exlist_core_old");
    bzero(ftmp2,sizeof(double)*ni*3);
    MR3calccoulomb_nlist_ij_emu(ni,xi,qi,ftmp2,nj,xj,qj,rscale,
				tblno,xmax,(periodicflag & 3)+(1-multiplyq)*4,
				numex,natex,-1.0);
    MR3calcvdw_nlist_ij_emu(ni,xi,atypei,ftmp2,nj,xj,atypej,nat,gscalesf,rscales,2,
			    xmax,(periodicflag & 3),numex,natex,-1.0);
    MR3_free_pointer(ftmp2,"MR3calccoulomb_vdw_ij_exlist_core_old");
  }
#endif
#if 0 // work before 090127
  MR3calccoulomb_nlist_ij_emu(ni,xi,qi,force,nj,xj,qj,rscale,
			      tblno,xmax,(periodicflag & 3)+(1-multiplyq)*4,
			      numex,natex,-1.0);
  MR3calcvdw_nlist_ij_emu(ni,xi,atypei,force,nj,xj,atypej,nat,gscalesf,rscales,2,
			  xmax,(periodicflag & 3),numex,natex,-1.0);
#endif
#if 0 // dummy new emulator
  {
    double *ftmp2=MR3_malloc_pointer(sizeof(double)*ni*3,"MR3calccoulomb_vdw_ij_exlist_core_old");
    bzero(ftmp2,sizeof(double)*ni*3);
    MR3calccoulomb_vdw_nlist_ij_emu(ni,xi,qi,atypei,ftmp2,nj,xj,qj,atypej,
				    nat,gscalesf,rscales,rscale,
				    tblno,xmax,(periodicflag & 3)+(1-multiplyq)*4,
				    numex,natex,-1.0);
    MR3_free_pointer(ftmp2,"MR3calccoulomb_vdw_ij_exlist_core_old");
  }
#endif
#if 1 // new routine
  MR3calccoulomb_vdw_nlist_ij_emu(ni,xi,qi,atypei,force,nj,xj,qj,atypej,
  				  nat,gscalesf,rscales,rscale,
  				  tblno,volume,(periodicflag & 3)+(1-multiplyq)*4,
  				  numex,natex,-1.0);
#endif
#ifdef MD_MEASURE_TIME
  vg_stop_and_accumulate_timer(56);
  vg_start_timer(57);
#endif
  if(overlapflag==0 || (overlapflag==1 && unit->gpuoverlapflag==0)){
    for(i=0;i<ni*3;i++) force[i]+=ftmp[i];
  }
  else{
    unit->gpuoverlapflag=2;
  }
  MR3_free_pointer(ftmp,"MR3calccoulomb_vdw_ij_exlist_core");
#endif
#ifdef MD_MEASURE_TIME
  vg_stop_and_accumulate_timer(57);
#endif
}


static void *MR3calccoulomb_vdw_ij_exlist_theread(void *thread_arg){
  thread_arg_coulombvdwijexlist *arg;

  arg=(thread_arg_coulombvdwijexlist *)thread_arg;
  MR3calccoulomb_vdw_ij_exlist_core(arg->ni,arg->xi,arg->qi,arg->atypei,arg->force, 
				    arg->nj,arg->xj,arg->qj,arg->atypej,
				    arg->nat,arg->gscalesf,arg->gscalesp,
				    arg->rscales,arg->rscale,
				    arg->tblno,arg->xmax,arg->potc,arg->potv,
				    arg->periodicflag,arg->potflag,arg->changeflag,
				    arg->numex,arg->natex);
  return NULL;
}


void MR3calccoulomb_vdw_ij_exlist(int ni, double xi[], double qi[], 
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
  /* tblno 0 ---------- coulomb 
           6 ---------- real space of coulomb (not supported now)
     periodicflag
           bit 0 : 0 -- non periodic
                   1 -- periodic
           bit 1 : 0 -- natex does not include duplicate list
                   1 -- natex includes duplicate list
           bit 12: 0 -- no overlap          
                   1 -- overlap with pthread 
  */
  double *ftmp;
  int i,j,overlapflag=0;
  static thread_arg_coulombvdwijexlist thread_arg;
  //  pthread_t thread;
  int ret,s;
  void *thread_ret;
  VG_UNIT *unit;
  static int *numexlocal=NULL,*natexlocal=NULL;
  //  double potclocal=0.0,potvlocal=0.0;
  static int *atypeilocal=NULL,*atypejlocal=NULL;

#ifdef MD_MEASURE_TIME
  vg_start_timer(50);
#endif
  //  printf("in MR3calccoulomb_vdw_ij_exlist : ni=%d nj=%d nat=%d tblno=%d xmax=%e periodicflag=%d potflag=%d changeflag=%d\n",ni,nj,nat,tblno,xmax,periodicflag,potflag,changeflag);
  if((periodicflag & (1<<12))!=0 && potflag==0) overlapflag=1;
  //  if(potflag==0) overlapflag=1;printf("overlapflag is set manually\n");
#ifdef MD_PRINT_WARN
  {
    static int ini=0;
    if(ini==0){
      printf("** warning : overlap is used only when potflag=0 **\n");
      ini=1;
    }
  }
#endif

  if(overlapflag==0){
#ifdef MD_MEASURE_TIME
    vg_start_timer(51);
#endif
    MR3calccoulomb_vdw_ij_exlist_core(ni,xi,qi,atypei,force,nj,xj,qj,atypej,
				      nat,gscalesf,gscalesp,rscales,rscale,
				      tblno,xmax,potc,potv,periodicflag, 
				      potflag,changeflag,numex,natex);
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(51);
#endif
  }
  else{
    static int ini=0;
#ifdef MD_MEASURE_TIME
    vg_start_timer(52);
#endif
    if(ini==0){
      printf("overlap is specified in MR3calccoulomb_vdw_ij_exlist\n");
      ini=1;
    }
#if 0
    if(M3_unit==NULL){
      if((M3_unit=(M3_UNIT *)MR3_malloc_pointer(sizeof(M3_UNIT),"M3_unit in MR3calccoulomb_vdw_ij_exlist"))==NULL){
	fprintf(stderr,"** error : can't malloc M3_unit **\n");
	MR3_exit(1);
      }
    }
    m3_initialize_unit_variables(M3_unit);
#endif
    unit=m3_get_unit();
    if(unit==NULL){
      fprintf(stderr,"** error : unit is NULL **\n");
      vg_exit(1);
    }
    MR3_free_pointer(unit->fthread,"MR3calccoulomb_vdw_ij_exlist");
    if((unit->fthread=(double *)MR3_malloc_pointer(sizeof(double)*ni*3,"MR3calccoulomb_vdw_ij_exlist"))==NULL){
      fprintf(stderr,"** error : can't malloc unit->fthread in MR3calccoulomb_vdw_ij_exlist **\n");
      MR3_exit(1);
    }
    /*    {
      static int ini=0;
      if(ini==0){
	printf("** warning : size of fthread is multiplied by 2 temporally to ensure large enought **\n");
	ini=1;
      }
      }*/
    unit->ni_overlap=ni;
    //    printf("after malloc of fthread, unit->fthread=%016llx, ni_overlap=%d\n",unit->fthread,unit->ni_overlap);
    for(i=0;i<ni*3;i++) unit->fthread[i]=force[i];

    MR3_free_pointer(numexlocal,"MR3calccoulomb_vdw_ij_exlist");
    MR3_free_pointer(natexlocal,"MR3calccoulomb_vdw_ij_exlist");
    if((numexlocal=(int *)MR3_malloc_pointer(sizeof(int)*ni,"MR3calccoulomb_vdw_ij_exlist"))==NULL){
      fprintf(stderr,"** error : can't malloc numexlocal **\n");
      MR3_exit(1);
    }
    for(i=s=0;i<ni;i++){
      numexlocal[i]=numex[i];
      s+=numex[i];
    }
    if((natexlocal=(int *)MR3_malloc_pointer(sizeof(int)*s,"MR3calccoulomb_vdw_ij_exlist"))==NULL){
      fprintf(stderr,"** error : can't malloc natexlocal **\n");
      MR3_exit(1);
    }
    for(i=0;i<s;i++) natexlocal[i]=natex[i];

    MR3_free_pointer(atypeilocal,"MR3calccoulomb_vdw_ij_exlist");
    MR3_free_pointer(atypejlocal,"MR3calccoulomb_vdw_ij_exlist");
    if((atypeilocal=(int *)MR3_malloc_pointer(sizeof(int)*ni,"MR3calccoulomb_vdw_ij_exlist"))==NULL){
      fprintf(stderr,"** error : can't malloc atypeilocal **\n");
      MR3_exit(1);
    }
    for(i=0;i<ni;i++) atypeilocal[i]=atypei[i];
    if((atypejlocal=(int *)MR3_malloc_pointer(sizeof(int)*nj,"MR3calccoulomb_vdw_ij_exlist"))==NULL){
      fprintf(stderr,"** error : can't malloc atypejlocal **\n");
      MR3_exit(1);
    }
    for(i=0;i<nj;i++) atypejlocal[i]=atypej[i];
    
    thread_arg.ni=ni;
    thread_arg.xi=xi;
    thread_arg.qi=qi;
    //    thread_arg.atypei=atypei;
    thread_arg.atypei=atypeilocal;
    thread_arg.force=unit->fthread;
    thread_arg.nj=nj;
    thread_arg.xj=xj;
    thread_arg.qj=qj;
    //    thread_arg.atypej=atypej;
    thread_arg.atypej=atypejlocal;
    thread_arg.nat=nat;
    thread_arg.gscalesf=gscalesf;
    thread_arg.gscalesp=gscalesp;
    thread_arg.rscales=rscales;
    thread_arg.rscale=rscale;
    thread_arg.tblno=tblno;
    thread_arg.xmax=xmax;
    //    thread_arg.potc=potc;
    //    thread_arg.potv=potv;
    thread_arg.potc=&(unit->potc);
    thread_arg.potv=&(unit->potv);
    thread_arg.periodicflag=periodicflag & 3;
    thread_arg.potflag=potflag;
    thread_arg.changeflag=changeflag;
    //    thread_arg.numex=numex;
    //    thread_arg.natex=natex;
    thread_arg.numex=numexlocal;
    thread_arg.natex=natexlocal;
    ret=pthread_create(&(unit->thread),NULL,MR3calccoulomb_vdw_ij_exlist_theread,
		       &thread_arg);
    if(ret){
      fprintf(stderr,"** can't cread therad **\n");
      MR3_exit(1);
    }
    {
      static int ini=0;
      if(ini==0){
	printf("pthread_create is called in MR3calccoulomb_vdw_ij_exlist\n");fflush(stdout);
	ini=1;
      }
    }

#if 0
    ret=pthread_join(unit->thread,&thread_ret);
    if(ret){
      fprintf(stderr,"** can't join therad **\n");
      MR3_exit(1);
    }
#endif
#if 0
    printf("MR3_get_forces_overlap and get_potc_and_potv is called in MR3calccoulomb_vdw_ij_exlist\n");
    MR3_get_forces_overlap(force);
    MR3_get_potc_and_potv_overlap(potc,potv);
#endif
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(52);
#endif
  }
#ifdef MD_MEASURE_TIME
  vg_stop_and_accumulate_timer(50);
#endif
}


void mr3calccoulomb_vdw_ij_exlist_(int *ni, double xi[], double qi[], 
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
  static int *atypej3=NULL;
  int i;

#if 0 // for debugging, works
  printf("Fortran library is used in mr3calccoulomb_vdw_ij_exlist_\n");
  {
    double cvf=Coulomb_vdw_factor,*ftmp;
    int modop,iperiod=1;
    Coulomb_vdw_factor=1.0;
    if((ftmp=(double *)MR3_malloc_pointer(sizeof(double)*(*ni)*3,"mr3calccoulomb_vdw_ij_exlist_"))==NULL){
      fprintf(stderr,"** error : can't malloc ftmp **\n");
      vg_exit(1);
    }

    if((*potflag) & 1){
      // coulomb potential
      modop=(*tblno)+1;
      bzero(ftmp,sizeof(double)*(*ni)*3);
      mr3calccoulomb_ij_exlist_(ni,xi,qi,ftmp,nj,xj,qj,rscale,
				&modop,xmax,&iperiod,numex,natex);
      for(i=0,*potc=0.0;i<*ni;i++) (*potc)+=cvf*ftmp[i*3];
    }
    // coulomb force
    modop=*tblno;
    bzero(ftmp,sizeof(double)*(*ni)*3);
    mr3calccoulomb_ij_exlist_(ni,xi,qi,ftmp,nj,xj,qj,rscale,
			      &modop,xmax,&iperiod,numex,natex);
    for(i=0;i<(*ni)*3;i++) force[i]+=cvf*ftmp[i];

    if((*potflag) & 1){
      // vdw potential
      modop=3;
      bzero(ftmp,sizeof(double)*(*ni)*3);
      mr3calcvdw_ij_exlist_(ni,xi,atypei,ftmp,nj,xj,atypej,
			    nat,gscalesp,rscales,
			    &modop,xmax,&iperiod,numex,natex);
      for(i=0,*potv=0.0;i<*ni;i++) (*potv)+=ftmp[i*3];
    }
    // vdw force
    modop=2;
    bzero(ftmp,sizeof(double)*(*ni)*3);
    mr3calcvdw_ij_exlist_(ni,xi,atypei,ftmp,nj,xj,atypej,
			  nat,gscalesf,rscales,
			  &modop,xmax,&iperiod,numex,natex);
    for(i=0;i<(*ni)*3;i++) force[i]+=ftmp[i];

    Coulomb_vdw_factor=cvf;
    MR3_free_pointer(ftmp,"mr3calccoulomb_vdw_ij_exlist_");
  }
  return;
#endif

#if 1 // malloc
  for(i=s=0;i<*ni;i++){
    s+=numex[i];
  }
  //  printf("allocating ni=%d, nj=%d, s=%d \n",*ni,*nj,s);
  if((atypei2=(int *)MR3_malloc_pointer(sizeof(int)*(*ni),"atypei2 in mr3calccoulomb_vdw_ij_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc atypei2 in mr3calccoulomb_vdw_ij_exlist_ **\n");
    MR3_exit(1);
  }
  if((atypej2=(int *)MR3_malloc_pointer(sizeof(int)*(*nj),"atypej2 in mr3calccoulomb_vdw_ij_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc atypej2 in mr3calccoulomb_vdw_ij_exlist_ **\n");
    MR3_exit(1);
  }
  if((natex2=(int *)MR3_malloc_pointer(sizeof(int)*s,"natex2 in mr3calccoulomb_vdw_ij_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc natex2 in mr3calccoulomb_ij_exlist_ **\n");
    MR3_exit(1);
  }
#endif
#if 1
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
      atypej2[i]=0;
    }
  }
  for(i=0;i<s;i++) natex2[i]=natex[i]-1;
#endif
#if 0
  if(atypej3==NULL){
    if((atypej3=(int *)MR3_malloc_pointer(sizeof(int)*(*nj),"mr3calccoulomb_vdw_ij_exlist_"))==NULL){
      fprintf(stderr,"** error : can't malloc atypej3 **\n");
      MR3_exit(1);
    }
  }
  for(i=0;i<*nj;i++) atypej3[i]=atypej2[i];
#endif

#if 1 // call c routine
  MR3calccoulomb_vdw_ij_exlist(*ni,xi,qi,atypei2,force, 
			       *nj,xj,qj,
			       atypej2,
			       //			       atypej3,
			       *nat,gscalesf,gscalesp,rscales,*rscale,
			       *tblno,*xmax,potc,potv,
			       *periodicflag,*potflag,*changeflag,
			       numex,natex2);
#else // use separate library
  printf("Separate library is used in mr3calccoulomb_vdw_ij_exlist_\n");
  {
    double cvf=Coulomb_vdw_factor,*ftmp;
    int modop,iperiod=1;
    Coulomb_vdw_factor=1.0;
    if((ftmp=(double *)MR3_malloc_pointer(sizeof(double)*(*ni)*3,"mr3calccoulomb_vdw_ij_exlist_"))==NULL){
      fprintf(stderr,"** error : can't malloc ftmp **\n");
      vg_exit(1);
    }

    if((*potflag) & 1){
      // coulomb potential
      modop=(*tblno)+1;
      bzero(ftmp,sizeof(double)*(*ni)*3);
#if 1
      mr3calccoulomb_ij_exlist_(ni,xi,qi,ftmp,nj,xj,qj,rscale,
				&modop,xmax,&iperiod,numex,natex);
#else
      MR3calccoulomb_ij_exlist(*ni,xi,qi,ftmp,*nj,xj,qj,*rscale,
			       modop,*xmax,iperiod,numex,natex2);
#endif
      for(i=0,*potc=0.0;i<*ni;i++) (*potc)+=cvf*ftmp[i*3];
    }
    // coulomb force
    modop=*tblno;
    bzero(ftmp,sizeof(double)*(*ni)*3);
#if 1 
    mr3calccoulomb_ij_exlist_(ni,xi,qi,ftmp,nj,xj,qj,rscale,
			      &modop,xmax,&iperiod,numex,natex);
#else
    MR3calccoulomb_ij_exlist(*ni,xi,qi,ftmp,*nj,xj,qj,*rscale,
			     modop,*xmax,iperiod,numex,natex2);
#endif
    for(i=0;i<(*ni)*3;i++) force[i]+=cvf*ftmp[i];
    for(i=50;i<55;i++) printf("fc[%d]=%e %e %e\n",i,ftmp[i*3],ftmp[i*3+1],ftmp[i*3+2]);

    if((*potflag) & 1){
      // vdw potential
      modop=3;
      bzero(ftmp,sizeof(double)*(*ni)*3);
#if 1
      mr3calcvdw_ij_exlist_(ni,xi,atypei,ftmp,nj,xj,atypej,
			    nat,gscalesp,rscales,
			    &modop,xmax,&iperiod,numex,natex);
#else
      MR3calcvdw_ij_exlist(*ni,xi,atypei2,ftmp,*nj,xj,atypej2,
			   *nat,gscalesp,rscales,
			   modop,*xmax,iperiod,numex,natex2);
#endif
      for(i=0,*potv=0.0;i<*ni;i++) (*potv)+=ftmp[i*3];
    }
    // vdw force
    modop=2;
    bzero(ftmp,sizeof(double)*(*ni)*3);
#if 1
    mr3calcvdw_ij_exlist_(ni,xi,atypei,ftmp,nj,xj,atypej,
			  nat,gscalesf,rscales,
			  &modop,xmax,&iperiod,numex,natex);
#else
    MR3calcvdw_ij_exlist(*ni,xi,atypei2,ftmp,*nj,xj,atypej2,
			 *nat,gscalesf,rscales,
			 modop,*xmax,iperiod,numex,natex2);
#endif
    for(i=0;i<(*ni)*3;i++) force[i]+=ftmp[i];
    for(i=50;i<55;i++) printf("fv[%d]=%e %e %e\n",i,ftmp[i*3],ftmp[i*3+1],ftmp[i*3+2]);

    Coulomb_vdw_factor=cvf;
    MR3_free_pointer(ftmp,"mr3calccoulomb_vdw_ij_exlist_");
  }
#endif // end of use separate library
#if 1
  MR3_free_pointer(atypei2,"atypei2 in mr3calccoulomb_vdw_ij_exlist_");
  MR3_free_pointer(natex2,"natex2 in mr3calccoulomb_vdw_ij_exlist_");
  MR3_free_pointer(atypej2,"atypej2 in mr3calccoulomb_vdw_ij_exlist_");
#endif
}


void mr3calccoulomb_vdw_ij_exlist_virial_(int *ni, double xi[], double qi[], 
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
  /*
    this routine does not work. 
   */
  int *atypei2,*atypej2,*index;
  int *natex2,s;
  static int *atypej3=NULL;
  int i,j;
  double *xi2,volume[3],xtest[3],*force2,*qi2;
  int ni2,ci[3],*numex2,*numex_offset;
  VG_UNIT *unit=m3_get_unit();
  double rcut=sqrt(unit->rcut2)+3.0;// margin is 3.0A

  for(i=0;i<3;i++) printf("** this routines does't work yet **\n");
#if 1
  {
    static int ini=0;
    if(ini==0){
      printf("** mr3calccoulomb_vdw_ij_exlist_virial is called instead of mr3calccoulomb_vdw_ij_exlist **\n");
      printf("** you need to modify two parts in grape.src to back to mr3calccoulomb_vdw_ij_exlist **\n");
      ini=1;
    }
  }
#endif
  setup_volume(*xmax,volume);
  // count the number of image atoms
  ni2=*ni;
  for(i=0;i<*ni;i++){
    for(ci[0]=-1;ci[0]<=1;ci[0]++){
      for(ci[1]=-1;ci[1]<=1;ci[1]++){
	for(ci[2]=-1;ci[2]<=1;ci[2]++){
	  if(ci[0]!=0 || ci[1]!=0 || ci[2]!=0){
	    for(j=0;j<3;j++){
	      xtest[j]=xi[i*3+j]+ci[j]*volume[j];
	    }
	    if(xtest[0]>=-rcut && xtest[0]<volume[0]+rcut &&
	       xtest[1]>=-rcut && xtest[1]<volume[1]+rcut &&
	       xtest[2]>=-rcut && xtest[2]<volume[2]+rcut){
	      ni2++;
	    }
	  }
	}
      }
    }
  }
  //  printf("ni2=%d\n",ni2);
  if((xi2=(double *)MR3_malloc_pointer(sizeof(double)*ni2*3,"xi2 in mr3calccoulomb_vdw_ij_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc xi2 in mr3calccoulomb_vdw_ij_exlist_ **\n");
    MR3_exit(1);
  }
  if((force2=(double *)MR3_malloc_pointer(sizeof(double)*ni2*3,"force2 in mr3calccoulomb_vdw_ij_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc force2 in mr3calccoulomb_vdw_ij_exlist_ **\n");
    MR3_exit(1);
  }
  if((qi2=(double *)MR3_malloc_pointer(sizeof(double)*ni2,"qi2 in mr3calccoulomb_vdw_ij_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc qi2 in mr3calccoulomb_vdw_ij_exlist_ **\n");
    MR3_exit(1);
  }
  if((atypei2=(int *)MR3_malloc_pointer(sizeof(int)*ni2,"atypei2 in mr3calccoulomb_vdw_ij_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc atypei2 in mr3calccoulomb_vdw_ij_exlist_ **\n");
    MR3_exit(1);
  }
  if((index=(int *)MR3_malloc_pointer(sizeof(int)*ni2,"index in mr3calccoulomb_vdw_ij_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc index in mr3calccoulomb_vdw_ij_exlist_ **\n");
    MR3_exit(1);
  }
  if((numex2=(int *)MR3_malloc_pointer(sizeof(int)*ni2,"numex2 in mr3calccoulomb_vdw_ij_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc numex2 in mr3calccoulomb_vdw_ij_exlist_ **\n");
    MR3_exit(1);
  }

  // copy the image atoms
  for(i=0;i<*ni;i++){
    for(j=0;j<3;j++) xi2[i*3+j]=xi[i*3+j];
    qi2[i]=qi[i];
    index[i]=i;
    if(atypei[i]>0){
      atypei2[i]=atypei[i]-1;
    }
    else{
      printf("  warning : atypei[%d]=%d should be positive\n",
	     i,atypei[i]);
      atypei2[i]=0;
    }
    numex2[i]=numex[i];
  }
  ni2=*ni;
  for(i=0;i<*ni;i++){
    for(ci[0]=-1;ci[0]<=1;ci[0]++){
      for(ci[1]=-1;ci[1]<=1;ci[1]++){
	for(ci[2]=-1;ci[2]<=1;ci[2]++){
	  if(ci[0]!=0 || ci[1]!=0 || ci[2]!=0){
	    for(j=0;j<3;j++){
	      xtest[j]=xi[i*3+j]+ci[j]*volume[j];
	    }
	    if(xtest[0]>=-rcut && xtest[0]<volume[0]+rcut &&
	       xtest[1]>=-rcut && xtest[1]<volume[1]+rcut &&
	       xtest[2]>=-rcut && xtest[2]<volume[2]+rcut){
	      for(j=0;j<3;j++) xi2[ni2*3+j]=xtest[j];
	      qi2[ni2]=qi[i];
	      index[ni2]=i;
	      atypei2[ni2]=atypei2[i];
	      numex2[ni2]=numex[i];
	      //	      if(ni2>=44950 && ni2<44960) printf("ni2=%d numex2=%d\n",ni2,numex[i]);
	      ni2++;
	    }
	  }
	}
      }
    }
  }

  if((numex_offset=(int *)MR3_malloc_pointer(sizeof(int)*(*ni),"natex_offset in mr3calccoulomb_vdw_ij_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc natex_offset in mr3calccoulomb_vdw_ij_exlist_ **\n");
    MR3_exit(1);
  }
  if((atypej2=(int *)MR3_malloc_pointer(sizeof(int)*(*nj),"atypej2 in mr3calccoulomb_vdw_ij_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc atypej2 in mr3calccoulomb_vdw_ij_exlist_ **\n");
    MR3_exit(1);
  }
  for(i=s=0;i<ni2;i++){
    s+=numex2[i];
  }
  if((natex2=(int *)MR3_malloc_pointer(sizeof(int)*s,"natex2 in mr3calccoulomb_vdw_ij_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc natex2 in mr3calccoulomb_ij_exlist_ **\n");
    MR3_exit(1);
  }

  for(i=0;i<*nj;i++){
    if(atypej[i]>0){
      atypej2[i]=atypej[i]-1;
    }
    else{
      printf("  warning : atypej[%d]=%d should be positive\n",
	     i,atypej[i]);
      atypej2[i]=0;
    }
  }

  s=0;
  numex_offset[0]=0;
  for(i=0;i<*ni;i++){
    if(i>=1) numex_offset[i]=numex_offset[i-1]+numex[i-1];
    for(j=0;j<numex[i];j++){
      natex2[s]=natex[s]-1;
      s++;
    }
  }
  for(i=*ni;i<ni2;i++){
    for(j=0;j<numex2[i];j++){
      natex2[s]=natex2[numex_offset[index[i]]+j];
      //      if(i<*ni+1000) printf("i=%d j=%d natex2[%d]=%d\n",i,j,s,natex2[s]);
      s++;
    }
  }
  for(i=0;i<ni2;i++) for(j=0;j<3;j++) force2[i*3+j]=0.0;
  MR3calccoulomb_vdw_ij_exlist(ni2,xi2,qi2,atypei2,force2, 
			       *nj,xj,qj,
			       atypej2,
			       //			       atypej3,
			       *nat,gscalesf,gscalesp,rscales,*rscale,
			       *tblno,*xmax,potc,potv,
			       *periodicflag,*potflag,*changeflag,
			       numex2,natex2);
  for(i=0;i<*ni;i++) for(j=0;j<3;j++) force[i*3+j]=0.0;
  for(i=0;i<ni2;i++){
    for(j=0;j<3;j++){
      force[index[i]*3+j]+=force2[i*3+j];
    }
  }

  MR3_free_pointer(atypei2,"atypei2 in mr3calccoulomb_vdw_ij_exlist_");
  MR3_free_pointer(qi2,"qi2 in mr3calccoulomb_vdw_ij_exlist_");
  MR3_free_pointer(force2,"qi2 in mr3calccoulomb_vdw_ij_exlist_");
  MR3_free_pointer(natex2,"natex2 in mr3calccoulomb_vdw_ij_exlist_");
  MR3_free_pointer(atypej2,"atypej2 in mr3calccoulomb_vdw_ij_exlist_");
  MR3_free_pointer(xi2,"xi2 in mr3calccoulomb_vdw_ij_exlist_");
  MR3_free_pointer(index,"index in mr3calccoulomb_vdw_ij_exlist_");
  MR3_free_pointer(numex2,"numex2 in mr3calccoulomb_vdw_ij_exlist_");
  MR3_free_pointer(numex_offset,"numex_offset in mr3calccoulomb_vdw_ij_exlist_");
}


void MR3calccoulomb_vdw_exlist(int n, double x[], double q[], 
			       int atype[], int nat, double force[],
			       double rscale, 
			       double gscalesf[], double gscalesp[],
			       double rscales[],
			       int tblno, int numex[], int natex[],
			       double volume1, 
			       double *potc, double *potv, 
			       int *convert_index, 
			       int potflag, int flag)
{
  /*
     flag         bit 0: 0 --- non periodic
                         1 --- periodic
                  bit 1: 0 --- natex does not include duplicate list
                         1 --- natex includes duplicate list
   */
#if 1
  /*
  printf("MR3calccoulomb_vdw_ij_exlist_host is called from MR3calccoulomb_vdw_exlist\n");
  printf("flag=%d potflag=%d\n",flag,potflag);
  MR3calccoulomb_vdw_ij_exlist_host(n,x,q,atype,force, 
				    n,x,q,atype,
				    nat,gscalesf,gscalesp,rscales,rscale,tblno,
				    volume1,potc,potv,0,1,0,
				    numex,natex);
  return;
  */
  MR3calccoulomb_vdw_ij_exlist(n,x,q,atype,force,n,x,q,atype,
			       nat,gscalesf,gscalesp,rscales,rscale,tblno,
			       volume1,potc,potv,flag,potflag,0,
			       numex,natex);
#else
  printf("** this routine has not been tested\n");
  MR3calccoulomb_ij(n,x,q,force,n,x,q,rscale,tblno,volume1,(flag & 1)+2);
  MR3calccoulomb_nlist_emu(x,n,q,rscale,tblno,volume1,(flag & 1),
			   numex,natex,-1.0,force);
  MR3calcvdw_ij(n,x,atype,force,n,x,atype,nat,gscalesf,rscales,
		2+(tblno % 2),volume1,(flag & 1));
  MR3calcvdw_nlist_emu2(x,n,atype,nat,gscalesf,rscales,2+(tblno % 2),
			volume1,(flag & 1),numex,natex,-1.0,force);
#endif
}


void MR3calccoulomb_vdw_exlist_from_cell(int n, double x[], double q[], 
			       int atype[], int nat, double force[],
			       double rscale, 
			       double gscalesf[], double gscalesp[],
			       double rscales[],
			       int tblno, int numex[], int natex[],
			       double volume1, 
			       double *potc, double *potv, 
			       int *convert_index, 
			       int potflag, int flag)
{
  /*
    flag bit 1 : 0 -- no overlap
                 1 -- overlap  
   */
  //  printf("MR3calccoulomb_vdw_exlist is called\n");
  double potc_nlist=0.0,potv_nlist=0.0;
  int i,j,nf_per=2;
  //  int changeflag=0;
  double (*ftmp)[3];
  int overlapflag=0;

#define CALC_NLIST_WITH_SPE

  {
    static int ini=0;
    static double volume1_bak=0.0;
    if(ini==0){
      double min=1e20,max=-1e20;
      for(i=0;i<n;i++){
	for(j=0;j<3;j++){
	  if(min>x[i*3+j]) min=x[i*3+j];
	  if(max<x[i*3+j]) max=x[i*3+j];
	}
      }
      printf("corrdinate min=%e max=%e ",min,max);
      printf("volume1=%f",volume1);
      //      volume1=200.0;
      volume1=(max-min)*1.5*2.0;
      printf(" is changed to %f\n",volume1);
      volume1_bak=volume1;
      ini=1;
    }
    else{
      volume1=volume1_bak;
    }
  }
    
  //  printf("rscales[0]=%e at the beggining of MR3calccoulomb_vdw_exlist\n",rscales[0]);
  if((flag & 2)!=0) overlapflag=2;

#if defined(CALC_NLIST_WITH_SPE) 
#if 0
  MR3calccoulomb_vdw_nlist_host((double (*)[3])x,n,q,0.0,
				atype,nat,gscalesf,gscalesp,rscales,
				tblno,volume1,nf_per,
				numex,natex,
				&potc_nlist,&potv_nlist,
				ftmp);
#endif
#if 0
  printf("after nlist_host potc_nlist=%e potv_nlist=%e\n",
	 potc_nlist,potv_nlist);
#endif
#else // else of CALC_NLIST_WITH_SPE
  if((ftmp=(double (*)[3])MR3_malloc_pointer(sizeof(double)*n*3,"ftmp in MR3calccoulomb_vdw_exlist"))==NULL){
    fprintf(stderr,"** error : can't malloc ftmp in MR3calccoulomb_vdw_exlist\n");
    MR3_exit(1);
  }
  for(i=0;i<n;i++) ftmp[i][0]=ftmp[i][1]=ftmp[i][2]=0.0;
  //  for(i=0;i<nat;i++) for(j=0;j<nat;j++) if((gscalesf[i*nat+j]==0.0 || gscalesp[i*nat+j]==0.0) && rscales[i*nat+j]!=1.0) printf("i=%d j=%d gscalesf=%e gscalesp=%e rscales=%e\n",i,j,gscalesf[i*nat+j],gscalesp[i*nat+j],rscales[i*nat+j]);

  MR3calccoulomb_vdw_nlist_host((double (*)[3])x,n,q,rscale,
				atype,nat,gscalesf,gscalesp,rscales,
				tblno,volume1,nf_per,
				numex,natex,
				&potc_nlist,&potv_nlist,
				ftmp);
  printf("after nlist_host potc_nlist=%e potv_nlist=%e\n",
	 potc_nlist,potv_nlist);
  //  printf("rscales[0]=%e\n",rscales[0]);
#endif // end of CALC_NLIST_WITH_SPE

  {
    static int ini=0,*numex2=NULL,*natex2=NULL; 
#ifdef PS3_MPI
    int ni,ni_offset,natex2_offset;
    double total_potc=0.0,total_potv=0.0;
    int ni3s[PS3_MPI_MAX_PROCS],ni_offset3s[PS3_MPI_MAX_PROCS];
#endif
    if(ini==0){
#ifdef PS3_MPI
      int argc=0;
      char **argv=NULL;
      MPI_Init(&argc,&argv);
      MPI_Comm_rank(MPI_COMM_WORLD,&Myrank);
      MPI_Comm_size(MPI_COMM_WORLD,&Nprocs);
      printf("rank=%d/%d> MPI initialized\n",Myrank,Nprocs);
#endif
      make_natex2(n,numex,natex,&numex2,&natex2,1);
      ini=1;
    }
    //    printf("natex2[0]=%d\n",natex2[0]);

    *potc=0.0;
    *potv=0.0;
#ifdef PS3_MPI
    ni_offset=0;
    for(i=0;i<=Myrank;i++){
      if(i>0) ni_offset+=ni;
      ni=n/Nprocs+(i<(n % Nprocs) ? 1:0);
    }
    natex2_offset=0;
    for(i=0;i<ni_offset;i++) natex2_offset+=numex2[i];
    //    printf("rank=%d> ni=%d ni_offset=%d natex2_offset=%d\n",Myrank,ni,ni_offset,natex2_offset);
    for(i=0;i<Nprocs;i++) ni3s[i]=n/Nprocs+(i<(n % Nprocs) ? 1:0);
    ni_offset3s[0]=0;
    for(i=1;i<Nprocs;i++) ni_offset3s[i]=ni_offset3s[i-1]+ni3s[i-1];
    for(i=0;i<Nprocs;i++){
      ni3s[i]=ni3s[i]*3;
      ni_offset3s[i]=ni_offset3s[i]*3;
    }
    //    for(i=0;i<Nprocs;i++) printf("rank=%d> ni3s[%d]/3=%d ni_offset3s/3=%d\n",Myrank,i,ni3s[i]/3,ni_offset3s[i]/3);

    MR3calccoulomb_vdw_ij_exlist_cell(ni,x+ni_offset*3,q+ni_offset,
				      atype+ni_offset,force+ni_offset*3,
				      n,x,q,atype,
				      nat,gscalesf,gscalesp,rscales,
				      rscale,tblno,
				      volume1,potc,potv,
				      numex2+ni_offset,natex2+natex2_offset,
				      nf_per,potflag,overlapflag);
    //    for(i=0;i<5;i++) printf("rank=%d> force[%d]=%e %e %e\n",Myrank,ni_offset+i,force[(ni_offset+i)*3],force[(ni_offset+i)*3+1],force[(ni_offset+i)*3+2]);
    //    printf("rank=%d> potc=%e potv=%e\n",Myrank,*potc,*potv);
    MPI_Allreduce(potc,&total_potc,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(potv,&total_potv,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    //    printf("rank=%d> total_potc=%e total_potv=%e\n",Myrank,total_potc,total_potv);
    *potc=total_potc;
    *potv=total_potv;

    MPI_Allgatherv(force+ni_offset*3,ni*3,MPI_DOUBLE,
		   force,ni3s,ni_offset3s,MPI_DOUBLE,MPI_COMM_WORLD);
    //    for(i=0;i<5;i++) printf("rank=%d> force[%d]=%e %e %e\n",Myrank,i,force[i*3],force[i*3+1],force[i*3+2]);
#else
    printf("this routine does not exist\n");
    /*
    MR3calccoulomb_vdw_exlist_cell(n,x,q,atype,nat,force,
				   rscale,gscalesf,gscalesp,rscales,
				   tblno,numex2,natex2,volume1,potc,potv,
				   NULL,potflag,nf_per);
    */
#endif
  }
#if 0
  printf("after ij_cell potc=%e(%e) potv=%e(%e)\n",
	 *potc,(*potc)*0.5,*potv,(*potv)*0.5);
#endif

  //  for(i=0;i<nat;i++) for(j=0;j<nat;j++) if((gscalesf[i*nat+j]==0.0 || gscalesp[i*nat+j]==0.0) && rscales[i*nat+j]!=1.0) printf("i=%d j=%d gscalesf=%e gscalesp=%e rscales=%e\n",i,j,gscalesf[i*nat+j],gscalesp[i*nat+j],rscales[i*nat+j]);

#ifndef CALC_NLIST_WITH_SPE
  for(i=0;i<n;i++){
    //    if(i<3) printf("i=%d force=%e ftmp=%e\n",i,force[i*3+0],ftmp[i][0]);
    force[i*3]-=ftmp[i][0];
    force[i*3+1]-=ftmp[i][1];
    force[i*3+2]-=ftmp[i][2];
  }
  (*potc)-=potc_nlist*2.0;
  (*potv)-=potv_nlist*2.0;
#else
  //  printf("**** nlist is not added to ij_cell\n");
#endif

#if 0
  MR3_free_pointer(ftmp,"MR3calccoulomb_vdw_exlist_from_cell");
  return;
#endif

#ifndef CALC_NLIST_WITH_SPE
  MR3_free_pointer(ftmp,"ftmp in MR3calccoulomb_vdw_exlist");
#endif

#if 0
  printf("********* forces are set zero\n");
  for(i=0;i<n;i++) force[i*3]=force[i*3+1]=force[i*3+2]=0.0;
#endif
}


static void MR3calccoulomb_vdw_exlist_host2(int *n, double x[], double q[], 
					    int atype[], int *nat, double force[],
					    double *rscale, 
					    double gscalesf[], double gscalesp[],
					    double rscales[],
					    int *tblno, int numex[], int natex[],
					    double *volume1,
					    double *potc, double *potv, 
					    int *convert_index, 
					    int *potflag, int *flag)
{
  int i,j,atij,k,jj,iexcl;
  double gsf,gsp,rs,r,dr[3],dtmp,rrs;
  
  *potc=*potv=0.0;
  iexcl=0;
  for(i=0;i<*n;i++){
    if(i % 1000==0) printf("calc nlist i=%d\n",i);
    for(jj=iexcl;jj<iexcl+numex[i];jj++){
      j=natex[jj]-1;
      //      printf("i=%d j=%d\n",i,j);
      if(j<0){
	//	printf("i=%d natex2[%d]=%d\n",i,jj,j);
	continue;
      }
      atij=(atype[i]-1)*(*nat)+atype[j]-1;
      rs=rscales[atij];
      gsf=gscalesf[atij];
      gsp=gscalesp[atij];
      r=0.0;
      for(k=0;k<3;k++){
	dr[k]=x[i*3+k]-x[j*3+k];
	r+=dr[k]*dr[k];
      }
      rrs=r*rs;
      r=sqrt(r);
      dtmp=gsf*(2.0*pow(rrs,-7.0)-pow(rrs,-4.0));
      for(k=0;k<3;k++){
	force[i*3+k]-=dtmp*dr[k];
	force[j*3+k]+=dtmp*dr[k];
      }
      dtmp=gsp*(pow(rrs,-6.0)-pow(rrs,-3.0));
      *potv-=dtmp*2.0;
      
      dtmp=q[i]*q[j]*pow(r,-3.0);
      for(k=0;k<3;k++){
	force[i*3+k]-=dtmp*dr[k];
	force[j*3+k]+=dtmp*dr[k];
      }
      *potc-=q[i]*q[j]*pow(r,-1.0)*2.0;
    }
    iexcl+=numex[i];
    if(i==0) printf("potc=%e potv=%e f0=%e %e %e\n",*potc,*potv,force[0],force[1],force[2]);
  }
  for(i=0;i<*n;i++){
    if(i % 1000==0) printf("calc i=%d\n",i);
    for(j=0;j<*n;j++){
      r=0.0;
      for(k=0;k<3;k++){
	dr[k]=x[i*3+k]-x[j*3+k];
	r+=dr[k]*dr[k];
      }
      if(r!=0.0){
	atij=(atype[i]-1)*(*nat)+atype[j]-1;
	rs=rscales[atij];
	gsf=gscalesf[atij];
	gsp=gscalesp[atij];
	rrs=r*rs;
	r=sqrt(r);
	dtmp=gsf*(2.0*pow(rrs,-7.0)-pow(rrs,-4.0));
	for(k=0;k<3;k++){
	  force[i*3+k]+=dtmp*dr[k];
	}
	dtmp=gsp*(pow(rrs,-6.0)-pow(rrs,-3.0));
	*potv+=dtmp;
	
	dtmp=q[i]*q[j]*pow(r,-3.0);
	for(k=0;k<3;k++){
	  force[i*3+k]+=dtmp*dr[k];
	}
	*potc+=q[i]*q[j]*pow(r,-1.0);
      }
    }
    if(i==0) printf("potc=%e potv=%e f0=%e %e %e\n",*potc,*potv,force[0],force[1],force[2]);
  }
  printf("potc=%e potv=%e f0=%e %e %e\n",*potc,*potv,force[0],force[1],force[2]);
  //    return;
}


static void compare_forces(double *potc, double potch, double *potv, double potvh,
			   double *force, double *f2, int *n)
{
  double r,sum,fsizeavr,diffmax1,diffmax2,fc[3];
  int i,j;
  
  printf("potc=%e potch=%e diff=%e %e %f\n",*potc,potch,*potc-potch,fabs((*potc-potch)/potch),log(fabs((*potc-potch)/potch))/log(10.0));
  printf("potv=%e potvh=%e diff=%e %e %f\n",*potv,potvh,*potv-potvh,fabs((*potv-potvh)/potvh),log(fabs((*potv-potvh)/potvh))/log(10.0));
  // average size of force vector
  sum=0;
  for(i=0;i<*n;i++){
    r=sqrt(f2[i*3]*f2[i*3]+f2[i*3+1]*f2[i*3+1]+f2[i*3+2]*f2[i*3+2]);
    sum+=r;
  }
  fsizeavr=sum/(*n);
  sum=diffmax1=diffmax2=0;
  for(i=0;i<*n;i++){
    double rg,rh,rdiff,rgh[3],rghsize,rdiff2;
    rg=sqrt(force[i*3]*force[i*3]+force[i*3+1]*force[i*3+1]+force[i*3+2]*force[i*3+2]);
    rh=sqrt(f2[i*3]*f2[i*3]+f2[i*3+1]*f2[i*3+1]+f2[i*3+2]*f2[i*3+2]);
    //      rdiff=fabs((rg-rh)/rh);
    if(i<3) printf("  force[%d]=%e %e %e f2=%e %e %e\n",i,force[i*3],force[i*3+1],force[i*3+2],f2[i*3],f2[i*3+1],f2[i*3+2]);
    for(j=0;j<3;j++) rgh[j]=force[i*3+j]-f2[i*3+j];
    rghsize=sqrt(rgh[0]*rgh[0]+rgh[1]*rgh[1]+rgh[2]*rgh[2]);
    rdiff=fabs(rghsize/rh);
    if(diffmax1<rdiff) diffmax1=rdiff;
    if(rdiff>1e-4){
      printf("i=%d rdiff=%e rg=%e rh=%e rghsize=%e rgh=%e %e %e\n",i,rdiff,rg,rh,rghsize,rgh[0],rgh[1],rgh[2]);
      printf("     force=%e %e %e f2=%e %e %e\n",force[i*3],force[i*3+1],force[i*3+2],f2[i*3],f2[i*3+1],f2[i*3+2]);
    }
    rdiff=fabs(rghsize/fsizeavr);
    if(diffmax2<rdiff) diffmax2=rdiff;
    if(rdiff>1e-4){
      printf("i=%d rdiff_against_average_force=%e rg=%e rh=%e rghsize=%e\n",i,rdiff,rg,rh,rghsize);
      printf("     force=%e %e %e f2=%e %e %e\n",force[i*3],force[i*3+1],force[i*3+2],f2[i*3],f2[i*3+1],f2[i*3+2]);
    }
    sum+=rdiff;
  }
  sum/=(*n);
  printf("avr. force size=%e diffmax1(relative error)=%e diffmax2(against avr. force)=%e\n",
	 fsizeavr,diffmax1,diffmax2);
  printf("     avr. force diff against avr. force=%e (log=%e)\n",sum,log(sum)/log(10.0));
  
  for(j=0;j<3;j++) fc[j]=0.0;
  for(i=0;i<*n;i++) for(j=0;j<3;j++) fc[j]+=force[i*3+j];
  for(j=0;j<3;j++) fc[j]/=(*n);
  for(j=0,sum=0;j<3;j++) sum+=fc[j]*fc[j];
  sum=sqrt(sum);
  printf("force drift grape=%e (%e %e %e)\n",sum,fc[0],fc[1],fc[2]);
  
  for(j=0;j<3;j++) fc[j]=0.0;
  for(i=0;i<*n;i++) for(j=0;j<3;j++) fc[j]+=f2[i*3+j];
  for(j=0;j<3;j++) fc[j]/=(*n);
  for(j=0,sum=0;j<3;j++) sum+=fc[j]*fc[j];
  sum=sqrt(sum);
  printf("force drift host=%e (%e %e %e)\n",sum,fc[0],fc[1],fc[2]);
}


void mr3calccoulomb_vdw_exlist_(int *n, double x[], double q[], 
				int atype[], int *nat, double force[],
				double *rscale, 
				double gscalesf[], double gscalesp[],
				double rscales[],
				int *tblno, int numex[], int natex[],
				double *volume1,
				double *potc, double *potv, 
				int *convert_index, 
				int *potflag, int *flag)
{
  int i,s,j,iexcl,iexcl2,jj,j2,jj2;
  static int *atype2=NULL,*natex2=NULL,*numex2=NULL,*natex3=NULL;
  static double *f2=NULL;
  double potch,potvh;

#ifdef MD_MEASURE_TIME
  vg_start_timer(44);
  vg_start_timer(45);
#endif
  //  printf("mr3calccoulomb_vdw_exlist_ is called\n");

  //  printf("rscales[0]=%e at the beggining of mr3calccoulomb_vdw_exlist_\n",rscales[0]);
  if(atype2==NULL){
    if((atype2=(int *)MR3_malloc_pointer(sizeof(int)*(*n),"atype2 in mr3calccoulomb_vdw_exlist_"))==NULL){
      fprintf(stderr,"** error : can't malloc atype2 **\n");
      MR3_exit(1);
    }
    for(i=0;i<*n;i++){
      if(atype[i]>0) atype2[i]=atype[i]-1;
      else           atype2[i]=0;
    }
    for(i=s=0;i<*n;i++){
      s+=numex[i];
    }
    if(s==0){
#if 1
      fprintf(stderr,"** warning : exlist must exist but skipped **\n");
#else
      fprintf(stderr,"** error : exlist must exist **\n");
      MR3_exit(1);
      return;
#endif
      /*
	MR3calccoulomb_vdw_ci(*n,x,q,atype2,*nat,force,
	*rscale,gscales,rscales,*tblno,
	*rcut,volume,ldim);
	*/
    }
    
    if((natex2=(int *)MR3_malloc_pointer(sizeof(int)*s,"natex2 in mr3calccoulomb_vdw_exlist_"))==NULL){
      fprintf(stderr,
	      "** error at malloc natex2 in mr3calccoulomb_ci_exlist_ **\n");
      MR3_exit(1);
    }
    for(i=0;i<s;i++)  natex2[i]=natex[i]-1;
#if 0 // for host calculation. set 0 for SPE
    //remove minus j-index and sort
    printf("** sorting of natex is performed for host calculation **\n");
    if((natex3=(int *)MR3_malloc_pointer(sizeof(int)*s*2,"numex3 in mr3calccoulomb_vdw_exlist_"))==NULL){
      fprintf(stderr,
	      "** error at malloc natex3 in mr3calccoulomb_ci_exlist_ **\n");
      MR3_exit(1);
    }
    if((numex2=(int *)MR3_malloc_pointer(sizeof(int)*(*n),"numex2 in mr3calccoulomb_vdw_exlist_"))==NULL){
      fprintf(stderr,
	      "** error at malloc numex2 in mr3calccoulomb_ci_exlist_ **\n");
      MR3_exit(1);
    }
    /* count non minus j-index */
    for(i=0;i<*n;i++) numex2[i]=0;
    iexcl=0;
    iexcl2=0;
    for(i=0;i<*n;i++){
      //      printf("natex[%d]=%d ",i,numex[i]);
      for(jj=iexcl;jj<iexcl+numex[i];jj++){
	j=natex[jj]-1;
	//	printf("%d ",j);
	if(j>=0){
	  numex2[i]++;
	  numex2[j]++;
	}
      }
      //      printf("\n");
      iexcl+=numex[i];
    }
    {
      int *natex3_offset;
      if((natex3_offset=(int *)MR3_malloc_pointer(sizeof(int)*(*n),"natex3_offset in mr3calccoulomb_vdw_exlist_"))==NULL){
	fprintf(stderr,
		"** error at malloc natex3 in mr3calccoulomb_ci_exlist_ **\n");
	MR3_exit(1);
      }
      natex3_offset[0]=0;
      for(i=1;i<*n;i++) natex3_offset[i]=natex3_offset[i-1]+numex2[i-1];
      for(i=0;i<*n;i++) numex2[i]=0;
      iexcl=0;
      for(i=0;i<*n;i++){
	for(jj=iexcl;jj<iexcl+numex[i];jj++){
	  j=natex[jj]-1;
	  if(j>=0 && i!=j){
	    natex3[natex3_offset[i]+numex2[i]]=j;
	    numex2[i]++;
	    natex3[natex3_offset[j]+numex2[j]]=i;
	    numex2[j]++;
	  }
	}
	iexcl+=numex[i];
      }
      MR3_free_pointer(natex3_offset,"mr3calccoulomb_vdw_exlist_");
    }

    /* sort j-index */
    iexcl=0;
    for(i=0;i<*n;i++){
      //      printf("natex[%d]=%d ",i,numex2[i]);
      for(jj=iexcl;jj<iexcl+numex2[i];jj++){
	j=natex3[jj];
	for(jj2=jj+1;jj2<iexcl+numex2[i];jj2++){
	  j2=natex3[jj2];
	  if(j>j2){
	    int tmp;
	    tmp=j2;
	    natex3[jj]=j2;
	    natex3[jj2]=j;
	    j=j2;
	  }
	}
      }
      for(jj=iexcl;jj<iexcl+numex2[i];jj++){
	j=natex3[jj];
	//	printf("%d ",j);
      }
      iexcl+=numex2[i];
      //      printf("\n");
    }
#endif    
  }
#if 0 // not work
  {
    int nf_per=0,nf_chg=0,ntbl;
    if((f2=(double *)MR3_malloc_pointer(sizeof(double)*(*n)*3,"f2 in mr3calccoulomb_vdw_exlist_"))==NULL){
      fprintf(stderr,"** error : can't malloc f2 in mr3calccoulomb_vdw_exlist_ **\n");
      MR3_exit(1);
    }
    potch=potvh=0.0;
    
    for(i=0;i<*n*3;i++) f2[i]=0.0;
    ntbl=3;
    mr1calcvdw_host_nf_(x,n,atype,nat,gscalesp,rscales,
			&ntbl,volume1,&nf_per,&nf_chg,
			(double (*)[3])f2);
    for(i=0;i<*n;i++) potvh-=f2[i*3]*2.0;
    
    for(i=0;i<*n*3;i++) f2[i]=0.0;
    ntbl=*tblno+1;
    mr1calccoulomb_host_nf_(x,n,q,rscale,&ntbl,volume1,
			    &nf_per,&nf_chg,(double (*)[3])f2);
    for(i=0;i<*n;i++) potch-=f2[i*3]*2.0;
    printf("potch=%e potvh=%e f2[0]=%e %e %e\n",potch,potvh,f2[0],f2[1],f2[2]);
    
    potch=potvh=0.0;printf("********* potch and potvh are set zero\n");
    
    for(i=0;i<*n*3;i++) f2[i]=0.0;
    MR3calccoulomb_vdw_nlist_host(x,*n,q,*rscale,
				  atype2,*nat,
				  gscalesf,gscalesp,rscales,
				  *tblno,*volume1,nf_per,
				  numex,natex2,
				  &potch,&potvh,
				  f2);
    for(i=0;i<*n*3;i++) f2[i]=-f2[i];
    potch=-potch;
    potvh=-potvh;
    
    ntbl=2;
    mr1calcvdw_host_nf_(x,n,atype,nat,gscalesf,rscales,
			&ntbl,volume1,&nf_per,&nf_chg,
			(double (*)[3])f2);
    mr1calccoulomb_host_nf_(x,n,q,rscale,tblno,volume1,
			    &nf_per,&nf_chg,(double (*)[3])f2);
    
    printf("potch=%e potvh=%e f2[0]=%e %e %e\n",potch,potvh,f2[0],f2[1],f2[2]);
  }
#endif
#ifdef MD_MEASURE_TIME
  vg_stop_and_accumulate_timer(45);
  vg_start_timer(46);
#endif
#if 1 // use PS3 or simple host calculation
  //  printf("flag is changed %d->%d before MR3calccoulomb_vdw_exlist\n",*flag,(*flag & 1));
  {
    int flagnew=0;
    static int ini=0,allowoverlapflag=1;
    char *s;
    //    if((*flag & (1<<6))!=0) flagnew|=1<<12; // overlap is converted to overlap by pthread
    flagnew|=(*flag & 0x1001); // allow only pthread overlap
#ifdef MD_PRINT_WARN
    if(ini==0){
      printf("** warning :  flag is changed %d->%d before MR3calccoulomb_vdw_exlist\n",*flag,flagnew);
      if((s=getenv("VG_NOOVERLAP_CVEXLIST"))!=NULL){
	printf("VG_NOOVERLAP_CVEXLIST is set and overlap is not allowed in MR3calccoulomb_vdw_exlist\n");
	allowoverlapflag=0;
      }
      ini=1;
    }
#endif
    if(allowoverlapflag) flagnew|=(*flag & 0x1041); // allow overlap by GPU
    MR3calccoulomb_vdw_exlist(*n,x,q,atype2,*nat,force,
			      *rscale,gscalesf,gscalesp,rscales,*tblno,
			      numex,natex2,
			      *volume1,
			      potc,potv,NULL,
			      *potflag,flagnew);
  }
  //  MR3_free_pointer(natex2,"natex2 in mr3calccoulomb_vdw_exlist_");
  //  MR3_free_pointer(atype2,"atype2 in mr3calccoulomb_vdw_exlist_");
  //  printf("rscales[0]=%e at the end of mr3calccoulomb_vdw_exlist_\n",rscales[0]);
#else
  printf("*** warning : non PS3 calculation **\n");
  printf("*** you need to change if 0 statement upper from here **\n");
  printf("*** which is labeled as 'for host calculation. set 0 for SPE' **\n");
#if 1
  printf("MR3calccoulomb_vdw_exlist_host2 is called\n");
  MR3calccoulomb_vdw_exlist_host2(n,x,q, 
				  atype,nat,force,
				  rscale,gscalesf,gscalesp,rscales,
				  tblno,numex,natex,
				  volume1,potc,potv, 
				  convert_index,potflag,flag);
#elif 0
  printf("MR3calccoulomb_vdw_exlist_host3 is called\n");
  MR3calccoulomb_vdw_exlist_host3(*n,x,q, 
				  atype,*nat,force,
				  *rscale,gscalesf,gscalesp,rscales,
				  *tblno,numex,natex,
				  *volume1,potc,potv, 
				  convert_index,*potflag,*flag);
#elif 0
  printf("MR3calccoulomb_vdw_exlist_host5 is called\n");
  MR3calccoulomb_vdw_exlist_host5(*n,x,q, 
				  atype2,*nat,force,
				  *rscale,gscalesf,gscalesp,rscales,
				  *tblno,numex2,natex3,
				  *volume1,potc,potv, 
				  convert_index,*potflag,*flag);
#else
  printf("MR3calccoulomb_vdw_exlist_host6 is called\n");
  MR3calccoulomb_vdw_exlist_host6(*n,x,q, 
				  atype2,*nat,force,
				  *rscale,gscalesf,gscalesp,rscales,
				  *tblno,numex2,natex3,
				  *volume1,potc,potv, 
				  convert_index,*potflag,*flag);
#endif
#endif
#ifdef MD_MEASURE_TIME
  vg_stop_and_accumulate_timer(46);
  vg_start_timer(47);
#endif
  
#if 0 // for comparing with host calculation
      // you need to change if 0 of other 'comparing with host calculation' line
      // in MR3_get_forces_and_potentials_overlap_cell()
  if(f2==NULL){
    if((f2=(double *)MR3_malloc_pointer(sizeof(double)*(*n)*3,"f2 in mr3calccoulomb_vdw_exlist_"))==NULL){
      fprintf(stderr,"** error : can't malloc f2 in mr3calccoulomb_vdw_exlist_ **\n");
      MR3_exit(1);
    }
  }
  F2compare=f2;
  Ncompare=*n;
  potch=potvh=0.0;
  for(i=0;i<*n;i++) f2[i*3]=f2[i*3+1]=f2[i*3+2]=0.0;
  MR3calccoulomb_vdw_ij_exlist_host(*n,x,q,atype2,f2, 
				    *n,x,q,atype2,
				    *nat,gscalesf,gscalesp,rscales,*rscale,*tblno,
				    *volume1,&potch,&potvh,0,1,0,
				    numex,natex2);
  Potccompare=potch;
  Potvcompare=potvh;

#if 1 // call subroutine for comparing force
  printf("***** calling compare_forces\n");
  compare_forces(potc,potch,potv,potvh,force,f2,n);
#elif 0 // else of call subroutine for comparing force 
  {
    double r,sum,fsizeavr,diffmax1,diffmax2,fc[3];

    printf("potc=%e potch=%e diff=%e %e %f\n",*potc,potch,*potc-potch,fabs((*potc-potch)/potch),log(fabs((*potc-potch)/potch))/log(10.0));
    printf("potv=%e potvh=%e diff=%e %e %f\n",*potv,potvh,*potv-potvh,fabs((*potv-potvh)/potvh),log(fabs((*potv-potvh)/potvh))/log(10.0));
    // average size of force vector
    sum=0;
    for(i=0;i<*n;i++){
      r=sqrt(f2[i*3]*f2[i*3]+f2[i*3+1]*f2[i*3+1]+f2[i*3+2]*f2[i*3+2]);
      sum+=r;
    }
    fsizeavr=sum/(*n);
    sum=diffmax1=diffmax2=0;
    for(i=0;i<*n;i++){
      double rg,rh,rdiff,rgh[3],rghsize,rdiff2;
      rg=sqrt(force[i*3]*force[i*3]+force[i*3+1]*force[i*3+1]+force[i*3+2]*force[i*3+2]);
      rh=sqrt(f2[i*3]*f2[i*3]+f2[i*3+1]*f2[i*3+1]+f2[i*3+2]*f2[i*3+2]);
      //      rdiff=fabs((rg-rh)/rh);
      for(j=0;j<3;j++) rgh[j]=force[i*3+j]-f2[i*3+j];
      rghsize=sqrt(rgh[0]*rgh[0]+rgh[1]*rgh[1]+rgh[2]*rgh[2]);
      rdiff=fabs(rghsize/rh);
      if(diffmax1<rdiff) diffmax1=rdiff;
      if(rdiff>1e-4){
	printf("i=%d rdiff=%e rg=%e rh=%e rghsize=%e rgh=%e %e %e\n",i,rdiff,rg,rh,rghsize,rgh[0],rgh[1],rgh[2]);
	printf("     force=%e %e %e f2=%e %e %e\n",force[i*3],force[i*3+1],force[i*3+2],f2[i*3],f2[i*3+1],f2[i*3+2]);
      }
      rdiff=fabs(rghsize/fsizeavr);
      if(diffmax2<rdiff) diffmax2=rdiff;
      if(rdiff>1e-4){
	printf("i=%d rdiff_against_average_force=%e rg=%e rh=%e rghsize=%e\n",i,rdiff,rg,rh,rghsize);
	printf("     force=%e %e %e f2=%e %e %e\n",force[i*3],force[i*3+1],force[i*3+2],f2[i*3],f2[i*3+1],f2[i*3+2]);
      }
      sum+=rdiff;
    }
    sum/=(*n);
    printf("avr. force size=%e diffmax1(relative error)=%e diffmax2(against avr. force)=%e\n",
	   fsizeavr,diffmax1,diffmax2);
    printf("     avr. force diff against avr. force=%e (log=%e)\n",sum,log(sum)/log(10.0));
    
    for(j=0;j<3;j++) fc[j]=0.0;
    for(i=0;i<*n;i++) for(j=0;j<3;j++) fc[j]+=force[i*3+j];
    for(j=0;j<3;j++) fc[j]/=(*n);
    for(j=0,sum=0;j<3;j++) sum+=fc[j]*fc[j];
    sum=sqrt(sum);
    printf("force drift grape=%e (%e %e %e)\n",sum,fc[0],fc[1],fc[2]);
    
    for(j=0;j<3;j++) fc[j]=0.0;
    for(i=0;i<*n;i++) for(j=0;j<3;j++) fc[j]+=f2[i*3+j];
    for(j=0;j<3;j++) fc[j]/=(*n);
    for(j=0,sum=0;j<3;j++) sum+=fc[j]*fc[j];
    sum=sqrt(sum);
    printf("force drift host=%e (%e %e %e)\n",sum,fc[0],fc[1],fc[2]);
    
  }
#endif // end of call subroutine for comparing force 
#endif // end of comparing with host calculation
  //  printf("potc=%e potv=%e force[0]=%e %e %e\n",*potc,*potv,force[0],force[1],force[2]);
#ifdef MD_MEASURE_TIME
  vg_stop_and_accumulate_timer(46);
  vg_stop_and_accumulate_timer(44);
#endif
}


void MR3_get_forces_overlap(double force[])
{
  VG_UNIT *unit;
  void *thread_ret;
  int idum=0,i,j,ret;

  unit=m3_get_unit();
  if(unit!=NULL){
    if(unit->gpuoverlapflag==2 && unit->fthread!=NULL){ // GPU overlap
      //      printf("in MR3_get_forces_overlap: function_index=%d **\n",unit->function_index);
      unit->calculate[unit->function_index]((void *)unit);
      //      printf("gpuoverlapflag!=0, copying fvec->fthread\n");
      vg_convert_forces(unit,unit->nf,unit->function_index,(double (*)[3])(unit->fthread));
      if(unit->function_index<40){
        MR3calccoulomb_vdw_ij_ci(0,NULL,NULL,NULL,force,0,NULL,NULL,NULL,0,NULL,NULL,0.0,0,0.0,0.0,NULL,NULL);
      }
      else if(unit->function_index<50){
        MR3calccoulomb_vdw_ij_ci_virial(0,NULL,NULL,NULL,force,0,NULL,NULL,NULL,0,NULL,NULL,0.0,0,0.0,0.0,NULL,NULL,NULL,0);
      }
      //      for(i=0;i<9;i++) printf("get_forces_overlap: f[%d]=%e %e %e\n",i,force[i*3],force[i*3+1],force[i*3+2]);
      unit->gpuoverlapflag=unit->function_index=0;
      unit->fthread=NULL;
    }
    if(unit->fthread!=NULL){                            // pthread overlap
      if(unit->thread!=0){
	ret=pthread_join(unit->thread,&thread_ret);
	if(ret){
	  fprintf(stderr,"** can't join therad in MR3_get_forces_and_virial_overlap **\n");
	  MR3_exit(1);
	}
      }
      //    printf("in MR3_get_forces_and_virial_overlap, forces are copied ni_overlap=%d\n",unit->ni_overlap);fflush(stdout);
      for(i=0;i<(unit->ni_overlap)*3;i++) force[i]+=unit->fthread[i];
      MR3_free_pointer(unit->fthread,"unit->fthread in MR3_get_forces_and_virial_overlap");
      unit->fthread=NULL;
      unit->thread=0;
      unit->ni_overlap=0;
      unit->gpuoverlapflag=unit->function_index=0;
#if 0
      for(i=0;i<3;i++){
	printf("get_forces_and_virial_overlap: force[%d]=%e %e %e\n",
	       i,force[i*3],force[i*3+1],force[i*3+2]);
      }
#endif
    }
#if 0
    else{
      static int ini=0;
      if(ini==0){
	fprintf(stderr,"** warning : MR3_get_forces_overlap supports only pthread overlap **\n");
	ini=1;
	//    vg_exit(1);
      }
    }
#endif
  }
}


void mr3_get_forces_overlap_(double force[])
{
  MR3_get_forces_overlap(force);
}


void MR3_get_potc_and_potv_overlap(double *potc, double *potv)
{
  VG_UNIT *unit;
  static int ini=0;

  if(ini==0){
    printf("MR3_get_potc_and_potv_overlap is called\n");
    printf("This routine should not be called when pthread is not used.\n");
    printf("You need to modify ew_force.f to skip this.\n");
    ini=1;
  }
  
  unit=m3_get_unit();
  *potc=unit->potc;
  *potv=unit->potv;
  unit->potc=unit->potv=0.0;
}


void mr3_get_potc_and_potv_overlap_(double *potc, double *potv)
{
  MR3_get_potc_and_potv_overlap(potc,potv);
}


void mr3_get_potc_and_potv_overlap__(double *potc, double *potv)
{
  mr3_get_potc_and_potv_overlap_(potc,potv);
}


void MR3_get_forces_and_virial_overlap(double force[], double virial[3][3])
{
  fprintf(stderr,"** error : MR3_get_forces_and_virial_overlap is not supported now **\n");
  vg_exit(1);
}


void mr3_get_forces_and_virial_overlap_(double force[], double virial[3][3])
{
  MR3_get_forces_and_virial_overlap(force,virial);
}


void MR3calccoulomb_ci_exlist(int n, double x[], double q[], double force[],
			      double rscale, int tblno,
			      int numex[], int natex[],
			      double rcut, double skinnb,
			      double volume[3], int ldim[3], 
			      int *convert_index,
			      int flag)
{
}


void mr3calccoulomb_ci_exlist_(int *n, double x[], double q[], double force[],
			       double *rscale, int *tblno,
			       int numex[], int natex[],
			       double *rcut, double *skinnb,
			       double volume[3], int ldim[3], 
			       int *convert_index,
			       int *flag)
{
  int i,*natex2,s,*convert_index2=NULL;

  for(i=s=0;i<*n;i++){
    s+=numex[i];
  }
  if(s==0){
    MR3calccoulomb_ci(*n,x,q,force,*rscale,*tblno,
		      *rcut,*skinnb,volume,ldim);
    return;
  }

  if((natex2=(int *)MR3_malloc_pointer(sizeof(int)*s,"natex2 in mr3calccoulomb_ci_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc natex2 in mr3calccoulomb_ci_exlist_ **\n");
    MR3_exit(1);
  }
  if(convert_index!=NULL){
    if((convert_index2=(int *)MR3_malloc_pointer(sizeof(int)*(*n),"convert_index2 in mr3calccoulomb_ci_exlist_"))==NULL){
      fprintf(stderr,
	      "** error at malloc convert_index2 in mr3calccoulomb_ci_exlist **\n");
      MR3_exit(1);
    }
    for(i=0;i<*n;i++) convert_index2[i]=convert_index[i]-1;
  }
  for(i=0;i<s;i++)  natex2[i]=natex[i]-1;
  MR3calccoulomb_ci_exlist(*n,x,q,force,*rscale,*tblno,
			   numex,natex2,*rcut,*skinnb,volume,ldim, 
			   convert_index2,*flag);
  MR3_free_pointer(natex2,"natex2 in mr3calccoulomb_ci_exlist_");
  if(convert_index!=NULL) MR3_free_pointer(convert_index2,"convert_index2 in mr3calccoulomb_ci_exlist_");
}


void MR3calcvdw_ci(int n, double x[], int atype[], int nat, double force[],
		   double gscales[], double rscales[], int tblno,
		   double rcut, double skinnb, double volume[3], int ldim[3])
{
}


void mr3calcvdw_ci_(int *n, double x[], int atype[], int *nat, double force[],
		    double gscales[], double rscales[], int *tblno,
		    double *rcut, double *skinnb, 
		    double volume[3], int ldim[3])
{
  int i,*atype2;

  if((atype2=(int *)MR3_malloc_pointer(sizeof(int)*(*n),"atype2 in mr3calcvdw_ci_"))==NULL){
    fprintf(stderr,"** error : can't malloc atype2 **\n");
    MR3_exit(1);
  }
  for(i=0;i<*n;i++){
    if(atype[i]>0) atype2[i]=atype[i]-1;
    else           atype2[i]=0;
  }
  MR3calcvdw_ci(*n,x,atype2,*nat,force,gscales,rscales,*tblno,
		*rcut,*skinnb,volume,ldim);
  MR3_free_pointer(atype2,"atype2 in mr3calcvdw_ci_");
}


void MR3calcvdw_ci_exlist(int n, double x[], int atype[], int nat, 
			  double force[],
			  double gscales[], double rscales[], int tblno,
			  int numex[], int natex[],
			  double rcut, double skinnb,
			  double volume[3], int ldim[3],
			  int *convert_index,
			  int flag)
{
}


void mr3calcvdw_ci_exlist_(int *n, double x[], int atype[], int *nat, 
			   double force[],
			   double gscales[], double rscales[], int *tblno,
			   int numex[], int natex[],
			   double *rcut, double *skinnb,
			   double volume[3], int ldim[3],
			   int *convert_index,
			   int *flag)
{
  int i,*atype2;
  int *natex2,s,*convert_index2=NULL;

#if 0
  {
    int flag1=1,flag2=0;
#if 0
    mr3calcvdw_ij_(n,x,atype,force,n,x,atype,
		  nat,gscales,rscales,tblno,volume,&flag1);
#endif
#if 0
    mr3calcvdw_ij_exlist_(n,x,atype,force,n,x,atype,
			  nat,gscales,rscales,tblno,volume,&flag1,
			  numex,natex);
#endif
#if 1
    mr3calcvdw_ci_(n,x,atype,nat,force,
		   gscales,rscales,tblno,rcut,volume,ldim);
#endif
    printf("returning in mr3calcvdw_ci_exlist\n");
    return;
  }
#endif  
  if((atype2=(int *)MR3_malloc_pointer(sizeof(int)*(*n),"atype2 in mr3calcvdw_ci_exlist_"))==NULL){
    fprintf(stderr,"** error : can't malloc atype2 **\n");
    MR3_exit(1);
  }
  for(i=0;i<*n;i++){
    if(atype[i]>0) atype2[i]=atype[i]-1;
    else           atype2[i]=0;
  }

  for(i=s=0;i<*n;i++){
    s+=numex[i];
  }
  if(s==0){
    MR3calcvdw_ci(*n,x,atype2,*nat,force,gscales,rscales,*tblno,
		  *rcut,*skinnb,volume,ldim);
    return;
  }

  if((natex2=(int *)MR3_malloc_pointer(sizeof(int)*s,"natex2 in mr3calcvdw_ci_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc natex2 in mr3calcvdw_ci_exlist_ **\n");
    MR3_exit(1);
  }
  if(convert_index!=NULL){
    if((convert_index2=(int *)MR3_malloc_pointer(sizeof(int)*(*n),"convert_index2 in mr3calcvdw_ci_exlist_"))==NULL){
      fprintf(stderr,
	      "** error at malloc convert_index2 in mr3calcvdw_ci_exlist **\n");
      MR3_exit(1);
    }
    for(i=0;i<*n;i++) convert_index2[i]=convert_index[i]-1;
  }
  for(i=0;i<s;i++)  natex2[i]=natex[i]-1;

  MR3calcvdw_ci_exlist(*n,x,atype2,*nat,force,gscales,rscales,*tblno,
		       numex,natex2,*rcut,*skinnb,volume,ldim,convert_index2,
		       *flag);
  MR3_free_pointer(atype2,"atype2 in mr3calcvdw_ci_exlist_");
  MR3_free_pointer(natex2,"natex2 in mr3calcvdw_ci_exlist_");
  if(convert_index!=NULL) MR3_free_pointer(convert_index2,"convert_index2 in mr3calcvdw_ci_exlist_");
}


void MR3calccoulomb_vdw_ci_exlist(int n, double x[], double q[], 
				  int atype[], int nat, double force[],
				  double rscale, 
				  double gscalesf[], double gscalesp[],
				  double rscales[],
				  int tblno, int numex[], int natex[],
				  double rcut, double skinnb,
				  double volume[3], int ldim[3],
				  double *potc, double *potv, 
				  int *convert_index, 
				  int potflag, int flag)
{
  double *ftmp,xmax;
  int periodicflag,i;
  VG_UNIT *unit=m3_get_unit();

  if(unit==NULL){
    fprintf(stderr,"** error : unit is NULL **\n");
    vg_exit(1);
  }
  vg_set_rcut2(unit,rcut*rcut);
  if(volume[0]!=volume[1] || volume[0]!=volume[2]){
    fprintf(stderr,"** error : only cubic cell is supported **\n");
    vg_exit(1);
  }
  xmax=volume[0];
  periodicflag=flag & 3;
  if((ftmp=(double *)MR3_malloc_pointer(sizeof(double)*n*3,"MR3calccoulomb_vdw_ci_exlist"))==NULL){
    fprintf(stderr,"** error : can't maloc ftmp **\n");
    vg_exit(1);
  }
  if((potflag & 1)==1){
    bzero(ftmp,sizeof(double)*n*3);
    *potc=0.0;
    MR3calccoulomb_ij_exlist(n,x,q,ftmp,
			     n,x,q,rscale,
			     tblno+1,xmax,periodicflag,
			     numex,natex);
    for(i=0;i<n;i++) (*potc)+=ftmp[i*3];
    bzero(ftmp,sizeof(double)*n*3);
    *potv=0.0;
    MR3calcvdw_ij_exlist(n,x,atype,ftmp,
			 n,x,atype,
			 nat,gscalesp,rscales,
			 3,xmax,periodicflag,numex,natex);
    for(i=0;i<n;i++) (*potv)+=ftmp[i*3];
  }
  bzero(ftmp,sizeof(double)*n*3);
  MR3calccoulomb_ij_exlist(n,x,q,ftmp,
			   n,x,q,rscale,
			   tblno,xmax,periodicflag,
			   numex,natex);
  for(i=0;i<n*3;i++) force[i]+=ftmp[i];
  bzero(ftmp,sizeof(double)*n*3);
  MR3calcvdw_ij_exlist(n,x,atype,ftmp,
		       n,x,atype,
		       nat,gscalesf,rscales,
		       2,xmax,periodicflag,numex,natex);
  for(i=0;i<n*3;i++) force[i]+=ftmp[i];
  MR3_free_pointer(ftmp,"MR3calccoulomb_vdw_ci_exlist");
}


void mr3calccoulomb_vdw_ci_exlist_(int *n, double x[], double q[], 
				   int atype[], int *nat, double force[],
				   double *rscale, 
				   double gscalesf[], double gscalesp[],
				   double rscales[],
				   int *tblno, int *numex, int natex[],
				   double *rcut, double *skinnb,
				   double volume[3], int ldim[3],
				   double *potc, double *potv, 
				   int *convert_index, 
				   int *potflag, int *flag)
{
  int i,*natex2,s,*convert_index2=NULL;
  int *atype2;

  if((atype2=(int *)MR3_malloc_pointer(sizeof(int)*(*n),"atype2 in mr3calccoulomb_vdw_ci_exlist_"))==NULL){
    fprintf(stderr,"** error : can't malloc atype2 **\n");
    MR3_exit(1);
  }
  for(i=0;i<*n;i++){
    if(atype[i]>0) atype2[i]=atype[i]-1;
    else           atype2[i]=0;
  }
  for(i=s=0;i<*n;i++){
    s+=numex[i];
  }
#if 1
  if(s==0){
    fprintf(stderr,"** error : exlist must exist **\n");
    MR3_exit(1);
    /*
    MR3calccoulomb_vdw_ci(*n,x,q,atype2,*nat,force,
			  *rscale,gscales,rscales,*tblno,
			  *rcut,volume,ldim);
			  */
    return;
  }
  if((natex2=(int *)MR3_malloc_pointer(sizeof(int)*s,"natex2 in mr3calccoulomb_vdw_ci_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc natex2 in mr3calccoulomb_ci_exlist_ **\n");
    MR3_exit(1);
  }
#else
  if(s==0){
    fprintf(stderr,"** warning : exlist is zero **\n");
    natex2=NULL;
    numex=NULL;
  }
  else{
    if((natex2=(int *)MR3_malloc_pointer(sizeof(int)*s,"natex2 in mr3calccoulomb_vdw_ci_exlist_"))==NULL){
      fprintf(stderr,
	      "** error at malloc natex2 in mr3calccoulomb_ci_exlist_ **\n");
      MR3_exit(1);
    }
  }
#endif

  if(convert_index!=NULL){
    if((convert_index2=(int *)MR3_malloc_pointer(sizeof(int)*(*n),"convert_index2 in mr3calccoulomb_vdw_ci_exlist_"))==NULL){
      fprintf(stderr,
	      "** error at malloc convert_index2 in mr3calccoulomb_ci_exlist **\n");
      MR3_exit(1);
    }
    for(i=0;i<*n;i++) convert_index2[i]=convert_index[i]-1;
  }
  for(i=0;i<s;i++)  natex2[i]=natex[i]-1;
  MR3calccoulomb_vdw_ci_exlist(*n,x,q,atype2,*nat,force,
			       *rscale,gscalesf,gscalesp,rscales,*tblno,
			       numex,natex2,
			       *rcut,*skinnb,volume,ldim,
			       potc,potv,convert_index2,
			       *potflag,*flag);
  MR3_free_pointer(natex2,"natex2 in mr3calccoulomb_vdw_ci_exlist_");
  if(convert_index!=NULL) MR3_free_pointer(convert_index2,"convert_index2 in mr3calccoulomb_vdw_ci_exlist_");
  MR3_free_pointer(atype2,"atype2 in mr3calccoulomb_vdw_ci_exlist_");
}
