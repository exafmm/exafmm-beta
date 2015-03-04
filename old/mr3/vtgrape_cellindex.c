static void check_rcut_and_skinnb(double rcut, double skinnb,
				  int ldim[3], double volume[3])
{
  if((rcut+skinnb)*ldim[0]>volume[0] || 
     (rcut+skinnb)*ldim[1]>volume[1] ||
     (rcut+skinnb)*ldim[2]>volume[2]){
    fprintf(stderr,"** error : rcut(%f)+skinnb(%f) must be smaller than a cell size (%f,%f,%f), ldim=(%d,%d,%d) **\n",
	    rcut,skinnb,volume[0]/ldim[0],volume[1]/ldim[1],volume[2]/ldim[2],
	    ldim[0],ldim[1],ldim[2]);
    MR3_exit(1);
  }
}


static void assign_cell(int n, double x_org[], 
			double xmin[3], double xmax[3], int ldimtotal[3],
			int blocksize,
			int (*cip)[4], double (*op)[3], M3_CELL *ci,
			int *cinthp)
{
  /* 
    input
     blocksize ---  1 means no dummy particle
                   >1 means some dummy particles in each cell
    output
     cip[i][0,1,2] -- index of a cell (cx,cy,cz) in which a paricle-i is.
     cip[i][3] ------ cellindex in which a particle-i is.
     op[i][0,1,2] --- if a particle-i is outsize of a simulation cell
                      (xmin-xmax), the offset of the coordinate is stored.
     ci[ci].base ---- index of the first particle in a cell-ci.
                      this becomes a multiple of a blocksize.
     ci[ci].size ---- number of particles in a cell-ci.
                      this does not become a multiple of a blocksize.
     cinthp[i] ------ offset of particles in the cell which a particle-i is.
                      particle-i is the cinthp[i]-th particle in a cell.
  */

  int i,j;
  double lv[3],volume[3];
  int cid[3],cit,ncell;
  int ipdb_last=-99999;
  double *x;

  ncell=ldimtotal[0]*ldimtotal[1]*ldimtotal[2];
  for(i=0;i<3;i++){
    volume[i]=xmax[i]-xmin[i];
    lv[i]=ldimtotal[i]/volume[i];
  }

  /* convert particle position to [0,volume) coordinate */
  if((x=(double *)MR3_malloc_pointer(sizeof(double)*n*3,"x in assign_cell"))==NULL){
    fprintf(stderr,"** error : can't malloc x in assign_cell **\n");
    vg_exit(1);
  }
  for(i=0;i<n;i++){
    for(j=0;j<3;j++){
      x[i*3+j]=x_org[i*3+j];
      while(x[i*3+j]<0.0)        x[i*3+j]+=volume[j];
      while(x[i*3+j]>=volume[j]) x[i*3+j]-=volume[j];
    }
  }

  /* assign a particle into a cell */
  for(i=0;i<ncell;i++) ci[i].size=0;
  for(i=0;i<n;i++){
    for(j=0;j<3;j++){
      cid[j]=(x[i*3+j]-xmin[j])*lv[j];
      op[i][j]=0.0;
#if 0
      if(cid[j]<0 || cid[j]>=ldimtotal[j]){
	printf("assign_cell: i=%d j=%d x=%e cid=%d\n",i,j,x[i*3+j],cid[j]);
      }
#endif
      while(cid[j]<0){
	cid[j]+=ldimtotal[j];
	op[i][j]+=volume[j];
      }
      while(cid[j]>=ldimtotal[j]){
	cid[j]-=ldimtotal[j];
	op[i][j]-=volume[j];
      }
    }
    cit=(cid[2]*ldimtotal[1]+cid[1])*ldimtotal[0]+cid[0];
    cip[i][3]=cit;
    for(j=0;j<3;j++) cip[i][j]=cid[j];
    cinthp[i]=ci[cit].size;
    ci[cit].size++;
  }

  /* make cell index list */
  ci[0].base=0;
  for(i=1;i<ncell;i++){
    //    ci[i].base=ci[i-1].base+ci[i-1].size;
    ci[i].base=ci[i-1].base+(ci[i-1].size+blocksize-1)/blocksize*blocksize;
  }

  /* free */
  MR3_free_pointer(x,"x in assign_cell");

#if 0
  for(i=0;i<ncell;i++){
    printf("ci[%d].base=%d size=%d\n",i,ci[i].base,ci[i].size);
  }
#endif  
}


static void cell_range(int n, double x[], 
		       double xmin[3], double xmax[3], int ldimtotal[3],
		       int *ncell, int cmin[3], int cmax[3],
		       int periodicflag)
{
  /*
    Returns the range of cells: cmin[3] and cmax[3]. Cellindex should
    be [ cmin[0], cmax[0] ], [ cmin[1], cmax[1] ], [ cmin[2], cmax[2] ].
    It means that the particle in the simulation cell has the cellindex
    in the range of [0, ldimtotal[x,y,z]-1 ].

    periodicflag -- 0 : non periodic
                    1 : periodic
                        cmax[0,1,2]-cmin[0,1,2]+1 should be the same
                        as ldimtotal[0,1,2]
   */
  int i,j,cid[3];
  double volume[3],lv[3];

  for(i=0;i<3;i++){
    cmin[i]=99999999;
    cmax[i]=-99999999;
  }
  for(i=0;i<3;i++){
    volume[i]=xmax[i]-xmin[i];
    lv[i]=ldimtotal[i]/volume[i];
  }
  for(i=0;i<n;i++){
    for(j=0;j<3;j++){
      cid[j]=floor((x[i*3+j]-xmin[j])*lv[j]);
      //      printf("i=%d j=%d cid=%d x=%e xmin=%e ldimtotal=%d volume=%e\n",i,j,cid[j],x[i*3+j],xmin[j],ldimtotal[j],volume[j]);
#if 0
      if(cid[j]<-100 || cid[j]>100){
	fprintf(stderr,"** error : i=%d j=%d x=%e xmin=%e lv=%e cid=%d is out of range **\n",i,j,x[i*3+j],xmin[j],lv[j],cid[j]);
	MR3_exit(1);
      }
#endif
      if(cid[j]<cmin[j]) cmin[j]=cid[j];
      if(cid[j]>cmax[j]) cmax[j]=cid[j];
    }
  }
  if(periodicflag){
    //    printf("in cell_range : cmax is modified because of periodicflag\n");
    for(j=0;j<3;j++) cmax[j]=cmin[j]+ldimtotal[j]-1;
  }
  *ncell=1;
  for(i=0;i<3;i++){
    *ncell*=(cmax[i]-cmin[i]+1);
  }
#if 1
  {
    static int ini=0;
    if(ini<2){
#ifdef MD_PRINT_WARN
      printf("in cell_range : ncell=%d cmin=(%d,%d,%d) cmax=(%d,%d,%d)\n",
	     *ncell,cmin[0],cmin[1],cmin[2],cmax[0],cmax[1],cmax[2]);
#endif
      ini++;
    }
  }
#endif
#if 1
  for(i=0;i<3;i++){
    if(cmin[i]<-100 || cmax[i]>100){
      fprintf(stderr,"** error : cmin=%d cmax=%d is out of range **\n",cmin[i],cmax[i]);
      MR3_exit(1);
    }
  }
  {
    int ncell=1;
    for(i=0;i<3;i++) ncell*=(cmax[i]-cmin[i]+1);
    if(ncell>MD_CELLINDEX_MAXIBLOCKS){
      fprintf(stderr,"** error : too many cells (%d) in cell_range **\n",ncell);
      MR3_exit(1);
    }
  }
#endif
}


static void assign_cell_noperiodic(int n, double *x,
				   double xmin[3], double xmax[3], int ldimtotal[3],
				   int blocksize, 
				   int ncell, int cmin[3], int cmax[3],
				   int (*cip)[4], M3_CELL *ci,
				   int *cinthp,
				   int periodicflag)
{
  /* 
    input
     blocksize ---  1 means no dummy particle
                   >1 means some dummy particles in each cell
     periodicflag - 0 non-periodic
                    1 periodic
                      cellindex is forced to be in the range of ldimtotal
    output
     ncell, cmin, cmax
                   -- output from cell_range()
     cip[i][0,1,2] -- index of a cell (cx,cy,cz) in which a paricle-i is.
                      actually, it is shifted by cmin[x,y,z]. Therefore
                      cip[i][0,1,2] starts from (0,0,0) and max value is
                      (cmax[0]-cmmin[0],cmax[1]-cmin[1],cmax[2]-cmin[2]).
     cip[i][3] ------ cellindex in which a particle-i is.
     ci[ci].base ---- index of the first particle in a cell-ci.
                      this becomes a multiple of a blocksize.
     ci[ci].size ---- number of particles in a cell-ci.
                      this does not become a multiple of a blocksize.
     cinthp[i] ------ offset of particles in the cell which a particle-i is.
                      particle-i is the cinthp[i]-th particle in a cell.
  */

  int i,j,cminmax[3];
  double lv[3],volume[3];
  int ipdb_last=-99999;
  int cid[3],cit;

  for(i=0;i<3;i++){
    volume[i]=xmax[i]-xmin[i];
    lv[i]=ldimtotal[i]/volume[i];
    cminmax[i]=cmax[i]-cmin[i]+1;
  }

  /* assign a particle into a cell */
  for(i=0;i<ncell;i++) ci[i].size=0;
  for(i=0;i<n;i++){
    for(j=0;j<3;j++){
      cid[j]=floor((x[i*3+j]-xmin[j])*lv[j]);
      cid[j]-=cmin[j];
      if(periodicflag){
	while(cid[j]>=cminmax[j]){
	  //	  printf("  i=%d j=%d x=%e cid=%d->",i,j,x[i*3+j],cid[j]);
	  cid[j]-=cminmax[j];
	  //	  printf("%d\n",cid[j]);
	}
      }
    }
    cit=(cid[2]*cminmax[1]+cid[1])*cminmax[0]+cid[0];
    cip[i][3]=cit;
    for(j=0;j<3;j++) cip[i][j]=cid[j];
    cinthp[i]=ci[cit].size;
    ci[cit].size++;
  }

  /* make cell index list */
  ci[0].base=0;
  for(i=1;i<ncell;i++){
    //    ci[i].base=ci[i-1].base+ci[i-1].size;
    ci[i].base=ci[i-1].base+(ci[i-1].size+blocksize-1)/blocksize*blocksize;
  }
#if 0
  for(i=0;i<ncell;i++){
    printf("in assign_cell_noperiodic : ci[%d].base=%d size=%d\n",i,ci[i].base,ci[i].size);
  }
#endif  
}


static void make_contiguous_xq(int n, double *x, double *q,
			       double xmin[3], double xmax[3],
			       int (*cip)[4], M3_CELL *ci, int *cinthp,
			       double *xcell, double *q2)
{
  int i,cellindex,xcellpointer,xyz;
  double tmp,volume[3];

  for(xyz=0;xyz<3;xyz++) volume[xyz]=xmax[xyz]-xmin[xyz];
  for(i=0;i<n;i++){
    cellindex=cip[i][3];
    xcellpointer=ci[cellindex].base+cinthp[i];
    for(xyz=0;xyz<3;xyz++){
      tmp=x[i*3+xyz]-xmin[xyz];
      while(tmp<0.0)          tmp+=volume[xyz];
      while(tmp>=volume[xyz]) tmp-=volume[xyz];
      xcell[xcellpointer*3+xyz]=tmp;
    }
    q2[xcellpointer]=q[i];
  }
}


static void make_contiguous_xa(int n, double *x, int *atype,
			       double xmin[3], double xmax[3],
			       int (*cip)[4], M3_CELL *ci, int *cinthp,
			       double *xcell, int *atype2)
{
  int i,cellindex,xcellpointer,xyz;
  double tmp,volume[3];

  for(xyz=0;xyz<3;xyz++) volume[xyz]=xmax[xyz]-xmin[xyz];
  for(i=0;i<n;i++){
    cellindex=cip[i][3];
    xcellpointer=ci[cellindex].base+cinthp[i];
    for(xyz=0;xyz<3;xyz++){
      tmp=x[i*3+xyz]-xmin[xyz];
      while(tmp<0.0)          tmp+=volume[xyz];
      while(tmp>=volume[xyz]) tmp-=volume[xyz];
      xcell[xcellpointer*3+xyz]=x[i*3+xyz];
    }
    atype2[xcellpointer]=atype[i]+1;
  }
}


static void make_contiguous_xqa(int n, double *x, double *q, int *atype,
				double xmin[3], double xmax[3],
				int (*cip)[4], M3_CELL *ci, int *cinthp,
				double *xcell, double *q2, int *atype2)
{
  int i,cellindex,xcellpointer,xyz;
  double tmp,volume[3];

  for(xyz=0;xyz<3;xyz++) volume[xyz]=xmax[xyz]-xmin[xyz];
  for(i=0;i<n;i++){
    cellindex=cip[i][3];
    xcellpointer=ci[cellindex].base+cinthp[i];
    for(xyz=0;xyz<3;xyz++){
      tmp=x[i*3+xyz]-xmin[xyz];
      while(tmp<0.0)          tmp+=volume[xyz];
      while(tmp>=volume[xyz]) tmp-=volume[xyz];
      xcell[xcellpointer*3+xyz]=x[i*3+xyz];
    }
    q2[xcellpointer]=q[i];
    atype2[xcellpointer]=atype[i]+1;
  }
}


void MR3calccoulomb_ci(int n, double x[], double q[], double force[],
		       double rscale, int tblno,
		       double rcut, double skinnb,
		       double volume[3], int ldim[3])
{
  double xmin[3]={0.0,0.0,0.0};
  int (*cip)[4],*cinthp,ncell,i,jcell,cid[3],cit,nj,ioffset,periodicflag=0;
  int rclist[27][3]={{-1,-1,-1},{ 0,-1,-1},{ 1,-1,-1},
		     {-1, 0,-1},{ 0, 0,-1},{ 1, 0,-1},
		     {-1, 1,-1},{ 0, 1,-1},{ 1, 1,-1},
		     {-1,-1, 0},{ 0,-1, 0},{ 1,-1, 0},
		     {-1, 0, 0},{ 0, 0, 0},{ 1, 0, 0},
		     {-1, 1, 0},{ 0, 1, 0},{ 1, 1, 0},
		     {-1,-1, 1},{ 0,-1, 1},{ 1,-1, 1},
		     {-1, 0, 1},{ 0, 0, 1},{ 1, 0, 1},
		     {-1, 1, 1},{ 0, 1, 1},{ 1, 1, 1}};
  double (*op)[3],*xi,*xj;
  M3_CELL *ci;
  VG_UNIT *unit=m3_get_unit();

  if(unit==NULL){
    fprintf(stderr,"** error : unit is NULL **\n");
    vg_exit(1);
  }
  vg_set_rcut2(unit,rcut*rcut);

#ifdef MD_PRINT_WARN
  printf("MR3calccoulomb_ci is called\n");
#endif

  if((rcut+skinnb)*ldim[0]>volume[0] || 
     (rcut+skinnb)*ldim[1]>volume[1] ||
     (rcut+skinnb)*ldim[2]>volume[2]){
    fprintf(stderr,"** error : rcut(%f)+skinnb(%f) must be smaller than a cell size (%f,%f,%f), ldim=(%d,%d,%d) **\n",
	    rcut,skinnb,volume[0]/ldim[0],volume[1]/ldim[1],volume[2]/ldim[2],
	    ldim[0],ldim[1],ldim[2]);
    MR3_exit(1);
  }

  ncell=ldim[0]*ldim[1]*ldim[2];

#if 0
  if((cip=(int (*)[4])MR3_malloc_pointer(sizeof(int)*4*n,"MR3calccoulomb_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc cip in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((op=(double (*)[3])MR3_malloc_pointer(sizeof(double)*3*n,"MR3calccoulomb_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc op in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((ci=(M3_CELL *)MR3_malloc_pointer(sizeof(M3_CELL)*ncell,"MR3calccoulomb_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc ci in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
#endif
#if 0
  if((cinthp=(int *)MR3_malloc_pointer(sizeof(int)*n,"MR3calccoulomb_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc cinthp in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
#endif
#if 0
  if((xi=(double *)MR3_malloc_pointer(sizeof(double)*n*3,"MR3calccoulomb_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc xi in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((xj=(double *)MR3_malloc_pointer(sizeof(double)*n*3,"MR3calccoulomb_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc xj in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  //  assign_cell(n,x,xmin,volume,ldim,cip,op,ci,cinthp);
  //  make_contiguous_x(n,x,cip,ci,cinthp,xi);
  //  make_contiguous_x(n,x,cip,ci,cinthp,xj);
#endif

#if 1
  // work
#ifdef MD_PRINT_WARN
  printf("periodicflag=%d in MR3calccoulomb_ci\n",periodicflag);
#endif
  MR3calccoulomb_ij(n,x,q,force,n,x,q,rscale,tblno,volume[0],periodicflag);
#endif

#if 0
  for(icell=0;icell<ncell;icell++){ // icell loop
    ni=ci[icell].size;
    ioffset=ci[icell].base;
    for(xyz=0;xyz<3;xyz++) cid[xyz]=cip[ioffset][xyz];
    cit=(cid[2]*ldim[1]+cid[1])*ldim[0]+cid[0];
#if 1
    if(cit!=icell){
      fprintf(stderr,"** error : cit must be same as icell **\n");
      MR3_exit(1);
    } 
#endif
    for(i=0;i<ni;i++){
      for(xyz=0;xyz<3;xyz++){
	xi[i*3+xyz]=
      }
    }
    for(jcell=0;jcell<27;jcell++){ // jcell loop
      nj=0;
      for(j=0;j<ci[cit].size;j++){
	for(xyz=0;xyz<3;xyz++){
	  cid[xyz]=cip[ioffset][xyz];
	}
	cit=(cid[2]*ldim[1]+cid[1])*ldim[0]+cid[0];
	for(xyz=0;xyz<3;xyz++){
	  xj[nj*3+xyz]=x[(ci[cit].base+j)*3+xyz];
	}
	nj++;
      }
      //      MR3calccoulomb_ij(,x,q,
    }
  }
#endif

#if 0
  MR3_free_pointer(cip,"MR3calccoulomb_ci");
  MR3_free_pointer(op,"MR3calccoulomb_ci");
  MR3_free_pointer(ci,"MR3calccoulomb_ci");
#endif
#if 0
  MR3_free_pointer(cinthp,"MR3calccoulomb_ci");
#endif
#if 0
  MR3_free_pointer(xi,"MR3calccoulomb_ci");
  MR3_free_pointer(xj,"MR3calccoulomb_ci");
#endif
}


void mr3calccoulomb_ci_(int *n, double x[], double q[], double force[],
			double *rscale, int *tblno,
			double *rcut, double *skinnb,
			double volume[3], int ldim[3])
{
  MR3calccoulomb_ci(*n,x,q,force,*rscale,*tblno,
		    *rcut,*skinnb,volume,ldim);
}


void MR3calccoulomb_ij_ci(int ni, double xi[], double qi[], double force[],
			  int nj, double xj[], double qj[],
			  double rscale, int tblno,
			  double rcut, double skinnb,
			  double volume[3], int ldim[3])
{
#ifdef MD_CELLINDEX
  double xmin[3]={0.0,0.0,0.0};
  int (*cipi)[4],*cinthpi,ncell;
  int (*cipj)[4],*cinthpj;
  int i,j,jcell,cidi[3],cidj[3],ioffset,joffset,cit,periodicflag=1;
  int nii,njj,icell,xyz,jci,ni2,nj2,ni2max,nj2max;
  int rclist[27][3]={{-1,-1,-1},{ 0,-1,-1},{ 1,-1,-1},
		     {-1, 0,-1},{ 0, 0,-1},{ 1, 0,-1},
		     {-1, 1,-1},{ 0, 1,-1},{ 1, 1,-1},
		     {-1,-1, 0},{ 0,-1, 0},{ 1,-1, 0},
		     {-1, 0, 0},{ 0, 0, 0},{ 1, 0, 0},
		     {-1, 1, 0},{ 0, 1, 0},{ 1, 1, 0},
		     {-1,-1, 1},{ 0,-1, 1},{ 1,-1, 1},
		     {-1, 0, 1},{ 0, 0, 1},{ 1, 0, 1},
		     {-1, 1, 1},{ 0, 1, 1},{ 1, 1, 1}};
  double (*opi)[3],(*opj)[3],*xi2,*xj2,*f2,*qi2,*qj2,*f3=NULL;
  M3_CELL *cii,*cij;
  int blocksize=VG_MINIMUM_PARTICLE_BLOCK_I;
  //  int blocksize=1;

  VG_UNIT *unit=m3_get_unit();

  if(unit==NULL){
    fprintf(stderr,"** error : unit is NULL **\n");
    vg_exit(1);
  }
  vg_set_rcut2(unit,rcut*rcut);

  {
    static int ini=0;
    if(ini==0){
#ifdef MD_PRINT_WARN
      printf("*** MR3calccoulomb_ci is called. periodicflag=%d, rscale=%f\n",periodicflag,rscale);
#endif
      ini=1;
    }
  }

  // check rcut and skinnb
  if((rcut+skinnb)*ldim[0]>volume[0] || 
     (rcut+skinnb)*ldim[1]>volume[1] ||
     (rcut+skinnb)*ldim[2]>volume[2]){
    fprintf(stderr,"** error : rcut(%f)+skinnb(%f) must be smaller than a cell size (%f,%f,%f), ldim=(%d,%d,%d) **\n",
	    rcut,skinnb,volume[0]/ldim[0],volume[1]/ldim[1],volume[2]/ldim[2],
	    ldim[0],ldim[1],ldim[2]);
    MR3_exit(1);
  }

  // assign cell
  ncell=ldim[0]*ldim[1]*ldim[2];
  if((cipi=(int (*)[4])MR3_malloc_pointer(sizeof(int)*4*ni,"MR3calccoulomb_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc cipi in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((opi=(double (*)[3])MR3_malloc_pointer(sizeof(double)*3*ni,"MR3calccoulomb_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc opi in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((cii=(M3_CELL *)MR3_malloc_pointer(sizeof(M3_CELL)*ncell,"MR3calccoulomb_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc cii in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((cinthpi=(int *)MR3_malloc_pointer(sizeof(int)*ni,"MR3calccoulomb_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc cinthpi in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((cipj=(int (*)[4])MR3_malloc_pointer(sizeof(int)*4*nj,"MR3calccoulomb_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc cipj in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((opj=(double (*)[3])MR3_malloc_pointer(sizeof(double)*3*nj,"MR3calccoulomb_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc opj in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((cij=(M3_CELL *)MR3_malloc_pointer(sizeof(M3_CELL)*ncell,"MR3calccoulomb_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc cij in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((cinthpj=(int *)MR3_malloc_pointer(sizeof(int)*nj,"MR3calccoulomb_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc cinthpj in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  assign_cell(ni,xi,xmin,volume,ldim,blocksize,cipi,opi,cii,cinthpi);
  assign_cell(nj,xj,xmin,volume,ldim,blocksize,cipj,opj,cij,cinthpj);

  // copy to contiguous array
  for(i=ni2=0;i<ncell;i++) ni2+=(cii[i].size+blocksize-1)/blocksize*blocksize;
  for(i=nj2=0;i<ncell;i++) nj2+=(cij[i].size+blocksize-1)/blocksize*blocksize;
  ni2max=ni2+ncell*blocksize;
  nj2max=nj2+ncell*blocksize;
  if((xi2=(double *)MR3_malloc_pointer(sizeof(double)*ni2*3,"MR3calccoulomb_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc xi2 in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((xj2=(double *)MR3_malloc_pointer(sizeof(double)*nj2*3,"MR3calccoulomb_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc xj2 in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((qi2=(double *)MR3_malloc_pointer(sizeof(double)*ni2,"MR3calccoulomb_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc qi2 in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((qj2=(double *)MR3_malloc_pointer(sizeof(double)*nj2,"MR3calccoulomb_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc qj2 in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((f2=(double *)MR3_malloc_pointer(sizeof(double)*ni2*3,"MR3calccoulomb_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc f2 in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  bzero(xi2,sizeof(double)*ni2*3);
  bzero(xj2,sizeof(double)*nj2*3);
  bzero(qi2,sizeof(double)*ni2);
  bzero(qj2,sizeof(double)*nj2);
  bzero(f2,sizeof(double)*ni2*3);
  make_contiguous_xq(ni,xi,qi,xmin,volume,cipi,cii,cinthpi,xi2,qi2);
  make_contiguous_xq(nj,xj,qj,xmin,volume,cipj,cij,cinthpj,xj2,qj2);

#if 0 // all-i all-j with contigous array, work
  if((f3=(double *)MR3_malloc_pointer(sizeof(double)*ni2*3,"MR3calccoulomb_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc f3 **\n");
    MR3_exit(1);
  }
  bzero(f3,sizeof(double)*ni2*3);
  MR3calccoulomb_ij(ni2,xi2,qi2,f3,nj2,xj2,qj2,rscale,tblno,volume[0],periodicflag);
  for(i=0;i<ni;i++){
    int cellindex;
    cellindex=cipi[i][3];
    ioffset=cii[cellindex].base+cinthpi[i];
    for(xyz=0;xyz<3;xyz++){
      force[i*3+xyz]+=f3[ioffset*3+xyz];
    }
  }
#endif

#if 0 // cellindex-i, all-j with contigous array, work
  tblno+=10;printf("tblno is changed to %d in MR3calccoulomb_ij_ci\n",tblno);
  m3_get_unit()->jsize=sizeof(VG_JVEC)*nj2max;
  m3_get_unit()->isize=sizeof(VG_IVEC)*ni2max;
  for(icell=0;icell<ncell;icell++){
    nii=(cii[icell].size+blocksize-1)/blocksize*blocksize;
    ioffset=cii[icell].base;
    MR3calccoulomb_ij(nii,xi2+ioffset*3,qi2+ioffset,f2+ioffset*3,nj2,xj2,qj2,rscale,tblno,volume[0],periodicflag);
  }
  for(i=0;i<ni;i++){
    int cellindex;
    cellindex=cipi[i][3];
    ioffset=cii[cellindex].base+cinthpi[i];
    for(xyz=0;xyz<3;xyz++){
      force[i*3+xyz]+=f2[ioffset*3+xyz];
    }
  }
#endif

#if 1 // cellindex-i, cellindex-j
  tblno+=10;
  {
    static int ini=0;
    if(ini==0){
#ifdef MD_PRINT_WARN
      printf("tblno is changed to %d in MR3calccoulomb_ij_ci\n",tblno);
#endif
      ini=1;
    }
  }
  m3_get_unit()->jsize=sizeof(VG_JVEC)*nj2max;
  m3_get_unit()->isize=sizeof(VG_IVEC)*ni2max;
  m3_get_unit()->fsize=sizeof(VG_FVEC)*ni2max;
  {
    int nicell,ii;
    VG_PSCALER *pscal;
    VG_UNIT *unit=m3_get_unit();

    if(unit->pscalers==NULL){
      if((unit->pscalers=(VG_PSCALER *)MR3_malloc_pointer(sizeof(VG_PSCALER),"pscalers in MR3calccoulomb_ij_ci"))==NULL){
	fprintf(stderr,"** error : can't malloc unit->pscalers in MR3calccoulomb_ij_ci **\n");
	vg_exit(1);
      }
    }
#if 0
    else{
      fprintf(stderr,"unit->pscal is not NULL\n");
    }
#endif
    pscal=(VG_PSCALER *)unit->pscalers;
    if(ni2>blocksize*MD_CELLINDEX_MAXIBLOCKS){
      fprintf(stderr,"** error : too few MD_CELLINDEX_MAXIBLOCKS=%d **\n",MD_CELLINDEX_MAXIBLOCKS);
      MR3_exit(1);
    }
    for(icell=nicell=0;icell<ncell;icell++){
      VG_CELL cijtmp[MD_NUM_JCELLS];
      cijtmp[27].size=27;
      nii=cii[icell].size;
      ioffset=cii[icell].base;
      cidi[0]=icell % ldim[0];
      cidi[1]=(icell/ldim[0]) % ldim[1];
      cidi[2]=(icell/ldim[0]/ldim[1]) % ldim[2];
      for(jci=0;jci<27;jci++){
	for(xyz=0;xyz<3;xyz++){
	  cidj[xyz]=cidi[xyz]+rclist[jci][xyz];
	  if(cidj[xyz]<0)          cidj[xyz]+=ldim[xyz];
	  if(cidj[xyz]>=ldim[xyz]) cidj[xyz]-=ldim[xyz];
	}
	jcell=(cidj[2]*ldim[1]+cidj[1])*ldim[0]+cidj[0];
	cijtmp[jci].base=cij[jcell].base;
	cijtmp[jci].size=cij[jcell].size;
      }
      for(ii=0;ii<nii;ii+=blocksize){
#if 1
	memcpy(pscal->ci[nicell],cijtmp,sizeof(VG_CELL)*MD_NUM_JCELLS);
#else
	for(jci=0;jci<27;jci++){
	  pscal->ci[nicell][jci].base=cijtmp[jci].base;
	  pscal->ci[nicell][jci].size=cijtmp[jci].size;
	}
#endif
	nicell++;
	if(nicell>MD_CELLINDEX_MAXIBLOCKS){
	  fprintf(stderr,"** error : nicell is too large **\n");
	  vg_exit(1);
	}
      }
    }
  }
  //  for(icell=0;icell<ncell;icell++){
  //    nii=(cii[icell].size+blocksize-1)/blocksize*blocksize;
  //    ioffset=cii[icell].base;
  //    MR3calccoulomb_ij(nii,xi2+ioffset*3,qi2+ioffset,f2+ioffset*3,nj2,xj2,qj2,rscale,tblno,volume[0],periodicflag);
  //  }
  MR3calccoulomb_ij(ni2,xi2,qi2,f2,nj2,xj2,qj2,rscale,tblno,volume[0],periodicflag);
  for(i=0;i<ni;i++){
    int cellindex;
    cellindex=cipi[i][3];
    ioffset=cii[cellindex].base+cinthpi[i];
    for(xyz=0;xyz<3;xyz++){
      force[i*3+xyz]+=f2[ioffset*3+xyz];
    }
  }
#endif


  MR3_free_pointer(cipi,"MR3calccoulomb_ij_ci");
  MR3_free_pointer(opi,"MR3calccoulomb_ij_ci");
  MR3_free_pointer(cii,"MR3calccoulomb_ij_ci");
  MR3_free_pointer(cinthpi,"MR3calccoulomb_ij_ci");
  MR3_free_pointer(cipj,"MR3calccoulomb_ij_ci");
  MR3_free_pointer(opj,"MR3calccoulomb_ij_ci");
  MR3_free_pointer(cij,"MR3calccoulomb_ij_ci");
  MR3_free_pointer(cinthpj,"MR3calccoulomb_ij_ci");
  MR3_free_pointer(xi2,"MR3calccoulomb_ij_ci");
  MR3_free_pointer(xj2,"MR3calccoulomb_ij_ci");
  MR3_free_pointer(qi2,"MR3calccoulomb_ij_ci");
  MR3_free_pointer(qj2,"MR3calccoulomb_ij_ci");
  MR3_free_pointer(f2,"MR3calccoulomb_ij_ci");
  if(f3!=NULL) MR3_free_pointer(f3,"MR3calccoulomb_ij_ci");
#else
  fprintf(stderr,"** MR3calccoulomb_ij_ci is not supported **\n");
  exit(1);
#endif
}


void MR3calccoulomb_ij_ci_old(int ni, double xi[], double qi[], double force[],
			      int nj, double xj[], double qj[],
			      double rscale, int tblno,
			      double rcut, double skinnb,
			      double volume[3], int ldim[3])
{
#ifdef MD_CELLINDEX
  double xmin[3]={0.0,0.0,0.0};
  int (*cipi)[4],*cinthpi,ncell;
  int (*cipj)[4],*cinthpj;
  int i,j,jcell,cidi[3],cidj[3],ioffset,joffset,cit,periodicflag=1;
  int nii,njj,icell,xyz,jci,ni2max,nj2max,ni2,nj2;
  int rclist[27][3]={{-1,-1,-1},{ 0,-1,-1},{ 1,-1,-1},
		     {-1, 0,-1},{ 0, 0,-1},{ 1, 0,-1},
		     {-1, 1,-1},{ 0, 1,-1},{ 1, 1,-1},
		     {-1,-1, 0},{ 0,-1, 0},{ 1,-1, 0},
		     {-1, 0, 0},{ 0, 0, 0},{ 1, 0, 0},
		     {-1, 1, 0},{ 0, 1, 0},{ 1, 1, 0},
		     {-1,-1, 1},{ 0,-1, 1},{ 1,-1, 1},
		     {-1, 0, 1},{ 0, 0, 1},{ 1, 0, 1},
		     {-1, 1, 1},{ 0, 1, 1},{ 1, 1, 1}};
  double (*opi)[3],(*opj)[3],*xi2,*xj2,*f2,*qi2,*qj2,*f3=NULL;
  M3_CELL *cii,*cij;
  //  int blocksize=VG_MINIMUM_PARTICLE_BLOCK_I;
  int blocksize=1;

  VG_UNIT *unit=m3_get_unit();

  if(unit==NULL){
    fprintf(stderr,"** error : unit is NULL **\n");
    vg_exit(1);
  }
  vg_set_rcut2(unit,rcut*rcut);

#ifdef MD_PRINT_WARN
  printf("*** MR3calccoulomb_ij_ci_old is called. periodicflag=%d, rscale=%f\n",periodicflag,rscale);
#endif

  if((rcut+skinnb)*ldim[0]>volume[0] || 
     (rcut+skinnb)*ldim[1]>volume[1] ||
     (rcut+skinnb)*ldim[2]>volume[2]){
    fprintf(stderr,"** error : rcut(%f)+skinnb(%f) must be smaller than a cell size (%f,%f,%f), ldim=(%d,%d,%d) **\n",
	    rcut,skinnb,volume[0]/ldim[0],volume[1]/ldim[1],volume[2]/ldim[2],
	    ldim[0],ldim[1],ldim[2]);
    MR3_exit(1);
  }

  ncell=ldim[0]*ldim[1]*ldim[2];
  ni2max=ni+ncell*blocksize;
  nj2max=nj+ncell*blocksize;
  if((cipi=(int (*)[4])MR3_malloc_pointer(sizeof(int)*4*ni,"MR3calccoulomb_ij_ci_old"))==NULL){
    fprintf(stderr,"** error : can't malloc cipi in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((opi=(double (*)[3])MR3_malloc_pointer(sizeof(double)*3*ni,"MR3calccoulomb_ij_ci_old"))==NULL){
    fprintf(stderr,"** error : can't malloc opi in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((cii=(M3_CELL *)MR3_malloc_pointer(sizeof(M3_CELL)*ncell,"MR3calccoulomb_ij_ci_old"))==NULL){
    fprintf(stderr,"** error : can't malloc cii in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((cinthpi=(int *)MR3_malloc_pointer(sizeof(int)*ni,"MR3calccoulomb_ij_ci_old"))==NULL){
    fprintf(stderr,"** error : can't malloc cinthpi in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((cipj=(int (*)[4])MR3_malloc_pointer(sizeof(int)*4*nj,"MR3calccoulomb_ij_ci_old"))==NULL){
    fprintf(stderr,"** error : can't malloc cipj in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((opj=(double (*)[3])MR3_malloc_pointer(sizeof(double)*3*nj,"MR3calccoulomb_ij_ci_old"))==NULL){
    fprintf(stderr,"** error : can't malloc opj in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((cij=(M3_CELL *)MR3_malloc_pointer(sizeof(M3_CELL)*ncell,"MR3calccoulomb_ij_ci_old"))==NULL){
    fprintf(stderr,"** error : can't malloc cij in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((cinthpj=(int *)MR3_malloc_pointer(sizeof(int)*nj,"MR3calccoulomb_ij_ci_old"))==NULL){
    fprintf(stderr,"** error : can't malloc cinthpj in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((xi2=(double *)MR3_malloc_pointer(sizeof(double)*ni2max*3,"MR3calccoulomb_ij_ci_old"))==NULL){
    fprintf(stderr,"** error : can't malloc xi2 in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((xj2=(double *)MR3_malloc_pointer(sizeof(double)*nj2max*3,"MR3calccoulomb_ij_ci_old"))==NULL){
    fprintf(stderr,"** error : can't malloc xj2 in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((qi2=(double *)MR3_malloc_pointer(sizeof(double)*ni2max,"MR3calccoulomb_ij_ci_old"))==NULL){
    fprintf(stderr,"** error : can't malloc qi2 in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((qj2=(double *)MR3_malloc_pointer(sizeof(double)*nj2max,"MR3calccoulomb_ij_ci_old"))==NULL){
    fprintf(stderr,"** error : can't malloc qj2 in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((f2=(double *)MR3_malloc_pointer(sizeof(double)*ni2max*3,"MR3calccoulomb_ij_ci_old"))==NULL){
    fprintf(stderr,"** error : can't malloc f2 in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  bzero(xi2,sizeof(double)*ni2max*3);
  bzero(xj2,sizeof(double)*nj2max*3);
  bzero(qi2,sizeof(double)*ni2max);
  bzero(qj2,sizeof(double)*nj2max);
  bzero(f2,sizeof(double)*ni2max*3);
  assign_cell(ni,xi,xmin,volume,ldim,blocksize,cipi,opi,cii,cinthpi);
  assign_cell(nj,xj,xmin,volume,ldim,blocksize,cipj,opj,cij,cinthpj);
  make_contiguous_xq(ni,xi,qi,xmin,volume,cipi,cii,cinthpi,xi2,qi2);
  make_contiguous_xq(nj,xj,qj,xmin,volume,cipj,cij,cinthpj,xj2,qj2);

#if 0 // all-i all-j, work
  if((f3=(double *)MR3_malloc_pointer(sizeof(double)*ni*3,"MR3calccoulomb_ij_ci_old"))==NULL){
    fprintf(stderr,"** error : can't malloc f3 **\n");
    MR3_exit(1);
  }
  bzero(f3,sizeof(double)*ni*3);

#if 0
  m3_get_unit()->debug_flag=1;fprintf(stderr,"debug_flag is set 1\n");
  //  MR3calccoulomb_ij(ni,xi,qi,f3,nj,xj,qj,rscale,tblno,volume[0],periodicflag);
  {
    int size=sizeof(VG_JVEC)*NJ_ROUNDUP(nj);
    fprintf(stderr,"unit->jvectors=%x before free\n",m3_get_unit()->jvectors);
    MR3_free_pointer(m3_get_unit()->jvectorMR3calccoulomb_ij_ci_olds);
    m3_get_unit()->jvectors=MR3_malloc_pinter(size,"MR3calccoulomb_ij_ci_old");
    fprintf(stderr,"unit->jvectors=%x after malloc, size=%d\n",m3_get_unit()->jvectors,size);
    //  vg_set_vectors_pos_charge(m3_get_unit(),nj,(double (*)[3])xj,qj);
    debug_mother(m3_get_unit());
  }
  m3_get_unit()->debug_flag=0;fprintf(stderr,"debug_flag is set 0\n");
#endif
#if 1 // adding to force (this should be used for correct value)
  for(i=0;i<ni;i++){
    for(xyz=0;xyz<3;xyz++){
      force[i*3+xyz]+=f3[i*3+xyz];
    }
  }
  for(i=0;i<3;i++){
    //    printf("force[%d]=%e %e %e\n",i,force[i*3],force[i*3+1],force[i*3+2]);
  }
#else // adding to f3 (this should be used for comparing cellindex)
  for(i=0;i<ni;i++){
    for(xyz=0;xyz<3;xyz++){
      f3[i*3+xyz]+=force[i*3+xyz];
    }
  }
  for(i=0;i<3;i++){
    //    printf("f3[%d]=%e %e %e\n",i,f3[i*3],f3[i*3+1],f3[i*3+2]);
  }
#endif
#endif

#if 0 // dummy command
  debug_mother(m3_get_unit());printf("debug_mother is called \n");
#endif
#if 0 // special-command to avoid the bug for cellindex-i:all-j
  {
    int size=sizeof(VG_JVEC)*NJ_ROUNDUP(nj);
    fprintf(stderr,"unit->jvectors=%x before free\n",m3_get_unit()->jvectors);
    MR3_free_pointer(m3_get_unit()->jvectorMR3calccoulomb_ij_ci_olds);
    m3_get_unit()->jvectors=MR3_malloc_pointer(size,"MR3calccoulomb_ij_ci_old");
    fprintf(stderr,"unit->jvectors=%x after malloc, size=%d\n",m3_get_unit()->jvectors,size);
    //  vg_set_vectors_pos_charge(m3_get_unit(),nj,(double (*)[3])xj,qj);
    debug_mother(m3_get_unit());
  }
#endif
#if 0 // cellindex-i, all-j, work
  for(icell=0;icell<ncell;icell++){
    nii=cii[icell].size;
    ioffset=cii[icell].base;
    //    printf("calculating icell=%d nii=%d ioffset=%d\n",icell,nii,ioffset);
    MR3calccoulomb_ij(nii,xi2+ioffset*3,qi2+ioffset,f2+ioffset*3,nj,xj,qj,rscale,tblno,volume[0],periodicflag);
    //    MR3calccoulomb_ij_host(nii,xi2+ioffset*3,qi2+ioffset,f2+ioffset*3,nj,xj,qj,rscale,tblno,volume[0],periodicflag);
  }
  for(i=0;i<ni;i++){
    int cellindex;
    cellindex=cipi[i][3];
    ioffset=cii[cellindex].base+cinthpi[i];
    for(xyz=0;xyz<3;xyz++){
      force[i*3+xyz]+=f2[ioffset*3+xyz];
    }
  }
#if 0
  for(i=0;i<ni;i++){
    printf("force[%d]=%e %e %e\n",i,force[i*3],force[i*3+1],force[i*3+2]);
  }
#endif
#if 0
  for(i=0;i<ni*3;i++){
    if(force[i]!=f3[i]){
      printf("i=%d xyz=%d force=%e f3=%e\n",i/3,i % 3,force[i],f3[i]);
    }
  }
#endif
#endif

#if 1 // cellindex-i, cellindex-j
  m3_get_unit()->jsize=sizeof(VG_JVEC)*NJ_ROUNDUP(nj);
  m3_get_unit()->isize=sizeof(VG_IVEC)*NJ_ROUNDUP(ni);
  bzero(f2,sizeof(double)*ni*3);
  for(icell=0;icell<ncell;icell++){
    nii=cii[icell].size;
    ioffset=cii[icell].base;
    cidi[0]=icell % ldim[0];
    cidi[1]=(icell/ldim[0]) % ldim[1];
    cidi[2]=(icell/ldim[0]/ldim[1]) % ldim[2];
    for(jci=0;jci<27;jci++){
      for(xyz=0;xyz<3;xyz++){
	cidj[xyz]=cidi[xyz]+rclist[jci][xyz];
	if(cidj[xyz]<0)          cidj[xyz]+=ldim[xyz];
	if(cidj[xyz]>=ldim[xyz]) cidj[xyz]-=ldim[xyz];
      }
      jcell=(cidj[2]*ldim[1]+cidj[1])*ldim[0]+cidj[0];
      njj=cij[jcell].size;
      joffset=cij[jcell].base;
      //      printf("calculating icell=%d nii=%d ioffset=%d jci=%d jcell=%d njj=%d joffset=%d\n",icell,nii,ioffset,jci,jcell,njj,joffset);
      //      MR3calccoulomb_ij(nii,xi2+ioffset*3,qi2+ioffset,f2+ioffset*3,njj,xj2+joffset*3,qj2+joffset,rscale,tblno,volume[0],periodicflag);
      MR3calccoulomb_ij_host(nii,xi2+ioffset*3,qi2+ioffset,f2+ioffset*3,njj,xj2+joffset*3,qj2+joffset,rscale,tblno,volume[0],periodicflag);
    }
  }
  for(i=0;i<ni;i++){
    int cellindex;
    cellindex=cipi[i][3];
    ioffset=cii[cellindex].base+cinthpi[i];
    for(xyz=0;xyz<3;xyz++){
      force[i*3+xyz]+=f2[ioffset*3+xyz];
    }
  }
  for(i=0;i<3;i++){
    //    printf("force[%d]=%e %e %e\n",i,force[i*3],force[i*3+1],force[i*3+2]);
  }
#endif

  MR3_free_pointer(cipi,"MR3calccoulomb_ij_ci_old");
  MR3_free_pointer(opi,"MR3calccoulomb_ij_ci_old");
  MR3_free_pointer(cii,"MR3calccoulomb_ij_ci_old");
  MR3_free_pointer(cinthpi,"MR3calccoulomb_ij_ci_old");
  MR3_free_pointer(cipj,"MR3calccoulomb_ij_ci_old");
  MR3_free_pointer(opj,"MR3calccoulomb_ij_ci_old");
  MR3_free_pointer(cij,"MR3calccoulomb_ij_ci_old");
  MR3_free_pointer(cinthpj,"MR3calccoulomb_ij_ci_old");
  MR3_free_pointer(xi2,"MR3calccoulomb_ij_ci_old");
  MR3_free_pointer(xj2,"MR3calccoulomb_ij_ci_old");
  MR3_free_pointer(qi2,"MR3calccoulomb_ij_ci_old");
  MR3_free_pointer(qj2,"MR3calccoulomb_ij_ci_old");
  MR3_free_pointer(f2,"MR3calccoulomb_ij_ci_old");
  if(f3!=NULL) MR3_free_pointer(f3,"MR3calccoulomb_ij_ci_old");
#else
  fprintf(stderr,"** MR3calccoulomb_ij_ci_old is not supported **\n");
  exit(1);
#endif
}


void MR3calcvdw_ij_ci(int ni, double xi[], int atypei[], double force[],
		      int nj, double xj[], int atypej[],
		      int nat, double gscale[], double rscale[],
		      int tblno,
		      double rcut, double skinnb,
		      double volume[3], int ldim[3])
{
#ifdef MD_CELLINDEX
  double xmin[3]={0.0,0.0,0.0};
  int (*cipi)[4],*cinthpi,ncell;
  int (*cipj)[4],*cinthpj;
  int i,j,jcell,cidi[3],cidj[3],ioffset,joffset,cit,periodicflag=1;
  int nii,njj,icell,xyz,jci,ni2,nj2,ni2max,nj2max;
  int rclist[27][3]={{-1,-1,-1},{ 0,-1,-1},{ 1,-1,-1},
		     {-1, 0,-1},{ 0, 0,-1},{ 1, 0,-1},
		     {-1, 1,-1},{ 0, 1,-1},{ 1, 1,-1},
		     {-1,-1, 0},{ 0,-1, 0},{ 1,-1, 0},
		     {-1, 0, 0},{ 0, 0, 0},{ 1, 0, 0},
		     {-1, 1, 0},{ 0, 1, 0},{ 1, 1, 0},
		     {-1,-1, 1},{ 0,-1, 1},{ 1,-1, 1},
		     {-1, 0, 1},{ 0, 0, 1},{ 1, 0, 1},
		     {-1, 1, 1},{ 0, 1, 1},{ 1, 1, 1}};
  double (*opi)[3],(*opj)[3],*xi2,*xj2,*f2,*f3=NULL;
  int *atypei2,*atypej2;
  M3_CELL *cii,*cij;
  double *gscale2,*rscale2;
  int blocksize=VG_MINIMUM_PARTICLE_BLOCK_I;
  //  int blocksize=1;

  VG_UNIT *unit=m3_get_unit();

  if(unit==NULL){
    fprintf(stderr,"** error : unit is NULL **\n");
    vg_exit(1);
  }
  vg_set_rcut2(unit,rcut*rcut);

#ifndef MD_CELLINDEX
  fprintf(stderr,"** MR3calcvdw_ij_ci is not supported **\n");
  exit(1);
#endif
  {
    static int ini=0;
    if(ini==0){
#ifdef MD_PRINT_WARN
      printf("*** MR3calcvdw_ij_ci is called. periodicflag=%d\n",periodicflag);
#endif
      ini=1;
    }
  }

#if 0
  MR3calcvdw_ij(ni,xi,atypei,force,nj,xj,atypej,nat,gscale,rscale,tblno,volume[0],periodicflag);
  return;
#endif

  // check rcut and skinnb
  if((rcut+skinnb)*ldim[0]>volume[0] || 
     (rcut+skinnb)*ldim[1]>volume[1] ||
     (rcut+skinnb)*ldim[2]>volume[2]){
    fprintf(stderr,"** error : rcut(%f)+skinnb(%f) must be smaller than a cell size (%f,%f,%f), ldim=(%d,%d,%d) **\n",
	    rcut,skinnb,volume[0]/ldim[0],volume[1]/ldim[1],volume[2]/ldim[2],
	    ldim[0],ldim[1],ldim[2]);
    MR3_exit(1);
  }

  // assign cell
  ncell=ldim[0]*ldim[1]*ldim[2];
  if((cipi=(int (*)[4])MR3_malloc_pointer(sizeof(int)*4*ni,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc cipi in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((opi=(double (*)[3])MR3_malloc_pointer(sizeof(double)*3*ni,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc opi in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((cii=(M3_CELL *)MR3_malloc_pointer(sizeof(M3_CELL)*ncell,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc cii in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((cinthpi=(int *)MR3_malloc_pointer(sizeof(int)*ni,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc cinthpi in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((cipj=(int (*)[4])MR3_malloc_pointer(sizeof(int)*4*nj,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc cipj in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((opj=(double (*)[3])MR3_malloc_pointer(sizeof(double)*3*nj,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc opj in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((cij=(M3_CELL *)MR3_malloc_pointer(sizeof(M3_CELL)*ncell,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc cij in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((cinthpj=(int *)MR3_malloc_pointer(sizeof(int)*nj,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc cinthpj in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  assign_cell(ni,xi,xmin,volume,ldim,blocksize,cipi,opi,cii,cinthpi);
  assign_cell(nj,xj,xmin,volume,ldim,blocksize,cipj,opj,cij,cinthpj);

  // copy to contiguous array
  for(i=ni2=0;i<ncell;i++) ni2+=(cii[i].size+blocksize-1)/blocksize*blocksize;
  for(i=nj2=0;i<ncell;i++) nj2+=(cij[i].size+blocksize-1)/blocksize*blocksize;
  ni2max=ni2+ncell*blocksize;
  nj2max=nj2+ncell*blocksize;
  if((xi2=(double *)MR3_malloc_pointer(sizeof(double)*ni2*3,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc xi2 in MR3calcvdw_ci **\n");
    MR3_exit(1);
  }
  if((xj2=(double *)MR3_malloc_pointer(sizeof(double)*nj2*3,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc xj2 in MR3calcvdw_ci **\n");
    MR3_exit(1);
  }
  if((atypei2=(int *)MR3_malloc_pointer(sizeof(int)*ni2,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc atypei2 in MR3calcvdw_ci **\n");
    MR3_exit(1);
  }
  if((atypej2=(int *)MR3_malloc_pointer(sizeof(int)*nj2,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc atypej2 in MR3calcvdw_ci **\n");
    MR3_exit(1);
  }
  if((f2=(double *)MR3_malloc_pointer(sizeof(double)*ni2*3,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc f2 in MR3calcvdw_ci **\n");
    MR3_exit(1);
  }
  bzero(xi2,sizeof(double)*ni2*3);
  bzero(xj2,sizeof(double)*nj2*3);
  bzero(atypei2,sizeof(int)*ni2);
  bzero(atypej2,sizeof(int)*nj2);
  bzero(f2,sizeof(double)*ni2*3);
  make_contiguous_xa(ni,xi,atypei,xmin,volume,cipi,cii,cinthpi,xi2,atypei2);
  make_contiguous_xa(nj,xj,atypej,xmin,volume,cipj,cij,cinthpj,xj2,atypej2);

  // new gscale, rscale
  nat++;
  if((gscale2=(double *)MR3_malloc_pointer(sizeof(double)*nat*nat,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc gscales in MR3calcvdw_ij_ci **\n");
    MR3_exit(1);
  }
  if((rscale2=(double *)MR3_malloc_pointer(sizeof(double)*nat*nat,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc rscales in MR3calcvdw_ij_ci **\n");
    MR3_exit(1);
  }
  for(i=0;i<nat;i++){
    for(j=0;j<nat;j++){
      if(i==0 || j==0){
	gscale2[i*nat+j]=0.0;
	rscale2[i*nat+j]=1.0;
      }
      else{
	gscale2[i*nat+j]=gscale[(i-1)*(nat-1)+j-1];
	rscale2[i*nat+j]=rscale[(i-1)*(nat-1)+j-1];
      }
    }
  }

#if 0 // all-i all-j with contigous array, work
  //  tblno+=10;printf("tblno is changed to %d in MR3calcvdw_ij_ci\n",tblno);
  if((f3=(double *)MR3_malloc_pointer(sizeof(double)*ni2*3,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc f3 **\n");
    MR3_exit(1);
  }
  bzero(f3,sizeof(double)*ni2*3);
  MR3calcvdw_ij_core(ni2,xi2,atypei2,f3,nj2,xj2,atypej2,nat,gscale2,rscale2,tblno,volume[0],periodicflag);
  for(i=0;i<ni;i++){
    int cellindex;
    cellindex=cipi[i][3];
    ioffset=cii[cellindex].base+cinthpi[i];
    for(xyz=0;xyz<3;xyz++){
      force[i*3+xyz]+=f3[ioffset*3+xyz];
    }
  }
#endif

#if 0 // cellindex-i, all-j with contigous array, not work
  //  tblno+=10;printf("tblno is changed to %d in MR3calccoulomb_ij_ci\n",tblno);
  m3_get_unit()->jsize=sizeof(VG_JVEC)*nj2max;
  m3_get_unit()->isize=sizeof(VG_IVEC)*ni2max;
  for(icell=0;icell<ncell;icell++){
    nii=(cii[icell].size+blocksize-1)/blocksize*blocksize;
    ioffset=cii[icell].base;
    //    printf("nii=%d(%d) ioffset=%d(%d)\n",nii,nii % 64,ioffset,ioffset % 64);
    MR3calcvdw_ij_core(nii,xi2+ioffset*3,atypei2+ioffset,f2+ioffset*3,nj2,xj2,atypej2,nat,gscale2,rscale2,tblno,volume[0],periodicflag);
  }
  for(i=0;i<ni;i++){
    int cellindex;
    cellindex=cipi[i][3];
    ioffset=cii[cellindex].base+cinthpi[i];
    for(xyz=0;xyz<3;xyz++){
      force[i*3+xyz]+=f2[ioffset*3+xyz];
    }
  }
#endif

#if 1 // cellindex-i, cellindex-j
  tblno+=10;
  {
    static int ini=0;
    if(ini==0){
#ifdef MD_PRINT_WARN
      printf("tblno is changed to %d in MR3calcvdw_ij_ci\n",tblno);
#endif
      ini=1;
    }
  }
  m3_get_unit()->jsize=sizeof(VG_JVEC)*nj2max;
  m3_get_unit()->isize=sizeof(VG_IVEC)*ni2max;
  m3_get_unit()->fsize=sizeof(VG_FVEC)*ni2max;
  {
    int nicell,ii;
    VG_PSCALER *pscal;
    VG_UNIT *unit=m3_get_unit();

    if(unit->pscalers==NULL){
      if((unit->pscalers=(VG_PSCALER *)MR3_malloc_pointer(sizeof(VG_PSCALER),"MR3calcvdw_ij_ci"))==NULL){
	fprintf(stderr,"** error : can't malloc unit->pscalers in MR3calccoulomb_ij_ci **\n");
	vg_exit(1);
      }
    }
#if 0
    else{
      fprintf(stderr,"unit->pscal is not NULL\n");
    }
#endif
    pscal=(VG_PSCALER *)unit->pscalers;
    if(ni2>blocksize*MD_CELLINDEX_MAXIBLOCKS){
      fprintf(stderr,"** error : too few MD_CELLINDEX_MAXIBLOCKS=%d **\n",MD_CELLINDEX_MAXIBLOCKS);
      MR3_exit(1);
    }
    for(icell=nicell=0;icell<ncell;icell++){
      VG_CELL cijtmp[MD_NUM_JCELLS];
      cijtmp[27].size=27;
      nii=cii[icell].size;
      ioffset=cii[icell].base;
      cidi[0]=icell % ldim[0];
      cidi[1]=(icell/ldim[0]) % ldim[1];
      cidi[2]=(icell/ldim[0]/ldim[1]) % ldim[2];
      for(jci=0;jci<27;jci++){
	for(xyz=0;xyz<3;xyz++){
	  cidj[xyz]=cidi[xyz]+rclist[jci][xyz];
	  if(cidj[xyz]<0)          cidj[xyz]+=ldim[xyz];
	  if(cidj[xyz]>=ldim[xyz]) cidj[xyz]-=ldim[xyz];
	}
	jcell=(cidj[2]*ldim[1]+cidj[1])*ldim[0]+cidj[0];
	cijtmp[jci].base=cij[jcell].base;
	cijtmp[jci].size=cij[jcell].size;
      }
      for(ii=0;ii<nii;ii+=blocksize){
#if 1
	memcpy(pscal->ci[nicell],cijtmp,sizeof(VG_CELL)*MD_NUM_JCELLS);
#else
	for(jci=0;jci<27;jci++){
	  pscal->ci[nicell][jci].base=cijtmp[jci].base;
	  pscal->ci[nicell][jci].size=cijtmp[jci].size;
	}
#endif
	nicell++;
	if(nicell>MD_CELLINDEX_MAXIBLOCKS){
	  fprintf(stderr,"** error : nicell is too large **\n");
	  vg_exit(1);
	}
      }
    }
  }
  //  for(icell=0;icell<ncell;icell++){
  //    nii=(cii[icell].size+blocksize-1)/blocksize*blocksize;
  //    ioffset=cii[icell].base;
  //    MR3calccoulomb_ij(nii,xi2+ioffset*3,qi2+ioffset,f2+ioffset*3,nj2,xj2,qj2,rscale,tblno,volume[0],periodicflag);
  //  }
  MR3calcvdw_ij_core(ni2,xi2,atypei2,f2,nj2,xj2,atypej2,nat,gscale2,rscale2,tblno,volume[0],periodicflag);
  for(i=0;i<ni;i++){
    int cellindex;
    cellindex=cipi[i][3];
    ioffset=cii[cellindex].base+cinthpi[i];
    for(xyz=0;xyz<3;xyz++){
      force[i*3+xyz]+=f2[ioffset*3+xyz];
    }
  }
#endif


  MR3_free_pointer(cipi,"MR3calcvdw_ij_ci");
  MR3_free_pointer(opi,"MR3calcvdw_ij_ci");
  MR3_free_pointer(cii,"MR3calcvdw_ij_ci");
  MR3_free_pointer(cinthpi,"MR3calcvdw_ij_ci");
  MR3_free_pointer(cipj,"MR3calcvdw_ij_ci");
  MR3_free_pointer(opj,"MR3calcvdw_ij_ci");
  MR3_free_pointer(cij,"MR3calcvdw_ij_ci");
  MR3_free_pointer(cinthpj,"MR3calcvdw_ij_ci");
  MR3_free_pointer(xi2,"MR3calcvdw_ij_ci");
  MR3_free_pointer(xj2,"MR3calcvdw_ij_ci");
  MR3_free_pointer(atypei2,"MR3calcvdw_ij_ci");
  MR3_free_pointer(atypej2,"MR3calcvdw_ij_ci");
  MR3_free_pointer(f2,"MR3calcvdw_ij_ci");
  if(f3!=NULL) MR3_free_pointer(f3,"MR3calcvdw_ij_ci");
  MR3_free_pointer(gscale2,"MR3calcvdw_ij_ci");
  MR3_free_pointer(rscale2,"MR3calcvdw_ij_ci");
#else
  fprintf(stderr,"** MR3calcvdw_ij_ci is not supported **\n");
  exit(1);
#endif
}


void MR3calccoulomb_vdw_ij_ci_virial(int ni, double xi[], double qi[], int atypei[], double force[],
				     int nj, double xj[], double qj[], int atypej[],
				     int nat, double gscales[], double rscales[], double rscale,
				     int tblnoc,
				     double rcut_org, double skinnb,
				     double volume[3], int ldim[3],
				     double virial[3][3],
                                     int flag)
{
  /*
    flag 0 -- non-periodic
         1 -- periodic
   */
  int tblnov=2;
  double xmin[3],xmax[3];
  static int (*cipi)[4],*cinthpi;
  int ncelli,ncellj,cmini[3],cmaxi[3],cminj[3],cmaxj[3],cminmaxi[3],cminmaxj[3];
  static int (*cipj)[4],*cinthpj;
  int i,j,jcell,cidi[3],cidj[3],ioffset,joffset,cit;
  int nii,njj,icell,xyz,jci,ni2,nj2,ni2max,nj2max;
  int rclist[27][3]={{-1,-1,-1},{ 0,-1,-1},{ 1,-1,-1},
		     {-1, 0,-1},{ 0, 0,-1},{ 1, 0,-1},
		     {-1, 1,-1},{ 0, 1,-1},{ 1, 1,-1},
		     {-1,-1, 0},{ 0,-1, 0},{ 1,-1, 0},
		     {-1, 0, 0},{ 0, 0, 0},{ 1, 0, 0},
		     {-1, 1, 0},{ 0, 1, 0},{ 1, 1, 0},
		     {-1,-1, 1},{ 0,-1, 1},{ 1,-1, 1},
		     {-1, 0, 1},{ 0, 0, 1},{ 1, 0, 1},
		     {-1, 1, 1},{ 0, 1, 1},{ 1, 1, 1}};
  static double (*opi)[3],(*opj)[3],*xi2,*xj2,*f2,*qi2,*qj2;
  static int *atypei2,*atypej2;
  static M3_CELL *cii,*cij;
  static double *gscale2,*rscale2;
  int blocksize=VG_MINIMUM_PARTICLE_BLOCK_I;
  VG_UNIT *unit=m3_get_unit();
  double rcut=rcut_org,potc,potv,alpha3;
  int multiplyq=0,tblnocv;
  int nicell,ii;
  VG_PSCALER *pscal;
  static int ni_bak=0;
  int cellperiodic=flag & 1; // 0 : non-periodic, 1 : periodic

#ifdef MD_MEASURE_TIME
  vg_start_timer(30);
#endif
#ifndef MD_CELLINDEX
  fprintf(stderr,"** MR3calccoulomb_vdw_ij_ci_virial is not supported **\n");
  exit(1);
#endif
#if !defined(MD_USE_QAUNION) || defined(MD_SIGMA_EPSION_IN_VGSTRUCT) || defined(MD_SORT_ATYPEI) || defined(MD_MATRIX_IN_SCALER)
  fprintf(stderr,"** error : MD_USE_QAUNION must be defined **\n");
  vg_exit(1);
#endif

  //  rcut_org=-rcut;printf("** rcut is made negative **\n");
#if defined(MD_PRINT_WARN) && 0
  printf("unit->gpuoverlapflag=%d in MR3calccoulomb_vdw_ij_ci_virial\n",unit->gpuoverlapflag);
#endif  

  if(unit->gpuoverlapflag==0 || unit->gpuoverlapflag==1){
#ifdef MD_MEASURE_TIME
    vg_start_timer(31);
#endif
#if 0
    if(ni>nj){
      static int ini=0;
      if(ini==0){
        printf("** cellperiodic is reset when ni>nj **\n");
        ini=1;
      }
      rcut_org=(rcut_org<0 ? rcut_org:-rcut_org);
    }
#endif
    if(rcut_org<0){
      rcut=(rcut_org<0 ? -rcut_org:rcut_org);
      cellperiodic=0;
    }
    
    if(unit==NULL){
      fprintf(stderr,"** error : unit is NULL **\n");
      vg_exit(1);
    }
    vg_set_rcut2(unit,rcut*rcut);

    {
      static int ini=0;
      if(ini<3){
#ifdef MD_PRINT_WARN
	printf("*** MR3calccoulomb_vdw_ij_ci_virial is called. cellperiodicflag=%d\n",cellperiodic);
#endif
	ini++;
      }
    }
    ni_bak=ni;

    for(i=0;i<3;i++){
      xmin[i]=-volume[i]*0.5;
      xmax[i]= volume[i]*0.5;
    }

    // check rcut and skinnb
    check_rcut_and_skinnb(rcut,skinnb,ldim,volume);
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(31);
    vg_start_timer(32);
#endif
    
    // assign cell
    cell_range(ni,xi,xmin,xmax,ldim,&ncelli,cmini,cmaxi,cellperiodic);
    cell_range(nj,xj,xmin,xmax,ldim,&ncellj,cminj,cmaxj,cellperiodic);
    if((cipi=(int (*)[4])MR3_malloc_pointer(sizeof(int)*4*ni,"MR3calccoulomb_vdw_ij_ci_virial"))==NULL){
      fprintf(stderr,"** error : can't malloc cipi in MR3calccoulomb_vdw_ij_ci_virial **\n");
      MR3_exit(1);
    }
    if((opi=(double (*)[3])MR3_malloc_pointer(sizeof(double)*3*ni,"MR3calccoulomb_vdw_ij_ci_virial"))==NULL){
      fprintf(stderr,"** error : can't malloc opi in MR3calccoulomb_vdw_ij_ci_virial **\n");
      MR3_exit(1);
    }
    if((cii=(M3_CELL *)MR3_malloc_pointer(sizeof(M3_CELL)*ncelli,"MR3calccoulomb_vdw_ij_ci_virial"))==NULL){
      fprintf(stderr,"** error : can't malloc cii in MR3calccoulomb_vdw_ij_ci_virial **\n");
      MR3_exit(1);
    }
    if((cinthpi=(int *)MR3_malloc_pointer(sizeof(int)*ni,"MR3calccoulomb_vdw_ij_ci_virial"))==NULL){
      fprintf(stderr,"** error : can't malloc cinthpi in MR3calccoulomb_vdw_ij_ci_virial **\n");
      MR3_exit(1);
    }
    if((cipj=(int (*)[4])MR3_malloc_pointer(sizeof(int)*4*nj,"MR3calccoulomb_vdw_ij_ci_virial"))==NULL){
      fprintf(stderr,"** error : can't malloc cipj in MR3calccoulomb_vdw_ij_ci_virial **\n");
      MR3_exit(1);
    }
    if((opj=(double (*)[3])MR3_malloc_pointer(sizeof(double)*3*nj,"MR3calccoulomb_vdw_ij_ci_virial"))==NULL){
      fprintf(stderr,"** error : can't malloc opj in MR3calccoulomb_vdw_ij_ci_virial **\n");
      MR3_exit(1);
    }
    if((cij=(M3_CELL *)MR3_malloc_pointer(sizeof(M3_CELL)*ncellj,"MR3calccoulomb_vdw_ij_ci_virial"))==NULL){
      fprintf(stderr,"** error : can't malloc cij in MR3calccoulomb_vdw_ij_ci_virial **\n");
      MR3_exit(1);
    }
    if((cinthpj=(int *)MR3_malloc_pointer(sizeof(int)*nj,"MR3calccoulomb_vdw_ij_ci_virial"))==NULL){
      fprintf(stderr,"** error : can't malloc cinthpj in MR3calccoulomb_vdw_ij_ci_virial **\n");
      MR3_exit(1);
    }
    assign_cell_noperiodic(ni,xi,xmin,xmax,ldim,blocksize,ncelli,cmini,cmaxi,cipi,cii,cinthpi,cellperiodic);
    assign_cell_noperiodic(nj,xj,xmin,xmax,ldim,blocksize,ncellj,cminj,cmaxj,cipj,cij,cinthpj,cellperiodic);
    for(i=0;i<3;i++){
      cminmaxi[i]=cmaxi[i]-cmini[i]+1;
      cminmaxj[i]=cmaxj[i]-cminj[i]+1;
    }
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(32);
    vg_start_timer(33);
#endif
    
    // copy to contiguous array
    for(i=ni2=0;i<ncelli;i++) ni2+=(cii[i].size+blocksize-1)/blocksize*blocksize;
    for(i=nj2=0;i<ncellj;i++) nj2+=(cij[i].size+blocksize-1)/blocksize*blocksize;
    ni2max=ni2+ncelli*blocksize;
    nj2max=nj2+ncellj*blocksize;
    if((xi2=(double *)MR3_malloc_pointer(sizeof(double)*ni2*3,"MR3calccoulomb_vdw_ij_ci_virial"))==NULL){
      fprintf(stderr,"** error : can't malloc xi2 in MR3calcvdw_ci **\n");
      MR3_exit(1);
    }
    if((xj2=(double *)MR3_malloc_pointer(sizeof(double)*nj2*3,"MR3calccoulomb_vdw_ij_ci_virial"))==NULL){
      fprintf(stderr,"** error : can't malloc xj2 in MR3calcvdw_ci **\n");
      MR3_exit(1);
    }
    if((qi2=(double *)MR3_malloc_pointer(sizeof(double)*ni2,"MR3calccoulomb_vdw_ij_ci"))==NULL){
      fprintf(stderr,"** error : can't malloc qi2 in MR3calccoulomb_vdw_ij_ci_virial **\n");
      MR3_exit(1);
    }
    if((qj2=(double *)MR3_malloc_pointer(sizeof(double)*nj2,"MR3calccoulomb_vdw_ij_ci"))==NULL){
      fprintf(stderr,"** error : can't malloc qj2 in MR3calccoulomb_vdw_ij_ci_virial **\n");
      MR3_exit(1);
    }
    if((atypei2=(int *)MR3_malloc_pointer(sizeof(int)*ni2,"MR3calccoulomb_vdw_ij_ci_virial"))==NULL){
      fprintf(stderr,"** error : can't malloc atypei2 in MR3calcvdw_ci **\n");
      MR3_exit(1);
    }
    if((atypej2=(int *)MR3_malloc_pointer(sizeof(int)*nj2,"MR3calccoulomb_vdw_ij_ci_virial"))==NULL){
      fprintf(stderr,"** error : can't malloc atypej2 in MR3calcvdw_ci **\n");
      MR3_exit(1);
    }
    if((f2=(double *)MR3_malloc_pointer(sizeof(double)*ni2*3,"MR3calccoulomb_vdw_ij_ci_virial"))==NULL){
      fprintf(stderr,"** error : can't malloc f2 in MR3calcvdw_ci **\n");
      MR3_exit(1);
    }
    bzero(xi2,sizeof(double)*ni2*3);
    bzero(xj2,sizeof(double)*nj2*3);
    bzero(qi2,sizeof(double)*ni2);
    bzero(qj2,sizeof(double)*nj2);
    bzero(atypei2,sizeof(int)*ni2);
    bzero(atypej2,sizeof(int)*nj2);
    bzero(f2,sizeof(double)*ni2*3);
    make_contiguous_xqa(ni,xi,qi,atypei,xmin,xmax,cipi,cii,cinthpi,xi2,qi2,atypei2);
    make_contiguous_xqa(nj,xj,qj,atypej,xmin,xmax,cipj,cij,cinthpj,xj2,qj2,atypej2);
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(33);
    vg_start_timer(34);
#endif
    
    // new gscale, rscale
    nat++;
    if((gscale2=(double *)MR3_malloc_pointer(sizeof(double)*nat*nat,"MR3calccoulomb_vdw_ij_ci_virial"))==NULL){
      fprintf(stderr,"** error : can't malloc gscales in MR3calccoulomb_vdw_ij_ci_virial **\n");
      MR3_exit(1);
    }
    if((rscale2=(double *)MR3_malloc_pointer(sizeof(double)*nat*nat,"MR3calccoulomb_vdw_ij_ci_virial"))==NULL){
      fprintf(stderr,"** error : can't malloc rscales in MR3calccoulomb_vdw_ij_ci_virial **\n");
      MR3_exit(1);
    }
    for(i=0;i<nat;i++){
      for(j=0;j<nat;j++){
	if(i==0 || j==0){
	  gscale2[i*nat+j]=0.0;
	  rscale2[i*nat+j]=1.0;
	}
	else{
	  gscale2[i*nat+j]=gscales[(i-1)*(nat-1)+j-1];
	  rscale2[i*nat+j]=rscales[(i-1)*(nat-1)+j-1];
	}
      }
    }
    
    // cellindex-i, cellindex-j
    tblnov+=10;
    tblnoc+=10;
    {
      static int ini=0;
      if(ini==0){
#ifdef MD_PRINT_WARN
	printf("tblno is changed to %d and %d in MR3calccoulomb_vdw_ij_ci_virial\n",tblnov,tblnoc);
#endif
	ini=1;
      }
    }
    unit->jsize=sizeof(VG_JVEC)*nj2max;
    unit->isize=sizeof(VG_IVEC)*ni2max;
    unit->fsize=sizeof(VG_FVEC)*ni2max;
    if(unit->pscalers==NULL){
      if((unit->pscalers=(VG_PSCALER *)MR3_malloc_pointer(sizeof(VG_PSCALER),"pscalers in MR3calccoulomb_vdw_ij_ci_virial"))==NULL){
	fprintf(stderr,"** error : can't malloc unit->pscalers in MR3calccoulomb_ij_ci_virial **\n");
	vg_exit(1);
      }
    }
    pscal=(VG_PSCALER *)unit->pscalers;
    if(ni2>blocksize*MD_CELLINDEX_MAXIBLOCKS){
      fprintf(stderr,"** error : too few MD_CELLINDEX_MAXIBLOCKS=%d **\n",MD_CELLINDEX_MAXIBLOCKS);
      MR3_exit(1);
    }
    
    for(icell=nicell=0;icell<ncelli;icell++){
      VG_CELL cijtmp[MD_NUM_JCELLS];
      int jcellnum;
      cijtmp[27].size=27;
      nii=cii[icell].size;
      ioffset=cii[icell].base;
      cidi[0]=icell % cminmaxi[0];
      cidi[1]=(icell/cminmaxi[0]) % cminmaxi[1];
      cidi[2]=(icell/cminmaxi[0]/cminmaxi[1]) % cminmaxi[2];
      //      printf("nicell=%d(%d,%d,%d) nii=%d\n",icell,cidi[0],cidi[1],cidi[2],nii);
      if(nii>0){
	for(xyz=0;xyz<3;xyz++) cidi[xyz]+=cmini[xyz];// convert to absolute cellindex
	for(jci=jcellnum=0;jci<27;jci++){
	  for(xyz=0;xyz<3;xyz++){
	    cidj[xyz]=cidi[xyz]+rclist[jci][xyz];
	    cidj[xyz]-=cminj[xyz];                   // convert to relative cellindex
	    //	  if(cidj[xyz]<0)          cidj[xyz]+=ldim[xyz];
	    //	  if(cidj[xyz]>=ldim[xyz]) cidj[xyz]-=ldim[xyz];
	  }
	  if(cidj[0]>=0 && cidj[0]<cminmaxj[0] &&
	     cidj[1]>=0 && cidj[1]<cminmaxj[1] &&
	     cidj[2]>=0 && cidj[2]<cminmaxj[2]){
	    jcell=(cidj[2]*cminmaxj[1]+cidj[1])*cminmaxj[0]+cidj[0];
	    if(cij[jcell].size>0){
	      cijtmp[jcellnum].base=cij[jcell].base;
	      cijtmp[jcellnum].size=cij[jcell].size;
              //              printf("  jci=%d jcellnum=%d cidj=(%d,%d,%d) base=%d size=%d\n",jci,jcellnum,cidj[0],cidj[1],cidj[2],cijtmp[jcellnum].base,cijtmp[jcellnum].size);
	      jcellnum++;
	    }
	  }
#if 1
	  else if(cellperiodic){ // periodic j-cell is used
#if 0
	    static int ini=0;
	    if(ini==0){
	      printf("*** periodic is assumed in _ci_virial **\n");
	      printf("*** last flag of cell_range() and assign_cell_noperiodic() must be set 0 for non-periodic **\n");
	      ini=1;
	    }
#endif
	    for(xyz=0;xyz<3;xyz++){
	      while(cidj[xyz]<0)              cidj[xyz]+=cminmaxj[xyz];
	      while(cidj[xyz]>=cminmaxj[xyz]) cidj[xyz]-=cminmaxj[xyz];
	    }
	    jcell=(cidj[2]*cminmaxj[1]+cidj[1])*cminmaxj[0]+cidj[0];
	    cijtmp[jcellnum].base=cij[jcell].base;
	    cijtmp[jcellnum].size=cij[jcell].size;
	    jcellnum++;
	  }
#endif
#if 0 // fill size-0 cells with jcellnum
	  else { // no j-cell exists
	    //	  printf("icell=%d jci=%d jcellnum=%d cidj=(%d,%d,%d) does not exist\n",icell,jci,jcellnum,cidj[0],cidj[1],cidj[2]);
	    cijtmp[jcellnum].base=0;
	    cijtmp[jcellnum].size=0;
	    jcellnum++;
	  }
#endif
	}
	//	if(jcellnum>0){
	if(1){
	  cijtmp[27].size=jcellnum;
          //          printf("icell=%d nii=%d cidi=(%d,%d,%d) jcellnum=%d\n",icell,nii,cidi[0],cidi[1],cidi[2],jcellnum);
	  for(ii=0;ii<nii;ii+=blocksize){
	    memcpy(pscal->ci[nicell],cijtmp,sizeof(VG_CELL)*MD_NUM_JCELLS);
#if VG_JDIV!=1
	    pscal->niblock[nicell]=(nii-ii<blocksize ? nii-ii:blocksize);
#endif
	    nicell++;
	    if(nicell>MD_CELLINDEX_MAXIBLOCKS){
	      fprintf(stderr,"** error : nicell=%d is too large **\n",nicell);
	      fprintf(stderr,"**         MD_CELLINDEX_MAXIBLOCKS=%d\n",
		      MD_CELLINDEX_MAXIBLOCKS);
	      vg_exit(1);
	    }
	  }
	}
      }
    }
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(34);
    vg_start_timer(35);
#endif

#if MD_PERIODIC_FIXED==1
    if(cellperiodic==0){
      static int ini=0;
      if(ini==0){
#ifdef MD_PRINT_WARN
	printf("** when cellperiodicflag is set 0,         **\n");
        printf("** xmax=volume[0] is not safe              **\n");
	printf("** so xmax*=4 in coulomb_vdw_ij_ci_virial  **\n");
#endif
	ini=1;
      }
      unit->xmax=volume[0]*4.0;
    }
    else{
      unit->xmax=volume[0];
    }
#endif
    tblnocv=(tblnoc % 10)+30;    
    //if(cellperiodic==0 && (tblnoc & 1)==1) tblnocv+=10;
    //    if(cellperiodic==0) tblnocv+=10;
    tblnocv+=10;
    //    printf("tblnocv=%d\n",tblnocv);
    
    vg_set_vectors_pos_charge_atype(unit,nj2,(double (*)[3])xj2,qj2,atypej2);
    vg_set_scalers_volume_alpha(unit,volume,rscale);
    vg_set_matrices_gscales_rscales(unit,nat,nat,gscale2,rscale2);
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(35);
    vg_start_timer(36);
#endif
    vg_set_pipeline_vectors_pos_charge_atype(unit,ni2,(double (*)[3])xi2,qi2,atypei2);
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(36);
    vg_start_timer(37);
#endif
#if 1
    vg_calculate_forces_potentials(unit,
				   tblnocv,
				   ni2,(double (*)[3])f2,&potc,&potv);
#endif
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(37);
    vg_start_timer(38);
#endif
  }
  if(unit->gpuoverlapflag==0 || unit->gpuoverlapflag==2){
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(38);
    vg_start_timer(39);
#endif
    for(i=0;i<ni_bak;i++){
      int cellindex;
      cellindex=cipi[i][3];
      ioffset=cii[cellindex].base+cinthpi[i];
      for(xyz=0;xyz<3;xyz++){
	force[i*3+xyz]+=f2[ioffset*3+xyz];
      }
#if 0
      printf("f2[%d]=(%e,%e,%e) is added to force[%d]=(%e,%e,%e)\n",
             ioffset,f2[ioffset*3],f2[ioffset*3+1],f2[ioffset*3+2],
             i,force[i*3],force[i*3+1],force[i*3+2]);
#endif
    }
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(39);
#endif

    MR3_free_pointer(cipi,"MR3calccoulomb_vdw_ij_ci_virial");
    MR3_free_pointer(opi,"MR3calccoulomb_vdw_ij_ci_virial");
    MR3_free_pointer(cii,"MR3calccoulomb_vdw_ij_ci_virial");
    MR3_free_pointer(cinthpi,"MR3calccoulomb_vdw_ij_ci_virial");
    MR3_free_pointer(cipj,"MR3calccoulomb_vdw_ij_ci_virial");
    MR3_free_pointer(opj,"MR3calccoulomb_vdw_ij_ci_virial");
    MR3_free_pointer(cij,"MR3calccoulomb_vdw_ij_ci_virial");
    MR3_free_pointer(cinthpj,"MR3calccoulomb_vdw_ij_ci_virial");
    MR3_free_pointer(xi2,"MR3calccoulomb_vdw_ij_ci_virial");
    MR3_free_pointer(xj2,"MR3calccoulomb_vdw_ij_ci_virial");
    MR3_free_pointer(qi2,"MR3calccoulomb_vdw_ij_ci_virial");
    MR3_free_pointer(qj2,"MR3calccoulomb_vdw_ij_ci_virial");
    MR3_free_pointer(atypei2,"MR3calccoulomb_vdw_ij_ci_virial");
    MR3_free_pointer(atypej2,"MR3calccoulomb_vdw_ij_ci_virial");
    MR3_free_pointer(f2,"MR3calccoulomb_vdw_ij_ci_virial");
    MR3_free_pointer(gscale2,"MR3calccoulomb_vdw_ij_ci_virial");
    MR3_free_pointer(rscale2,"MR3calccoulomb_vdw_ij_ci_virial");
  }
#ifdef MD_MEASURE_TIME
  vg_stop_and_accumulate_timer(30);
#endif
}


void MR3calccoulomb_vdw_ij_ci(int ni, double xi[], double qi[], int atypei[], double force[],
			      int nj, double xj[], double qj[], int atypej[],
			      int nat, double gscales[], double rscales[], double rscale,
			      int tblnoc,
			      double rcut, double skinnb,
			      double volume[3], int ldim[3])
{
#ifdef MD_MEASURE_TIME
  vg_start_timer(48);
#endif
#if 1 // work in 2009
  MR3calccoulomb_vdw_ij_ci_old_09(ni,xi,qi,atypei,force,
				  nj,xj,qj,atypej,
				  nat,gscales,rscales,rscale,
				  tblnoc,
				  rcut,skinnb,
				  volume,ldim);
#else
  double virial[3][3];
  int flag=1; // periodic is assumed
#ifdef MD_PRINT_WARN
  {
    static int ini=0;
    if(ini==0){
      printf("** warning : MR3calccoulomb_vdw_ij_ci_virial is used instead of MR3calccoulomb_vdw_ij_ci_old_09 **\n");
      ini=1;
    }
  }
#endif
  MR3calccoulomb_vdw_ij_ci_virial(ni,xi,qi,atypei,force,
				  nj,xj,qj,atypej,
				  nat,gscales,rscales,rscale,
				  tblnoc,
				  rcut,skinnb,
				  volume,ldim,
				  virial,flag);
#endif
#ifdef MD_MEASURE_TIME
  vg_stop_and_accumulate_timer(48);
#endif
}


void mr3calccoulomb_vdw_ij_ci_(int *ni, double xi[], double qi[], int atypei[], double force[],
                               int *nj, double xj[], double qj[], int atypej[],
                               int *nat, double gscales[], double rscales[], double *rscale,
                               int *tblnoc,
                               double *rcut, double *skinnb,
                               double volume[3], int ldim[3])
/*
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
                               int numex[], int natex[])*/
{
  int *atypei2,*atypej2;
  static int *atypej3=NULL;
  int i;

  if((atypei2=(int *)MR3_malloc_pointer(sizeof(int)*(*ni),"atypei2 in mr3calccoulomb_vdw_ij_ci_"))==NULL){
    fprintf(stderr,
	    "** error at malloc atypei2 in mr3calccoulomb_vdw_ij_ci_ **\n");
    MR3_exit(1);
  }
  if((atypej2=(int *)MR3_malloc_pointer(sizeof(int)*(*nj),"atypej2 in mr3calccoulomb_vdw_ij_ci_"))==NULL){
    fprintf(stderr,
	    "** error at malloc atypej2 in mr3calccoulomb_vdw_ij_ci_ **\n");
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
      atypej2[i]=0;
    }
  }

  MR3calccoulomb_vdw_ij_ci(*ni,xi,qi,atypei2,force, 
                           *nj,xj,qj,atypej2,
                           *nat,gscales,rscales,*rscale,
                           *tblnoc,
                           *rcut,*skinnb,
                           volume,ldim);
  MR3_free_pointer(atypei2,"atypei2 in mr3calccoulomb_vdw_ij_ci_");
  MR3_free_pointer(atypej2,"atypej2 in mr3calccoulomb_vdw_ij_ci_");
}


void MR3calccoulomb_vdw_ij_ci_old_09(int ni, double xi[], double qi[], int atypei[], double force[],
				     int nj, double xj[], double qj[], int atypej[],
				     int nat, double gscales[], double rscales[], double rscale,
				     int tblnoc,
				     double rcut, double skinnb,
				     double volume[3], int ldim[3])
{
#ifdef MD_CELLINDEX
  int tblnov=2;
  double xmin[3]={0.0,0.0,0.0};
  static int (*cipi)[4],*cinthpi;
  int ncell;
  static int (*cipj)[4],*cinthpj;
  int i,j,jcell,cidi[3],cidj[3],ioffset,joffset,cit,periodicflag=1;
  int nii,njj,icell,xyz,jci,ni2,nj2,ni2max,nj2max;
  int rclist[27][3]={{-1,-1,-1},{ 0,-1,-1},{ 1,-1,-1},
		     {-1, 0,-1},{ 0, 0,-1},{ 1, 0,-1},
		     {-1, 1,-1},{ 0, 1,-1},{ 1, 1,-1},
		     {-1,-1, 0},{ 0,-1, 0},{ 1,-1, 0},
		     {-1, 0, 0},{ 0, 0, 0},{ 1, 0, 0},
		     {-1, 1, 0},{ 0, 1, 0},{ 1, 1, 0},
		     {-1,-1, 1},{ 0,-1, 1},{ 1,-1, 1},
		     {-1, 0, 1},{ 0, 0, 1},{ 1, 0, 1},
		     {-1, 1, 1},{ 0, 1, 1},{ 1, 1, 1}};
  static double (*opi)[3],(*opj)[3],*xi2,*xj2,*f2,*qi2,*qj2;
  static int *atypei2,*atypej2;
  static M3_CELL *cii,*cij;
  static double *gscale2,*rscale2;
  int blocksize=VG_MINIMUM_PARTICLE_BLOCK_I;
  VG_UNIT *unit=m3_get_unit();
  double potc,potv,alpha3,xmax;
  int multiplyq=0,tblnocv;
  int nicell,ii;
  VG_PSCALER *pscal;
  static int ni_bak=0;

#ifdef MD_MEASURE_TIME
  vg_start_timer(60);
#endif
  //  nj=6;printf("** nj is changed to 6 **\n");
#ifndef MD_CELLINDEX
  fprintf(stderr,"** MR3calccoulomb_vdw_ij_ci is not supported **\n");
  exit(1);
#endif
#if !defined(MD_USE_QAUNION) || defined(MD_SIGMA_EPSION_IN_VGSTRUCT) || defined(MD_SORT_ATYPEI) || defined(MD_MATRIX_IN_SCALER)
  fprintf(stderr,"** error : MD_USE_QAUNION must be defined **\n");
  vg_exit(1);
#endif
  if(unit==NULL){
    fprintf(stderr,"** error : unit is NULL **\n");
    vg_exit(1);
  }

  if(unit->gpuoverlapflag==0 || unit->gpuoverlapflag==1){
#ifdef MD_MEASURE_TIME
    vg_start_timer(61);
#endif
    {
      static int ini=0;
      if(ini==0){
#ifdef MD_PRINT_WARN
	printf("*** MR3calccoulomb_vdw_ij_ci is called. periodicflag=%d\n",periodicflag);
#endif
	ini=1;
      }
    }
    ni_bak=ni;

    // set rcut2
    vg_set_rcut2(unit,rcut*rcut);

    // check rcut and skinnb
    check_rcut_and_skinnb(rcut,skinnb,ldim,volume);

    // assign cell
    ncell=ldim[0]*ldim[1]*ldim[2];
    if((cipi=(int (*)[4])MR3_malloc_pointer(sizeof(int)*4*ni,"MR3calcvdw_ij_ci"))==NULL){
      fprintf(stderr,"** error : can't malloc cipi in MR3calccoulomb_vdw_ij_ci **\n");
      MR3_exit(1);
    }
    if((opi=(double (*)[3])MR3_malloc_pointer(sizeof(double)*3*ni,"MR3calcvdw_ij_ci"))==NULL){
      fprintf(stderr,"** error : can't malloc opi in MR3calccoulomb_vdw_ij_ci **\n");
      MR3_exit(1);
    }
    if((cii=(M3_CELL *)MR3_malloc_pointer(sizeof(M3_CELL)*ncell,"MR3calcvdw_ij_ci"))==NULL){
      fprintf(stderr,"** error : can't malloc cii in MR3calccoulomb_vdw_ij_ci **\n");
      MR3_exit(1);
    }
    if((cinthpi=(int *)MR3_malloc_pointer(sizeof(int)*ni,"MR3calcvdw_ij_ci"))==NULL){
      fprintf(stderr,"** error : can't malloc cinthpi in MR3calccoulomb_ci **\n");
      MR3_exit(1);
    }
    if((cipj=(int (*)[4])MR3_malloc_pointer(sizeof(int)*4*nj,"MR3calcvdw_ij_ci"))==NULL){
      fprintf(stderr,"** error : can't malloc cipj in MR3calccoulomb_ci **\n");
      MR3_exit(1);
    }
    if((opj=(double (*)[3])MR3_malloc_pointer(sizeof(double)*3*nj,"MR3calcvdw_ij_ci"))==NULL){
      fprintf(stderr,"** error : can't malloc opj in MR3calccoulomb_ci **\n");
      MR3_exit(1);
    }
    if((cij=(M3_CELL *)MR3_malloc_pointer(sizeof(M3_CELL)*ncell,"MR3calcvdw_ij_ci"))==NULL){
      fprintf(stderr,"** error : can't malloc cij in MR3calccoulomb_ci **\n");
      MR3_exit(1);
    }
    if((cinthpj=(int *)MR3_malloc_pointer(sizeof(int)*nj,"MR3calcvdw_ij_ci"))==NULL){
      fprintf(stderr,"** error : can't malloc cinthpj in MR3calccoulomb_ci **\n");
      MR3_exit(1);
    }
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(61);
    vg_start_timer(62);
#endif
    assign_cell(ni,xi,xmin,volume,ldim,blocksize,cipi,opi,cii,cinthpi);
    assign_cell(nj,xj,xmin,volume,ldim,blocksize,cipj,opj,cij,cinthpj);
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(62);
    vg_start_timer(63);
#endif
    
    // copy to contiguous array
    for(i=ni2=0;i<ncell;i++) ni2+=(cii[i].size+blocksize-1)/blocksize*blocksize;
    for(i=nj2=0;i<ncell;i++) nj2+=(cij[i].size+blocksize-1)/blocksize*blocksize;
    ni2max=ni2+ncell*blocksize;
    nj2max=nj2+ncell*blocksize;
    if((xi2=(double *)MR3_malloc_pointer(sizeof(double)*ni2*3,"MR3calcvdw_ij_ci"))==NULL){
      fprintf(stderr,"** error : can't malloc xi2 in MR3calcvdw_ci **\n");
      MR3_exit(1);
    }
    if((xj2=(double *)MR3_malloc_pointer(sizeof(double)*nj2*3,"MR3calcvdw_ij_ci"))==NULL){
      fprintf(stderr,"** error : can't malloc xj2 in MR3calcvdw_ci **\n");
      MR3_exit(1);
    }
    if((qi2=(double *)MR3_malloc_pointer(sizeof(double)*ni2,"MR3calccoulomb_vdw_ij_ci"))==NULL){
      fprintf(stderr,"** error : can't malloc qi2 in MR3calccoulomb_ci **\n");
      MR3_exit(1);
    }
    if((qj2=(double *)MR3_malloc_pointer(sizeof(double)*nj2,"MR3calccoulomb_vdw_ij_ci"))==NULL){
      fprintf(stderr,"** error : can't malloc qj2 in MR3calccoulomb_ci **\n");
      MR3_exit(1);
    }
    if((atypei2=(int *)MR3_malloc_pointer(sizeof(int)*ni2,"MR3calcvdw_ij_ci"))==NULL){
      fprintf(stderr,"** error : can't malloc atypei2 in MR3calcvdw_ci **\n");
      MR3_exit(1);
    }
    if((atypej2=(int *)MR3_malloc_pointer(sizeof(int)*nj2,"MR3calcvdw_ij_ci"))==NULL){
      fprintf(stderr,"** error : can't malloc atypej2 in MR3calcvdw_ci **\n");
      MR3_exit(1);
    }
    if((f2=(double *)MR3_malloc_pointer(sizeof(double)*ni2*3,"MR3calcvdw_ij_ci"))==NULL){
      fprintf(stderr,"** error : can't malloc f2 in MR3calcvdw_ci **\n");
      MR3_exit(1);
    }
    bzero(xi2,sizeof(double)*ni2*3);
    bzero(xj2,sizeof(double)*nj2*3);
    bzero(qi2,sizeof(double)*ni2);
    bzero(qj2,sizeof(double)*nj2);
    bzero(atypei2,sizeof(int)*ni2);
    bzero(atypej2,sizeof(int)*nj2);
    bzero(f2,sizeof(double)*ni2*3);
    make_contiguous_xqa(ni,xi,qi,atypei,xmin,volume,cipi,cii,cinthpi,xi2,qi2,atypei2);
    make_contiguous_xqa(nj,xj,qj,atypej,xmin,volume,cipj,cij,cinthpj,xj2,qj2,atypej2);
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(63);
    vg_start_timer(64);
#endif
    
    // new gscale, rscale
    nat++;
    if((gscale2=(double *)MR3_malloc_pointer(sizeof(double)*nat*nat,"MR3calcvdw_ij_ci"))==NULL){
      fprintf(stderr,"** error : can't malloc gscales in MR3calcvdw_ij_ci **\n");
      MR3_exit(1);
    }
    if((rscale2=(double *)MR3_malloc_pointer(sizeof(double)*nat*nat,"MR3calcvdw_ij_ci"))==NULL){
      fprintf(stderr,"** error : can't malloc rscales in MR3calcvdw_ij_ci **\n");
      MR3_exit(1);
    }
    for(i=0;i<nat;i++){
      for(j=0;j<nat;j++){
	if(i==0 || j==0){
	  gscale2[i*nat+j]=0.0;
	  rscale2[i*nat+j]=1.0;
	}
	else{
	  gscale2[i*nat+j]=gscales[(i-1)*(nat-1)+j-1];
	  rscale2[i*nat+j]=rscales[(i-1)*(nat-1)+j-1];
	}
      }
    }
    
    // cellindex-i, cellindex-j
    tblnov+=10;
    tblnoc+=10;
    {
      static int ini=0;
      if(ini==0){
#ifdef MD_PRINT_WARN
	printf("tblno is changed to %d and %d in MR3calccoulomb_vdw_ij_ci\n",tblnov,tblnoc);
#endif
	ini=1;
      }
    }
    unit->jsize=sizeof(VG_JVEC)*nj2max;
    unit->isize=sizeof(VG_IVEC)*ni2max;
    unit->fsize=sizeof(VG_FVEC)*ni2max;
    if(unit->pscalers==NULL){
      if((unit->pscalers=(VG_PSCALER *)MR3_malloc_pointer(sizeof(VG_PSCALER),"MR3calccoulomb_vdw_ij_ci_old_09"))==NULL){
	fprintf(stderr,"** error : can't malloc unit->pscalers in MR3calccoulomb_vdw_ij_ci_old_09 **\n");
	vg_exit(1);
      }
    }
    pscal=(VG_PSCALER *)unit->pscalers;
    if(ni2>blocksize*MD_CELLINDEX_MAXIBLOCKS){
      fprintf(stderr,"** error : too few MD_CELLINDEX_MAXIBLOCKS=%d **\n",MD_CELLINDEX_MAXIBLOCKS);
      MR3_exit(1);
    }
    for(icell=nicell=0;icell<ncell;icell++){
      VG_CELL cijtmp[MD_NUM_JCELLS];
      cijtmp[27].size=27;
      nii=cii[icell].size;
      ioffset=cii[icell].base;
      cidi[0]=icell % ldim[0];
      cidi[1]=(icell/ldim[0]) % ldim[1];
      cidi[2]=(icell/ldim[0]/ldim[1]) % ldim[2];
      for(jci=0;jci<27;jci++){
	for(xyz=0;xyz<3;xyz++){
	  cidj[xyz]=cidi[xyz]+rclist[jci][xyz];
	  if(cidj[xyz]<0)          cidj[xyz]+=ldim[xyz];
	  if(cidj[xyz]>=ldim[xyz]) cidj[xyz]-=ldim[xyz];
	}
	jcell=(cidj[2]*ldim[1]+cidj[1])*ldim[0]+cidj[0];
	cijtmp[jci].base=cij[jcell].base;
	cijtmp[jci].size=cij[jcell].size;
      }
      for(ii=0;ii<nii;ii+=blocksize){
	memcpy(pscal->ci[nicell],cijtmp,sizeof(VG_CELL)*MD_NUM_JCELLS);
#if VG_JDIV!=1
	pscal->niblock[nicell]=(nii-ii<blocksize ? nii-ii:blocksize);
#endif
	nicell++;
	if(nicell>MD_CELLINDEX_MAXIBLOCKS){
	  fprintf(stderr,"** error : nicell is too large **\n");
	  vg_exit(1);
	}
      }
    }

    xmax=volume[0];
    if((periodicflag & 2)!=0) multiplyq=1;
    if((periodicflag & 1)==0) xmax*=2;
#if MD_PERIODIC_FIXED==1
    unit->xmax=xmax;
#endif
    tblnocv=(tblnoc % 10)+30;    
    
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(64);
    vg_start_timer(65);
#endif
    vg_set_vectors_pos_charge_atype(unit,nj2,(double (*)[3])xj2,qj2,atypej2);
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(65);
    vg_start_timer(66);
#endif
    vg_set_scalers_volume_alpha(unit,volume,rscale);
    vg_set_matrices_gscales_rscales(unit,nat,nat,gscale2,rscale2);
    vg_set_pipeline_vectors_pos_charge_atype(unit,ni2,(double (*)[3])xi2,qi2,atypei2);
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(66);
    vg_start_timer(67);
#endif
    vg_calculate_forces_potentials(unit,
				   tblnocv,
				   ni2,(double (*)[3])f2,&potc,&potv);
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(67);
    vg_start_timer(68);
#endif
  }
  if(unit->gpuoverlapflag==0 || unit->gpuoverlapflag==2){
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(68);
    vg_start_timer(69);
#endif
    for(i=0;i<ni_bak;i++){
      int cellindex;
      cellindex=cipi[i][3];
      ioffset=cii[cellindex].base+cinthpi[i];
      for(xyz=0;xyz<3;xyz++){
	force[i*3+xyz]+=f2[ioffset*3+xyz];
      }
    }

    MR3_free_pointer(cipi,"MR3calcvdw_ij_ci");
    MR3_free_pointer(opi,"MR3calcvdw_ij_ci");
    MR3_free_pointer(cii,"MR3calcvdw_ij_ci");
    MR3_free_pointer(cinthpi,"MR3calcvdw_ij_ci");
    MR3_free_pointer(cipj,"MR3calcvdw_ij_ci");
    MR3_free_pointer(opj,"MR3calcvdw_ij_ci");
    MR3_free_pointer(cij,"MR3calcvdw_ij_ci");
    MR3_free_pointer(cinthpj,"MR3calcvdw_ij_ci");
    MR3_free_pointer(xi2,"MR3calcvdw_ij_ci");
    MR3_free_pointer(xj2,"MR3calcvdw_ij_ci");
    MR3_free_pointer(qi2,"MR3calccoulomb_vdw_ij_ci");
    MR3_free_pointer(qj2,"MR3calccoulomb_vdw_ij_ci");
    MR3_free_pointer(atypei2,"MR3calcvdw_ij_ci");
    MR3_free_pointer(atypej2,"MR3calcvdw_ij_ci");
    MR3_free_pointer(f2,"MR3calcvdw_ij_ci");
    MR3_free_pointer(gscale2,"MR3calcvdw_ij_ci");
    MR3_free_pointer(rscale2,"MR3calcvdw_ij_ci");
#ifdef MD_MEASURE_TIME
    vg_stop_and_accumulate_timer(69);
#endif
  }
#else
  fprintf(stderr,"** MR3calccoulomb_vdw_ij_ci is not supported **\n");
  exit(1);
#endif
#ifdef MD_MEASURE_TIME
  vg_stop_and_accumulate_timer(60);
#endif
}


void MR3calccoulomb_vdw_ij_ci_old_090121(int ni, double xi[], double qi[], int atypei[], double force[],
					 int nj, double xj[], double qj[], int atypej[],
					 int nat, double gscales[], double rscales[], double rscale,
					 int tblnoc,
					 double rcut, double skinnb,
					 double volume[3], int ldim[3])
{
#ifdef MD_CELLINDEX
  int tblnov=2;
  double xmin[3]={0.0,0.0,0.0};
  int (*cipi)[4],*cinthpi,ncell;
  int (*cipj)[4],*cinthpj;
  int i,j,jcell,cidi[3],cidj[3],ioffset,joffset,cit,periodicflag=1;
  int nii,njj,icell,xyz,jci,ni2,nj2,ni2max,nj2max;
  int rclist[27][3]={{-1,-1,-1},{ 0,-1,-1},{ 1,-1,-1},
		     {-1, 0,-1},{ 0, 0,-1},{ 1, 0,-1},
		     {-1, 1,-1},{ 0, 1,-1},{ 1, 1,-1},
		     {-1,-1, 0},{ 0,-1, 0},{ 1,-1, 0},
		     {-1, 0, 0},{ 0, 0, 0},{ 1, 0, 0},
		     {-1, 1, 0},{ 0, 1, 0},{ 1, 1, 0},
		     {-1,-1, 1},{ 0,-1, 1},{ 1,-1, 1},
		     {-1, 0, 1},{ 0, 0, 1},{ 1, 0, 1},
		     {-1, 1, 1},{ 0, 1, 1},{ 1, 1, 1}};
  double (*opi)[3],(*opj)[3],*xi2,*xj2,*f2,*f3=NULL,*qi2,*qj2;
  int *atypei2,*atypej2;
  M3_CELL *cii,*cij;
  double *gscale2,*rscale2;
  int blocksize=VG_MINIMUM_PARTICLE_BLOCK_I;
  VG_UNIT *unit=m3_get_unit();

  vg_set_rcut2(unit,rcut*rcut);
#if !defined(MD_USE_QAUNION) || defined(MD_SIGMA_EPSION_IN_VGSTRUCT) || defined(MD_SORT_ATYPEI) || defined(MD_MATRIX_IN_SCALER)
  fprintf(stderr,"** error : MD_USE_QAUNION must be defined **\n");
  vg_exit(1);
#endif
  {
    static int ini=0;
    if(ini==0){
#ifdef MD_PRINT_WARN
      printf("*** MR3calccoulomb_vdw_ij_ci is called. periodicflag=%d\n",periodicflag);
#endif
      ini=1;
    }
  }

  // check rcut and skinnb
  if((rcut+skinnb)*ldim[0]>volume[0] || 
     (rcut+skinnb)*ldim[1]>volume[1] ||
     (rcut+skinnb)*ldim[2]>volume[2]){
    fprintf(stderr,"** error : rcut(%f)+skinnb(%f) must be smaller than a cell size (%f,%f,%f), ldim=(%d,%d,%d) **\n",
	    rcut,skinnb,volume[0]/ldim[0],volume[1]/ldim[1],volume[2]/ldim[2],
	    ldim[0],ldim[1],ldim[2]);
    MR3_exit(1);
  }

  // assign cell
  ncell=ldim[0]*ldim[1]*ldim[2];
  if((cipi=(int (*)[4])MR3_malloc_pointer(sizeof(int)*4*ni,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc cipi in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((opi=(double (*)[3])MR3_malloc_pointer(sizeof(double)*3*ni,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc opi in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((cii=(M3_CELL *)MR3_malloc_pointer(sizeof(M3_CELL)*ncell,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc cii in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((cinthpi=(int *)MR3_malloc_pointer(sizeof(int)*ni,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc cinthpi in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((cipj=(int (*)[4])MR3_malloc_pointer(sizeof(int)*4*nj,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc cipj in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((opj=(double (*)[3])MR3_malloc_pointer(sizeof(double)*3*nj,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc opj in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((cij=(M3_CELL *)MR3_malloc_pointer(sizeof(M3_CELL)*ncell,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc cij in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((cinthpj=(int *)MR3_malloc_pointer(sizeof(int)*nj,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc cinthpj in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  assign_cell(ni,xi,xmin,volume,ldim,blocksize,cipi,opi,cii,cinthpi);
  assign_cell(nj,xj,xmin,volume,ldim,blocksize,cipj,opj,cij,cinthpj);

  // copy to contiguous array
  for(i=ni2=0;i<ncell;i++) ni2+=(cii[i].size+blocksize-1)/blocksize*blocksize;
  for(i=nj2=0;i<ncell;i++) nj2+=(cij[i].size+blocksize-1)/blocksize*blocksize;
  ni2max=ni2+ncell*blocksize;
  nj2max=nj2+ncell*blocksize;
  if((xi2=(double *)MR3_malloc_pointer(sizeof(double)*ni2*3,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc xi2 in MR3calcvdw_ci **\n");
    MR3_exit(1);
  }
  if((xj2=(double *)MR3_malloc_pointer(sizeof(double)*nj2*3,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc xj2 in MR3calcvdw_ci **\n");
    MR3_exit(1);
  }
  if((qi2=(double *)MR3_malloc_pointer(sizeof(double)*ni2,"MR3calccoulomb_vdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc qi2 in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((qj2=(double *)MR3_malloc_pointer(sizeof(double)*nj2,"MR3calccoulomb_vdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc qj2 in MR3calccoulomb_ci **\n");
    MR3_exit(1);
  }
  if((atypei2=(int *)MR3_malloc_pointer(sizeof(int)*ni2,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc atypei2 in MR3calcvdw_ci **\n");
    MR3_exit(1);
  }
  if((atypej2=(int *)MR3_malloc_pointer(sizeof(int)*nj2,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc atypej2 in MR3calcvdw_ci **\n");
    MR3_exit(1);
  }
  if((f2=(double *)MR3_malloc_pointer(sizeof(double)*ni2*3,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc f2 in MR3calcvdw_ci **\n");
    MR3_exit(1);
  }
  bzero(xi2,sizeof(double)*ni2*3);
  bzero(xj2,sizeof(double)*nj2*3);
  bzero(qi2,sizeof(double)*ni2);
  bzero(qj2,sizeof(double)*nj2);
  bzero(atypei2,sizeof(int)*ni2);
  bzero(atypej2,sizeof(int)*nj2);
  bzero(f2,sizeof(double)*ni2*3);
  make_contiguous_xqa(ni,xi,qi,atypei,xmin,volume,cipi,cii,cinthpi,xi2,qi2,atypei2);
  make_contiguous_xqa(nj,xj,qj,atypej,xmin,volume,cipj,cij,cinthpj,xj2,qj2,atypej2);

  // new gscale, rscale
  nat++;
  if((gscale2=(double *)MR3_malloc_pointer(sizeof(double)*nat*nat,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc gscales in MR3calcvdw_ij_ci **\n");
    MR3_exit(1);
  }
  if((rscale2=(double *)MR3_malloc_pointer(sizeof(double)*nat*nat,"MR3calcvdw_ij_ci"))==NULL){
    fprintf(stderr,"** error : can't malloc rscales in MR3calcvdw_ij_ci **\n");
    MR3_exit(1);
  }
  for(i=0;i<nat;i++){
    for(j=0;j<nat;j++){
      if(i==0 || j==0){
	gscale2[i*nat+j]=0.0;
	rscale2[i*nat+j]=1.0;
      }
      else{
	gscale2[i*nat+j]=gscales[(i-1)*(nat-1)+j-1];
	rscale2[i*nat+j]=rscales[(i-1)*(nat-1)+j-1];
      }
    }
  }

  // cellindex-i, cellindex-j
  tblnov+=10;
  tblnoc+=10;
  {
    static int ini=0;
    if(ini==0){
#ifdef MD_PRINT_WARN
      printf("tblno is changed to %d and %d in MR3calccoulomb_vdw_ij_ci\n",tblnov,tblnoc);
#endif
      ini=1;
    }
  }
  unit->jsize=sizeof(VG_JVEC)*nj2max;
  unit->isize=sizeof(VG_IVEC)*ni2max;
  unit->fsize=sizeof(VG_FVEC)*ni2max;
  {
    int nicell,ii;
    VG_PSCALER *pscal;

    if(unit->pscalers==NULL){
      if((unit->pscalers=(VG_PSCALER *)MR3_malloc_pointer(sizeof(VG_PSCALER),"pscalers in MR3calccoulomb_vdw_ij_ci_old_090121"))==NULL){
	fprintf(stderr,"** error : can't malloc unit->pscalers in MR3calccoulomb_vdw_ij_ci_old_090121 **\n");
	vg_exit(1);
      }
    }
    pscal=(VG_PSCALER *)unit->pscalers;
    if(ni2>blocksize*MD_CELLINDEX_MAXIBLOCKS){
      fprintf(stderr,"** error : too few MD_CELLINDEX_MAXIBLOCKS=%d **\n",MD_CELLINDEX_MAXIBLOCKS);
      MR3_exit(1);
    }
    for(icell=nicell=0;icell<ncell;icell++){
      VG_CELL cijtmp[MD_NUM_JCELLS];
      cijtmp[27].size=27;
      nii=cii[icell].size;
      ioffset=cii[icell].base;
      cidi[0]=icell % ldim[0];
      cidi[1]=(icell/ldim[0]) % ldim[1];
      cidi[2]=(icell/ldim[0]/ldim[1]) % ldim[2];
      for(jci=0;jci<27;jci++){
	for(xyz=0;xyz<3;xyz++){
	  cidj[xyz]=cidi[xyz]+rclist[jci][xyz];
	  if(cidj[xyz]<0)          cidj[xyz]+=ldim[xyz];
	  if(cidj[xyz]>=ldim[xyz]) cidj[xyz]-=ldim[xyz];
	}
	jcell=(cidj[2]*ldim[1]+cidj[1])*ldim[0]+cidj[0];
	cijtmp[jci].base=cij[jcell].base;
	cijtmp[jci].size=cij[jcell].size;
      }
      for(ii=0;ii<nii;ii+=blocksize){
	memcpy(pscal->ci[nicell],cijtmp,sizeof(VG_CELL)*MD_NUM_JCELLS);
	nicell++;
	if(nicell>MD_CELLINDEX_MAXIBLOCKS){
	  fprintf(stderr,"** error : nicell is too large **\n");
	  vg_exit(1);
	}
      }
    }
  }

#if 0
  MR3calccoulomb_ij(ni2,xi2,qi2,f2,nj2,xj2,qj2,rscale,tblnoc,volume[0],periodicflag);
  MR3calcvdw_ij_core(ni2,xi2,atypei2,f2,nj2,xj2,atypej2,nat,gscale2,rscale2,tblnov,volume[0],periodicflag);
#else
  {
    double potc,potv,alpha3,xmax=volume[0],(*ftmp)[3];
    int multiplyq=0,tblno=tblnoc,tblnocv=(tblnoc % 10)+30;
    if((periodicflag & 2)!=0) multiplyq=1;
    if((periodicflag & 1)==0) xmax*=2;
#if MD_PERIODIC_FIXED==1
    unit->xmax=xmax;
#endif
    
    vg_set_vectors_pos_charge_atype(unit,nj2,(double (*)[3])xj2,qj2,atypej2);
    vg_set_scalers_volume_alpha(unit,volume,rscale);
    vg_set_matrices_gscales_rscales(unit,nat,nat,gscale2,rscale2);
    vg_set_pipeline_vectors_pos_charge_atype(unit,ni2,(double (*)[3])xi2,qi2,atypei2);
#if 1 // work for real and vdw
    vg_calculate_forces_potentials(unit,
				   //				   tblno,
				   //36,
				   tblnocv,
				   ni2,(double (*)[3])f2,&potc,&potv);
    //    for(i=0;i<ni2;i++) printf("tblno=%d f2[%d]=%e %e %e\n",tblnocv,i,f2[i*3],f2[i*3+1],f2[i*3+2]);
#endif
#if 0 // work for real
    vg_calculate_forces_potentials(unit,
				   tblno,
				   ni2,(double (*)[3])f2,&potc,&potv);
    if(tblno>=100) tblno=(tblno-100) % 10;
    switch(tblno){
    case 0:
    case 1:
      if(multiplyq){
	for(i=0;i<ni2;i++) for(j=0;j<3;j++) f2[i*3+j]*=qi2[i]*Coulomb_vdw_factor;
      }
      break;
    case 6:
    case 16:
      alpha3=rscale*rscale*rscale;
      for(i=0;i<ni2;i++) for(j=0;j<3;j++) f2[i*3+j]*=qi2[i]*alpha3*Coulomb_vdw_factor;
      break;
    case 7:
    case 17:
      for(i=0;i<ni2;i++) f2[i*3]*=qi2[i]*rscale*Coulomb_vdw_factor;
      break;
    default:
      fprintf(stderr,"** error : not supported mode=%d **\n",tblno);
      vg_exit(1);
      break;
    }
    //    for(i=0;i<3;i++) printf("OK  f2[%d]=%e %e %e\n",i,f2[i*3],f2[i*3+1],f2[i*3+2]);
#endif
#if 0 // work for vdw
    if((ftmp=(double (*)[3])MR3_malloc_pointer(sizeof(double)*ni2*3,"MR3calcvdw_ij_core"))==NULL){
      fprintf(stderr,"** error : can't malloc ftmp in MR3calcvdw **\n");
      vg_exit(1);
    }
    
    vg_calculate_forces_potentials(unit,
				   tblnov,
				   //				   36,
				   ni2,ftmp,&potc,&potv);
    for(i=0;i<ni2;i++) for(j=0;j<3;j++) f2[i*3+j]+=ftmp[i][j];
    
    MR3_free_pointer(ftmp,"MR3calcvdw_ij_core");
#endif
  }
#endif


  for(i=0;i<ni;i++){
    int cellindex;
    cellindex=cipi[i][3];
    ioffset=cii[cellindex].base+cinthpi[i];
    //    printf("ioffset=%d i=%d f2=%e %e %e\n",ioffset,i,f2[ioffset*3],f2[ioffset*3+1],f2[ioffset*3+2]);
    for(xyz=0;xyz<3;xyz++){
      force[i*3+xyz]+=f2[ioffset*3+xyz];
    }
  }

  MR3_free_pointer(cipi,"MR3calcvdw_ij_ci");
  MR3_free_pointer(opi,"MR3calcvdw_ij_ci");
  MR3_free_pointer(cii,"MR3calcvdw_ij_ci");
  MR3_free_pointer(cinthpi,"MR3calcvdw_ij_ci");
  MR3_free_pointer(cipj,"MR3calcvdw_ij_ci");
  MR3_free_pointer(opj,"MR3calcvdw_ij_ci");
  MR3_free_pointer(cij,"MR3calcvdw_ij_ci");
  MR3_free_pointer(cinthpj,"MR3calcvdw_ij_ci");
  MR3_free_pointer(xi2,"MR3calcvdw_ij_ci");
  MR3_free_pointer(xj2,"MR3calcvdw_ij_ci");
  MR3_free_pointer(qi2,"MR3calccoulomb_vdw_ij_ci");
  MR3_free_pointer(qj2,"MR3calccoulomb_vdw_ij_ci");
  MR3_free_pointer(atypei2,"MR3calcvdw_ij_ci");
  MR3_free_pointer(atypej2,"MR3calcvdw_ij_ci");
  MR3_free_pointer(f2,"MR3calcvdw_ij_ci");
  if(f3!=NULL) MR3_free_pointer(f3,"MR3calcvdw_ij_ci");
  MR3_free_pointer(gscale2,"MR3calcvdw_ij_ci");
  MR3_free_pointer(rscale2,"MR3calcvdw_ij_ci");
#else
  fprintf(stderr,"** MR3calccoulomb_ij_ci_old_090121 is not supported **\n");
  exit(1);
#endif

}


void MR3calccoulomb_vdw_ij_ci_exlist(int ni, double xi[], double qi[], int atypei[], double force[],
                                     int nj, double xj[], double qj[], int atypej[],
                                     int nat, double gscales[], double rscales[], double rscale,
                                     int tblno,
                                     double rcut, double skinnb,
                                     double volume[3], int ldim[3],
                                     int numex[], int natex[],
                                     int flag)
{
  /*
    flag bit 0 -- 0 : non periodic cell (currently ni>nj makes it default)
                  1 : periodic cell
         bit 1 -- 0 : natex does not include duplicate list
                  1 : natex includes duplicate list
         bit 6 -- 0 : non-overlap
                  1 : overlap
    */
#if !defined(MD_CELLINDEX) || !defined(MD_QAUNION_ATYPEBIT) || \
  !defined(__USE_COULOMBVDW)
  fprintf(stderr,"** this routine is not supported : compile again **\n");
  MR3_exit(1);
#else
  double *ftmp,*forcep;
  int i,j,multiplyq=1,overlapflag=0;
  VG_UNIT *unit=m3_get_unit();
  //  double volume[3]={0.0,0.0,0.0};
  double stress[3][3],volume2[3],sc,sv;

#ifdef MD_MEASURE_TIME
  vg_start_timer(41);
#endif
  if((flag & 1)==0){
    for(i=0;i<3;i++) volume2[i]=volume[i]*4.0;
    {
      static int ini=0;
      if(ini==0){
#ifdef MD_PRINT_WARN
        printf("** volume2 is set 4 times larger than volume **\n");
#endif
        ini=1;
      }
    }
  }
  else{
    for(i=0;i<3;i++) volume2[i]=volume[i];
  }
#if 0
  {
    int offset=0;
    for(i=0;i<ni;i++){
      printf("numex[%d]=%d natex=",i,numex[i]);
      for(j=0;j<numex[i];j++){
        printf("%d ",natex[offset+j]);
      }
      printf("\n");
      offset+=numex[i];
    }
  }
#endif  
  if(ldim[0]<=1 && ldim[1]<=1 && ldim[2]<=1){
    fprintf(stderr,"** ldim<=1 is not allowed in MR3calccoulomb_vdw_ij_ci_exlist (ldim=%d %d %d)**\n",ldim[0],ldim[1],ldim[2]);
    MR3_exit(1);
  }
  if((flag & (1<<6))!=0){
#ifdef MD_PRINT_WARN
    printf("overlapflag=1 is set in MR3calccoulomb_vdw_ij_ci_exlist\n");
#endif
    overlapflag=1;
  }

  if((ftmp=(double *)MR3_malloc_pointer(sizeof(double)*ni*3,"MR3calccoulomb_vdw_ij_ci_exlist"))==NULL){
    fprintf(stderr,"** error : can't malloc ftmp **\n");
    vg_exit(1);
  }

  if(tblno==6 || tblno==7) multiplyq=0;
  if(overlapflag!=0){
    unit->gpuoverlapflag=1;
    unit->ni_overlap=ni;
  }
  else{
    unit->gpuoverlapflag=0;
  }
#ifdef MD_MEASURE_TIME
  vg_start_timer(42);
#endif
  bzero(ftmp,sizeof(double)*ni*3);
  MR3calccoulomb_vdw_ij_ci_virial(ni,xi,qi,atypei,ftmp,
                                  nj,xj,qj,atypej,
                                  nat,gscales,rscales,rscale,
                                  tblno,rcut,skinnb,volume,ldim,stress,flag & 1);
#ifdef MD_MEASURE_TIME
  vg_stop_and_accumulate_timer(42);
  vg_start_timer(43);
#endif
  //  for(i=0;i<ni;i++) printf("grape3 after ci_virial : ftmp[%d]=%e %e %e\n",i,ftmp[i*3],ftmp[i*3+1],ftmp[i*3+2]);
  if(tblno & 1){ // potential calculation
    if(overlapflag!=0){
      fprintf(stderr,"** overlap is now allowed for potential calculation **\n");
      MR3_exit(1);
    }

    //    for(i=0;i<10;i++) printf("grape3 after pot ci_virial : ftmp[%d]=%e %e %e\n",i,ftmp[i*3],ftmp[i*3+1],ftmp[i*3+2]);
    //    for(i=0;i<10;i++) printf("after vdw sum   : ftmp[%d]=%e %e %e\n",i,ftmp[i*3],ftmp[i*3+1],ftmp[i*3+2]);
    MR3calccoulomb_vdw_nlist_ij_emu(ni,xi,qi,atypei,ftmp,nj,xj,qj,atypej,
                                    nat,gscales,rscales,rscale,
                                    tblno,volume2,(flag & 3)+(1-multiplyq)*4,
                                    numex,natex,-1.0);
    //    for(i=0;i<10;i++) printf("after summing up cv emu1 : ftmp[%d]=%e %e %e\n",i,ftmp[i*3],ftmp[i*3+1],ftmp[i*3+2]);
    //    for(i=0,sc=sv=0.0;i<ni;i++){ if(i<10) printf("after summing up cv emu2 : ftmp[%d]=%e %e %e\n",i,ftmp[i*3],ftmp[i*3+1],ftmp[i*3+2]);sc+=ftmp[i*3];sv+=ftmp[i*3+1];}
    //    printf("sc=%e sv=%e\n",sc,sv);
  }
  else{
    MR3calccoulomb_vdw_nlist_ij_emu(ni,xi,qi,atypei,force,nj,xj,qj,atypej,
                                    nat,gscales,rscales,rscale,
                                    tblno,volume2,(flag & 3)+(1-multiplyq)*4,
                                    numex,natex,-1.0);
  }
#ifdef MD_MEASURE_TIME
  vg_stop_and_accumulate_timer(43);
  vg_start_timer(44);
#endif
  if(overlapflag==0 || (overlapflag==1 && unit->gpuoverlapflag==0)){
    for(i=0;i<ni*3;i++) force[i]+=ftmp[i];
  }
  else{
    unit->gpuoverlapflag=2;
  }
  MR3_free_pointer(ftmp,"MR3calccoulomb_vdw_ij_ci_exlist");
#endif
#ifdef MD_MEASURE_TIME
  vg_stop_and_accumulate_timer(41);
#endif
}


void mr3calccoulomb_vdw_ij_ci_exlist_(int *ni, double xi[], double qi[], int atypei[], double force[],
                                      int *nj, double xj[], double qj[], int atypej[],
                                      int *nat, double gscales[], double rscales[], double *rscale,
                                      int *tblnoc,
                                      double *rcut, double *skinnb,
                                      double volume[3], int ldim[3],
                                      int numex[], int natex[], int *flag)
{
  int *atypei2,*atypej2;
  int *natex2,s;
  static int *atypej3=NULL;
  int i;

#ifdef MD_MEASURE_TIME
  vg_start_timer(40);
#endif
  for(i=s=0;i<*ni;i++){
    s+=numex[i];
  }
  if((atypei2=(int *)MR3_malloc_pointer(sizeof(int)*(*ni),"atypei2 in mr3calccoulomb_vdw_ij_ci_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc atypei2 in mr3calccoulomb_vdw_ij_ci_exlist_ **\n");
    MR3_exit(1);
  }
  if((atypej2=(int *)MR3_malloc_pointer(sizeof(int)*(*nj),"atypej2 in mr3calccoulomb_vdw_ij_ci_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc atypej2 in mr3calccoulomb_vdw_ij_ci_exlist_ **\n");
    MR3_exit(1);
  }
  if((natex2=(int *)MR3_malloc_pointer(sizeof(int)*s,"natex2 in mr3calccoulomb_vdw_ij_ci_exlist_"))==NULL){
    fprintf(stderr,
	    "** error at malloc natex2 in mr3calccoulomb_ij_ci_exlist_ **\n");
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
      atypej2[i]=0;
    }
  }
  for(i=0;i<s;i++) natex2[i]=natex[i]-1;

  MR3calccoulomb_vdw_ij_ci_exlist(*ni,xi,qi,atypei2,force, 
                                  *nj,xj,qj,atypej2,
                                  *nat,gscales,rscales,*rscale,
                                  *tblnoc,
                                  *rcut,*skinnb,
                                  volume,ldim,
                                  numex,natex2,*flag);
  MR3_free_pointer(atypei2,"atypei2 in mr3calccoulomb_vdw_ij_ci_exlist_");
  MR3_free_pointer(natex2,"natex2 in mr3calccoulomb_vdw_ij_ci_exlist_");
  MR3_free_pointer(atypej2,"atypej2 in mr3calccoulomb_vdw_ij_ci_exlist_");
#ifdef MD_MEASURE_TIME
  vg_stop_and_accumulate_timer(40);
#endif
}






