#define kernel(func) \
  int i,j,k,at,njj,jj;						\
  float dr[3],rr2,qj,vol[3],vol1[3],ftmp,alpha2,dn2,dn6,sqdn;\
  float api2=(float)(2.e0/sqrt(M_PI));\
  VG_JVEC *jvec;\
  \
  alpha2=(scalers->alpha)*(scalers->alpha);\
  for(k=0;k<3;k++){\
    vol[k]=scalers->volume[k];\
    vol1[k]=1.0f/vol[k];\
  }\
  for(i=0;i<nii;i++){\
    for(j=0;j<nj;j+=VG_MINIMUM_PARTICLE_BLOCK_J){\
      njj=(j+VG_MINIMUM_PARTICLE_BLOCK_J<nj) ? \
	VG_MINIMUM_PARTICLE_BLOCK_J:nj-j;     \
      jvec=jvectors+j;\
      for(jj=0;jj<njj;jj++){\
        rr2=0.0f;\
        for(k=0;k<3;k++){\
	  dr[k]=ivec[i].r[k]-jvec[jj].r[k];\
	  dr[k]-=rintf(dr[k]*vol1[k])*vol[k];\
	  rr2+=dr[k]*dr[k];\
        }\
        if(rr2!=0.0){\
	  func;\
        }\
      }\
    }\
  }

#define mother(kernel_func) \
  int i,j,natj,nii,nj;	    \
  VG_UNIT *unit;\
  VG_IVEC *ivec;\
  VG_JVEC *jvectors;\
  VG_MATRIX *mat;\
  VG_FVEC *fvec;\
  VG_SCALER *scal;\
  VG_PSCALER *pscal;\
  \
  unit=(VG_UNIT *)unit_org;\
  mat=unit->matrices;\
  scal=unit->scalers;\
  pscal=unit->pscalers;\
  natj=unit->natj;\
  nj=unit->nj;\
  jvectors=unit->jvectors;\
  for(i=0;i<unit->ni;i+=VG_MINIMUM_PARTICLE_BLOCK_I){\
    nii=(i+VG_MINIMUM_PARTICLE_BLOCK_I<(unit->ni) ? \
	 VG_MINIMUM_PARTICLE_BLOCK_I:(unit->ni)-i);\
    ivec=(VG_IVEC *)(unit->ivectors)+i;\
    fvec=(VG_FVEC *)(unit->fvectors)+i;\
    kernel_func(nii,nj,natj,ivec,jvectors,mat,scal,fvec,pscal);	\
  }


static void grav_kernel(int nii, int nj, int natj, 
			VG_IVEC *ivec, VG_JVEC *jvectors, 
			VG_MATRIX *matrix, VG_SCALER *scalers,
			VG_FVEC *fvec, VG_PSCALER *pscal)
{
  kernel(grav);
}


static void gravpot_kernel(int nii, int nj, int natj, 
			   VG_IVEC *ivec, VG_JVEC *jvectors, 
			   VG_MATRIX *matrix, VG_SCALER *scalers,
			   VG_FVEC *fvec, VG_PSCALER *pscal)
{
  kernel(gravpot);
}


static void lj_kernel(int nii, int nj, int natj, 
		      VG_IVEC *ivec, VG_JVEC *jvectors, 
		      VG_MATRIX *matrix, VG_SCALER *scalers,
		      VG_FVEC *fvec, VG_PSCALER *pscal)
{
  kernel(lj);
}


static void ljpot_kernel(int nii, int nj, int natj, 
			 VG_IVEC *ivec, VG_JVEC *jvectors,
			 VG_MATRIX *matrix, VG_SCALER *scalers,
			 VG_FVEC *fvec, VG_PSCALER *pscal)
{
  kernel(ljpot);
}


static void real_kernel(int nii, int nj, int natj, 
			VG_IVEC *ivec, VG_JVEC *jvectors,
			VG_MATRIX *matrix, VG_SCALER *scalers,
			VG_FVEC *fvec, VG_PSCALER *pscal)
{
  kernel(real);
}


static void realpot_kernel(int nii, int nj, int natj, 
			   VG_IVEC *ivec, VG_JVEC *jvectors,
			   VG_MATRIX *matrix, VG_SCALER *scalers,
			   VG_FVEC *fvec, VG_PSCALER *pscal)
{
  kernel(realpot);
}


static void wave_kernel(int nii, int nj, int natj, 
			VG_IVEC *ivec, VG_JVEC *jvectors,
			VG_MATRIX *matrix, VG_SCALER *scalers,
			VG_FVEC *fvec, VG_PSCALER *pscal)
{
  kernel(wave);
}


static void wavepot_kernel(int nii, int nj, int natj, 
			   VG_IVEC *ivec, VG_JVEC *jvectors,
			   VG_MATRIX *matrix, VG_SCALER *scalers,
			   VG_FVEC *fvec, VG_PSCALER *pscal)
{
  kernel(wavepot);
}


static void md_grav(void *unit_org)
{
  mother(grav_kernel);
}


static void md_gravpot(void *unit_org)
{
  mother(gravpot_kernel);
}


static void md_lj(void *unit_org)
{
  mother(lj_kernel);
}


static void md_ljpot(void *unit_org)
{
  mother(ljpot_kernel);
}


static void md_real(void *unit_org)
{
  mother(real_kernel);
}


static void md_realpot(void *unit_org)
{
  mother(realpot_kernel);
}


static void md_wave(void *unit_org)
{
  mother(wave_kernel);
}


static void md_wavepot(void *unit_org)
{
  mother(wavepot_kernel);
}


