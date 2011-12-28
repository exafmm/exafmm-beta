#define NMAX      65536
#define KMAX      16384
#define NTHRE      128
#define NLOAD       32
#define ATYPE            64
#define ATYPE2         (ATYPE * ATYPE)

#define MD_REAL_R2MIN 0.0001f
#define MD_REAL_R2MAX 10.0f
#define MD_LJ_R2MIN   0.25
#define MD_LJ_R2MAX   64.0

typedef union {
  float q;
  int atype;
} VG_QATYPE;

typedef struct {
  int r[3];
  VG_QATYPE qatype;
} VG_XVEC;

typedef struct {
  float k[3];
  float factor1;
} VG_KVEC;

typedef union {
  int i;
  float f;
} FI;

typedef union {
  struct{
    FI fi0;
    FI fi1;
  } fi2;
  double d;
} DI2;

#ifdef MR3_MALLOC
#define MR3_malloc_pointer(x,y) MR3_my_malloc2(x,y)
#define MR3_free_pointer(x,y)   MR3_my_free2((void **)(&(x)),y)
#else
#define MR3_malloc_pointer(x,y) malloc(x)
#define MR3_free_pointer(x,y)   free(x)
#endif

#define MR3_exit(x) exit(x)

void MR3init(void);

void MR3free(void);

void MR3SetTable(char *filename, int tblno, int flag);

void MR3calccoulomb_ij(int ni, double xi[], double qi[], double force[],
                       int nj, double xj[], double qj[],
                       double rscale,
                       int tblno, double xmax, int periodicflag);

void MR3calcewald(int *k, int knum_org, double *x, int n, double *q,
                  double alpha, double epsilon, double cell[3][3],
                  double *force, double *tpot, double stress[3][3]);

/* mr3_host.c */
void MR3calccoulomb_ij_host(int ni, double xi[], double qi[], double force[], int nj, double xj[], double qj[], double rscale, int tblno, double xmax, int periodicflag);

void MR3calccoulomb_host(double x[], int n, double q[], double rscale, int tblno, double xmax, int periodicflag, int natchangeflag, double force[]);

void MR3calcewald_dft_host(int k[], int knum, double x[], int n, double chg[], double cellsize[3], double bs[], double bc[]);

void MR3calcewald_idft_eng_host(int k[], double bs[], double bc[], int knum, double x[], int n, double cellsize[3], double force[]);

void MR3calcewald_idft_force_host(int k[], double bs[], double bc[], int knum, double x[], int n, double cellsize[3], double force[]);

void MR3calcewald_host(int *k, int knum_org, double *x, int n, double *chg, double alpha, double epsilon, double cell[3][3], double *force, double *tpot, double stress[3][3]);

int get_knum(double ksize);

void init_kvec(double ksize, int *kvec);
