#define NMAX   65536
#define KMAX   16384
#define NTHRE  128
#define NLOAD  32
#define ATYPE  64
#define ATYPE2 (ATYPE * ATYPE)

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

typedef struct {
  float gscale;
  float rscale;
} VG_MATRIX;

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

void MR3calccoulomb_ij(int ni, double xi[], double qi[], double force[],
                       int nj, double xj[], double qj[],
                       double rscale,
                       int tblno, double xmax, int periodicflag);

void MR3calcvdw_ij(int ni, double xi[], int atypei[], double force[],
                   int nj, double xj[], int atypej[],
                   int nat, double gscale[], double rscale[],
                   int tblno, double xmax, int periodicflag);

void MR3calcewald(int *k, int knum_org, double *x, int n, double *q,
                  double alpha, double epsilon, double cell[3][3],
                  double *force, double *tpot, double stress[3][3]);

int get_knum(double ksize);

void init_kvec(double ksize, int *kvec);
