void getCoef(double *C, const double *dX, double invR2, const double invR) {
  C[0] = invR;
  invR2 = -invR2;
  double x = dX[0], y = dX[1], z = dX[2];
  double invR3 = invR * invR2;
  C[1] = x * invR3;
  C[2] = y * invR3;
  C[3] = z * invR3;

  double invR5 = 3 * invR3 * invR2;
  double t = x * invR5;
  C[4] = x * t + invR3;
  C[5] = y * t;
  C[6] = z * t;
  t = y * invR5;
  C[7] = y * t + invR3;
  C[8] = z * t;
  C[9] = z * z * invR5 + invR3;

  double invR7 = 5 * invR5 * invR2;
  t = x * x * invR7;
  C[10] = x * (t + 3 * invR5);
  C[11] = y * (t +     invR5);
  C[12] = z * (t +     invR5);
  t = y * y * invR7;
  C[13] = x * (t +     invR5);
  C[16] = y * (t + 3 * invR5);
  C[17] = z * (t +     invR5);
  t = z * z * invR7;
  C[15] = x * (t +     invR5);
  C[18] = y * (t +     invR5);
  C[19] = z * (t + 3 * invR5);
  C[14] = x * y * z * invR7;
}

void M2LSum(double *L, const double *C, const double*M) {
  for (int l=0; l<LTERM; l++) L[l] += M[0] * C[l];

  L[0] += M[1]*C[1]+M[2]*C[2]+M[3]*C[3];
  L[1] += M[1]*C[4]+M[2]*C[5]+M[3]*C[6];
  L[2] += M[1]*C[5]+M[2]*C[7]+M[3]*C[8];
  L[3] += M[1]*C[6]+M[2]*C[8]+M[3]*C[9];

  for (int l=4; l<10; l++) L[0] += M[l] * C[l];
  L[1] += M[4]*C[10]+M[5]*C[11]+M[6]*C[12]+M[7]*C[13]+M[8]*C[14]+M[9]*C[15];
  L[2] += M[4]*C[11]+M[5]*C[13]+M[6]*C[14]+M[7]*C[16]+M[8]*C[17]+M[9]*C[18];
  L[3] += M[4]*C[12]+M[5]*C[14]+M[6]*C[15]+M[7]*C[17]+M[8]*C[18]+M[9]*C[19];
  L[4] += M[1]*C[10]+M[2]*C[11]+M[3]*C[12];
  L[5] += M[1]*C[11]+M[2]*C[13]+M[3]*C[14];
  L[6] += M[1]*C[12]+M[2]*C[14]+M[3]*C[15];
  L[7] += M[1]*C[13]+M[2]*C[16]+M[3]*C[17];
  L[8] += M[1]*C[14]+M[2]*C[17]+M[3]*C[18];
  L[9] += M[1]*C[15]+M[2]*C[18]+M[3]*C[19];
}

void powerM(double *C, const double *dX) {
  C[1] = C[0] * dX[0];
  C[2] = C[0] * dX[1];
  C[3] = C[0] * dX[2];

  C[4] = C[1] * dX[0] / 2;
  C[5] = C[2] * dX[0];
  C[6] = C[3] * dX[0];
  C[7] = C[2] * dX[1] / 2;
  C[8] = C[3] * dX[1];
  C[9] = C[3] * dX[2] / 2;
}

void powerL(double *C, const double *dX) {
  C[1] = C[0] * dX[0];
  C[2] = C[0] * dX[1];
  C[3] = C[0] * dX[2];

  C[4] = C[1] * dX[0] / 2;
  C[5] = C[2] * dX[0];
  C[6] = C[3] * dX[0];
  C[7] = C[2] * dX[1] / 2;
  C[8] = C[3] * dX[1];
  C[9] = C[3] * dX[2] / 2;
}

void M2MSum(double *MI, const double *C, const double *MJ) {
  for (int i=1; i<MTERM; i++) MI[i] += MJ[i];

  MI[4] += C[1]*MJ[1];
  MI[5] += C[1]*MJ[2]+C[2]*MJ[1];
  MI[6] += C[1]*MJ[3]+C[3]*MJ[1];
  MI[7] += C[2]*MJ[2];
  MI[8] += C[2]*MJ[3]+C[3]*MJ[2];
  MI[9] += C[3]*MJ[3];
}

void L2LSum(double *LI, const double *C, const double *LJ) {
  LI[1] += C[1]*LJ[4]+C[2]*LJ[5]+C[3]*LJ[6];
  LI[2] += C[1]*LJ[5]+C[2]*LJ[7]+C[3]*LJ[8];
  LI[3] += C[1]*LJ[6]+C[2]*LJ[8]+C[3]*LJ[9];

  LI[1] += C[4]*LJ[10]+C[5]*LJ[11]+C[6]*LJ[12]+C[7]*LJ[13]+C[8]*LJ[14]+C[9]*LJ[15];
  LI[2] += C[4]*LJ[11]+C[5]*LJ[13]+C[6]*LJ[14]+C[7]*LJ[16]+C[8]*LJ[17]+C[9]*LJ[18];
  LI[3] += C[4]*LJ[12]+C[5]*LJ[14]+C[6]*LJ[15]+C[7]*LJ[17]+C[8]*LJ[18]+C[9]*LJ[19];
  LI[4] += C[1]*LJ[10]+C[2]*LJ[11]+C[3]*LJ[12];
  LI[5] += C[1]*LJ[11]+C[2]*LJ[13]+C[3]*LJ[14];
  LI[6] += C[1]*LJ[12]+C[2]*LJ[14]+C[3]*LJ[15];
  LI[7] += C[1]*LJ[13]+C[2]*LJ[16]+C[3]*LJ[17];
  LI[8] += C[1]*LJ[14]+C[2]*LJ[17]+C[3]*LJ[18];
  LI[9] += C[1]*LJ[15]+C[2]*LJ[18]+C[3]*LJ[19];
}

void L2PSum(float *TRG, const double *C, const double *L) {
  TRG[1] += C[1]*L[4]+C[2]*L[5]+C[3]*L[6];
  TRG[2] += C[1]*L[5]+C[2]*L[7]+C[3]*L[8];
  TRG[3] += C[1]*L[6]+C[2]*L[8]+C[3]*L[9];

  TRG[1] += C[4]*L[10]+C[5]*L[11]+C[6]*L[12]+C[7]*L[13]+C[8]*L[14]+C[9]*L[15];
  TRG[2] += C[4]*L[11]+C[5]*L[13]+C[6]*L[14]+C[7]*L[16]+C[8]*L[17]+C[9]*L[18];
  TRG[3] += C[4]*L[12]+C[5]*L[14]+C[6]*L[15]+C[7]*L[17]+C[8]*L[18]+C[9]*L[19];
}
