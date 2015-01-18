void getAnm(real_t Anm1[][P+1], real_t Anm2[][P+1]) {
  Anm1[0][0] = 1;
  Anm2[0][0] = 1;
  for (int m=0; m<=P; m++) {
    if(m>0) Anm1[m][m] = sqrt((2 * m - 1.0) / (2 * m));
    if(m<P) Anm1[m+1][m] = sqrt(2 * m + 1.0);
    for (int n=m+2; n<=P; n++) {
      Anm1[n][m] = 2 * n - 1;
      Anm2[n][m] = sqrt((n + m - 1.0) * (n - m - 1.0));
      Anm1[n][m] /= sqrt(real_t(n - m) * (n + m));
      Anm2[n][m] /= sqrt(real_t(n - m) * (n + m));
    }
  }
}
