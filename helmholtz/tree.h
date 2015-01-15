void getBounds(vec3 * Xj, int numBodies, vec3 & X0, real_t & R0) {
  vec3 Xmin = Xj[1];
  vec3 Xmax = Xj[1];
  for (int i=0; i<numBodies; i++) {
    Xmin = min(Xj[i],Xmin);
    Xmax = max(Xj[i],Xmax);
  }
  real_t diameter = 0;
  for (int d=0; d<3; d++) {
    real_t dX = Xmax[d] - Xmin[d];
    diameter = dX > diameter ? dX : diameter;
    X0[d] = (Xmax[d] + Xmin[d]) * 0.5;
  }
  R0 = diameter * 0.5;
}
