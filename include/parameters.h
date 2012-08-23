#ifndef parameters_h
#define parameters_h

// Compile-time parameters
static const int P = 3;                                         //!< Order of expansions

// Runtime parameters
struct Parameters {
  int NCRIT;                                                    //!< Number of bodies per leaf cell
  int NSPAWN;                                                   //!< Threshold of NDLEAF for spawning new threads
  int IMAGES;                                                   //!< Number of periodic image sublevels
  float THETA;                                                  //!< Multipole acceptance criteria
  float EPS2;                                                   //!< Softening parameter (squared)
  Parameters() : NCRIT(10), NSPAWN(1000), IMAGES(0), THETA(.6), EPS2(.0) {}
};

#endif
