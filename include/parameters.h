#ifndef parameters_h
#define parameters_h

// Compile-time parameters
const int P = 3;                                                //!< Order of expansions
const float EPS2 = .0;                                          //!< Softening parameter (squared)

// Runtime parameters
struct Parameters {
  int NCRIT;                                                    //!< Number of bodies per leaf cell
  int NSPAWN;                                                   //!< Threshold of NDLEAF for spawning new threads
  int IMAGES;                                                   //!< Number of periodic image sublevels
  float THETA;                                                  //!< Multipole acceptance criteria
  Parameters() : NCRIT(10), NSPAWN(1000), IMAGES(0), THETA(.6) {}
};

#endif
