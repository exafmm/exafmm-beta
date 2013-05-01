#ifndef matern_h
#define matern_h
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include "types.h"
using boost::math::cyl_bessel_k;
using boost::math::tgamma;

const real_t NU = 1.5;                                        //!< Order of modified Bessel function of the second kind
const real_t SIGMA = 10;                                      //!< Scaling of Matern covariance function

#endif
