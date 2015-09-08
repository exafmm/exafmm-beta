#ifndef verify_h
#define verify_h
#include "logger.h"

#if Laplace
typedef double target_t;                                        // Target type is double
#elif Helmholtz
typedef std::complex<double> target_t;                          // Target type is complex double
#endif 

//! Verify results
class Verify {
public:
  //! Get sum of scalar component of a vector of target bodies
  target_t getSumScalar(Bodies & bodies) {
    target_t v = 0;                                             // Initialize difference
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      v += B->TRG[0] * B->SRC;                                  //  Sum of scalar component
    }                                                           // End loop over bodies
    return v;                                                   // Return difference
  }

  //! Get norm of scalar component of a vector of target bodies
  target_t getNrmScalar(Bodies & bodies) {
    target_t v = 0;                                             // Initialize norm
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      v += B->TRG[0] * B->TRG[0];                               //  Norm of scalar component
    }                                                           // End loop over bodies
    return v;                                                   // Return norm
  }

  //! Get difference between scalar component of two vectors of target bodies
  target_t getDifScalar(Bodies & bodies, Bodies & bodies2) {
    target_t v = 0;                                             // Initialize difference
    B_iter B2 = bodies2.begin();                                // Set iterator of bodies2
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++, B2++) { // Loop over bodies & bodies2
      v += (B->TRG[0] - B2->TRG[0]) * (B->TRG[0] - B2->TRG[0]); //  Difference of scalar component
    }                                                           // End loop over bodies & bodies2
    return v;                                                   // Return difference
  }

  //! Get difference between scalar component of two vectors of target bodies
  target_t getRelScalar(Bodies & bodies, Bodies & bodies2) {
    target_t v = 0;                                             // Initialize difference
    B_iter B2 = bodies2.begin();                                // Set iterator of bodies2
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++, B2++) { // Loop over bodies & bodies2
      v += ((B->TRG[0] - B2->TRG[0]) * (B->TRG[0] - B2->TRG[0]))
	/ (B2->TRG[0] * B2->TRG[0]);                            //  Difference of scalar component
    }                                                           // End loop over bodies & bodies2
    return v;                                                   // Return difference
  }

  //! Get norm of scalar component of a vector of target bodies
  target_t getNrmVector(Bodies & bodies) {
    target_t v = 0;                                             // Initialize norm
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      v += B->TRG[1] * B->TRG[1]                                //  Norm of vector x component
	+  B->TRG[2] * B->TRG[2]                                //  Norm of vector y component
	+  B->TRG[3] * B->TRG[3];                               //  Norm of vector z component
    }                                                           // End loop over bodies
    return v;                                                   // Return norm
  }

  //! Get difference between scalar component of two vectors of target bodies
  target_t getDifVector(Bodies & bodies, Bodies & bodies2) {
    target_t v = 0;                                             // Initialize difference
    B_iter B2 = bodies2.begin();                                // Set iterator of bodies2
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++, B2++) { // Loop over bodies & bodies2
      v += (B->TRG[1] - B2->TRG[1]) * (B->TRG[1] - B2->TRG[1])  //  Difference of vector x component
	+  (B->TRG[2] - B2->TRG[2]) * (B->TRG[2] - B2->TRG[2])  //  Difference of vector y component
	+  (B->TRG[3] - B2->TRG[3]) * (B->TRG[3] - B2->TRG[3]); //  Difference of vector z component
    }                                                           // End loop over bodies & bodies2
    return v;                                                   // Return difference
  }

  //! Get difference between scalar component of two vectors of target bodies
  target_t getRelVector(Bodies & bodies, Bodies & bodies2) {
    target_t v = 0;                                             // Initialize difference
    B_iter B2 = bodies2.begin();                                // Set iterator of bodies2
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++, B2++) { // Loop over bodies & bodies2
      v += ((B->TRG[1] - B2->TRG[1]) * (B->TRG[1] - B2->TRG[1]) +//  Difference of vector x component
	    (B->TRG[2] - B2->TRG[2]) * (B->TRG[2] - B2->TRG[2]) +//  Difference of vector y component
	    (B->TRG[3] - B2->TRG[3]) * (B->TRG[3] - B2->TRG[3]))//  Difference of vector z component
	/ (B2->TRG[1] * B2->TRG[1] +                            //  Norm of vector x component
	   B2->TRG[2] * B2->TRG[2] +                            //  Norm of vector y component
	   B2->TRG[3] * B2->TRG[3]);                            //  Norm of vector z component
    }                                                           // End loop over bodies & bodies2
    return v;                                                   // Return difference
  }

  //! Print relative L2 norm scalar error
  void print(std::string title, target_t v) {
    if (logger::verbose) {                                      // If verbose flag is true
      std::cout << std::setw(logger::stringLength) << std::left //  Set format
		<< title << " : " << std::setprecision(logger::decimal) << std::scientific // Set title
                << v << std::endl;                              //  Print potential error
    }                                                           // End if for verbose flag
  }
};

#endif
