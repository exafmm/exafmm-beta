#include "types.h"
#include "vtk.h"

int main()
{
  int const Nbody=5000;
  Bodies bodies(Nbody);
  B_iter B;

  for( B=bodies.begin(); B!=bodies.end(); ++B ) {               // Loop over all bodies
    for( int d=0; d!=3; ++d )                                   //  Loop over each dimension
      B->pos[d] = rand()/(1.+RAND_MAX)*2-1;                     //   Initialize positions
    real r = sqrt(B->pos[0]*B->pos[0]+B->pos[1]*B->pos[1]+B->pos[2]*B->pos[2]);
    for( int d=0; d!=3; ++d )                                   //  Loop over each dimension
      B->pos[d] /= r;                                           //   Normalize positions
    B->scal = 1./bodies.size();                                 //  Initialize source value
  }                                                             // End loop over all bodies

  real r0(0);                                                   // Root radius
  vect xmin,xmax,x0(0);                                         // Min,Max,Center of domain
  B = bodies.begin();                                           // Reset bodies counter
  xmin = xmax = B->pos;                                         // Initialize xmin,xmax
  for( B=bodies.begin(); B!=bodies.end(); ++B ) {               // Loop over all bodies
    for( int d=0; d!=3; ++d ) {                                 //  Loop over each dimension
      if     (B->pos[d] < xmin[d]) xmin[d] = B->pos[d];         //   Determine xmin
      else if(B->pos[d] < xmax[d]) xmax[d] = B->pos[d];         //   Determine xmax
    }                                                           //  End loop over each dimension
    x0 += B->pos;                                               //  Sum positions
  }                                                             // End loop over all bodies
  x0 /= bodies.size();                                          // Calculate average position
  for( int d=0; d!=3; ++d ) {                                   // Loop over each dimension
    x0[d] = int(x0[d]+0.5);                                     //  Shift center to nearest integer
    r0 = std::max(xmax[d]-x0[d],r0);                            //  Calculate max distance from center
    r0 = std::max(x0[d]-xmin[d],r0);                            //  Calculate max distance from center
  }                                                             // End loop over each dimension
  r0 = pow(2.0,int(1.0+log(r0)/M_LN2));                         // Add some leeway to root radius

  vtkPlot vtk;
  vtk.setDomain(r0,x0);
  vtk.setGroup(0,bodies.size()/2);
  for( B=bodies.begin(); B!=bodies.begin()+bodies.size()/2; ++B ) {
    vtk.setPoints(0,B->pos);
  }
  vtk.setGroup(1,bodies.size()/2);
  for( B=bodies.begin()+bodies.size()/2; B!=bodies.end(); ++B ) {
    vtk.setPoints(1,B->pos);
  }
  vtk.plot(2);
}
