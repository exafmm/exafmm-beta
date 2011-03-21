#include "dataset.h"
#include "construct.h"
#include "kernel.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  assert( IMAGES == 0 );
  int numBodies = 1000;
  Bodies bodies(numBodies);
  Cells cells;
  Dataset D;
  Evaluator E;
  TreeConstructor T;
  std::ofstream file;
  file.open("data",std::ofstream::binary);

  for( int it=0; it!=25; ++it ) {
    numBodies = int(pow(10,(it+24)/8.0));
    std::cout << "N             : " << numBodies << std::endl;
    bodies.resize(numBodies);
    D.random(bodies,1,1);
    T.startTimer("FMM          ");
    T.setDomain(bodies);
#ifdef TOPDOWN
    T.topdown(bodies,cells);
#else
    T.bottomup(bodies,cells);
#endif

    T.downward(cells,cells,1);
    T.stopTimer("FMM          ",true);
    T.eraseTimer("FMM          ");

    T.startTimer("Direct sum   ");
    T.buffer = bodies;
#if 1
    for( B_iter B=T.buffer.begin(); B!=T.buffer.end(); ++B ) {
      B->pot = -B->Q / std::sqrt(EPS2);
      B->acc = 0;
    }
    E.evalP2P(T.buffer,T.buffer);
    for( B_iter B=T.buffer.begin(); B!=T.buffer.end(); ++B ) {
      file << B->pot << std::endl;
      file << B->acc[0] << std::endl;
      file << B->acc[1] << std::endl;
      file << B->acc[2] << std::endl;
    }
#else
    for( B_iter B=T.buffer.begin(); B!=T.buffer.end(); ++B ) {
      file >> B->pot;
      file >> B->acc[0];
      file >> B->acc[1];
      file >> B->acc[2];
    }
#endif
    T.stopTimer("Direct sum   ");
    T.printAllTime();
    T.eraseTimer("Direct sum   ");
    T.writeTime();
    T.resetTimer();

    B_iter B  = bodies.begin();
    B_iter B2 = T.buffer.begin();
    real potDiff = 0, potNorm = 0, accDiff = 0, accNorm = 0;
    for( int i=0; i!=numBodies; ++i,++B,++B2 ) {
      B->pot -= B->Q / std::sqrt(EPS2);
#ifdef DEBUG
      std::cout << B->I << " " << B->acc[0] << " " << B2->acc[0] << std::endl;
#endif
      potDiff += (B->pot - B2->pot) * (B->pot - B2->pot);
      potNorm += B2->pot * B2->pot;
      accDiff += norm(B->acc - B2->acc);
      accNorm += norm(B2->acc);
    }
    std::cout << "Error (pot)   : " << std::sqrt(potDiff/potNorm) << std::endl;
    std::cout << "Error (acc)   : " << std::sqrt(accDiff/accNorm) << std::endl;
    cells.clear();
  }
  file.close();
#ifdef VTK
  int Ncell = 0;
  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroupOfPoints(bodies,Ncell);
  vtk.plot(Ncell);
#endif
}
