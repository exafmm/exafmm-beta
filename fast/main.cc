//-----------------------------------------------------------------------------+
// Copyright (C) 2000-2004, 2008, 2010 Walter Dehnen                           |
// GNU General Public License                                                  |
//-----------------------------------------------------------------------------+
#include <tree.h>
#include <fstream>

int main(int argc, const char* argv[])
{
#ifdef MANY
  for ( int it=0; it<25; it++ ) {
#else
  for ( int it=8; it<9; it++ ) {
#endif
  unsigned N=int(pow(10,(it+24)/8.0));
  double   d,dx,dy,dz,tic,toc,tree,approx;
  if(argc>1) N = atoi(argv[1]);
  vect     *ACC = new vect[N];
  real     *POT = new real[N];

  bodies body(N);
  for(body=0; body!=body.end(); ++body) {
    body.pos()[0]   = rand()/(1.+double(RAND_MAX));
    body.pos()[1]   = rand()/(1.+double(RAND_MAX));
    body.pos()[2]   = rand()/(1.+double(RAND_MAX));
    body.scal()     = 1.0/real(N);
  }

  tic = get_time();
  TreeBuilder TB(body);
  TB.build();
  Evaluator *FMM = new Evaluator(body,TB.RAD,TB.LEVEL,TB.NLEAF,TB.NCELL);
  TB.link(FMM->C0,FMM->L0);
  toc = get_time();
  tree = toc-tic;
  tic = get_time();
  FMM->approximate();
  toc = get_time();
  approx = toc-tic;

  std::cout<<std::setprecision(3);
#ifdef DIRECT
  for(body=0; body!=body.end(); ++body) {
    ACC[body.index()] = body.acc();
    POT[body.index()] = body.pot();
  }
  tic = get_time();
  FMM->exact();
  toc = get_time();
  std::cout<<"direct time      = "<< toc-tic <<"\n";
#endif
  real dacc=0.,dnac=0., dpot=0., dnpt=0.;
  for(body=0; body!=body.end(); ++body) {
    vect& Acc(ACC[body.index()]);
    real& Pot(POT[body.index()]);
//    std::cout << body << " " << Pot << " " << body.pot() << std::endl;
    d  = Pot    - body.pot();
    dx = Acc[0] - body.acc()[0];
    dy = Acc[1] - body.acc()[1];
    dz = Acc[2] - body.acc()[2];
    dpot += d*d;
    dnpt += Pot*Pot;
    dacc += dx*dx+dy*dy+dz*dz;
    dnac += norm(Acc);
  }
#ifdef DIRECT
  std::cout<<"FMM time         = "<< tree+approx <<"\n";
  std::cout<<"Acc error        = "<< std::sqrt(dacc/dnac) <<"\n";
  std::cout<<"Pot error        = "<< std::sqrt(dpot/dnpt) <<"\n";
#else
  std::cout<< tree+approx <<std::endl;
#endif
  std::ofstream file("time",std::ios::out | std::ios::app);
  file << tree+approx << std::endl;
  file.close();
  delete[] ACC;
  delete[] POT;
  delete FMM;
  }
}
