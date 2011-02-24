#include <iostream>
#include <assert.h>

int main() {
  int p = 20, nc = 0;
  for( int n=0; n!=p; ++n ) {
    for( int nx=n; nx>=0; --nx ) {
      for( int ny=n-nx; ny>=0; --ny ) {
        int nz = n-nx-ny;
        int i = n*(n+1)*(n+2)/6+(n-nx)*(n-nx+1)/2+n-nx-ny;
        assert( i == nc );
        nc++;
      }
    }
  }

  int pc = 0;
  for( int n=0; n!=p; ++n ) {
    for( int px=n; px>=0; --px ) {
      for( int py=n-px; py>=0; --py ) {
        for( int pz=n-px-py; pz>=0; --pz ) {
          int i = n*(n+1)*(n+2)*(n+3)/24+(n-px)*(n-px+1)*(n-px+2)/6+(n-px-py)*(n-px-py+1)/2+n-px-py-pz;
          assert( i == pc );
          pc++;
        }
      }
    }
  }
  std::cout << p << " " << nc << " " << pc << std::endl;

}
