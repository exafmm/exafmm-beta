#pragma once
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>

class Dataset {
public:
  std::vector<kvec4> pos;
  Dataset(unsigned long n,
	  const char *filename = "plummer.dat") : pos(n) {
    std::fstream file;
    file.open(filename,std::ios::in);
    if (!file.fail()) {
      unsigned long ntmp;
      file.read((char *)&ntmp, sizeof(unsigned long));
      if (n == ntmp) {
	file.read((char *)&pos[0], n*sizeof(double4));
	return;
      }
    }
    file.close();
    unsigned long i = 0;
#if 0
    const float scale = 3.0 * M_PI / 16.0;
    while (i < n) {
      float R = 1.0 / sqrt( (pow(drand48(), -2.0 / 3.0) - 1.0) );
      if (R < 100.0) {
	float Z = (1.0 - 2.0 * drand48()) * R;
        float theta = 2.0 * M_PI * drand48();
	float X = sqrt(R * R - Z * Z) * cos(theta);
	float Y = sqrt(R * R - Z * Z) * sin(theta);
	X *= scale; Y *= scale; Z *= scale;
	pos[i][0] = X;
	pos[i][1] = Y;
	pos[i][2] = Z;
	pos[i][3] = drand48() / n;
	ldiv_t tmp_i = ldiv(i, n/33);
	if(tmp_i.rem == 0) {
	  printf(".");
	  fflush(stdout);
	}
	i++;
      }
    }
#else
    for (i=0; i<n; i++) {
      pos[i][0] = drand48();
      pos[i][1] = drand48();
      pos[i][2] = drand48();
      pos[i][3] = drand48() / n;
    }
#endif
    kvec4 com(0.0);
    for (i=0; i<n; i++) {
      com[0] += abs(pos[i][3]) * pos[i][0];
      com[1] += abs(pos[i][3]) * pos[i][1];
      com[2] += abs(pos[i][3]) * pos[i][2];
      com[3] += abs(pos[i][3]);
    }
    com[0] /= com[3];
    com[1] /= com[3];
    com[2] /= com[3];

    for(i=0; i<n; i++) {
      pos[i][0] -= com[0];
      pos[i][1] -= com[1];
      pos[i][2] -= com[2];
    }
    printf("\n");
    file.open(filename,std::ios::out);
    file.write((char *)&n, sizeof(unsigned long));
    file.write((char *)&pos[0], n*sizeof(double4));
    file.close();
  }
};
