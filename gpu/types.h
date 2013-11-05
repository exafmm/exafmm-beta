#ifndef _TYPES_H_
#define _TYPES_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cuda_runtime.h>
#include <fstream>
#include <iostream>
#include <stdint.h>
#include <string>
#include <sys/time.h>
#include <vector>
#include "cudavec.h"
#include "macros.h"
#include "vec.h"

typedef float real;

const int  P     = 3;
const real EPS   = 1e-6;
const real EPS2  = 0.0001;
const real THETA = .75;

const int MTERM = P*(P+1)*(P+2)/6;

#endif
