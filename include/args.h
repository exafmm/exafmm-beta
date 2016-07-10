#ifndef args_h
#define args_h
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <stdint.h>
#ifndef _SX
#include <getopt.h>
#endif

namespace exafmm {
#ifndef _SX
  static struct option long_options[] = {
    {"ncrit",        required_argument, 0, 'c'},
    {"cutoff",       required_argument, 0, 'C'},
    {"distribution", required_argument, 0, 'd'},
    {"dual",         no_argument,       0, 'D'},
    {"graft",        no_argument,       0, 'g'},
    {"getMatrix",    no_argument,       0, 'G'},
    {"help",         no_argument,       0, 'h'},
    {"images",       required_argument, 0, 'i'},
    {"IneJ",         no_argument,       0, 'j'},
    {"mutual",       no_argument,       0, 'm'},
    {"numBodies",    required_argument, 0, 'n'},
    {"useRopt",      no_argument,       0, 'o'},
    {"PP",           required_argument, 0, 'P'},
    {"repeat",       required_argument, 0, 'r'},
    {"nspawn",       required_argument, 0, 's'},
    {"theta",        required_argument, 0, 't'},
    {"threads",      required_argument, 0, 'T'},
    {"verbose",      no_argument,       0, 'v'},
    {"write",        no_argument,       0, 'w'},
    {"useRmax",      no_argument,       0, 'x'},
    {0, 0, 0, 0}
  };
#endif

  class Args {
  public:
    int ncrit;
    double cutoff;
    const char * distribution;
    int dual;
    int graft;
    int getMatrix;
    int images;
    int IneJ;
    int mutual;
    int numBodies;
    int useRopt;
    int PP;
    int repeat;
    int nspawn;
    double theta;
    int threads;
    int verbose;
    int write;
    int useRmax;

  private:
    void usage(char * name) {
      fprintf(stderr,
	      "Usage: %s [options]\n"
	      "Long option (short option)     : Description (Default value)\n"
	      " --ncrit (-c)                  : Number of bodies per leaf cell (%d)\n"
	      " --cutoff (-C)                 : Cutoff distance of interaction (%f)\n"
	      " --distribution (-d) [l/c/s/p] : lattice, cube, sphere, octant, plummer (%s)\n"
	      " --dual (-D)                   : Use dual tree traversal (%d)\n"
	      " --graft (-g)                  : Graft remote trees to global tree (%d)\n"
	      " --getMatrix (-G)              : Write G matrix to file (%d)\n"
	      " --help (-h)                   : Show this help document\n"
	      " --images (-i)                 : Number of periodic image levels (%d)\n"
	      " --IneJ (-j)                   : Use different sources & targets (%d)\n"
	      " --mutual (-m)                 : Use mutual interaction (%d)\n"
	      " --numBodies (-n)              : Number of bodies (%d)\n"
	      " --useRopt (-o)                : Use error optimized theta for MAC (%d)\n"
	      " --P (-P) not working          : Order of expansion (%d)\n"
	      " --repeat (-r)                 : Number of iteration loops (%d)\n"
	      " --nspawn (-s)                 : Threshold for stopping task creation during recursion (%d)\n"
	      " --theta (-t)                  : Multipole acceptance criterion (%f)\n"
	      " --threads (-T)                : Number of threads (%d)\n"
	      " --verbose (-v)                : Print information to screen (%d)\n"
	      " --write (-w)                  : Write timings to file (%d)\n"
	      " --useRmax (-x)                : Use maximum distance for MAC (%d)\n",
	      name,
	      ncrit,
	      cutoff,
	      distribution,
	      dual,
	      graft,
	      getMatrix,
	      images,
	      IneJ,
	      mutual,
	      numBodies,
	      useRopt,
	      PP,
	      repeat,
	      nspawn,
	      theta,
	      threads,
	      verbose,
	      write,
	      useRmax);
    }

    const char * parse(const char * arg) {
      switch (arg[0]) {
      case 'c':
	return "cube";
      case 'l':
	return "lattice";
      case 'o':
	return "octant";
      case 'p':
	return "plummer";
      case 's':
	return "sphere";
      default:
	fprintf(stderr, "invalid distribution %s\n", arg);
	abort();
      }
      return "";
    }

    uint64_t getDistNum(const char * distribution) {
      if (distribution == "cube")
        return 0;
      else if (distribution == "lattice")
        return 1;
      else if (distribution == "octant")
        return 2;
      else if (distribution == "plummer")
        return 3;
      else if (distribution == "sphere")
        return 4;
      else {
        fprintf(stderr, "invalid distribution %s\n", distribution);
        abort();
      }
    }

    uint64_t getKernelNum() {
      uint64_t key = 0;
#if EXAFMM_LAPLACE
      key |= 0;
#elif EXAFMM_HELMHOLTZ
      key |= 1;
#elif EXAFMM_STOKES
      key |= 2;
#elif EXAFMM_BIOTSAVART
      key |= 3;
#endif
#if EXAFMM_CARTESIAN
      key |= 0 << 3;
#elif EXAFMM_SPHERICAL
      key |= 1 << 3;
#elif EXAFMM_PLANEWAVE
      key |= 2 << 3;
#elif EXAFMM_EQUIVALENT
      key |= 3 << 3;
#endif
      return key;
    }

  public:
    Args(int argc=0, char ** argv=NULL) :
      ncrit(64),
      cutoff(.0),
      distribution("cube"),
      dual(0),
      graft(0),
      getMatrix(0),
      images(0),
      IneJ(0),
      mutual(0),
      numBodies(1000000),
      useRopt(0),
      PP(4),
      repeat(1),
      nspawn(5000),
      theta(.4),
      threads(16),
      verbose(0),
      write(0),
      useRmax(0) {
      while (1) {
#if _SX
#warning SX does not have getopt_long
	int c = getopt(argc, argv, "c:d:DgGhi:jmn:oP:r:s:t:T:vwx");
#else
	int option_index;
	int c = getopt_long(argc, argv, "c:d:DgGhi:jmn:oP:r:s:t:T:vwx", long_options, &option_index);
#endif
	if (c == -1) break;
	switch (c) {
	case 'c':
	  ncrit = atoi(optarg);
	  break;
	case 'C':
	  cutoff = atof(optarg);
	  break;
	case 'd':
	  distribution = parse(optarg);
	  break;
	case 'D':
	  dual = 1;
	  break;
	case 'g':
	  graft = 1;
	  break;
	case 'G':
	  getMatrix = 1;
	  break;
	case 'h':
	  usage(argv[0]);
	  abort();
	case 'i':
	  images = atoi(optarg);
	  break;
	case 'j':
	  IneJ = 1;
	  break;
	case 'm':
	  mutual = 1;
	  break;
	case 'n':
	  numBodies = atoi(optarg);
	  break;
	case 'o':
	  useRopt = 1;
	  break;
	case 'P':
	  PP = atoi(optarg);
	  break;
	case 'r':
	  repeat = atoi(optarg);
	  break;
	case 's':
	  nspawn = atoi(optarg);
	  break;
	case 't':
	  theta = atof(optarg);
	  break;
	case 'T':
	  threads = atoi(optarg);
	  break;
	case 'v':
	  verbose = 1;
	  break;
	case 'w':
	  write = 1;
	  break;
	case 'x':
	  useRmax = 1;
	  break;
	default:
	  usage(argv[0]);
	  abort();
	}
      }
    }

    uint64_t getKey(int mpisize=1) {
      uint64_t key = 0;
      key |= uint64_t(log(ncrit)/log(2)); assert(key <= 1 << 4);
      key |= getDistNum(distribution) << 4; assert(key <= 1 << 7);
      key |= dual << 7; assert(key <= 1 << 8);
      key |= graft << 8; assert(key <= 1 << 9);
      key |= images << 9; assert(key <= 1 << 11);
      key |= IneJ << 11; assert(key <= 1 << 12);
      key |= mutual << 12; assert(key <= 1 << 13);
      key |= uint64_t(log(numBodies)/log(10)) << 13; assert(key <= 1 << 17);
      key |= useRopt << 17; assert(key <= 1 << 8);
      key |= PP << 18; assert(key <= 1 << 24);
      key |= uint64_t(log(nspawn)/log(10)) << 24; assert(key <= 1 << 27);
      key |= uint64_t(theta*14) << 27; assert(key <= 1 << 30);
      key |= uint64_t(log(threads)/log(2)) << 30; assert(key <= uint64_t(1) << 33);
      key |= uint64_t(useRmax) << 33; assert(key <= uint64_t(1) << 34);
      key |= getKernelNum() << 34; assert(key <= uint64_t(1) << 40);
      key |= uint64_t(log(mpisize)/log(2)) << 40;
      return key;
    }

    void print(int stringLength, int PP) {
      if (verbose) {
	std::cout << std::setw(stringLength) << std::fixed << std::left
		  << "ncrit" << " : " << ncrit << std::endl
		  << std::setw(stringLength)
		  << "cutoff" << " : " << cutoff << std::endl
		  << std::setw(stringLength)
		  << "distribution" << " : " << distribution << std::endl
		  << std::setw(stringLength)
		  << "dual" << " : " << dual << std::endl
		  << std::setw(stringLength)
		  << "graft" << " : " << graft << std::endl
		  << std::setw(stringLength)
		  << "getMatrix" << " : " << getMatrix << std::endl
		  << std::setw(stringLength)
		  << "images" << " : " << images << std::endl
		  << std::setw(stringLength)
		  << "IneJ" << " : " << IneJ << std::endl
		  << std::setw(stringLength)
		  << "mutual" << " : " << mutual << std::endl
		  << std::setw(stringLength)
		  << "numBodies" << " : " << numBodies << std::endl
		  << std::setw(stringLength)
		  << "useRopt" << " : " << useRopt << std::endl
		  << std::setw(stringLength)
		  << "P" << " : " << PP << std::endl
		  << std::setw(stringLength)
		  << "repeat" << " : " << repeat << std::endl
		  << std::setw(stringLength)
		  << "nspawn" << " : " << nspawn << std::endl
		  << std::setw(stringLength)
		  << "theta" << " : " << theta << std::endl
		  << std::setw(stringLength)
		  << "threads" << " : " << threads << std::endl
		  << std::setw(stringLength)
		  << "verbose" << " : " << verbose << std::endl
		  << std::setw(stringLength)
		  << "write" << " : " << write << std::endl
		  << std::setw(stringLength)
		  << "useRmax" << " : " << useRmax << std::endl;
      }
    }
  };
}
#endif
