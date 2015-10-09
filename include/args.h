// Argument parser based on GNU getopt
#ifndef args_h
#define args_h
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <iomanip>

// In alphabetical order of the short option
static struct option long_options[] = {
  {"ncrit",        required_argument, 0, 'c'},
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
  {"repeat",       required_argument, 0, 'r'},
  {"nspawn",       required_argument, 0, 's'},
  {"theta",        required_argument, 0, 't'},
  {"threads",      required_argument, 0, 'T'},
  {"verbose",      no_argument,       0, 'v'},
  {"write",        no_argument,       0, 'w'},
  {"useRmax",      no_argument,       0, 'x'},
  {0, 0, 0, 0}
};

class Args {
public:
  int ncrit;
  const char * distribution;
  int dual;
  int graft;
  int getMatrix;
  int images;
  int IneJ;
  int mutual;
  int numBodies;
  int useRopt;
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
            " --repeat (-r)                 : Number of iteration loops (%d)\n"
            " --nspawn (-s)                 : Threshold for stopping task creation during recursion (%d)\n"
            " --theta (-t)                  : Multipole acceptance criterion (%f)\n"
            " --threads (-T)                : Number of threads (%d)\n"
	    " --verbose (-v)                : Print information to screen (%d)\n"
            " --write (-w)                  : Write timings to file (%d)\n"
	    " --useRmax (-x)                : Use maximum distance for MAC (%d)\n",
            name,
            ncrit,
            distribution,
	    dual,
	    graft,
	    getMatrix,
            images,
	    IneJ,
            mutual,
            numBodies,
	    useRopt,
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
      exit(0);
    }
    return "";
  }

public:
  Args(int argc=0, char ** argv=NULL) :
#if Helmholtz
					ncrit(1000),
#else
					ncrit(64),
#endif
					distribution("cube"),
					dual(0),
					graft(0),
					getMatrix(0),
					images(0),
					IneJ(0),
					mutual(0),
					numBodies(1000000),
					useRopt(0),
					repeat(1),
					nspawn(5000),
#if Helmholtz
					theta(1.),
#else
					theta(.4),
#endif
					threads(16),
					verbose(0),
					write(0),
					useRmax(0) {
    while (1) {
      int option_index;
      int c = getopt_long(argc, argv, "c:d:DgGhi:jmn:or:s:t:T:vwx", long_options, &option_index);
      if (c == -1) break;
      switch (c) {
      case 'c':
        ncrit = atoi(optarg);
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
        exit(0);
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
        exit(0);
      }
    }
  }

  void print(int stringLength, int PP) {
    if (verbose) {                                              // If verbose flag is true
      std::cout << std::setw(stringLength) << std::fixed << std::left// Set format
		<< "ncrit" << " : " << ncrit << std::endl       //  Print ncrit
		<< std::setw(stringLength)                      //  Set format
		<< "distribution" << " : " << distribution << std::endl// Print distribution
		<< std::setw(stringLength)                      //  Set format
		<< "dual" << " : " << dual << std::endl         //  Print images
		<< std::setw(stringLength)                      //  Set format
		<< "graft" << " : " << graft << std::endl       //  Print graft
		<< std::setw(stringLength)                      //  Set format
		<< "getMatrix" << " : " << getMatrix << std::endl// Print getMatrix
		<< std::setw(stringLength)                      //  Set format
		<< "images" << " : " << images << std::endl     //  Print images
		<< std::setw(stringLength)                      //  Set format
		<< "IneJ" << " : " << IneJ << std::endl         //  Print IneJ
		<< std::setw(stringLength)                      //  Set format
		<< "mutual" << " : " << mutual << std::endl     //  Print mutual
		<< std::setw(stringLength)                      //  Set format
		<< "numBodies" << " : " << numBodies << std::endl// Print numBodies  
		<< std::setw(stringLength)                      //  Set format
		<< "useRopt" << " : " << useRopt << std::endl   //  Print useRopt
		<< std::setw(stringLength)                      //  Set format
		<< "P" << " : " << PP << std::endl              //  Print P
		<< std::setw(stringLength)                      //  Set format
		<< "repeat" << " : " << repeat << std::endl     //  Print repeat
		<< std::setw(stringLength)                      //  Set format
		<< "nspawn" << " : " << nspawn << std::endl     //  Print nspawn
		<< std::setw(stringLength)                      //  Set format
		<< "theta" << " : " << theta << std::endl       //  Print theta
		<< std::setw(stringLength)                      //  Set format
		<< "threads" << " : " << threads << std::endl   //  Print threads
		<< std::setw(stringLength)                      //  Set format
		<< "verbose" << " : " << verbose << std::endl   //  Print verbose
		<< std::setw(stringLength)                      //  Set format
		<< "write" << " : " << write << std::endl       //  Print write
		<< std::setw(stringLength)                      //  Set format
		<< "useRmax" << " : " << useRmax << std::endl;  //  Print useRmax
    }                                                           // End if for verbose flag
  }
};
#endif
