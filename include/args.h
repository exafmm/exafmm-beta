#ifndef args_h
#define args_h
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <iomanip>

static struct option long_options[] = {
  {"numBodies",    1, 0, 'n'},
  {"ncrit",        1, 0, 'c'},
  {"nspawn",       1, 0, 's'},
  {"threads",      1, 0, 'T'},
  {"images",       1, 0, 'i'},
  {"theta",        1, 0, 't'},
  {"useRmax",      1, 0, 'x'},
  {"useRopt",      1, 0, 'o'},
  {"mutual",       1, 0, 'm'},
  {"graft",        1, 0, 'g'},
  {"verbose",      1, 0, 'v'},
  {"distribution", 1, 0, 'd'},
  {"repeat",       1, 0, 'r'},
  {"help",         0, 0, 'h'},
  {0, 0, 0, 0}
};

class Args {
public:
  int numBodies;
  int ncrit;
  int nspawn;
  int threads;
  int images;
  double theta;
  int useRmax;
  int useRopt;
  int mutual;
  int graft;
  int verbose;
  const char * distribution;
  int repeat;

private:
  void usage(char * name) {
    fprintf(stderr,
            "Usage: %s [options]\n"
            "Long option (short option)     : Description (Default value)\n"
            " --numBodies (-n)              : Number of bodies (%d)\n"
            " --ncrit (-c)                  : Number of bodies per leaf cell (%d)\n"
            " --nspawn (-s)                 : Threshold for stopping task creation during recursion (%d)\n"
            " --threads (-T)                : Number of threads (%d)\n"
            " --images (-i)                 : Number of periodic image levels (%d)\n"
            " --theta (-t)                  : Multipole acceptance criterion (%f)\n"
	    " --useRmax (-x) [0/1]          : Use maximum distance for MAC (%d)\n"
	    " --useRopt (-o) [0/1]          : Use error optimized theta for MAC (%d)\n"
            " --mutual (-m) [0/1]           : Use mutual interaction (%d)\n"
	    " --graft (-g) [0/1]            : Graft remote trees to global tree (%d)\n"
	    " --verbose (-v) [0/1]          : Print information to screen (%d)\n"
            " --distribution (-d) [l/c/s/p] : lattice, cube, sphere, plummer (%s)\n"
            " --repeat (-r)                 : Number of iteration loops (%d)\n"
            " --help (-h)                   : Show this help document\n",
            name,
            numBodies,
            ncrit,
            nspawn,
	    threads,
            images,
            theta,
	    useRmax,
	    useRopt,
            mutual,
	    graft,
	    verbose,
            distribution,
	    repeat);
  }

  const char * parse(const char * arg) {
    switch (arg[0]) {
    case 'l':
      return "lattice";
    case 'c':
      return "cube";
    case 's':
      return "sphere";
    case 'p':
      return "plummer";
    default:
      fprintf(stderr, "invalid distribution %s\n", arg);
      exit(0);
    }
    return "";
  }

public:
  Args(int argc=0, char ** argv=NULL) : numBodies(1000000), ncrit(16), nspawn(1000), threads(16), images(0),
					theta(.4), useRmax(1), useRopt(1), mutual(1), graft(1),
					verbose(1), distribution("cube"), repeat(1) {
    while (1) {
      int option_index;
      int c = getopt_long(argc, argv, "n:c:s:T:i:t:x:o:m:g:v:d:r:h", long_options, &option_index);
      if (c == -1) break;
      switch (c) {
      case 'n':
        numBodies = atoi(optarg);
        break;
      case 'c':
        ncrit = atoi(optarg);
        break;
      case 's':
        nspawn = atoi(optarg);
        break;
      case 'T':
        threads = atoi(optarg);
        break;
      case 'i':
        images = atoi(optarg);
        break;
      case 't':
        theta = atof(optarg);
        break;
      case 'x':
        useRmax = atof(optarg);
        break;
      case 'o':
        useRopt = atof(optarg);
        break;
      case 'm':
        mutual = atoi(optarg);
        break;
      case 'g':
	graft = atoi(optarg);
	break;
      case 'v':
	verbose= atoi(optarg);
	break;
      case 'd':
        distribution = parse(optarg);
        break;
      case 'r':
        repeat = atoi(optarg);
        break;
      case 'h':
        usage(argv[0]);
        exit(0);
      default:
        usage(argv[0]);
        exit(0);
      }
    }
  }

  void print(int stringLength, int P) {
    if (verbose) {                                              // If verbose flag is true
      std::cout << std::setw(stringLength) << std::fixed << std::left// Set format
		<< "numBodies" << " : " << numBodies << std::endl // Print numBodies  
		<< std::setw(stringLength)                      //  Set format
		<< "P" << " : " << P << std::endl               //  Print P
		<< std::setw(stringLength)                      //  Set format
		<< "theta" << " : " << theta << std::endl       //  Print theta
		<< std::setw(stringLength)                      //  Set format
		<< "ncrit" << " : " << ncrit << std::endl       //  Print ncrit
		<< std::setw(stringLength)                      //  Set format
		<< "nspawn" << " : " << nspawn << std::endl     //  Print nspawn
		<< std::setw(stringLength)                      //  Set format
		<< "threads" << " : " << threads << std::endl   //  Print threads
		<< std::setw(stringLength)                      //  Set format
		<< "images" << " : " << images << std::endl     //  Print images
		<< std::setw(stringLength)                      //  Set format
		<< "useRmax" << " : " << useRmax << std::endl   //  Print useRmax
		<< std::setw(stringLength)                      //  Set format
		<< "useRopt" << " : " << useRopt << std::endl   //  Print useRopt
		<< std::setw(stringLength)                      //  Set format
		<< "mutual" << " : " << mutual << std::endl     //  Print mutual
		<< std::setw(stringLength)                      //  Set format
		<< "graft" << " : " << graft << std::endl       //  Print graft
		<< std::setw(stringLength)                      //  Set format
		<< "verbose" << " : " << verbose << std::endl   //  Print verbose
		<< std::setw(stringLength)                      //  Set format
		<< "distribution" << " : " << distribution << std::endl// Print distribution
		<< std::setw(stringLength)                      //  Set format
		<< "repeat" << " : " << repeat << std::endl;    //  Print distribution
    }                                                           // End if for verbose flag
  }
};
#endif
