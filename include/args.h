#ifndef args_h
#define args_h
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <iomanip>

struct option long_options[] = {
  {"numBodies",    1, 0, 'n'},
  {"numTarget",    1, 0, 't'},
  {"ncrit",        1, 0, 'c'},
  {"nspawn",       1, 0, 's'},
  {"images",       1, 0, 'i'},
  {"theta",        1, 0, 'o'},
  {"mutual",       1, 0, 'm'},
  {"verbose",      1, 0, 'v'},
  {"distribution", 1, 0, 'd'},
  {"help",         0, 0, 'h'},
  {0, 0, 0, 0}
};

class Args {
 public:
  int numBodies;
  int numTarget;
  int NCRIT;
  int NSPAWN;
  int IMAGES;
  double THETA;
  int mutual;
  int verbose;
  const char * distribution;

 private:
  void usage(char * name) {
    fprintf(stderr,
            "Usage: %s [options]\n"
            "Option : Description (Default value):\n"
            " --numBodies : Number of bodies (%d)\n"
            " --numTarget : Number of targets for error checking (%d)\n"
            " --ncrit : Number of bodies per leaf cell (%d)\n"
            " --nspawn : Threshold for stopping thread spawning during recursion (%d)\n"
            " --images : Number of periodic image levels (%d)\n"
            " --theta : Multipole acceptance criterion (%f)\n"
            " --mutual [0/1] : Use mutual interaction (%d)\n"
	    " --verbose [0/1] : Print information to screen (%d)\n"
            " --distribution [l/c/s/p] : lattice, cube, sphere, plummer (%s)\n"
            " --help : Show this help document\n",
            name,
            numBodies,
            numTarget,
            NCRIT,
            NSPAWN,
            IMAGES,
            THETA,
            mutual,
	    verbose,
            distribution);
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
  Args(int argc=0, char ** argv=NULL) : numBodies(1000000), numTarget(100), NCRIT(16), NSPAWN(1000), IMAGES(0),
    THETA(.6), mutual(1), verbose(1), distribution("cube") {
    while (1) {
      int option_index;
      int c = getopt_long(argc, argv, "", long_options, &option_index);
      if (c == -1) break;
      switch (c) {
      case 'n':
        numBodies = atoi(optarg);
        break;
      case 't':
        numTarget = atoi(optarg);
        break;
      case 'c':
        NCRIT = atoi(optarg);
        break;
      case 's':
        NSPAWN = atoi(optarg);
        break;
      case 'i':
        IMAGES = atoi(optarg);
        break;
      case 'o':
        THETA = atof(optarg);
        break;
      case 'm':
        mutual = atoi(optarg);
        break;
      case 'v':
	verbose= atoi(optarg);
	break;
      case 'd':
        distribution = parse(optarg);
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
    if (verbose) {
      std::cout << std::setw(stringLength) << std::left         // Set format
		<< "numBodies" << " : " << numBodies << std::endl // Print numBodies  
		<< std::setw(stringLength)                      // Set format
		<< "P" << " : " << P << std::endl               // Print P
		<< std::setw(stringLength)                      // Set format
		<< "THETA" << " : " << THETA << std::endl       // Print THETA
		<< std::setw(stringLength)                      // Set format
		<< "NCRIT" << " : " << NCRIT << std::endl       // Print NCRIT
		<< std::setw(stringLength)                      // Set format
		<< "NSPAWN" << " : " << NSPAWN << std::endl     // Print NSPAWN
		<< std::setw(stringLength)                      // Set format
		<< "IMAGES" << " : " << IMAGES << std::endl     // Print IMAGES
		<< std::setw(stringLength)                      // Set format
		<< "mutual" << " : " << mutual << std::endl     // Print mutual
		<< std::setw(stringLength)                      // Set format
		<< "verbose" << " : " << verbose << std::endl   // Print verbose
		<< std::setw(stringLength)                      // Set format
		<< "distribution" << " : " << distribution << std::endl;// Print distribution
    } else {
      std::cout << std::setw(stringLength) << std::left         // Set format
		<< "numBodies" << " : " << numBodies << std::endl; // Print numBodies  
    }
  }
};
#endif
