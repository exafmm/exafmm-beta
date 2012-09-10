#ifndef options_h
#define options_h
#include <cstdio>
#include <cstdlib>
#include <getopt.h>

struct option long_options[] = {
  {"numBodies",    1, 0, 'n'},
  {"numTarget",    1, 0, 't'},
  {"ncrit",        1, 0, 'c'},
  {"nspawn",       1, 0, 's'},
  {"images",       1, 0, 'i'},
  {"theta",        1, 0, 'o'},
  {"mutual",       1, 0, 'm'},
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
  const char * distribution;

private:
  void usage(char * name) {
    fprintf(stderr,
            "Usage: %s [options]\n"
            "Options:\n"
            " --numBodies : Number of bodies (%d)\n"
            " --numTarget : Number of targets for error checking (%d)\n"
            " --ncrit : Number of bodies per leaf cell (%d)\n"
            " --nspawn : Threshold for splitting both cells during recursion (%d)\n"
            " --images : Number of periodic image levels (%d)\n"
            " --theta : Multipole acceptance criterion (%f)\n"
            " --mutual [0/1] :  use mutual interaction (%d)\n"
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
  }

public:
  Args(int argc, char ** argv) : numBodies(1000000), numTarget(100), NCRIT(10), NSPAWN(1000), IMAGES(0),
                                 THETA(.6), mutual(1), distribution("cube") {
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
};

#endif
