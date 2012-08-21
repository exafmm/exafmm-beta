#pragma once
#include <getopt.h>

struct Args {
  int numBodies;
  int numTarget;
  int ncrit;
  int nspawn;
  int images;
  double theta;
  int buildOnly;
  int mutual;
  const char * distribution;
};

static struct option long_options[] = {
  {"numBodies",    1, 0, 'n'},
  {"numTarget",    1, 0, 't'},
  {"ncrit",        1, 0, 'c'},
  {"nspawn",       1, 0, 's'},
  {"images",       1, 0, 'i'},
  {"theta",        1, 0, 'h'},
  {"buildOnly",    1, 0, 'b'},
  {"mutual",       1, 0, 'm'},
  {"distribution", 1, 0, 'd'},
  {0, 0, 0, 0}
};

static void showArgs(Args * args) {
  printf("numBodies: %d\n", args->numBodies);
  printf("numTarget: %d\n", args->numTarget);
  printf("ncrit: %d\n", args->ncrit);
  printf("nspawn: %d\n", args->nspawn);
  printf("images: %d\n", args->images);
  printf("theta: %f\n", args->theta);
  printf("buildOnly: %d\n", args->buildOnly);
  printf("mutual: %d\n", args->mutual);
  printf("distribution: %s\n", args->distribution);
}

static Args defaultArgs() {
  Args args;
  args.numBodies = 1000000;
  args.numTarget = 100;
  args.ncrit = 10;
  args.nspawn = 1000;
  args.images = 0;
  args.theta = 0.6;
  args.buildOnly = 0;
  args.mutual = 1;
  args.distribution = "cube";
  return args;
}

static void usage(char * name) {
  Args args = defaultArgs();
  fprintf(stderr, 
          "Usage: %s [options]\n"
          "Options:\n"
          " --numBodies : Number of bodies (%d)\n"
          " --numTarget : Number of targets for error checking (%d)\n"
          " --ncrit : Number of bodies per leaf cell (%d)\n"
          " --nspawn : Threshold for splitting both cells during recursion (%d)\n"
          " --images : Number of periodic image levels (%d)\n"
          " --theta : Multipole acceptance criterion (%f)\n"
          " --buildOnly [0/1] : build tree and do not evaluate force (%d)\n"
          " --mutual [0/1] :  use mutual interaction (%d)\n"
          " --distribution [l/c/s/p] : lattice, cube, sphere, plummer (%s)\n",
          name,
          args.numBodies,
          args.numTarget,
          args.ncrit,
          args.nspawn,
          args.images,
          args.theta,
          args.buildOnly,
          args.mutual,
          args.distribution);
}

static const char * parse_distribution(const char * arg) {
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
    return NULL;
  }
}

static Args * parse_cmdline_args(int argc, char ** argv, Args * args) {
  *args = defaultArgs();
  while (1) {
    int option_index;
    int c = getopt_long(argc, argv, "", long_options, &option_index);
    if (c == -1) break;
    switch (c) {
    case 'n':
      args->numBodies = atoi(optarg);
      break;
    case 't':
      args->numTarget = atoi(optarg);
      break;
    case 'c':
      args->ncrit = atoi(optarg);
      break;
    case 's':
      args->nspawn = atoi(optarg);
      break;
    case 'i':
      args->images = atoi(optarg);
      break;
    case 'h':
      args->theta = atof(optarg);
      break;
    case 'b':
      args->buildOnly = atoi(optarg);
      break;
    case 'm':
      args->mutual = atoi(optarg);
      break;
    case 'd':
      args->distribution = parse_distribution(optarg);
      if (args->distribution == NULL) return NULL;
      break;
    default:
      usage(argv[0]);
      return NULL;
    }
  }
  return args;
}

static int gendata(Dataset& data, Bodies& bodies, const char * distribution) {
  switch (distribution[0]) {
  case 'l':
    data.lattice(bodies);
    return 1;
  case 'c':
    data.cube(bodies);
    return 1;
  case 's':
    data.sphere(bodies);
    return 1;
  case 'p':
    data.plummer(bodies);
    return 1;                        
  default:
    fprintf(stderr, "unknown data distribution %s\n", distribution);
    return 0;                        
  }
}
