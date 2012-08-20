#ifndef exafmm_config_h
#define exafmm_config_h

#include <unistd.h>
#include <getopt.h>
#include <errno.h>

struct exafmm_config {
  int numBodies;
  const char * distribution;	// cube, lattice, sphere, plummer
#if IMPL_MUTUAL
  int mutual;
#endif
#if PARALLEL_EVERYTHING
  int parallelEverything;
#endif
  int splitBothThreshold;
  int splitParallelThreshold;
  int steps;
  int images;
  double theta;
  int ncrit;
#if SIMDIZATION
  simdize_option simdize;
#endif
  unsigned int evalError;
  int buildOnly;
  int dumpcells;
};

static struct option exafmm_options[] = {
  {"numBodies",              1, 0, 'n'},
  {"distribution",           1, 0, 0},
#if IMPL_MUTUAL
  {"mutual",                 1, 0, 0},
#endif
#if PARALLEL_EVERYTHING
  {"parallelEverything",     1, 0, 0},
#endif
  {"splitBothThreshold",     1, 0, 0},
  {"splitParallelThreshold", 1, 0, 0},
  {"steps",                  1, 0, 0},
  {"images",                 1, 0, 0},
  {"theta",                  1, 0, 0},
  {"ncrit",                  1, 0, 0},
#if SIMDIZATION
  {"simdize",                1, 0, 0},
#endif
  {"evalError",              1, 0, 0},
  {"buildOnly",              1, 0, 0},
  {"dumpcells",              1, 0, 0},
  {0, 0, 0, 0}
};

static void show_exafmm_config(exafmm_config * o) {
  printf("numBodies: %d\n", o->numBodies);
  printf("distribution: %s\n", o->distribution);
#if IMPL_MUTUAL
  printf("mutual: %d\n", o->mutual);
#endif
#if PARALLEL_EVERYTHING
  printf("parallelEverything: %d\n", o->parallelEverything);
#endif
  printf("splitBothThreshold: %d\n", o->splitBothThreshold);
  printf("splitParallelThreshold: %d\n", o->splitParallelThreshold);
  printf("steps: %d\n", o->steps);
  printf("images: %d\n", o->images);
  printf("theta: %f\n", o->theta);
  printf("ncrit: %d\n", o->ncrit);
#if SIMDIZATION
  printf("simdize: %d\n", o->simdize);
#endif
  printf("evalError: %u\n", o->evalError);
  printf("buildOnly: %d\n", o->buildOnly);
  printf("dumpcells: %d\n", o->dumpcells);
}

static exafmm_config mk_default_exafmm_config() {
  exafmm_config o;
  o.numBodies = 1000000;
  o.distribution = "plummer";
#if IMPL_MUTUAL
  o.mutual = 1;
#endif
#if PARALLEL_EVERYTHING
  o.parallelEverything = true;
#endif
  o.splitBothThreshold = 1000;
  o.splitParallelThreshold = 0;
  o.steps = 2;
  o.images = 0;
  o.theta = 0.6;
  o.ncrit = 32;
#if SIMDIZATION
  o.simdize = simdize_sse;
#endif
  o.evalError = 100;
  o.buildOnly = 0;
  o.dumpcells = 0;
  return o;
}

static void exafmm_usage(char * progname) {
  exafmm_config o = mk_default_exafmm_config();
  fprintf(stderr, 
	  "Usage: %s [options]\n"
	  "Options:\n"
	  " --numBodies N :\n"
	  "  simulate N particles (%d)\n"
	  " --distribution [l/c/s/p] :\n"
	  "  lattice, cube, sphere, plummer (%s)\n"
#if IMPL_MUTUAL
	  " --mutual [0/1] :\n"
	  "  use mutual interaction (%d)\n"
#endif
#if PARALLEL_EVERYTHING
	  " --parallelEverything [0/1] :\n"
	  "  parallelize tree building, upward/downward passes (%d)\n"
#endif
	  " --splitBothThreshold N :\n"
	  "  split both cells on recursion when the two cells contain >= N particles (%d)\n"
	  " --splitParallelThreshold N :\n"
	  "  parallel recursion for >= N particles (%d)\n"
	  " --steps N :\n"
	  "  set number of steps to N (%d)\n"
	  " --images [0/1] :\n"
	  "  when set, periodic (%d)\n"
	  " --theta T :\n"
	  "  set theta to T (%f)\n"
	  " --ncrit N :\n"
	  "  leaf nodes have <= N particles (%d)\n"
#if SIMDIZATION
	  " --simdize [none 0/sse 1/avx 2] :\n"
	  "  use simd instructions (%d)\n"
#endif
	  " --evalError D :\n"
	  "  check and show error of D particles (%d)\n"
	  " --buildOnly [0/1] :\n"
	  "  build tree and do not evaluate force (%d)\n"
#if 0
	  " --dumpcells [0/1] :\n"
	  "  dump cells after building tree (%d)\n"
#endif
	  ,
	  progname,
	  o.numBodies,
	  o.distribution,
#if IMPL_MUTUAL
	  o.mutual,
#endif
#if PARALLEL_EVERYTHING
	  o.parallelEverything,
#endif
	  o.splitBothThreshold,
	  o.splitParallelThreshold,
	  o.steps,
	  o.images,
	  o.theta,
	  o.ncrit,
#if SIMDIZATION
	  o.simdize,
#endif
	  o.evalError,
	  o.buildOnly
#if 0
	  , o.dumpcells
#endif
	  );
}

static int safe_atoi(const char * arg, int * x) {
  char * endptr;
  errno = 0;
  *x = strtol(arg, &endptr, 10);
  if (errno) {
    perror("strtol");
    return 0;			// NG
  }
  if (endptr == arg) {
    fprintf(stderr, "No digits were found\n");
    return 0;			// NG
  }
  return 1;			// OK
}

static int safe_atou(const char * arg, unsigned int * x) {
  char * endptr;
  errno = 0;
  *x = strtoul(arg, &endptr, 10);
  if (errno) {
    perror("strtoul");
    return 0;			// NG
  }
  if (endptr == arg) {
    fprintf(stderr, "No digits were found\n");
    return 0;			// NG
  }
  return 1;			// OK
}

static int safe_atof(const char * arg, double * x) {
  char * endptr;
  errno = 0;
  *x = strtod(arg, &endptr);
  if (errno) {
    perror("strtod");
    return 0;			// NG
  }
  if (endptr == arg) {
    fprintf(stderr, "No digits were found\n");
    return 0;			// NG
  }
  return 1;			// OK
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

#if SIMDIZATION
static simdize_option parse_simdize_option(const char * arg) {
  switch (arg[0]) {
  case 'n':
    return simdize_none;
  case 's':
    return simdize_sse;
  case 'a':
    return simdize_avx;
  default:
    {
      int x;
      if (safe_atoi(arg, &x)) {
	return (simdize_option)x;
      } else {
	fprintf(stderr, "invalid simdize argument %s, no simdization\n", arg);
	return simdize_none;
      } 
    }
  }
}
#endif

static exafmm_config * parse_cmdline_args(int argc, char ** argv, exafmm_config * o) {
  *o = mk_default_exafmm_config();
  while (1) {
    int option_index;
    int c = getopt_long(argc, argv, "", exafmm_options, &option_index);
    if (c == -1) break;
    switch (c) {
    case 0:
      {
	const char * name = exafmm_options[option_index].name;
	if (strcmp(name, "distribution") == 0) {
	  o->distribution = parse_distribution(optarg);
	  if (o->distribution == NULL) return NULL;
#if IMPL_MUTUAL
	} else if (strcmp(name, "mutual") == 0) {
	  if (!safe_atoi(optarg, &o->mutual)) return NULL;
#endif
#if PARALLEL_EVERYTHING
	} else if (strcmp(name, "parallelEverything") == 0) {
	  if (!safe_atoi(optarg, &o->parallelEverything)) return NULL;
#endif
	} else if (strcmp(name, "splitBothThreshold") == 0) {
	  if (!safe_atoi(optarg, &o->splitBothThreshold)) return NULL;
	} else if (strcmp(name, "splitParallelThreshold") == 0) {
	  if (!safe_atoi(optarg, &o->splitParallelThreshold)) return NULL;
	} else if (strcmp(name, "steps") == 0) {
	  if (!safe_atoi(optarg, &o->steps)) return NULL;
	} else if (strcmp(name, "theta") == 0) {
	  if (!safe_atof(optarg, &o->theta)) return NULL;
	} else if (strcmp(name, "ncrit") == 0) {
	  if (!safe_atoi(optarg, &o->ncrit)) return NULL;
#if SIMDIZATION
	} else if (strcmp(name, "simdize") == 0) {
	  o->simdize = parse_simdize_option(optarg);
#endif
	} else if (strcmp(name, "evalError") == 0) {
	  if (!safe_atou(optarg, &o->evalError)) return NULL;
	} else if (strcmp(name, "buildOnly") == 0) {
	  if (!safe_atoi(optarg, &o->buildOnly)) return NULL;
	} else if (strcmp(name, "dumpcells") == 0) {
	  if (!safe_atoi(optarg, &o->dumpcells)) return NULL;
	} else {
	  fprintf(stderr, "?? getopt returned long opt %s\n", name);
	  exafmm_usage(argv[0]);
	  return NULL;
	}
	break;
      }
    case 'n':
      if (!safe_atoi(optarg, &o->numBodies)) return NULL;
      break;
    default: /* '?' */
      exafmm_usage(argv[0]);
      return NULL;
    }
  }
  return o;
}

#endif
