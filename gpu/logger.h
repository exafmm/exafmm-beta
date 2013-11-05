#pragma once
#include <sys/time.h>

static inline double get_time() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return double(tv.tv_sec + tv.tv_usec * 1e-6);
}
