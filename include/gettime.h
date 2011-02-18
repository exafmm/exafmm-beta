#ifndef gettime_h
#define gettime_h
#include <sys/time.h>

double get_time() {                                             // Timer function
  struct timeval tv;                                            // Time value
  gettimeofday(&tv, NULL);                                      // Get time of day in seconds and microseconds
  return double(tv.tv_sec+tv.tv_usec*1e-6);                     // Combine seconds and microseconds and return
}

#endif
