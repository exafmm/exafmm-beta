#ifndef timer_h
#define timer_h
#include <sys/time.h>

namespace {
  double get_time(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return double(tv.tv_sec+tv.tv_usec*1e-6);
  }

  double timer;

  void start_timer(){
    timer = get_time();
  }

  void stop_timer(char *event_name){
    printf("%-20s: %8.4lf s\n",event_name, get_time()-timer);
  }
}

#endif
