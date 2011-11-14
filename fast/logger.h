/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#ifndef logger_h
#define logger_h
#include <sys/time.h>
#include "types.h"

class Logger {
private:
  std::ofstream file;
  Event         tic;
  Event         memory;

  double get_time() const {                                     // Timer function
    struct timeval tv;                                          // Time value
    gettimeofday(&tv, NULL);                                    // Get time of day in seconds and microseconds
    return double(tv.tv_sec+tv.tv_usec*1e-6);                   // Combine seconds and microseconds and return
  }

public:
  bool  printNow;
  Event timer;

  Logger() {
    file.open("time.dat");
    printNow = false;
  }
  ~Logger() {
    file.close();
  }

  void startTimer(const std::string &event) {
    tic[event] = get_time();
  }

  double stopTimer(const std::string &event, bool print=false) {
    double toc = get_time();
    timer[event] += toc - tic[event];
    if(print) std::cout << event << " : " << timer[event] << std::endl;
    return toc - tic[event];
  }

  void eraseTimer(std::string event) {
    timer.erase(event);
  }

  void resetTimer() {
    timer.clear();
  }

  void allocMemory(std::string event, double bytes) {
    memory[event] += bytes;
  }

  void freeMemory(std::string event, double bytes) {
    memory[event] -= bytes;
  }

  void printTime(std::string event) {
    std::cout << event << " : " << timer[event] << std::endl;
  }

  void printMemory(std::string event) {
    std::cout << event << " : " << memory[event] << std::endl;
  }

  void printAllTime() {
    for( E_iter E=timer.begin(); E!=timer.end(); ++E ) {
      std::cout << E->first << " : " << E->second << std::endl;
    }
  }

  void writeTime() {
    for( E_iter E=timer.begin(); E!=timer.end(); ++E ) {
      file <<  E->first << " " << E->second << std::endl;
    }
  }

  void copyTime(Event &timer2) {
    for( E_iter E2=timer2.begin(); E2!=timer2.end(); ++E2 ) {
      timer[E2->first] = E2->second;
    }
  }
};

#endif
