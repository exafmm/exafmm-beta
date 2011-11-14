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

//! Timer and memory logger
class Logger {
private:
  std::ofstream file;                                           //!< File ID to store log
  Event         tic;                                            //!< Timer (map from string to double)
  Event         memory;                                         //!< Memory logger (map from string to double)

//! Timer function
  double get_time() {
    struct timeval tv;                                          // Time value
    gettimeofday(&tv, NULL);                                    // Get time of day in seconds and microseconds
    return double(tv.tv_sec+tv.tv_usec*1e-6);                   // Combine seconds and microseconds and return
  }

public:
  bool  printNow;                                               //!< Switch to print timings
  Event timer;                                                  //!< Stores timings for all events

//! Constructor
  Logger() {
    file.open("time.dat");                                      // Open timer log file
    printNow = false;                                           // Don't print by default
  }
//! Destructor
  ~Logger() {
    file.close();                                               // Close timer log file
  }

//! Start timer for given event
  void startTimer(std::string event) {
    tic[event] = get_time();                                    // Get time of day and store in tic
  }

//! Stop timer for given event
  double stopTimer(std::string event, bool print=false) {
    double toc = get_time();                                    // Get time of day and store in toc
    timer[event] += toc - tic[event];                           // Accumulate event time to timer
    if(print) std::cout << event << " : " << timer[event] << std::endl;// Print event and timer to screen
    return toc - tic[event];                                    // Return the event time
  }

//! Erase entry in timer
  void eraseTimer(std::string event) {
    timer.erase(event);                                         // Erase event from timer
  }

//! Erase all events in timer
  void resetTimer() {
    timer.clear();                                              // Clear timer
  }

//! Record memory allocation
  void allocMemory(std::string event, double bytes) {
    memory[event] += bytes;                                     // Add allocation event to memory
  }

//! Record memory free
  void freeMemory(std::string event, double bytes) {
    memory[event] -= bytes;                                     // Delete allocation event from memory
  }

//! Print timings of a specific event
  void printTime(std::string event) {
    std::cout << event << " : " << timer[event] << std::endl;   // Print event and timer
  }

//! Print memory usage of a specific event
  void printMemory(std::string event) {
    std::cout << event << " : " << memory[event] << std::endl;  // Print event and memory
  }

//! Print timings of all events
  void printAllTime() {
    for( E_iter E=timer.begin(); E!=timer.end(); ++E ) {        // Loop over all events
      std::cout << E->first << " : " << E->second << std::endl; //  Print event and timer
    }                                                           // End loop over all events
  }

//! Write timings of all events
  void writeTime() {
    for( E_iter E=timer.begin(); E!=timer.end(); ++E ) {        // Loop over all events
      file <<  E->first << " " << E->second << std::endl;       //  Print event and timer
    }                                                           // End loop over all events
  }

//! Copy timer to another timer
  void copyTime(Event &timer2) {
    for( E_iter E2=timer2.begin(); E2!=timer2.end(); ++E2 ) {   // Loop over all events
      timer[E2->first] = E2->second;                            //  Copy timer
    }                                                           // End loop over all events
  }
};

#endif
