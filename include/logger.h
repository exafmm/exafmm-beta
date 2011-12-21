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

//! Timer and Trace logger
class Logger {
private:
  std::ofstream   file;                                         //!< File ID to store log
  Timer           beginTimer;                                   //!< Timer base value
  Timer           timer;                                        //!< Stores timings for all events
  Traces          traces;                                       //!< Stores traces for all events
  pthread_mutex_t mutex;                                        //!< Pthread communicator

//! Timer function
  double get_time() const {
    struct timeval tv;                                          // Time value
    gettimeofday(&tv, NULL);                                    // Get time of day in seconds and microseconds
    return double(tv.tv_sec+tv.tv_usec*1e-6);                   // Combine seconds and microseconds and return
  }

public:
  bool printNow;                                                //!< Switch to print timings

//! Constructor
  Logger() {
    file.open("time.dat");                                      // Open timer log file
    printNow = false;                                           // Don't print by default
    pthread_mutex_init(&mutex,NULL);                            // Initialize pthread communicator
  }
//! Destructor
  ~Logger() {
    file.close();                                               // Close timer log file
  }

//! Start timer for given event
  inline void startTimer(std::string event) {
    beginTimer[event] = get_time();                             // Get time of day and store in beginTimer
  }

//! Stop timer for given event
  inline double stopTimer(std::string event, bool print=false) {
    double endTimer = get_time();                               // Get time of day and store in endTimer
    timer[event] += endTimer - beginTimer[event];               // Accumulate event time to timer
    if(print) std::cout << event << " : " << timer[event] << std::endl;// Print event and timer to screen
    return endTimer - beginTimer[event];                        // Return the event time
  }

//! Erase entry in timer
  inline void eraseTimer(std::string event) {
    timer.erase(event);                                         // Erase event from timer
  }

//! Erase all events in timer
  inline void resetTimer() {
    timer.clear();                                              // Clear timer
  }

//! Print timings of a specific event
  inline void printTime(std::string event) {
    std::cout << event << " : " << timer[event] << std::endl;   // Print event and timer
  }

//! Print timings of all events
  inline void printAllTime() {
    for( TI_iter E=timer.begin(); E!=timer.end(); ++E ) {       // Loop over all events
      std::cout << E->first << " : " << E->second << std::endl; //  Print event and timer
    }                                                           // End loop over all events
  }

//! Write timings of all events
  inline void writeTime() {
    for( TI_iter E=timer.begin(); E!=timer.end(); ++E ) {       // Loop over all events
      file << E->first << " " << E->second << std::endl;        //  Print event and timer
    }                                                           // End loop over all events
  }

//! Start tracer for given event
  inline void startTracer(ThreadTrace &beginTrace) {
    pthread_mutex_lock(&mutex);                                 // Lock shared variable access
    beginTrace[pthread_self()] = get_time();                    // Get time of day and store in beginTrace
    pthread_mutex_unlock(&mutex);                               // Unlock shared variable access
  }

//! Stop tracer for given event
  inline void stopTracer(ThreadTrace &beginTrace, int color) {
    pthread_mutex_lock(&mutex);                                 // Lock shared variable access
    Trace trace;                                                // Define trace structure
    trace.thread = pthread_self();                              // Store pthread id
    trace.begin  = beginTrace[pthread_self()];                  // Store tic
    trace.end    = get_time();                                  // Store toc
    trace.color  = color;                                       // Store color of event
    traces.push(trace);                                         // Push trace to queue of traces
    pthread_mutex_unlock(&mutex);                               // Unlock shared variable access
  }

//! Write traces of all events
  inline void writeTrace() {
    FILE *fid = fopen("trace.svg", "w");
    double scale = 30000.0;
    fprintf(fid,
      "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
      "<!DOCTYPE svg PUBLIC \"-_W3C_DTD SVG 1.0_EN\" \"http://www.w3.org/TR/SVG/DTD/svg10.dtd\">\n"
      "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n"
      "  width=\"200mm\" height=\"40mm\" viewBox=\"0 0 20000 4000\">\n"
      "  <g>\n" );
    int num_thread = 0;
    ThreadMap threadMap;
    double base = traces.front().begin;
    while( !traces.empty() ) {
      Trace trace = traces.front();
      traces.pop();
      pthread_t thread = trace.thread;
      double begin  = trace.begin;
      double end    = trace.end;
      int    color  = trace.color;
      if( threadMap[thread] == 0 ) {
        threadMap[thread] = ++num_thread;
      }
      begin -= base;
      end   -= base;
      fprintf(fid,
        "    "
        "<rect x=\"%.2lf\" y=\"%.0lf\" width=\"%.2lf\" height=\"%.0lf\" "
        "fill=\"#%06x\" stroke=\"#000000\" stroke-width=\"1\"/>\n",
        begin * scale,
        (threadMap[thread] - 1) * 100.0,
        (end - begin) * scale,
        90.0,
        color);
    }
    fprintf(fid,"  </g>\n" "</svg>\n");
    fclose(fid);
  }
};

#endif
