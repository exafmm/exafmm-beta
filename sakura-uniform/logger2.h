#ifndef logger_h
#define logger_h
#include <iomanip>
#include <iostream>
#include <map>
#include <sys/time.h>

//! Timer and Tracer logger
namespace logger {
  typedef std::map<std::string,double> Timer;                   //!< Map of timer event name to timed value

  Timer beginTimer;                                             //!< Timer base value
  Timer timer;                                                  //!< Timings of all events
  int stringLength = 20;                                        //!< Max length of event name
  int decimal = 7;                                              //!< Decimal precision
  bool verbose = false;                                         //!< Print to screen

  //! Timer function
  double get_time() {
    struct timeval tv;                                          // Time value
    gettimeofday(&tv, NULL);                                    // Get time of day in seconds and microseconds
    return double(tv.tv_sec+tv.tv_usec*1e-6);                   // Combine seconds and microseconds and return
  }

  //! Print message to standard output
  inline void printTitle(std::string title) {
    if (verbose) {                                              // If verbose flag is true
      title += " ";                                             //  Append space to end of title
      std::cout << "--- " << std::setw(stringLength)            //  Align string length
                << std::left                                    //  Left shift
                << std::setfill('-')                            //  Set to fill with '-'
                << title << std::setw(10) << "-"                //  Fill until end of line
                << std::setfill(' ') << std::endl;              //  Set back to fill with ' '
    }                                                           // End if for verbose flag
  }

  //! Start timer for given event
  inline void startTimer(std::string event) {
    beginTimer[event] = get_time();                             // Get time of day and store in beginTimer
  }

  //! Print timings of a specific event
  inline void printTime(std::string event) {
    if (verbose) {                                              // If verbose flag is true
      std::cout << std::setw(stringLength) << std::left         //  Set format
		<< event << " : " << std::setprecision(decimal) << std::fixed
		<< timer[event] << " s" << std::endl;           //  Print event and timer
    }                                                           // End if for verbose flag
  }

  //! Stop timer for given event
  double stopTimer(std::string event, int print=1) {
    double endTimer = get_time();                               // Get time of day and store in endTimer
    timer[event] += endTimer - beginTimer[event];               // Accumulate event time to timer
    if (verbose && print) printTime(event);                     // Print event and timer to screen
    return endTimer - beginTimer[event];                        // Return the event time
  }

  //! Print error between FMM and direct
  inline void printError(std::string title, double v) {
    if (logger::verbose) {                                      // If verbose flag is true
      std::cout << std::setw(logger::stringLength) << std::left //  Set format
                << title << " : " << std::setprecision(logger::decimal) << std::scientific // Set title
                << v << std::endl;                              //  Print potential error
    }                                                           // End if for verbose flag
  }
};
#endif
