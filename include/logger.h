#ifndef logger_h
#define logger_h
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <pthread.h>
#include <queue>
#include <stdint.h>
#include <string>
#include <sstream>
#include <sys/time.h>
#include "thread.h"
#include <vector>

#if PAPI
#include <cstring>
#include <papi.h>
#endif

//! Structure for pthread based tracer
struct Tracer {
  pthread_t thread;                                             //!< pthread id
  double    begin;                                              //!< Begin timer of tracer
  double    end;                                                //!< End timer of tracer
  Tracer() {}                                                   //!< Constructor
};

//! Timer and Tracer logger
namespace logger {
  typedef std::map<std::string,double> Timer;                   //!< Map of timer event name to timed value
  typedef Timer::iterator              T_iter;                  //!< Iterator of timer event name map
  typedef std::queue<Tracer>           Tracers;                 //!< Queue of tracers
  typedef std::map<pthread_t,int>      ThreadMap;               //!< Map of pthread id to thread id

  Timer           beginTimer;                                   //!< Timer base value
  Timer           timer;                                        //!< Timings of all events
  Tracers         tracers;                                      //!< Tracers for all events
  pthread_mutex_t mutex;                                        //!< Pthread communicator
#if PAPI
  int                    PAPIEventSet = PAPI_NULL;              //!< PAPI event set
  std::vector<char*>     PAPIEventNames;                        //!< Vector of PAPI event names
  std::vector<int>       PAPIEventCodes;                        //!< Vector of PAPI event codes
  std::vector<int64_t>   PAPIEventValues;                       //!< Vector of PAPI event values
#endif

  int stringLength = 20;                                        //!< Max length of event name
  int decimal = 7;                                              //!< Decimal precision
  bool verbose = false;                                         //!< Print to screen

  //! Timer function
  double get_time() {
    struct timeval tv;                                          // Time value
    gettimeofday(&tv, NULL);                                    // Get time of day in seconds and microseconds
    return double(tv.tv_sec+tv.tv_usec*1e-6);                   // Combine seconds and microseconds and return
  }

  //! Cycle counter
  inline uint64_t get_cycle() {
    uint32_t low = 0, high = 0;                                 // Define low and high 32 bits of cycle counter
#ifndef __FUJITSU
    asm volatile ("rdtsc" : "=a" (low), "=d" (high));           // Call rdtsc
#endif
    return (uint64_t(high) << 32) | uint64_t(low);              // Return 64 bit cycle counter
  }

  //! Cycle counter with thread ID
  inline uint64_t get_cycle(uint32_t * id) {
    uint32_t low = 0, high = 0;                                 // Define low and high 32 bits of cycle counter
    if (!id) return 0;                                          // Count only for valid thread ID
#ifndef __FUJITSU
    asm volatile ("rdtscp" : "=a" (low), "=d" (high), "=c" (*id));// Call rdtscp
#endif
    return (uint64_t(high) << 32) | uint64_t(low);              // Return 64 bit cycle counter
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

  //! Write timings of all events
  inline void writeTime(int mpirank=0) {
    std::stringstream name;                                     // File name
    name << "time" << std::setfill('0') << std::setw(6)         // Set format
         << mpirank << ".dat";                                  // Create file name for timer
    std::ofstream timerFile(name.str().c_str(), std::ios::app); // Open timer log file
    for (T_iter E=timer.begin(); E!=timer.end(); E++) {         // Loop over all events
      timerFile << std::setw(stringLength) << std::left         //  Set format
		<< E->first << " " << E->second << std::endl;   //  Print event and timer
    }                                                           // End loop over all events
    timerFile.close();                                          // Close timer log file
    timer.clear();                                              // Clear timer
  }

  //! Erase single event in timer
  inline void resetTimer(std::string event) {
    timer.erase(event);                                         // Erase event from timer
  }

  //! Erase all events in timer
  inline void resetTimer() {
    timer.clear();                                              // Clear timer
  }

  //! Start PAPI event
  inline void startPAPI() {
#if PAPI
    PAPI_library_init(PAPI_VER_CURRENT);                        // Initialize PAPI library
    char * allEvents = getenv("EXAFMM_PAPI_EVENTS");            // Get all PAPI event strings
    char eventName[256];                                        // PAPI event name
    while (allEvents) {                                         // While event string is not empty
      char * event = strchr(allEvents, ',');                    //  Get single event string
      int n = (event == NULL ? (int)strlen(allEvents) : event - allEvents);// Count string length
      int eventCode;                                            //  PAPI event code
      snprintf(eventName, n+1, "%s", allEvents);                //  Get PAPI event name
      if (PAPI_event_name_to_code(eventName, &eventCode) == PAPI_OK) {// Event name to event code
        PAPIEventNames.push_back(strdup(eventName));            //   Push event name to vector
        PAPIEventCodes.push_back(eventCode);                    //   Push event code to vector
      }                                                         //  End if for event name to event code
      if (event == NULL) break;                                 //  Stop if event string is empty
      else allEvents = event + 1;                               //  Else move to next event string
    };                                                          // End while loop for event string
    if (!PAPIEventCodes.empty()) {                              // If PAPI events are set
      PAPI_create_eventset(&PAPIEventSet);                      // Create PAPI event set
      for (int i=0; i<int(PAPIEventCodes.size()); i++) {        // Loop over PAPI events
	PAPI_add_event(PAPIEventSet, PAPIEventCodes[i]);        //  Add PAPI event
      }                                                         // End loop over PAPI events
      PAPI_start(PAPIEventSet);                                 // Start PAPI counter
    }
#endif
  }

  //! Stop PAPI event
  inline void stopPAPI() {
#if PAPI
    if (!PAPIEventCodes.empty()) {                              // If PAPI events are set
      PAPIEventValues.resize(PAPIEventCodes.size());            //  Resize PAPI event value vector
      PAPI_stop(PAPIEventSet, &PAPIEventValues[0]);             //  Stop PAPI counter
    }                                                           // End if for PAPI events
#endif
  }

  //! Print PAPI event
  inline void printPAPI() {
#if PAPI
    if (!PAPIEventCodes.empty() && verbose) {                   // If PAPI events are set and verbose is true
      printTitle("PAPI stats ");
      for (int i=0; i<int(PAPIEventCodes.size()); i++) {        //  Loop over PAPI events
        std::cout << std::setw(stringLength) << std::left       //   Set format
		  << PAPIEventNames[i] << " : " << std::setprecision(decimal) << std::fixed
		  << PAPIEventValues[i] << std::endl;           //   Print event and timer
      }                                                         //  End loop over PAPI events
    }                                                           // End if for PAPI events
#endif
  }

#if TRACE
  //! Initialize tracer
  inline void initTracer() {
    pthread_mutex_init(&mutex,NULL);                            // Initialize pthread communicator
  }

  //! Start tracer for given event
  inline void startTracer(Tracer &tracer) {
    pthread_mutex_lock(&mutex);                                 // Lock shared variable access
    tracer.thread = pthread_self();                             // Store pthread id
    tracer.begin  = get_time();                                 // Start timer
    pthread_mutex_unlock(&mutex);                               // Unlock shared variable access
  }

  //! Stop tracer for given event
  inline void stopTracer(Tracer & tracer) {
    pthread_mutex_lock(&mutex);                                 // Lock shared variable access
    tracer.end = get_time();                                    // Stop timer
    tracers.push(tracer);                                       // Push tracer to queue of tracers
    pthread_mutex_unlock(&mutex);                               // Unlock shared variable access
  }

  //! Write tracers of all events
  inline void writeTracer(int mpirank=0) {
    startTimer("Write tracer");                                 // Start timer
    std::stringstream name;                                     // File name
    name << "trace" << std::setfill('0') << std::setw(6)        // Set format
         << mpirank << ".svg";                                  // Create file name for tracer
    std::ofstream traceFile(name.str().c_str());                // Open tracer log file
    traceFile << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n" // Header statements for tracer log file
	      << "<!DOCTYPE svg PUBLIC \"-_W3C_DTD SVG 1.0_EN\" \"http://www.w3.org/TR/SVG/DTD/svg10.dtd\">\n"
	      << "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n"
	      << "  width=\"200mm\" height=\"40mm\" viewBox=\"0 0 20000 4000\">\n"
	      << "  <g>\n";
    int num_thread = 0;                                         // Counter for number of threads to trace
    ThreadMap threadMap;                                        // Map pthread ID to thread ID
    double base = tracers.front().begin;                        // Base time
    double scale = 30000.0;                                     // Scale the length of bar plots
    while (!tracers.empty()) {                                  // While queue of traces is not empty
      Tracer tracer = tracers.front();                          //  Get tracer at front of the queue
      tracers.pop();                                            //  Pop tracer at front
      pthread_t thread = tracer.thread;                         //  Set pthread ID of tracer
      double begin  = tracer.begin;                             //  Set begin time of tracer
      double end    = tracer.end;                               //  Set end time of tracer
      int    color  = 0x0000ff;                                 //  Set color of tracer
      if (threadMap[thread] == 0) {                             //  If it's a new pthread ID
        threadMap[thread] = num_thread++;                       //   Map it to an incremented thread ID
      }                                                         //  End if for new pthread ID
      begin -= base;                                            //  Subtract base time from begin time
      end   -= base;                                            //  Subtract base time from end time
      traceFile << "    <rect x=\"" << begin * scale            //  x position of bar plot
		<< "\" y=\"" << threadMap[thread] * 100.0       //  y position of bar plot
		<< "\" width=\"" << (end - begin) * scale       //  width of bar
		<< "\" height=\"90.0\" fill=\"#"<< std::setfill('0')
		<< std::setw(6) << std::hex << color            // height of bar
		<< "\" stroke=\"#000000\" stroke-width=\"1\"/>\n";//  stroke color and width
    }                                                           // End while loop for queue of tracers
    traceFile << "  </g>\n" "</svg>\n";                         // Footer for tracer log file
    traceFile.close();                                          // Close tracer log file
    stopTimer("Write tracer",verbose);                          // Stop timer
  }
#else
  inline void initTracer() {}
  inline void startTracer(Tracer) {}
  inline void stopTracer(Tracer) {}
  inline void writeTracer() {}
  inline void writeTracer(int) {}
#endif

#if DAG_RECORDER == 2
  //! Start DAG recorder
  inline void startDAG() {
    dr_start(0);                                                // Start DAG recorder
  }

  //! Stop DAG recorder
  inline void stopDAG() {
    dr_stop();                                                  // Stop DAG recorder
  }

  //! Write DAG to file
  inline void writeDAG() {
    dr_dump();                                                  // Write DAG to file
  }
#else
  inline void startDAG() {}
  inline void stopDAG() {}
  inline void writeDAG() {}
#endif
};
#endif
