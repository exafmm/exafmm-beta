#ifndef logger_h
#define logger_h
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <pthread.h>
#include <queue>
#include <string>
#include <sstream>
#include <sys/time.h>
#include <vector>

#if PAPI
#include <cstring>
#include <papi.h>
#endif

//! Structure for pthread based trace
struct Trace {
  pthread_t thread;                                             //!< pthread id
  double    begin;                                              //!< Begin timer of trace
  double    end;                                                //!< End timer of trace
  Trace() {}                                                    //!< Constructor
};

//! Timer and Trace logger
class Logger {
 typedef std::map<std::string,double>           Timer;          //!< Map of timer event name to timed value
 typedef std::map<std::string,double>::iterator T_iter;         //!< Iterator of timer event name map
 typedef std::queue<Trace>                      Traces;         //!< Queue of traces
 typedef std::map<pthread_t,int>                ThreadMap;      //!< Map of pthread id to thread id

 private:
  Timer           beginTimer;                                   //!< Timer base value
  Timer           timer;                                        //!< Timings of all events
  Traces          traces;                                       //!< Traces for all events
  pthread_mutex_t mutex;                                        //!< Pthread communicator
#if PAPI
  int                    PAPIEventSet;                          //!< PAPI event set
  std::vector<char*>     PAPIEventNames;                        //!< Vector of PAPI event names
  std::vector<int>       PAPIEventCodes;                        //!< Vector of PAPI event codes
  std::vector<long long> PAPIEventValues;                       //!< Vector of PAPI event values
#endif

 public:
  int stringLength;                                             //!< Max length of event name
  int decimal;                                                  //!< Decimal precision
  bool verbose;                                                 //!< Print to screen

 private:
//! Timer function
  double get_time() const {
    struct timeval tv;                                          // Time value
    gettimeofday(&tv, NULL);                                    // Get time of day in seconds and microseconds
    return double(tv.tv_sec+tv.tv_usec*1e-6);                   // Combine seconds and microseconds and return
  }

 public:
//! Constructor
  Logger() : beginTimer(), timer(), traces(), mutex(),          // Initializing class variables (empty)
#if PAPI
    PAPIEventSet(PAPI_NULL),                                    // Initializing PAPI event set
#endif
    stringLength(20),                                           // Max length of event name
    decimal(7),                                                 // Decimal precision
    verbose(false) {                                            // Don't print timings by default
    pthread_mutex_init(&mutex,NULL);                            // Initialize pthread communicator
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

//! Stop timer for given event
  double stopTimer(std::string event) {
    double endTimer = get_time();                               // Get time of day and store in endTimer
    timer[event] += endTimer - beginTimer[event];               // Accumulate event time to timer
    if (verbose) printTime(event);                              // Print event and timer to screen
    return endTimer - beginTimer[event];                        // Return the event time
  }

//! Print timings of a specific event
  inline void printTime(std::string event) {
    if (verbose) {                                              // If verbose flag is true
      std::cout << std::setw(stringLength) << std::left         //  Set format
        << event << " : " << std::setprecision(decimal) << std::fixed
        << timer[event] << " s" << std::endl;                   //  Print event and timer
    }                                                           // End if for verbose flag
  }

//! Write timings of all events
  inline void writeTime(int mpirank=0) {
    std::stringstream name;                                     // File name
    name << "time" << std::setfill('0') << std::setw(6)         // Set format
         << mpirank << ".dat";                                  // Create file name for timer
    std::ofstream timerFile(name.str().c_str(), std::ios::app); // Open timer log file
    for (T_iter E=timer.begin(); E!=timer.end(); E++) {         // Loop over all events
      timerFile << std::setw(stringLength) << std::left         //  Set format
        << E->first << " " << E->second << std::endl;           //  Print event and timer
    }                                                           // End loop over all events
    timerFile.close();                                          // Close timer log file
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
//! Start tracer for given event
  inline void startTracer(Trace &trace) {
    pthread_mutex_lock(&mutex);                                 // Lock shared variable access
    trace.thread = pthread_self();                              // Store pthread id
    trace.begin  = get_time();                                  // Start timer
    pthread_mutex_unlock(&mutex);                               // Unlock shared variable access
  }

//! Stop tracer for given event
  inline void stopTracer(Trace &trace) {
    pthread_mutex_lock(&mutex);                                 // Lock shared variable access
    trace.end    = get_time();                                  // Stop timer
    traces.push(trace);                                         // Push trace to queue of traces
    pthread_mutex_unlock(&mutex);                               // Unlock shared variable access
  }

//! Write traces of all events
  inline void writeTrace(int mpirank=0) {
    startTimer("Write trace");                                  // Start timer
    std::stringstream name;                                     // File name
    name << "trace" << std::setfill('0') << std::setw(6)        // Set format
         << mpirank << ".svg";                                  // Create file name for trace
    std::ofstream traceFile(name.str().c_str());                // Open trace log file
    traceFile << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n" // Header statements for trace log file
      << "<!DOCTYPE svg PUBLIC \"-_W3C_DTD SVG 1.0_EN\" \"http://www.w3.org/TR/SVG/DTD/svg10.dtd\">\n"
      << "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n"
      << "  width=\"200mm\" height=\"40mm\" viewBox=\"0 0 20000 4000\">\n"
      << "  <g>\n";
    int num_thread = 0;                                         // Counter for number of threads to trace
    ThreadMap threadMap;                                        // Map pthread ID to thread ID
    double base = traces.front().begin;                         // Base time
    double scale = 30000.0;                                     // Scale the length of bar plots
    while (!traces.empty()) {                                   // While queue of traces is not empty
      Trace trace = traces.front();                             //  Get trace at front of the queue
      traces.pop();                                             //  Pop trace at front
      pthread_t thread = trace.thread;                          //  Get pthread ID of trace
      double begin  = trace.begin;                              //  Get begin time of trace
      double end    = trace.end;                                //  Get end time of trace
      int    color  = 0x0000ff;                                 //  Get color of trace
      if (threadMap[thread] == 0) {                             //  If it's a new pthread ID
        threadMap[thread] = num_thread++;                       //   Map it to an incremented thread ID
      }                                                         //  End if for new pthread ID
      begin -= base;                                            //  Subtract base time from begin time
      end   -= base;                                            //  Subtract base time from end time
      traceFile << "    <rect x=\"" << begin * scale            //  x position of bar plot
        << "\" y=\"" << threadMap[thread] * 100.0               //  y position of bar plot
        << "\" width=\"" << (end - begin) * scale               //  width of bar
        << "\" height=\"90.0\" fill=\"#"<< std::setfill('0') << std::setw(6) << std::hex << color// height of bar
        << "\" stroke=\"#000000\" stroke-width=\"1\"/>\n";      //  stroke color and width
    }                                                           // End while loop for queue of traces
    traceFile << "  </g>\n" "</svg>\n";                         // Footer for trace log file
    traceFile.close();                                          // Close trace log file
    stopTimer("Write trace",verbose);                           // Stop timer
  }
#else
  inline void startTracer(Trace) {}
  inline void stopTracer(Trace) {}
  inline void writeTrace() {}
  inline void writeTrace(int) {}
#endif

  //! Print relative L2 norm error
  void printError(double diff1, double norm1) {
    if (verbose) {                                              // If verbose flag is true
      std::cout << std::setw(stringLength) << std::left << std::scientific//  Set format
                << "Rel. L2 Error (pot)" << " : " << std::sqrt(diff1/norm1) << std::endl;// Print potential error
    }                                                           // End if for verbose flag
  }
};
#endif
