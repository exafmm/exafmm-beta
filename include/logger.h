#ifndef logger_h
#define logger_h
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <pthread.h>
#include <queue>
#include <string>
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
};

//! Timer and Trace logger
class Logger {
 typedef std::map<std::string,double>           Timer;          //!< Map of timer event name to timed value
 typedef std::map<std::string,double>::iterator T_iter;         //!< Iterator of timer event name map
 typedef std::queue<Trace>                      Traces;         //!< Queue of traces
 typedef std::map<pthread_t,int>                ThreadMap;      //!< Map of pthread id to thread id

 private:
  std::ofstream   timerFile;                                    //!< File ID to store log
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
  bool printNow;                                                //!< Switch to print timings

 private:
//! Timer function
  double get_time() const {
    struct timeval tv;                                          // Time value
    gettimeofday(&tv, NULL);                                    // Get time of day in seconds and microseconds
    return double(tv.tv_sec+tv.tv_usec*1e-6);                   // Combine seconds and microseconds and return
  }

 public:
//! Constructor
  Logger() : timerFile("time.dat"),                             // Open timer log file
    beginTimer(), timer(), traces(), mutex(),                   // Initializing class variables (empty)
#if PAPI
    PAPIEventSet(PAPI_NULL),                                    // Initializing PAPI event set
#endif
    stringLength(20),                                           // Max length of event name
    decimal(7),                                                 // Decimal precision
    printNow(false) {                                           // Don't print timings by default
    pthread_mutex_init(&mutex,NULL);                            // Initialize pthread communicator
  }
//! Destructor
  ~Logger() {
    timerFile.close();                                          // Close timer log file
  }

//! Start timer for given event
  inline void startTimer(std::string event) {
    beginTimer[event] = get_time();                             // Get time of day and store in beginTimer
  }

//! Stop timer for given event
  double stopTimer(std::string event, bool print=false) {
    double endTimer = get_time();                               // Get time of day and store in endTimer
    timer[event] += endTimer - beginTimer[event];               // Accumulate event time to timer
    if (print) printTime(event);                                // Print event and timer to screen
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
    std::cout << std::setw(stringLength) << std::left           // Set format
      << event << " : " << std::setprecision(decimal) << std::fixed
      << timer[event] << " s" << std::endl;                     // Print event and timer
  }

//! Print timings of all events
  inline void printAllTime() {
    for (T_iter E=timer.begin(); E!=timer.end(); E++) {         // Loop over all events
      std::cout << std::setw(stringLength) << std::left         //  Set format
        << E->first << " : " << E->second << std::endl;         //  Print event and timer
    }                                                           // End loop over all events
  }

//! Write timings of all events
  inline void writeTime() {
    for (T_iter E=timer.begin(); E!=timer.end(); E++) {         // Loop over all events
      timerFile << std::setw(stringLength) << std::left         //  Set format
        << E->first << " " << E->second << std::endl;           //  Print event and timer
    }                                                           // End loop over all events
  }

//! Start PAPI event
  inline void startPAPI() {
#if PAPI
    PAPI_library_init(PAPI_VER_CURRENT);                        // Initialize PAPI library
#if _OPENMP
    PAPI_thread_init(pthread_self);                             // Initialize PAPI thread
#endif
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
    if (!PAPIEventCodes.empty()) {                              // If PAPI events are set
      std::cout << "--- PAPI stats -------------------" << std::endl;
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
    char fname[256];                                            // File name
    sprintf(fname,"trace%4.4d.svg",mpirank);                    // Create file name for trace
    std::ofstream traceFile(fname);                             // Open trace log file
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
  }
#else
  inline void startTracer(Trace) {}
  inline void stopTracer(Trace) {}
  inline void writeTrace() {}
  inline void writeTrace(int) {}
#endif

};
#endif
