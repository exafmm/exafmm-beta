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

#if PAPI_EX
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <papi.h>
#elif PAPI
#include <papi.h>
#endif

//! Structure for pthread based trace
struct Trace {
  pthread_t thread;                                             //!< pthread id
  double    begin;                                              //!< Begin timer of trace
  double    end;                                                //!< End timer of trace
  int       color;                                              //!< Color of trace
};
typedef std::queue<Trace>                      Traces;          //!< Queue of traces
typedef std::map<pthread_t,double>             ThreadTrace;     //!< Map of pthread id to traced value
typedef std::map<pthread_t,int>                ThreadMap;       //!< Map of pthread id to thread id
typedef std::map<std::string,double>           Timer;           //!< Map of timer event name to timed value
typedef std::map<std::string,double>::iterator T_iter;          //!< Iterator of timer event name map

//! Timer and Trace logger
class Logger {
 private:
  std::ofstream   timerFile;                                    //!< File ID to store log
  Timer           beginTimer;                                   //!< Timer base value
  Timer           timer;                                        //!< Stores timings for all events
  Traces          traces;                                       //!< Stores traces for all events
  pthread_mutex_t mutex;                                        //!< Pthread communicator
#if PAPI_EX
  int             nPAPIEvents;
  char         ** PAPIEventNames;
  int           * PAPIEventCodes;
  long long     * PAPIEventValues;
  int             PAPIEventSet;
#elif PAPI
  int PAPIEVENT;                                                //!< PAPI event handle
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

#if PAPI_EX
  /* papi related aux functions */
  void showPAPIError(int retval, const char * msg) {
    fprintf(stderr, "error in PAPI call %s (retval %d): %s\n", 
	    msg, retval, PAPI_strerror(retval));
  }

  /* events : a string like "L1_DCM,L2_DCM,L3_DCM"
     PAPIEventCodes : an array that can have max_events
     parse events into an array of event codes.
     return the number of events successfully parsed */
  
#define MAX_PAPI_EVENT_NAME 256
  int parsePAPIEvents(char * events_string) {
    char event_name[MAX_PAPI_EVENT_NAME];
    
    /* split events_string by ',' and check if each event
       name is supported. count the number of valid events
       into nPAPIEvents */
    char * p = events_string;
    int ne = 0;
    while (p) {
      char * q = strchr(p, ',');
      int n = (q == NULL ? (int)strlen(p) : q - p);
      unsigned int to_write = n + 1;
      if (to_write <= sizeof(event_name)) {
	int event_code;
	snprintf(event_name, to_write, "%s", p);
	int retval = PAPI_event_name_to_code(event_name, &event_code);
	if (retval == PAPI_OK) {
	  ne++;		/* OK, this is a vaild event */
	} else {
	  fprintf(stderr, 
		  "warning: could not find event %s; the papi error follows:\n", 
		  p);
	  showPAPIError(retval, "PAPI_event_name_to_code");
	}
      } else {
	char * x = strndup(p, n);
	fprintf(stderr, "warning: event name '%s' too long (ignored)\n", x);
	free(x);
      }
      if (q == NULL) break;
      else p = q + 1;
    };
    
    if (ne == 0) return 0;

    nPAPIEvents = ne;
    PAPIEventNames = new char * [nPAPIEvents];
    PAPIEventCodes = new int [nPAPIEvents];
    PAPIEventValues = new long long [nPAPIEvents];
    PAPIEventSet = 0;

    for (int k = 0; k < nPAPIEvents; k++) {
      PAPIEventNames[k] = NULL;
      PAPIEventCodes[k] = 0;
      PAPIEventValues[k] = 0;
    }
    
    p = events_string;
    int i = 0;
    while (p) {
      char * q = strchr(p, ',');
      int n = (q == NULL ? (int)strlen(p) : q - p);
      unsigned int to_write = n + 1;
      if (to_write <= sizeof(event_name)) {
	int event_code;
	snprintf(event_name, to_write, "%s", p);
	int retval = PAPI_event_name_to_code(event_name, &event_code);
	if (retval == PAPI_OK) {
	  assert(i < nPAPIEvents);
	  PAPIEventNames[i] = strdup(event_name);
	  PAPIEventCodes[i] = event_code;
	  i++;
	}
      }
      if (q == NULL) break;
      else p = q + 1;
    };
    assert(i == nPAPIEvents);
    return nPAPIEvents;		/* event_count */
  }

  int PAPIEventCodesToEventSet() {
    int es = PAPI_NULL;
    int retval = PAPI_create_eventset(&es);
    if (retval != PAPI_OK) {
      showPAPIError(retval, "PAPI_create_eventset");
      return PAPI_NULL;
    }
    for (int i = 0; i < nPAPIEvents; i++) {
      int retval0 = PAPI_add_event(es, PAPIEventCodes[i]);
      if (retval0 != PAPI_OK) {
	fprintf(stderr, "couldn't add event %s (code=%d)\n", 
		PAPIEventNames[i], PAPIEventCodes[i]);
	showPAPIError(retval, "PAPI_add_event");
	return PAPI_NULL;
      }
    }
    return es; 			/* OK */
  }

  void showPAPICounts() {
    for (int i = 0; i < nPAPIEvents; i++) {
      printf("%s : %lld\n", PAPIEventNames[i], PAPIEventValues[i]);
    }
  }

  int initPAPI() {
    int ver = PAPI_library_init(PAPI_VER_CURRENT);
    if (ver != PAPI_VER_CURRENT) { 
      fprintf(stderr, "error: could not init papi\n");
      return 0;			// NG
    }
#if _OPENMP
    if (PAPI_thread_init(pthread_self) != PAPI_OK) {
      fprintf(stderr, "error could not init papi\n");
      return 0;			// NG
    }
#endif
    return 1;
  }
#endif

 public:
//! Constructor
  Logger() : timerFile("time.dat"),                             // Open timer log file
    beginTimer(), timer(), traces(), mutex(),                   // Initializing class variables (empty)
#if PAPI
    PAPIEVENT(PAPI_NULL),                                       // PAPI event handle
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
      << timer[event] << " s" << std::endl;   // Print event and timer
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
    for (T_iter E=timer.begin(); E!=timer.end(); E++) {        // Loop over all events
      timerFile << std::setw(stringLength) << std::left         //  Set format
        << E->first << " " << E->second << std::endl;           //  Print event and timer
    }                                                           // End loop over all events
  }

//! Start PAPI event
  inline void startPAPI() {
#if PAPI_EX
    /* set nPAPIEvents, PAPIEventNames, and PAPIEventCodes  */
    if (initPAPI()) {
      nPAPIEvents = parsePAPIEvents(getenv("EXAFMM_PAPI_EVENTS"));
      if (nPAPIEvents) {
	PAPIEventSet = PAPIEventCodesToEventSet();
	PAPI_start(PAPIEventSet);
      }
    }
#elif PAPI
    int events[3] = { PAPI_L2_DCM, PAPI_L2_DCA, PAPI_TLB_DM };  // PAPI event type
    PAPI_library_init(PAPI_VER_CURRENT);                        // PAPI initialize
    PAPI_create_eventset(&PAPIEVENT);                           // PAPI create event set
    PAPI_add_events(PAPIEVENT, events, 3);                      // PAPI add events
    PAPI_start(PAPIEVENT);                                      // PAPI start
#endif
  }

//! Stop PAPI event
  inline void stopPAPI() {
#if PAPI_EX
    if (nPAPIEvents) {
      PAPI_stop(PAPIEventSet, PAPIEventValues);
      showPAPICounts();
    }
#elif PAPI
    long long values[3] = {0,0,0};                              // Values for each event
    PAPI_stop(PAPIEVENT,values);                                // PAPI stop
    std::cout << "--- PAPI stats -------------------" << std::endl
      << std::setw(stringLength) << std::left                   // Set format
      << "L2 Miss"    << " : " << values[0] << std::endl        // Print L2 Misses
      << std::setw(stringLength) << std::left                   // Set format
      << "L2 Access"  << " : " << values[1] << std::endl        // Print L2 Access
      << std::setw(stringLength) << std::left                   // Set format
      << "TLB Miss"   << " : " << values[2] << std::endl        // Print TLB Misses
      << "--- PAPI stats -------------------" << std::endl;
#endif
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
  inline void writeTrace(int mpirank) {
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
      int    color  = trace.color;                              //  Get color of trace
      if (threadMap[thread] == 0) {                             //  If it's a new pthread ID
        threadMap[thread] = ++num_thread;                       //   Map it to an incremented thread ID
      }                                                         //  End if for new pthread ID
      begin -= base;                                            //  Subtract base time from begin time
      end   -= base;                                            //  Subtract base time from end time
      traceFile << "    <rect x=\"" << begin * scale            //  x position of bar plot
        << "\" y=\"" << (threadMap[thread] - 1) * 100.0         //  y position of bar plot
        << "\" width=\"" << (end - begin) * scale               //  width of bar
        << "\" height=\"90.0\" fill=\"#"<< std::setfill('0') << std::setw(6) << std::hex << color// height of bar
        << "\" stroke=\"#000000\" stroke-width=\"1\"/>\n";      //  stroke color and width
    }                                                           // End while loop for queue of traces
    traceFile << "  </g>\n" "</svg>\n";                         // Footer for trace log file 
    traceFile.close();                                          // Close trace log file
  }
};
#endif
