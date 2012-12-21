/*
    Modified by T.Narumi RIKEN  Sep 5, 2005
             by T.Koishi RIKEN
*/
/*
 *  This file is provided for use with the unix-socket-faq.  It is public
 *  domain, and may be copied freely.  There is no copyright on it.  The
 *  original work was by Vic Metcalfe (vic@brutus.tlug.org), and any
 *  modifications made to that work were made with the understanding that
 *  the finished work would be in the public domain.
 *
 *  If you have found a bug, please pass it on to me at the above address
 *  acknowledging that there will be no copyright on your work.
 *
 *  The most recent version of this file, and the unix-socket-faq can be
 *  found at http://www.interlog.com/~vic/sock-faq/.
 */
#ifndef WIN32
#include "sockhelp.h"

#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <ctype.h>
#include <string.h>

#define SOCKET int
#define INVALID_SOCKET -1
#define CLOSE(x) close(x)
//#include "cfortran.h"
#else
#include <windows.h>
#include <stdio.h>
#define CLOSE(x) closesocket(x)
#endif

//#include "mr3.h"

//  int listensock = -1; /* So that we can close sockets on ctrl-c */
SOCKET listensock = INVALID_SOCKET;

int STRLEN=0;
  int new_orig = 0;

int MR3_close_ch(SOCKET s_sock)
{
    if (listensock!=INVALID_SOCKET) CLOSE(listensock);
    if (s_sock!=INVALID_SOCKET)     CLOSE(s_sock);
    return 0;
}
int mr3_close_ch_(SOCKET *s_sock)
{
    return MR3_close_ch(*s_sock);
}
#ifndef WIN32
/* This waits for all children, so that they don't become zombies. */
static void sig_chld(int signal_type)
{
  int pid;
  int status;

  while ( (pid = wait3(&status, WNOHANG, NULL)) > 0);
}


/* This ignores the SIGPIPE signal.  This is usually a good idea, since
   the default behaviour is to terminate the application.  SIGPIPE is
   sent when you try to write to an unconnected socket.  You should
   check your return codes to make sure you catch this error! */
static void ignore_pipe(void)
{
  struct sigaction sig;

  sig.sa_handler = SIG_IGN;
  sig.sa_flags = 0;
  sigemptyset(&sig.sa_mask);
  sigaction(SIGPIPE,&sig,NULL);
}
#endif // WIN32


/* Take a service name, and a service type, and return a port number.  If the
   service name is not found, it tries it as a decimal number.  The number
   returned is byte ordered for the network. */
static int atoport(service, proto)
char *service;
char *proto;
{
  int port;
  long int lport;
  struct servent *serv;
  char *errpos;

  /* First try to read it from /etc/services */
  serv = getservbyname(service, proto);
  if (serv != NULL)
    port = serv->s_port;
  else { /* Not in services, maybe a number? */
    lport = strtol(service,&errpos,0);
    if ( (errpos[0] != 0) || (lport < 1) || (lport > 65535) )
      return -1; /* Invalid port address */
    port = htons(lport);
  }
  return port;
}


/* This function listens on a port, and returns connections.  It forks
   returns off internally, so your main function doesn't have to worry
   about that.  This can be confusing if you don't know what is going on.
   The function will create a new process for every incoming connection,
   so in the listening process, it will never return.  Only when a connection
   comes in, and we create a new process for it will the function return.
   This means that your code that calls it should _not_ loop.

   The parameters are as follows:
     socket_type: SOCK_STREAM or SOCK_DGRAM (TCP or UDP sockets)
     port: The port to listen on.  Remember that ports < 1024 are
       reserved for the root user.  Must be passed in network byte
       order (see "man htons").
     listener: This is a pointer to a variable for holding the file
       descriptor of the socket which is being used to listen.  It
       is provided so that you can write a signal handler to close
       it in the event of program termination.  If you aren't interested,
       just pass NULL.  Note that all modern unixes will close file
       descriptors for you on exit, so this is not required. */
static SOCKET get_connection(socket_type, port, listener)
int socket_type;
u_short port;
int *listener;
{
  struct sockaddr_in address;
  SOCKET listening_socket;
  SOCKET connected_socket = INVALID_SOCKET;
  int new_process;
  int reuse_addr = 1;

  /* Setup internet address information.  
     This is used with the bind() call */
  memset((char *) &address, 0, sizeof(address));
  address.sin_family = AF_INET;
  address.sin_port = port;
  address.sin_addr.s_addr = htonl(INADDR_ANY);

  listening_socket = socket(AF_INET, socket_type, 0);
  if (listening_socket==INVALID_SOCKET) {
#ifndef WIN32
    perror("get_connection: socket");
#else
    fprintf(stderr, "get_connection::socket: %d\n", WSAGetLastError());
#endif
    MR3_exit(EXIT_FAILURE);
  }
  
  if (listener != NULL)
    *listener = listening_socket;

  setsockopt(listening_socket, SOL_SOCKET, SO_REUSEADDR, (char *)&reuse_addr, 
    sizeof(reuse_addr));

  if (bind(listening_socket, (struct sockaddr *) &address, 
    sizeof(address))) {
#ifndef WIN32
    perror("get_connection: bind");
#else
    fprintf(stderr, "get_connection::bind: %d\n", WSAGetLastError());
#endif
    CLOSE(listening_socket);
    MR3_exit(EXIT_FAILURE);
  }

  if (socket_type == SOCK_STREAM) {
    fprintf(stderr,"wait connection\n");
    listen(listening_socket, 1); /* Queue up to five connections before
                                  having them automatically rejected. */

    while(connected_socket == INVALID_SOCKET ) {
      connected_socket = accept(listening_socket, NULL, NULL);
      if (connected_socket == INVALID_SOCKET ) {
        /* Either a real error occured, or blocking was interrupted for
           some reason.  Only abort execution if a real error occured. */
#ifndef WIN32
        if (errno != EINTR) {
#else
        if (WSAGetLastError() != WSAEINTR) {
#endif
          perror("accept");
          CLOSE(listening_socket);
          MR3_exit(EXIT_FAILURE);
        } else {
          continue;    /* don't fork - do the accept again */
        }
      }

#ifndef WIN32
#ifdef FORK
      new_process = fork();
      if (new_process < 0) {
        perror("fork");
        close(connected_socket);
        connected_socket = -1;
      }
      else { /* We have a new process... */
        if (new_process == 0) {
          /* This is the new process. */
          close(listening_socket); /* Close our copy of this socket */
	  if (listener != NULL){
	          *listener = -1; /* Closed in this process.  We are not 
	 			     responsible for it. */
		  new_orig=1;
	  }
        }
        else {
          /* This is the main loop.  Close copy of connected socket, and
             continue loop. */
          close(connected_socket);
# ifdef LOOP_LISTEN
          connected_socket = -1;
# endif
	  new_orig=0;
	  {
	    int pid;
	    int status;
	    fprintf(stderr,"wait frok process\n");
	    while ( (pid = wait3(&status, WNOHANG, NULL)) == 0);
	    while ( (pid = wait3(&status, WNOHANG, NULL)) > 0){
	    }
	    fprintf(stderr,"frok process done\n");
	  }
# ifndef LOOP_LISTEN
	  MR3_exit(0);
# endif
        }
      }
#endif //FORK
#endif //WIN32
    }
    return connected_socket;
  }
  else
    return listening_socket;
}


SOCKET MR3_open_port(int port_number)
{
  char s_port_num[10];
  int s_port = -1;
#ifndef WIN32
  struct sigaction s_act, s_oldact;
#endif
  
  int listensock = -1; /* So that we can close sockets on ctrl-c */

#ifdef WIN32
  // Initialize Winsock
 {
   WORD winsockVersionRequired = MAKEWORD(2,0);
   WSADATA wsadata;
   if (WSAStartup(winsockVersionRequired, &wsadata)) {
     fprintf(stderr, "MR3_open_port: Failed to initialize winsock dll.\n");
     MR3_exit(EXIT_FAILURE);
   }
 }
#endif // WIN32

#if 1   /* added by Tetsu Narumi */
  /*  printf("%d\n",port_number);*/
  if(getenv("MDVIS_PORT_NUM")!=NULL){
    port_number=atoi(getenv("MDVIS_PORT_NUM"));
  }
  /*printf("%d\n",port_number);*/
#endif

  sprintf(s_port_num,"%d",port_number);

#ifndef WIN32
  ignore_pipe();
  sigemptyset(&s_act.sa_mask);
  s_act.sa_flags = 0;
  s_act.sa_handler = sig_chld;
  sigaction(SIGCHLD, &s_act, &s_oldact);
#endif

  s_port = atoport(s_port_num, "tcp");
  if (s_port == -1) {
    fprintf(stderr,"Unable to find service: %s\n",s_port_num);
    MR3_exit(EXIT_FAILURE);
  }
  return get_connection(SOCK_STREAM, s_port, &listensock);
}

/* Converts ascii text to in_addr struct.  NULL is returned if the address
   can not be found. */
static struct in_addr *atoaddr(address)
char *address;
{
  struct hostent *host;
  static struct in_addr saddr;

  /* First try it as aaa.bbb.ccc.ddd. */
  saddr.s_addr = inet_addr(address);
  if (saddr.s_addr != -1) {
    return &saddr;
  }
  host = gethostbyname(address);
  if (host != NULL) {
    return (struct in_addr *) *host->h_addr_list;
  }
  return NULL;
}


/* This is a generic function to make a connection to a given server/port.
   service is the port name/number,
   type is either SOCK_STREAM or SOCK_DGRAM, and
   netaddress is the host name to connect to.
   The function returns the socket, ready for action.*/
static SOCKET make_connection(service, type, netaddress)
char *service;
int type;
char *netaddress;
{
  /* First convert service from a string, to a number... */
  int port = -1;
  struct in_addr *addr;
  int sock, connected;
  struct sockaddr_in address;

  if (type == SOCK_STREAM) 
    port = atoport(service, "tcp");
  if (type == SOCK_DGRAM)
    port = atoport(service, "udp");
  if (port == -1) {
    fprintf(stderr,"make_connection:  Invalid socket type.\n");
    return INVALID_SOCKET;
  }
  addr = atoaddr(netaddress);
  if (addr == NULL) {
    fprintf(stderr,"make_connection:  Invalid network address.\n");
    return INVALID_SOCKET;
  }
 
  memset((char *) &address, 0, sizeof(address));
  address.sin_family = AF_INET;
  address.sin_port = (port);
  address.sin_addr.s_addr = addr->s_addr;

  sock = socket(AF_INET, type, 0);

  printf("Connecting to %s on port %d.\n",inet_ntoa(*addr),htons(port));

  if (type == SOCK_STREAM) {
    connected = connect(sock, (struct sockaddr *) &address, 
      sizeof(address));
    if (connected < 0) {
      perror("connect");
      return INVALID_SOCKET;
    }
    return sock;
  }
  /* Otherwise, must be for udp, so bind to address. */
  if (bind(sock, (struct sockaddr *) &address, sizeof(address))) {
    perror("bind");
    return -1;
  }
  return sock;
}

/* This is just like the read() system call, accept that it will make
   sure that all your data goes through the socket. */
int MR3_sock_read(SOCKET sockfd, char *buf, size_t count)
{
  size_t bytes_read = 0;
  int this_read;

  while (bytes_read < count) {
#ifndef WIN32
    do
      this_read = read(sockfd, buf, count - bytes_read);
    while ( (this_read < 0) && (errno == EINTR) );
#else
    do
      this_read = recv(sockfd, buf, count - bytes_read, 0);
    while ( (this_read < 0) && (WSAGetLastError() == WSAEINTR) );
#endif
    if (this_read < 0)
      return this_read;
    else if (this_read == 0)
      return bytes_read;
    bytes_read += this_read;
    buf += this_read;
  }
  return count;
}

int mr3_sock_read_(SOCKET *sockfd, char* buf, int* count)
{
  return MR3_sock_read(*sockfd, buf, *count);
}
int mr3_sock_read__(SOCKET *sockfd, char* buf, int* count)
{
  return MR3_sock_read(*sockfd, buf, *count);
}

int mr3_sock_read_i_(SOCKET* sockfd, int* buf)
{
  return MR3_sock_read(*sockfd,(char*)buf,(size_t)4);
}
int mr3_sock_read_i__(SOCKET* sockfd, int* buf)
{
  return MR3_sock_read(*sockfd,(char*)buf,(size_t)4);
}

int MR3_sock_read_iv(SOCKET sockfd, int* buf, int count)
{
  return MR3_sock_read(sockfd,(char*)buf,(size_t)count*4);
}

int mr3_sock_read_iv_(SOCKET* sockfd, int* buf, int* count)
{
  return MR3_sock_read(*sockfd,(char*)buf,(size_t)*count*4);
}
int mr3_sock_read_iv__(SOCKET* sockfd, int* buf, int* count)
{
  return MR3_sock_read(*sockfd,(char*)buf,(size_t)*count*4);
}

int MR3_sock_read_d(SOCKET sockfd, double* buf)
{
  return MR3_sock_read(sockfd,(char*)buf,(size_t)8);
}

int mr3_sock_read_d__(SOCKET* sockfd, double* buf)
{
  return MR3_sock_read(*sockfd,(char*)buf,(size_t)8);
}
int mr3_sock_read_d_(SOCKET* sockfd, double* buf)
{
  return MR3_sock_read(*sockfd,(char*)buf,(size_t)8);
}

int MR3_sock_read_dv(SOCKET sockfd, double* buf, int count)
{
  return MR3_sock_read(sockfd,(char*)buf,(size_t)count*8);
}
int mr3_sock_read_dv__(SOCKET* sockfd, double* buf, int* count)
{
  return MR3_sock_read(*sockfd,(char*)buf,(size_t)*count*8);
}
int mr3_sock_read_dv_(SOCKET* sockfd, double* buf, int* count)
{
  return MR3_sock_read(*sockfd,(char*)buf,(size_t)*count*8);
}

/* This is just like the write() system call, accept that it will
   make sure that all data is transmitted. */
int MR3_sock_write(SOCKET sockfd, char *buf, size_t count)
{
  size_t bytes_sent = 0;
  int this_write;

  while (bytes_sent < count) {
#ifndef WIN32
    do
      this_write = write(sockfd, buf, count - bytes_sent);
    while ( (this_write < 0) && (errno == EINTR) );
#else
    do
      this_write = send(sockfd, buf, count - bytes_sent,0);
    while ( (this_write < 0) && (WSAGetLastError() == WSAEINTR) );
#endif
    if (this_write <= 0)
      return this_write;
    bytes_sent += this_write;
    buf += this_write;
  }
  return count;
}


int MR3_sock_write_c(SOCKET sockfd, char* buf)
{
  int count;

  count = strlen(buf);
  return MR3_sock_write(sockfd,buf,(size_t)count);
}

int mr3_sock_write_c_(SOCKET *sockfd, char* buf)
{
  return MR3_sock_write_c(*sockfd,buf);
}

int mr3_sock_write_c__(SOCKET *sockfd, char* buf)
{
  return MR3_sock_write_c(*sockfd,buf);
}

int MR3_sock_write_i(SOCKET sockfd, int buf)
{
  return MR3_sock_write(sockfd,(char*)&buf,(size_t)4);
}

int mr3_sock_write_i_(SOCKET* sockfd, int* buf)
{
  return MR3_sock_write(*sockfd,(char*)buf,(size_t)4);
}
int mr3_sock_write_i__(SOCKET* sockfd, int* buf)
{
  return MR3_sock_write(*sockfd,(char*)buf,(size_t)4);
}

int MR3_sock_write_iv(SOCKET sockfd, int* buf, int count)
{
  return MR3_sock_write(sockfd,(char*)buf,(size_t)count*4);
}

int mr3_sock_write_iv_(SOCKET* sockfd, int* buf, int* count)
{
  return MR3_sock_write(*sockfd,(char*)buf,(size_t)*count*4);
}
int mr3_sock_write_iv__(SOCKET* sockfd, int* buf, int* count)
{
  return MR3_sock_write(*sockfd,(char*)buf,(size_t)*count*4);
}

int MR3_sock_write_dv(SOCKET sockfd, double* buf, int count)
{
  return MR3_sock_write(sockfd,(char*)buf,(size_t)count*8);
}

int mr3_sock_write_dv__(SOCKET* sockfd, double* buf, int* count)
{
  return MR3_sock_write(*sockfd,(char*)buf,(size_t)*count*8);
}
int mr3_sock_write_dv_(SOCKET* sockfd, double* buf, int* count)
{
  return MR3_sock_write(*sockfd,(char*)buf,(size_t)*count*8);
}

/* This function reads from a socket, until it recieves a linefeed
   character.  It fills the buffer "str" up to the maximum size "count".

   This function will return -1 if the socket is closed during the read
   operation.

   Note that if a single line exceeds the length of count, the extra data
   will be read and discarded!  You have been warned. */
int MR3_sock_gets(sockfd, str, count)
SOCKET sockfd;
char *str;
size_t count;
{
  int bytes_read;
  int total_count = 0;
  char *current_position;
  char last_read = 0;

  current_position = str;
  while (last_read != 10) {
#ifndef WIN32
    bytes_read = read(sockfd, &last_read, 1);
#else
    bytes_read = recv(sockfd, &last_read, 1, 0);
#endif
    if (bytes_read <= 0) {
      /* The other side may have closed unexpectedly */
      return -1; /* Is this effective on other platforms than linux? */
    }
    if ( (total_count < count) && (last_read != 10) && (last_read !=13) ) {
      current_position[0] = last_read;
      current_position++;
      total_count++;
    }
  }
  if (count > 0)
    current_position[0] = 0;
  return total_count;
}

/* This function writes a character string out to a socket.  It will 
   return -1 if the connection is closed while it is trying to write. */
int MR3_sock_puts(sockfd, str)
SOCKET sockfd;
char *str;
{
  return MR3_sock_write(sockfd, str, strlen(str));
}

void MR3_sock_init(SOCKET *s_sock_org, int port_num, double x[], int igraph[],
		   int natom, int str_len, int lbres[],
		   int ipres[], int nres)
{
  int i0,i;
  SOCKET s_sock;
  double d0,d1,d2,d3;
  double real_buf[100];
  
  s_sock=*s_sock_org=MR3_open_port(port_num);
  printf("open port %d\n",port_num);
    
  MR3_sock_read_d(s_sock,&d0);
  printf("%f\n",d0);

  // amber mode on    
  i = 1;
  MR3_sock_write_i(s_sock,i);

  // set atom type number
  i0 = 12;
  MR3_sock_write_i(s_sock,i0);

  // set atom radius
  for(i = 0;i<i0;i++){
    real_buf[i] = 1.0;
  }

  d0 = 2.0;
  d1 = 0.7;
  real_buf[0] = 1.20*d0;
  real_buf[1] = 1.52*d0;
  real_buf[2] = 1.55*d0;
  real_buf[3] = 1.70*d0;
  real_buf[4] = 1.80*d0;
  real_buf[5] = 1.80*d0;
  real_buf[6] = 1.73*d0;
  real_buf[10] = 1.20*d1;
  real_buf[11] = 1.52*d1;
  for(i = 0;i<i0;i++){
    printf("%f\n",real_buf[i]);
  }
  MR3_sock_write_dv(s_sock,real_buf,i0);

  //     set color of atom 1 H
  i = 0;
  real_buf[i]   = 1.0;
  real_buf[i+1] = 1.0;
  real_buf[i+2] = 1.0;

  //     set color of atom 2 O
  i = i+3;
  real_buf[i]   = 1.0;
  real_buf[i+1] = 0.0;
  real_buf[i+2] = 0.0;

  //     set color of atom 3 N
  i = i+3;
  real_buf[i]   = 0.0;
  real_buf[i+1] = 1.0;
  real_buf[i+2] = 0.0;

  //     set color of atom 4 C
  i = i+3;
  real_buf[i]   = 0.5;
  real_buf[i+1] = 0.5;
  real_buf[i+2] = 0.5;

  //     set color of atom 5 S
  i = i+3;
  real_buf[i]   = 1.0;
  real_buf[i+1] = 1.0;
  real_buf[i+2] = 0.0;

  //     set color of atom 6 P
  i = i+3;
  real_buf[i]   = 1.0;
  real_buf[i+1] = 0.6;
  real_buf[i+2] = 0.0;

  //     set color of atom 5 Mg
  i = i+3;
  real_buf[i]   = 1.0;
  real_buf[i+1] = 0.0;
  real_buf[i+2] = 1.0;

  //     set color of atom 11 OW
  i = 10*3;
  real_buf[i]   = 0.5;
  real_buf[i+1] = 1.0;
  real_buf[i+2] = 1.0;

  //     set color of atom 12 HW
  i = i+3;
  real_buf[i]   = 1.0;
  real_buf[i+1] = 0.5;
  real_buf[i+2] = 1.0;

  MR3_sock_write_dv(s_sock,real_buf,i0*3);

  //     send atom number and type

  MR3_sock_write_i(s_sock,natom);

  MR3_sock_write_i(s_sock,nres);
  MR3_sock_write_iv(s_sock,lbres,nres);
  MR3_sock_write_iv(s_sock,ipres,nres);

  MR3_sock_write_iv(s_sock,igraph,natom);
  MR3_sock_write_dv(s_sock,x,natom*3);

  MR3_sock_write_i(s_sock,str_len);
  STRLEN=str_len;
}


void mr3_sock_init_(SOCKET *s_sock, int *port_num, double x[], int igraph[],
		    int *natom, int *str_len, int lbres[],
		    int ipres[], int *nres)
{
  MR3_sock_init(s_sock,*port_num,x,igraph,
		*natom,*str_len,lbres,
		ipres,*nres);
}


void mr3_sock_init__(SOCKET *s_sock, int *port_num, double x[], int igraph[],
		    int *natom, int *str_len, int lbres[],
		    int ipres[], int *nres)
{
  mr3_sock_init_(s_sock,port_num,x,igraph,
		 natom,str_len,lbres,
		 ipres,nres);
}


void MR3_sock_1step(SOCKET s_sock, double t1, double time, double t2, 
		    int recv_buf[],
		    int natom, double x[], double f[], double f_temp[])
{
  int n;
  char s[1000];
  sprintf(s,"T=%6.1f(K) N=%6d|Temp:%6.1f(K) time:%6.2f(ps)",
	  t2,natom,t1,time);
  // fprintf(stderr,"%s\n",s);
  //MR3_sock_write_c(s_sock,s);
  MR3_sock_write(s_sock,s, STRLEN);
  MR3_sock_read_iv(s_sock,recv_buf,3);
  //printf("recv_buf=%d %d %d\n",recv_buf[0],recv_buf[1],recv_buf[2]);
  MR3_sock_write_dv(s_sock,x,natom*3);
  MR3_sock_write_dv(s_sock,f+recv_buf[1]*3,recv_buf[2]*3);
  MR3_sock_read_dv(s_sock,f_temp,recv_buf[2]*3);
}


void mr3_sock_1step_(SOCKET *s_sock, double *t1, double *time, double *t2,
		     int recv_buf[],
		     int *natom, double x[], double f[], double f_tmp[])
{
  MR3_sock_1step(*s_sock,*t1,*time,*t2,recv_buf,
		 *natom,x,f,f_tmp);
}


void mr3_sock_1step__(SOCKET *s_sock, double *t1, double *time, double *t2,
		     int recv_buf[],
		     int *natom, double x[], double f[], double f_tmp[])
{
  MR3_sock_1step(*s_sock,*t1,*time,*t2,recv_buf,
		 *natom,x,f,f_tmp);
}
