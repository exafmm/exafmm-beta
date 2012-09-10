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
#if VTK
#define VTK_EXCLUDE_STRSTREAM_HEADERS
#include <vtkPoints.h>
#include <vtkGraphLayoutView.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkGraphWriter.h>
#include <vtkViewTheme.h>
#endif
#include "parallelfmm.h"

//! Structure of source vertices
struct JVertex {
  int       Ista;                                               //!< Start index of connection list
  int       Iend;                                               //!< End index of connection list
  vect      X;                                                  //!< Position vector
};

//! Structure of vertices
struct Vertex : public JVertex {
#if VTK
  vtkIdType Id;                                                 //!< VTK vertex ID
#endif
  vect      F;                                                  //!< Force vector
};

const std::string  INPUT_PATH = "../../archive/fdgl_in/";       //!< Input file name
const std::string OUTPUT_PATH = "../../archive/fdgl_out/";      //!< Output file name

std::vector<int> edges;                                         // Vector of edges
std::vector<Vertex> vertices;                                   // Vector of vertices
std::vector<Vertex> verticesG;                                  // Global vector of vertices
std::vector<JVertex> jvertices;                                 // Vector of source vertices
typedef std::vector<Vertex>::iterator V_iter;                   // Iterator for vertex vector
typedef std::vector<JVertex>::iterator JV_iter;                 // Iterator for source vertex vector
#if VTK
vtkMutableUndirectedGraph *graph = vtkMutableUndirectedGraph::New();// VTK graph object
vtkGraphLayoutView *view = vtkGraphLayoutView::New();           // VTK view object
vtkGraphWriter *writer = vtkGraphWriter::New();                 // VTK writer object
#endif
int numVertices, numEdges, localVertices, maxEdges=0;           // Define global variables
int *offset = new int [MPISIZE];                                // Offset for Allgatherv

//! Timer function
double get_time() {
  struct timeval tv;                                            // Time value
  gettimeofday(&tv,NULL);                                       // Get time of day in seconds and microseconds
  return double(tv.tv_sec+tv.tv_usec*1e-6);                     // Combine seconds and microseconds and return
}

//! Split array into numSplit ranges
void splitRange(int &begin, int &end, int iSplit, int numSplit) {
  int size = end - begin;                                       // Total size of array
  int increment = size / numSplit;                              // Approximate size of range
  int remainder = size % numSplit;                              // Remainder of split range
  begin += iSplit * increment + std::min(iSplit,remainder);     // Begin index of range
  end = begin + increment;                                      // End index of range
  if( remainder > iSplit ) end++;                               // Account for remainders
}

//! Initialize vertices
void initVertices() {
  int begin = 0;                                                // Global Begin index
  int end = numVertices;                                        // Global End index
  splitRange(begin,end,MPIRANK,MPISIZE);                        // Get local begin & end index
  localVertices = end - begin;                                  // Get local vertex size
  srand48(0);                                                   // Seed random number generator
  for( int i=0; i<numVertices; ++i ) {                          // Loop over global vertices
    Vertex vertex;                                              // Instantiate vertex type
#if VTK
    vertex.Id = graph->AddVertex();                             // Set VTK vertex ID
#endif
    vertex.Ista = vertex.Iend = 0;                              // Initialize start & end index of connection list
    vertex.X[0] = 2 * drand48() - 1;                            // Initialize x positions
    vertex.X[1] = 2 * drand48() - 1;                            // Initialize y positions
    vertex.X[2] = 2 * drand48() - 1;                            // Initialize z positions
    vertex.F = 0;                                               // Initialize force vector
    if( begin <= i && i < end ) vertices.push_back(vertex);     // Push vertex into vector
  }
  verticesG.resize(numVertices);                                // Resize global vertex vector
}

//! Read edge info from file
void readEdges(std::ifstream &fid) {
  std::string line0, line1;                                     // Define temporary strings for reading data
  int begin = 0;                                                // Global Begin index
  int end = numVertices;                                        // Global End index
  splitRange(begin,end,MPIRANK,MPISIZE);                        // Get local begin & end index
  fid >> line0 >> line1;                                        // Read first data
  while( atoi(line0.c_str()) < begin ) {                        // While data is before beginning of range
    fid >> line0 >> line1;                                      //  Read data (not used)
  }                                                             // End while loop for range
  int targetOld=0, target=-1;                                   // Initialize previous and current target vertex
  while( !fid.eof() && atoi(line0.c_str()) < end ) {            // While not end of file and within range
    target = atoi(line0.c_str()) - begin;                       //  Current target vertex
    if (target != targetOld) {                                  //  If it's a new target vertex
      vertices[targetOld].Iend = edges.size();                  //   Set end index of previous range
      maxEdges = std::max(vertices[targetOld].Iend - vertices[targetOld].Ista, maxEdges);// Keep track of max edges
      vertices[target].Ista = edges.size();                     //   Set begin index of current range
      targetOld = target;                                       //   Save target in targetOld
    }                                                           //  Endif for new target
    edges.push_back(atoi(line1.c_str()));                       //  Push edge into vector
    fid >> line0 >> line1;                                      //  Read new line
  }                                                             // End while loop for file read
  vertices[targetOld].Iend = edges.size();                      // Set end index of last range
  maxEdges = std::max(vertices[targetOld].Iend - vertices[targetOld].Ista, maxEdges);// Keep track of max edges
}

//! Interface to initVertices and readEdges
void readGraph(std::string fileid) {
  std::string line0, line1;                                     // Define temporary strings for reading data
  std::string file=INPUT_PATH+fileid+".graph";                  // Set file name
  std::ifstream fid(file.c_str(),std::ios::in);                 // Open file
  fid >> line0 >> line1;                                        // Read first line
  numVertices = atoi(line0.c_str());                            // Set global number of vertices
  numEdges    = atoi(line1.c_str());                            // Set global number of edges
  initVertices();                                               // Initialize vertices
  readEdges(fid);                                               // Read edge info from file
  jvertices.resize(numEdges);                                   // Resize source vertex vector
  fid.close();                                                  // Close file
}

//! Communicate vertices round-robin
void shiftVertices() {
  int newSize = 0;                                              // Initialize new number of vertices
  int oldSize = vertices.size();                                // Get old number of vertices
  const int bytes = sizeof(Vertex) / 8;                         // Get word count of Vertex type
  const int isend = (MPIRANK + 1          ) % MPISIZE;          // Rank to send to
  const int irecv = (MPIRANK - 1 + MPISIZE) % MPISIZE;          // Rank to receive from
  MPI_Request sreq,rreq;                                        // MPI request flags

  MPI_Isend(&oldSize,1,MPI_INT,irecv,0,MPI_COMM_WORLD,&sreq);   // Send number of vertices to communicate
  MPI_Irecv(&newSize,1,MPI_INT,isend,0,MPI_COMM_WORLD,&rreq);   // Recv number of vertices to communicate
  MPI_Wait(&sreq,MPI_STATUS_IGNORE);                            // Wait for send to complete
  MPI_Wait(&rreq,MPI_STATUS_IGNORE);                            // Wait for receive to complete

  std::vector<Vertex> buffer(newSize);                          // Allocate send/receive buffer
  MPI_Isend(&vertices[0],oldSize*bytes,MPI_DOUBLE,irecv,        // Send vertices
            1,MPI_COMM_WORLD,&sreq);
  MPI_Irecv(&buffer[0]  ,newSize*bytes,MPI_DOUBLE,isend,        // Receive vertices
            1,MPI_COMM_WORLD,&rreq);
  MPI_Wait(&sreq,MPI_STATUS_IGNORE);                            // Wait for send to complete
  MPI_Wait(&rreq,MPI_STATUS_IGNORE);                            // Wait for receive to complete
  vertices = buffer;                                            // Copy buffer to vertex vector
}

//! Communicate edges round-robin
void shiftEdges() {
  int newSize = 0;                                              // Initialize new number of edges
  int oldSize = edges.size();                                   // Get old number of edges
  const int isend = (MPIRANK + 1          ) % MPISIZE;          // Rank to send to
  const int irecv = (MPIRANK - 1 + MPISIZE) % MPISIZE;          // Rank to receive from
  MPI_Request sreq,rreq;                                        // MPI request flags

  MPI_Isend(&oldSize,1,MPI_INT,irecv,0,MPI_COMM_WORLD,&sreq);   // Send number of edges to communicate
  MPI_Irecv(&newSize,1,MPI_INT,isend,0,MPI_COMM_WORLD,&rreq);   // Send number of edges to communicate
  MPI_Wait(&sreq,MPI_STATUS_IGNORE);                            // Wait for send to complete
  MPI_Wait(&rreq,MPI_STATUS_IGNORE);                            // Wait for receive to complete

  std::vector<int> buffer(newSize);                             // Allocate send/receive buffer
  MPI_Isend(&edges[0] ,oldSize,MPI_INT,irecv,                   // Send edges
            1,MPI_COMM_WORLD,&sreq);
  MPI_Irecv(&buffer[0],newSize,MPI_INT,isend,                   // Receive edges
            1,MPI_COMM_WORLD,&rreq);
  MPI_Wait(&sreq,MPI_STATUS_IGNORE);                            // Wait for send to complete
  MPI_Wait(&rreq,MPI_STATUS_IGNORE);                            // Wait for receive to complete
  edges = buffer;                                               // Copy buffer to edge vector
}

//! Gather vertices from all ranks
void gatherVertices() {
  int sendCnt = localVertices * sizeof(Vertex) / 8;             // Send count : number of local vertices in 8 byte units
  int *recvCnt = new int [MPISIZE];                             // Array of receive counts from all procs
  MPI_Allgather(&sendCnt,1,MPI_INT,recvCnt,1,MPI_INT,MPI_COMM_WORLD);// Gather receive counts
  offset[0] = 0;                                                // Initialize receive count offset
  for( int i=1; i!=MPISIZE; ++i ) {                             // Loop over MPI ranks
    offset[i] = offset[i-1] + recvCnt[i-1];                     //  Determine receive count offsets
  }                                                             // End loop over MPI ranks
  MPI_Allgatherv(&vertices[0],sendCnt,MPI_DOUBLE,&verticesG[0], // Gather vertices
                              recvCnt,offset,MPI_DOUBLE,MPI_COMM_WORLD);
  delete[] recvCnt;                                             // Delete array of receive counts
}

#if VTK
//! Set edges for VTK
void setEdges() {
  gatherVertices();                                             // Gather vertices from all ranks
  for( int irank=0; irank!=MPISIZE; ++irank ) {                 // Loop over MPI ranks
    for( V_iter V=vertices.begin(); V!=vertices.end(); ++V ) {  //  Loop over vertices
      for( int i=V->Ista; i<V->Iend; ++i ) {                    //   Loop over edges
        if (V->Id < verticesG[edges[i]].Id) {                   //    If target ID <= source ID
          graph->AddEdge(V->Id, verticesG[edges[i]].Id);        //     Add edge to VTK graph object
        }                                                       //    Endif for target ID <= source ID
      }                                                         //   End loop over edges
    }                                                           //  End loop over vertices
    shiftVertices();                                            //  Communicate vertices round-robin
    shiftEdges();                                               //  Communicate edges round-robin
  }                                                             // Eend loop over MPI ranks
}

//! Set vertices for VTK
void setVertices() {
  vtkPoints *points = vtkPoints::New();                         // Define VTK points object
  for( int irank=0; irank!=MPISIZE; ++irank ) {                 // Loop over MPI ranks
    for( V_iter V=vertices.begin(); V!=vertices.end(); ++V ) {  //  Loop over vertices
      points->InsertNextPoint(V->X[0], V->X[1], V->X[2]);       //   Add vertex to VTK points object
    }                                                           //  End loop over vertices
    shiftVertices();                                            //  Communicate vertices round-robin
  }                                                             // End loop over MPI ranks
  graph->SetPoints(points);                                     // Add points to VTK graph object
}
#endif

//! Repulsion force calculation (interface to FMM)
void repulsion(ParallelFMM<Laplace> &FMM) {
  Bodies bodies(localVertices);                                 // Define vector of bodies
  Bodies jbodies;                                               // Define vector of source bodies
  Cells cells, jcells;                                          // Define vector of cells
  B_iter B = bodies.begin();                                    // Iterator of first body
  vect Xmin = B->X;                                             // Initialize minimum of X
  vect Xmax = B->X;                                             // Initialize maximum of X
  for( V_iter V=vertices.begin(); V!=vertices.end(); ++V, ++B ) {// Loop over vertices
    B->X = V->X;                                                //  Copy vertex position to body position
    B->SRC = -100;                                              //  Set source value
    B->TRG = 0;                                                 //  Initialize target value
    B->IBODY = V-vertices.begin();                              //  Tag body with initial index
    B->IPROC = MPIRANK;                                         //  Tag body with initial MPI rank
    if( Xmin[0] < B->X[0] ) Xmin[0] = B->X[0];                  //  Determine xmin
    if( Xmin[1] < B->X[1] ) Xmin[1] = B->X[1];                  //  Determine ymin
    if( Xmin[2] < B->X[2] ) Xmin[2] = B->X[2];                  //  Determine zmin
    if( Xmax[0] > B->X[0] ) Xmax[0] = B->X[0];                  //  Determine xmax
    if( Xmax[1] > B->X[1] ) Xmax[1] = B->X[1];                  //  Determine ymax
    if( Xmax[2] > B->X[2] ) Xmax[2] = B->X[2];                  //  Determine zmax
  }                                                             // End loop over vertices
  FMM.setGlobDomain(bodies);                                    // Set global domain size of FMM
  FMM.octsection(bodies);                                       // Partition domain and redistribute bodies
  FMM.bottomup(bodies,cells);                                   // Tree construction (bottom up) & upward sweep
  FMM.commBodies(cells);                                        // Send bodies (not receiving yet)
  jbodies = bodies;                                             // Vector of source bodies
  jcells = cells;                                               // Vector of source cells
  FMM.commCells(jbodies,jcells);                                // Communicate cells (receive bodies here)
  FMM.downward(cells,jcells);                                   // Downward sweep
  FMM.unpartition(bodies);                                      // Send bodies back to where they came from
  for( B=bodies.begin(); B!=bodies.end(); ++B ) {               // Loop over bodies
    vertices[B->IBODY].F[0] = B->TRG[1];                        //  Copy body x force to vertex x force
    vertices[B->IBODY].F[1] = B->TRG[2];                        //  Copy body y force to vertex y force
    vertices[B->IBODY].F[2] = B->TRG[3];                        //  Copy body z force to vertex z force
  }                                                             // End loop over bodies
}

//! Do a gather to update source vertex values
void commVertices() {
  gatherVertices();                                             // Gather vertices from all ranks
  for( V_iter VI=vertices.begin(); VI!=vertices.end(); ++VI ) { // Loop over target vertices
    for( int i=VI->Ista; i<VI->Iend; ++i ) {                    //  Loop over edges
      V_iter VJ = verticesG.begin()+edges[i];                   //   Get source vertex
      jvertices[i].Ista = VJ->Ista;                             //   Update start index of connection list
      jvertices[i].Iend = VJ->Iend;                             //   Update end index of connection list
      jvertices[i].X    = VJ->X;                                //   Update position vector
    }                                                           //  End loop over edges
  }                                                             // End loop over target vertices
}

//! Spring force calculation
void spring() {
  float l = 2 / pow(numVertices,1./3);                          // Set natural length of spring
  for( V_iter VI=vertices.begin(); VI!=vertices.end(); ++VI ) { // Loop over target vertices
    for( int i=VI->Ista; i<VI->Iend; ++i ) {                    //  Loop over edges
      JV_iter VJ = jvertices.begin()+i;                         //   Get source vertex
      vect dist = VI->X - VJ->X;                                //   Distance vector from source to target
      float R = sqrtf(norm(dist) + EPS2);                       //   Scalar distance from source to target
      float weight = (VI->Iend-VI->Ista) * (VJ->Iend-VJ->Ista); //   Weight based on number of edges
      VI->F -= dist / R * (R - l / weight);                     //   Spring force vector
    }                                                           //  End loop over edges
  }                                                             // End loop over vertices
}

//! Move vertices
void moveVertices() {
  for( V_iter V=vertices.begin(); V!=vertices.end(); ++V ) {    // Loop over vertices
    vect dX;                                                    //  Position increment
    if( norm(V->F) < EPS ) dX = 0;                              //  Filter noisy movement
    else dX = V->F / std::sqrt(norm(V->F));                     //  Always move at constant speed
    V->X += dX;                                                 //  Update position
  }                                                             // End loop over vertices
}

#if VTK
//! draw VTK graph
void drawGraph() {
  vtkViewTheme* theme = vtkViewTheme::New();                    // VTK theme object
  theme->SetBackgroundColor(1.,1.,1.);                          // Set lower background to white
//  theme->SetBackgroundColor2(1.,1.,1.);                         // Set upper background to white
  theme->SetPointColor(.2,.2,.2);                               // Set point color to gray
  theme->SetCellColor(.2,.2,.2);                                // Set cell color to gray
  view->ApplyViewTheme(theme);                                  // Apply theme to VTK view object
  theme->Delete();                                              // Delete theme object
  view->AddRepresentationFromInput(graph);                      // Add graph to view
  view->SetLayoutStrategy("Pass Through");                      // Set layout strategy
  view->GetInteractor()->GetRenderWindow()->SetSize(700,700);   // Set window size
  view->ResetCamera();                                          // Reset camera
  view->Render();                                               // Render
}

//! Make plot interactive
void finalizeGraph() {
  vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();// VTK style object
  view->GetInteractor()->SetInteractorStyle(style);             // Set style to VTK view object
  view->GetInteractor()->Initialize();                          // Initialzie interactor
  view->GetInteractor()->Start();                               // Start interactor
}

//! VTK formatted I/O
void writeGraph(std::string fileid, int step) {
  std::stringstream fname;                                      // Stringstream for file name
  fname << OUTPUT_PATH << fileid << '_' << step << ".vtk";      // Create file name
  writer->SetFileName(fname.str().c_str());                     // Set file name
  writer->SetInput(graph);                                      // Input graph object
  writer->Write();                                              // Write to file
}
#endif

int main() {
  double t0, t[7] = {0,0,0,0,0,0,0};                            // Timer array
  IMAGES = 0;                                                   // Level of periodic image tree (0 for non-periodic)
  THETA = 0.6;                                                  // Multipole acceptance criteria
  ParallelFMM<Laplace> FMM;                                     // Instantiate ParallelFMM class
  FMM.initialize();                                             // Initialize FMM
  std::string fileid="jagmesh1_v936_e2664";                     // Set data file (choose from below)
//  std::string fileid="add32_v4960_e9462";
//  std::string fileid="4elt_v15606_e45878";
//  std::string fileid="finan512_v74752_e261120";
//  std::string fileid="uni-socialgraph-anonymized_v59216214_e92522012";
//  std::string fileid="mhrw-socialgraph-anonymized_v72261577_e159189465";
//  std::string fileid="uplrg_v1000_e2000";
//  std::string fileid="uplrg_v10000_e20000";
//  std::string fileid="uplrg_v100000_e200000";
//  std::string fileid="uplrg_v1000000_e2000000";
//  std::string fileid="uplrg_v10000000_e20000000";
//  std::string fileid="uplrg_v100000000_e200000000";

  t0 = get_time();                                              // Start timer
  readGraph(fileid);                                            // Read graph data
#if VTK
  setEdges();                                                   // Set VTK edges
#endif
  t[0] += get_time() - t0;                                      // Stop timer

  for( int step=0; step<1000; ++step ) {                        // Loop over time steps
    if( MPIRANK == 0 ) std::cout << step << std::endl;          //  Print step
    t0 = get_time();                                            //  Start timer
    repulsion(FMM);                                             //  Repulsion force calculation
    t[1] += get_time() - t0;                                    //  Stop timer

    t0 = get_time();                                            //  Start timer
    commVertices();                                             //  Update source vertex values
    spring();                                                   //  Spring force calculation
    t[2] += get_time() - t0;                                    //  Stop timer

    t0 = get_time();                                            //  Start timer
    moveVertices();                                             //  Move vertices
    t[3] += get_time() - t0;                                    //  Stop timer
#if VTK
    t0 = get_time();                                            //  Start timer
    setVertices();                                              //  Set VTK vertices
    t[4] += get_time() - t0;                                    //  Stop timer

    t0 = get_time();                                            //  Start timer
    if( MPIRANK == 0 && step % 10 == 0 ) drawGraph();           //  Draw VTK graph
    t[5] += get_time() - t0;                                    //  Stop timer

    t0 = get_time();                                            //  Start timer
//    if( MPIRANK == 0 ) writeGraph(fileid,step);                 //  VTK formatted I/O
    t[6] += get_time() - t0;                                    //  Stop timer
#endif
  }

  if( MPIRANK == 0 ) {
    std::cout << "V          : " << numVertices << std::endl;   // Print number of vertices
    std::cout << "E          : " << numEdges / 2 << std::endl;  // Print number of edges
    std::cout << "initialize : " << t[0] << std::endl;          // Print initialization time
    std::cout << "repulsion  : " << t[1] << std::endl;          // Print repulsion time
    std::cout << "spring     : " << t[2] << std::endl;          // Print spring time
    std::cout << "move       : " << t[3] << std::endl;          // Print move time
#if VTK
    std::cout << "setVer     : " << t[4] << std::endl;          // Print set vertex time
    std::cout << "draw       : " << t[5] << std::endl;          // Print draw time
    std::cout << "writeGraph : " << t[6] << std::endl;          // Print I/O time
#endif
  }

#if VTK
  if( MPIRANK == 0 ) finalizeGraph();                           // Finalize VTK graph (make it interactive)
#endif
  FMM.finalize();                                               // Finalize FMM
}
