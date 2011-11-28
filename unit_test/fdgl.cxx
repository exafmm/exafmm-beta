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
#define VTK_EXCLUDE_STRSTREAM_HEADERS
#include <vtkPoints.h>
#include <vtkGraphLayoutView.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkGraphWriter.h>
#include "parallelfmm.h"

struct JVertex {
  int       Ista;
  int       Iend;
  vect      X;
};

struct Vertex : JVertex {
  vtkIdType Id;
  vect      F;
};

const std::string  INPUT_PATH = "../../fdgl_in/";
const std::string OUTPUT_PATH = "../../fdgl_out/";

std::vector<int> edges;
std::vector<Vertex> vertices;
std::vector<Vertex> verticesG;
std::vector<JVertex> jvertices;
typedef std::vector<Vertex>::iterator V_iter;
typedef std::vector<JVertex>::iterator JV_iter;
vtkMutableUndirectedGraph *graph = vtkMutableUndirectedGraph::New();
vtkGraphLayoutView *view = vtkGraphLayoutView::New();
vtkGraphWriter *writer = vtkGraphWriter::New();
int numVertices, numEdges, localVertices, maxEdges=0;
int *offset = new int [MPISIZE];

double get_time() {
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return double(tv.tv_sec+tv.tv_usec*1e-6);
}

void splitRange(int &begin, int &end, int iSplit, int numSplit) {
  int size = end - begin;
  int increment = size / numSplit;
  int remainder = size % numSplit;
  begin += iSplit * increment + std::min(iSplit,remainder);
  end = begin + increment;
  if( remainder > iSplit ) end++;
}

void initVertices() {
  int begin = 0;
  int end = numVertices;
  splitRange(begin,end,MPIRANK,MPISIZE);
  localVertices = end - begin;
  for( int i=0; i<numVertices; ++i ) {
    Vertex vertex;
    vertex.Id = graph->AddVertex();
    vertex.Ista = vertex.Iend = 0;
    vertex.X[0] = 2 * drand48() - 1;
    vertex.X[1] = 2 * drand48() - 1;
    vertex.X[2] = 2 * drand48() - 1;
    vertex.F = 0;
    if( begin <= i && i < end ) vertices.push_back(vertex);
  }
  verticesG.resize(numVertices);
}

void readEdges(std::ifstream &fid) {
  std::string line0, line1;
  int begin = 0;
  int end = numVertices;
  splitRange(begin,end,MPIRANK,MPISIZE);
  fid >> line0 >> line1;
  while( atoi(line0.c_str()) < begin ) {
    fid >> line0 >> line1;
  }
  int prevSrc=0, crntSrc=-1;
  while( !fid.eof() ) {
    crntSrc = atoi(line0.c_str()) - begin;
    if (crntSrc != prevSrc) {
      vertices[prevSrc].Iend = edges.size();
      maxEdges = std::max(vertices[prevSrc].Iend - vertices[prevSrc].Ista, maxEdges);
      vertices[crntSrc].Ista = edges.size();
      prevSrc = crntSrc;
    }
    edges.push_back(atoi(line1.c_str()));
    fid >> line0 >> line1;
  }
  vertices[prevSrc].Iend = edges.size();
  maxEdges = std::max(vertices[prevSrc].Iend - vertices[prevSrc].Ista, maxEdges);
}

void readGraph(std::string fileid) {
  std::string line0, line1;
  std::string file=INPUT_PATH+fileid+".graph";
  std::ifstream fid(file.c_str(),std::ios::in);
  // First line in file contains <num nodes> <num edges>
  fid >> line0 >> line1;
  numVertices = atoi(line0.c_str());
  numEdges    = atoi(line1.c_str());
  if (numVertices > 0) { // handling bad input files
    initVertices();
    readEdges(fid);
  } else {
    std::cerr << "Bad input file: " << file << std::endl;
  }
  jvertices.resize(numEdges);
  fid.close();
}

void shiftVertices() {
  int newSize = 0;
  int oldSize = vertices.size();
  const int bytes = sizeof(Vertex);
  const int isend = (MPIRANK + 1          ) % MPISIZE;
  const int irecv = (MPIRANK - 1 + MPISIZE) % MPISIZE;
  MPI_Request sreq,rreq;

  MPI_Isend(&oldSize,1,MPI_INT,irecv,0,MPI_COMM_WORLD,&sreq);
  MPI_Irecv(&newSize,1,MPI_INT,isend,0,MPI_COMM_WORLD,&rreq);
  MPI_Wait(&sreq,MPI_STATUS_IGNORE);
  MPI_Wait(&rreq,MPI_STATUS_IGNORE);

  std::vector<Vertex> buffer(newSize);
  MPI_Isend(&vertices[0],oldSize*bytes,MPI_BYTE,irecv,
            1,MPI_COMM_WORLD,&sreq);
  MPI_Irecv(&buffer[0]  ,newSize*bytes,MPI_BYTE,isend,
            1,MPI_COMM_WORLD,&rreq);
  MPI_Wait(&sreq,MPI_STATUS_IGNORE);
  MPI_Wait(&rreq,MPI_STATUS_IGNORE);
  vertices = buffer;
}

void shiftEdges() {
  int newSize = 0;
  int oldSize = edges.size();
  const int bytes = sizeof(int);
  const int isend = (MPIRANK + 1          ) % MPISIZE;
  const int irecv = (MPIRANK - 1 + MPISIZE) % MPISIZE;
  MPI_Request sreq,rreq;

  MPI_Isend(&oldSize,1,MPI_INT,irecv,0,MPI_COMM_WORLD,&sreq);
  MPI_Irecv(&newSize,1,MPI_INT,isend,0,MPI_COMM_WORLD,&rreq);
  MPI_Wait(&sreq,MPI_STATUS_IGNORE);
  MPI_Wait(&rreq,MPI_STATUS_IGNORE);

  std::vector<int> buffer(newSize);
  MPI_Isend(&edges[0] ,oldSize*bytes,MPI_BYTE,irecv,
            1,MPI_COMM_WORLD,&sreq);
  MPI_Irecv(&buffer[0],newSize*bytes,MPI_BYTE,isend,
            1,MPI_COMM_WORLD,&rreq);
  MPI_Wait(&sreq,MPI_STATUS_IGNORE);
  MPI_Wait(&rreq,MPI_STATUS_IGNORE);
  edges = buffer;
}

void gatherVertices() {
  int sendCnt = localVertices * sizeof(Vertex);
  int *recvCnt = new int [MPISIZE];
  MPI_Allgather(&sendCnt,1,MPI_INT,recvCnt,1,MPI_INT,MPI_COMM_WORLD);
  offset[0] = 0;
  for( int i=1; i!=MPISIZE; ++i ) {
    offset[i] = offset[i-1] + recvCnt[i-1];
  }
  MPI_Allgatherv(&vertices[0],sendCnt,MPI_BYTE,&verticesG[0],recvCnt,offset,MPI_BYTE,MPI_COMM_WORLD);
  delete[] recvCnt;
}

void setEdges() {
  gatherVertices();
  for( int irank=0; irank!=MPISIZE; ++irank ) {
    for( V_iter V=vertices.begin(); V!=vertices.end(); ++V ) {
      for( int i=V->Ista; i<V->Iend; ++i ) {
        // Data comes in duplicates - only add src <= dest
        if (V->Id < verticesG[edges[i]].Id)
          graph->AddEdge(V->Id, verticesG[edges[i]].Id);
      }
    }
    shiftVertices();
    shiftEdges();
  }
}

void setVertices() {
  vtkPoints *points = vtkPoints::New();
  for( int irank=0; irank!=MPISIZE; ++irank ) {
    for( V_iter V=vertices.begin(); V!=vertices.end(); ++V ) {
      points->InsertNextPoint(V->X[0], V->X[1], V->X[2]);
    }
    shiftVertices();
  }
  graph->SetPoints(points);
}

void repulsion(ParallelFMM<Laplace> &FMM) {
  Bodies bodies(localVertices);
  Bodies jbodies;
  Cells cells, jcells;
  B_iter B = bodies.begin();
  for( V_iter V=vertices.begin(); V!=vertices.end(); ++V, ++B ) {
    B->X = V->X;
    B->SRC = -0.001;
    B->TRG = 0;
    B->IBODY = V-vertices.begin();
    B->IPROC = MPIRANK;
    B->ICELL = 0;
  }
  FMM.setGlobDomain(bodies);
  FMM.octsection(bodies);
  FMM.bottomup(bodies,cells);
  FMM.commBodies(cells);
  jbodies = bodies;
  jcells = cells;
  FMM.commCells(jbodies,jcells);
  FMM.downward(cells,jcells,1);
  FMM.unpartition(bodies);
  for( B=bodies.begin(); B!=bodies.end(); ++B ) {
    vertices[B->IBODY].F[0] = B->TRG[1];
    vertices[B->IBODY].F[1] = B->TRG[2];
    vertices[B->IBODY].F[2] = B->TRG[3];
  }
}

void commVertices() {
  gatherVertices();
  for( V_iter VI=vertices.begin(); VI!=vertices.end(); ++VI ) {
    for( int i=VI->Ista; i<VI->Iend; ++i ) {
      V_iter VJ = verticesG.begin()+edges[i];
      jvertices[i].Ista = VJ->Ista;
      jvertices[i].Iend = VJ->Iend;
      jvertices[i].X    = VJ->X;
    }
  }
}

void spring() {
  float l = 2 / pow(numVertices,1./3);
  for( V_iter VI=vertices.begin(); VI!=vertices.end(); ++VI ) {
    for( int i=VI->Ista; i<VI->Iend; ++i ) {
      JV_iter VJ = jvertices.begin()+i;
      vect dist = VI->X - VJ->X;
      float R = sqrtf(norm(dist) + EPS2);
      VI->F -= dist / R * (R - l);
    }
  }
}

void springWeighted() {
  float l = 2 / pow(numVertices,1./3);
  for( V_iter VI=vertices.begin(); VI!=vertices.end(); ++VI ) {
    for( int i=VI->Ista; i<VI->Iend; ++i ) {
      JV_iter VJ = jvertices.begin()+i;
      vect dist = VI->X - VJ->X;
      int vSrcDeg = VI->Iend - VI->Ista;
      int vTrgDeg = VJ->Iend - VJ->Ista;
      float weight = (vSrcDeg==vTrgDeg)?(1.0/vTrgDeg):
                  abs(vSrcDeg-vTrgDeg)/std::max(vSrcDeg,vTrgDeg);
      float R = sqrtf(norm(dist) + EPS2);
      VI->F -= dist / R * (R - l) * weight;
    }
  }
}

void central() {
  for( V_iter V=vertices.begin(); V!=vertices.end(); ++V ) {
    vect dist = V->X;
    float R2 = norm(dist) + EPS2;
    float R3 = std::sqrt(R2) * R2;
    V->F -= dist / R3 * ( (V->Iend - V->Ista) / maxEdges) * 0;
  }
}

void moveVertices(int step) {
  for( V_iter V=vertices.begin(); V!=vertices.end(); ++V ) {
    if( sqrtf(norm(V->F))*step < 100 && norm(V->X) < 100 ) V->X += V->F * step * 0.001;
  }
}

void drawGraph() {
  view->AddRepresentationFromInput(graph);
  view->SetLayoutStrategy("Pass Through");
  view->GetInteractor()->GetRenderWindow()->SetSize(700,700);
  view->ResetCamera();
  view->Render();
}

void finalizeGraph() {
  vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
  view->GetInteractor()->SetInteractorStyle(style);
  view->GetInteractor()->Initialize();
  view->GetInteractor()->Start();
}

void writeGraph(std::string fileid, int step) {
  std::stringstream ss;
  ss << OUTPUT_PATH << fileid << '_' << step << ".vtk";
  writer->SetFileName(ss.str().c_str());
  writer->SetInput(graph);
  writer->Write();
}

int main() {
  double t0, t[8] = {0,0,0,0,0,0,0,0};
  IMAGES = 0;
  THETA = 0.6;
  ParallelFMM<Laplace> FMM;
  FMM.initialize();
  std::string fileid="test3_v331_e361";
  bool useEW=false, useVW=false;

  t0 = get_time();
  readGraph(fileid);
  setEdges();
  t[0] += get_time() - t0;

  for( int step=0; step<100; ++step ) {
    if( MPIRANK == 0 ) std::cout << step << std::endl;

    t0 = get_time();
    repulsion(FMM);
    t[1] += get_time() - t0;

    t0 = get_time();
    commVertices();
    if (useEW) {
      springWeighted();
    } else {
      spring();
    }
    t[2] += get_time() - t0;

    t0 = get_time();
    if(useVW) {
      central();
    }
    t[3] += get_time() - t0;

    t0 = get_time();
    moveVertices(step);
    t[4] += get_time() - t0;

    t0 = get_time();
    setVertices();
    t[5] += get_time() - t0;

    t0 = get_time();
    if( MPIRANK == 0 ) drawGraph();
    t[6] += get_time() - t0;

    t0 = get_time();
//    if( MPIRANK == 0 ) writeGraph(fileid,step);
    t[7] += get_time() - t0;

    usleep(1000);
  }

  if( MPIRANK == 0 ) {
    std::cout << "V          : " << numVertices << std::endl;
    std::cout << "E          : " << numEdges / 2 << std::endl;
    std::cout << "initialize : " << t[0] << std::endl;
    std::cout << "repulsion  : " << t[1] << std::endl;
    std::cout << "spring     : " << t[2] << std::endl;
    if (useVW) std::cout << "central    : " << t[3] << std::endl;
    std::cout << "move       : " << t[4] << std::endl;
    std::cout << "setVer     : " << t[5] << std::endl;
    std::cout << "draw       : " << t[6] << std::endl;
    std::cout << "writeGraph : " << t[7] << std::endl;
  }

  FMM.finalize();
  if( MPIRANK == 0 ) finalizeGraph();
}
