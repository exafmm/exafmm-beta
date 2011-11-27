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
 
struct Vertex {
  vtkIdType Id;
  int       Ista;
  int       Iend;
  vect      X;
  vect      F;
};

const std::string  INPUT_PATH = "../../fdgl_in/";
const std::string OUTPUT_PATH = "../../fdgl_out/";

typedef std::vector<Vertex>::iterator V_iter;
std::vector<Vertex> vertices;
std::vector<int> edges;
vtkMutableUndirectedGraph *graph = vtkMutableUndirectedGraph::New();
vtkGraphLayoutView *view = vtkGraphLayoutView::New();
vtkGraphWriter *writer = vtkGraphWriter::New();
int eC, vC, maxDegree=0;

double get_time() {
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return double(tv.tv_sec+tv.tv_usec*1e-6);
}

void initVertices() {
  for (int i = 0; i < vC; i++) {
    Vertex vertex;
    vertex.Id = graph->AddVertex();   
    vertex.Ista = vertex.Iend = 0;
    vertex.X[0] = 2 * drand48() - 1;
    vertex.X[1] = 2 * drand48() - 1;
    vertex.X[2] = 2 * drand48() - 1;
    vertex.F = 0; 
    vertices.push_back(vertex);   
  }
}

void readEdges(std::ifstream &fid) {
  std::string line0, line1;
  
  fid >> line0 >> line1;
  int prevSrc=0, crntSrc=-1;
  while( !fid.eof() ) { 
    crntSrc = atoi(line0.c_str());
    if (crntSrc != prevSrc) {
      vertices[prevSrc].Iend = edges.size();
      vertices[crntSrc].Ista = edges.size();
      maxDegree = ((vertices[prevSrc].Iend - vertices[prevSrc].Ista) > maxDegree)
                   ?(vertices[prevSrc].Iend - vertices[prevSrc].Ista):maxDegree;
      
      prevSrc = crntSrc;
    }

    edges.push_back(atoi(line1.c_str()));
    fid >> line0 >> line1;
  }
  vertices[prevSrc].Iend = edges.size();
}

void readGraph(std::string fileid) {
  std::string line0, line1;
  
  std::string file=INPUT_PATH+fileid+".graph";
  
  std::ifstream fid(file.c_str(),std::ios::in);
  
  // First line in file contains <num nodes> <num edges>
  fid >> line0 >> line1;
  
  vC = atoi(line0.c_str());
  eC = atoi(line1.c_str());

  if (vC > 0) { // handling bad input files
    initVertices();
    readEdges(fid);
  } else {
    std::cerr << "Bad input file: " << file << endl;
  }
  fid.close();
}

void setEdges() {
  for( V_iter V=vertices.begin(); V!=vertices.end(); ++V ) {
    for( int i=V->Ista; i<V->Iend; ++i ) {
      // Data comes in duplicates - only add src <= dest
      if (V->Id <= vertices[edges[i]].Id) 
        graph->AddEdge(V->Id, vertices[edges[i]].Id);
    }
  }
}

void repulsion(ParallelFMM<Laplace> &FMM) {
  int begin = 0;
  int end = vC;
  FMM.splitRange(begin,end,MPIRANK,MPISIZE);
  Bodies bodies(end-begin);
  Bodies jbodies;
  Cells cells, jcells;
  B_iter B = bodies.begin();
  for( V_iter V=vertices.begin()+begin; V!=vertices.begin()+end; ++V, ++B ) {
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
  float *send0 = new float [vC];
  float *send1 = new float [vC];
  float *send2 = new float [vC];
  float *recv0 = new float [vC];
  float *recv1 = new float [vC];
  float *recv2 = new float [vC];
  for( int i=0; i!=vC; ++i ) {
    send0[i] = send1[i] = send2[i] = 0;
  }
  for( B=bodies.begin(); B!=bodies.end(); ++B ) {
    send0[B->IBODY] = B->TRG[1];
    send1[B->IBODY] = B->TRG[2];
    send2[B->IBODY] = B->TRG[3];
  }
  MPI_Allreduce(send0,recv0,vC,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(send1,recv1,vC,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(send2,recv2,vC,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
  for( int i=0; i!=vC; ++i ) {
    vertices[i].F[0] = recv0[i];
    vertices[i].F[1] = recv1[i];
    vertices[i].F[2] = recv2[i];
  }
  delete[] send0;
  delete[] send1;
  delete[] send2;
  delete[] recv0;
  delete[] recv1;
  delete[] recv2;
}

void spring() {
  float l = 2 / pow(vC,1./3);
  for( V_iter VI=vertices.begin(); VI!=vertices.end(); ++VI ) {
    for( int i=VI->Ista; i<VI->Iend; ++i ) {
      V_iter VJ = vertices.begin()+edges[i];
      vect dist = VI->X - VJ->X;
      float R = sqrtf(norm(dist) + EPS2);
      VI->F -= dist / R * (R - l);
    }
  }
}

void springEW() {
  float l = 2 / pow(vC,1./3);
  for( V_iter VI=vertices.begin(); VI!=vertices.end(); ++VI ) {
    for( int i=VI->Ista; i<VI->Iend; ++i ) {
      V_iter VJ = vertices.begin()+edges[i];
      vect dist = VI->X - VJ->X;
      
      int vSrcDeg  = VI->Iend - VI->Ista;
      int vDestDeg = VJ->Iend - VJ->Ista;
      
      float eW = (vSrcDeg==vDestDeg)?(1.0/vDestDeg):
                 abs(vSrcDeg - vDestDeg)/std::max(vSrcDeg,vDestDeg);
      
      float R = sqrtf(norm(dist) + EPS2);
      VI->F -= dist / R * (R - l) * eW;
    }
  }
}

void central() {
  for( V_iter V=vertices.begin(); V!=vertices.end(); ++V ) {
    vect dist = V->X;
    float R2 = norm(dist) + EPS2;
    float R3 = std::sqrt(R2) * R2;
    V->F -= dist / R3 * ( (V->Iend - V->Ista) / maxDegree) * 0;
  }
}

void moveVertices(int step) {
  for( V_iter V=vertices.begin(); V!=vertices.end(); ++V ) {
    if( sqrtf(norm(V->F))*step < 100 && norm(V->X) < 100 ) V->X += V->F * step * 0.001;
  } 
}

void setVertices() {
  vtkPoints *points = vtkPoints::New();
  for( V_iter V=vertices.begin(); V!=vertices.end(); ++V ) {
    points->InsertNextPoint(V->X[0], V->X[1], V->X[2]);
  }
  graph->SetPoints(points);
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
  bool useEW=false,useVW=false;
  
  t0 = get_time();
  readGraph(fileid);
  setEdges();
  t[0] += get_time() - t0;
  
  for( int step=0; step<300; ++step ) {
    std::cout << step << std::endl;
    
    t0 = get_time();
    repulsion(FMM);
    t[1] += get_time() - t0;
    
    if (useEW) {
      t0 = get_time();
      springEW();
      t[2] += get_time() - t0;
    } else {
      t0 = get_time();
      spring();
      t[2] += get_time() - t0;
    }
    
    if(useVW) {
      t0 = get_time();
      central();
      t[3] += get_time() - t0;
    }
    
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
//    writeGraph(fileid,step);
    t[7] += get_time() - t0;

    usleep(1000);
  }

  std::cout << "V          : " << vC << std::endl;
  std::cout << "E          : " << (edges.size()/2) << std::endl;
  std::cout << "initialize : " << t[0] << std::endl;
  std::cout << "repulsion  : " << t[1] << std::endl;
  if (useVW) std::cout << "springEW   : " << t[2] << std::endl;
  else       std::cout << "spring     : " << t[2] << std::endl;
  if (useVW) std::cout << "centralVW  : " << t[3] << std::endl;
  std::cout << "move       : " << t[4] << std::endl;
  std::cout << "setVer     : " << t[5] << std::endl;
  std::cout << "draw       : " << t[6] << std::endl;
  std::cout << "writeGraph : " << t[7] << std::endl;
  
  FMM.finalize();
  if( MPIRANK == 0 ) finalizeGraph();
}
