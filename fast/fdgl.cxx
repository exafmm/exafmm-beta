#define VTK_EXCLUDE_STRSTREAM_HEADERS
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <vtkPoints.h>
#include <vtkGraphLayoutView.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <sys/time.h>
#include "tree.h"
 
struct Vertex {
  vtkIdType Id;
  int       Ista;
  int       Iend;
  vect      X;
  vect      F;
};

typedef std::vector<Body>::iterator B_iter;
typedef std::vector<Vertex>::iterator V_iter;
std::vector<Body> bodies;
std::vector<Cell> cells;
std::vector<Vertex> vertices;
std::vector<int> edges;
vtkMutableUndirectedGraph *graph = vtkMutableUndirectedGraph::New();
vtkGraphLayoutView *view = vtkGraphLayoutView::New();
SerialFMM FMM;

double get_time() {
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return double(tv.tv_sec+tv.tv_usec*1e-6);
}

void readNodes(std::string fileid) {
  std::string line, file="../fdgl/"+fileid+"_oNodes.csv";
  std::ifstream fid(file.c_str(),std::ios::in);
  fid >> line >> line >> line >> line >> line >> line;
  while( !fid.eof() ) {
    Vertex vertex;
    vertex.Id = graph->AddVertex();
    vertex.Ista = vertex.Iend = 0;
    vertex.X[0] = 2 * drand48() - 1;
    vertex.X[1] = 2 * drand48() - 1;
    vertex.X[2] = 2 * drand48() - 1;
    vertex.F = 0;
    vertices.push_back(vertex);
    fid >> line >> line >> line >> line >> line >> line;
  }
  fid.close();
  bodies.resize(vertices.size());
}

void readEdges(std::string fileid) {
  std::string line0, line1, line;
  std::string file="../fdgl/"+fileid+"_oEdges.csv";
  std::ifstream fid(file.c_str(),std::ios::in);
  fid >> line >> line0 >> line1 >> line >> line;
  while( !fid.eof() ) {
    edges.push_back(atoi(line1.c_str()));
    vertices[atoi(line0.c_str())].Iend = edges.size();
    fid >> line >> line0 >> line1 >> line >> line;
  }
  fid.close();
}

void edgeRange() {
  vertices.begin()->Ista = 0;
  for( V_iter V=vertices.begin()+1; V!=vertices.end(); ++V ) {
    if( (V-1)->Iend == 0 ) {
      if( V-vertices.begin() == 1 ) {
        (V-1)->Iend = 0;
      } else {
        (V-1)->Iend = (V-2)->Iend;
      }
    }
    V->Ista = (V-1)->Iend;
  }
}

void setEdges() {
  for( V_iter V=vertices.begin(); V!=vertices.end(); ++V ) {
    for( int i=V->Ista; i<V->Iend; ++i ) {
      graph->AddEdge(V->Id, vertices[edges[i]].Id);
    }
  }
}

void repulsion() {
  B_iter B = bodies.begin();
  for( V_iter V=vertices.begin(); V!=vertices.end(); ++V, ++B ) {
    B->X = V->X;
    B->SRC = -0.001;
    B->TRG = 0;
    B->IBODY = 0;
    B->IPROC = 0;
    B->ICELL = V-vertices.begin();
  }
  FMM.topdown(bodies,cells);
  FMM.approximate(cells);
  for( B=bodies.begin(); B!=bodies.end(); ++B ) {
    vertices[B->ICELL].F[0] = B->TRG[1];
    vertices[B->ICELL].F[1] = B->TRG[2];
    vertices[B->ICELL].F[2] = B->TRG[3];
  }
}

void spring() {
  float l = 2 / pow(vertices.size(),1./3);
  for( V_iter VI=vertices.begin(); VI!=vertices.end(); ++VI ) {
    for( int i=VI->Ista; i<VI->Iend; ++i ) {
      V_iter VJ = vertices.begin()+edges[i];
      vec<3,float> dist = VI->X - VJ->X;
      float R = sqrtf(norm(dist) + EPS2);
      VI->F -= dist / R * (R - l);
      VJ->F += dist / R * (R - l);
    }
  }
}

void central() {
  for( V_iter V=vertices.begin(); V!=vertices.end(); ++V ) {
    vec<3,float> dist = V->X;
    float R2 = norm(dist) + EPS2;
    float R3 = std::sqrt(R2) * R2;
    V->F -= dist / R3 * 0;
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

int main() {
  double t0, t[5] = {0,0,0,0,0};
  IMAGES = 0;
  THETA = 0.6;
  t0 = get_time();
  readNodes("test3");
  readEdges("test3");
  edgeRange();
  setEdges();
  t[0] += get_time() - t0;
  for( int step=0; step<1000; ++step ) {
    std::cout << step << std::endl;
    t0 = get_time();
    repulsion();
    t[1] += get_time() - t0;
    t0 = get_time();
    spring();
    t[2] += get_time() - t0;
    t0 = get_time();
    central();
    t[3] += get_time() - t0;
    t0 = get_time();
    moveVertices(step);
    t[4] += get_time() - t0;
    setVertices();
    drawGraph();
    usleep(1000);
  }
  std::cout << "N          : " << vertices.size() << std::endl;
  std::cout << "initialize : " << t[0] << std::endl;
  std::cout << "repulsion  : " << t[1] << std::endl;
  std::cout << "spring     : " << t[2] << std::endl;
  std::cout << "central    : " << t[3] << std::endl;
  std::cout << "move       : " << t[4] << std::endl;
  finalizeGraph();
}
