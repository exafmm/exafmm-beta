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
#include "vec.h"
 
struct Body {
  vec<3,float> X;
  vec<3,float> V;
  vec<3,float> F;
  vec<2,int>   I;
};
typedef std::vector<Body>::iterator B_iter;
std::vector<Body> bodies;
std::vector<vtkIdType> vertices;
std::vector<int> edges;
vtkMutableUndirectedGraph *graph = vtkMutableUndirectedGraph::New();
vtkGraphLayoutView *view = vtkGraphLayoutView::New();

double get_time() {
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return double(tv.tv_sec+tv.tv_usec*1e-6);
}

void readNodes(std::string fileid) {
  std::string line, file="../fdgl/"+fileid+"_oNodes.csv";
  std::ifstream fid(file.c_str(),std::ios::in);
  while( !fid.eof() ) {
    fid >> line >> line >> line >> line >> line >> line;
    vertices.push_back(graph->AddVertex());
  }
  fid.close();
}

void initBodies() {
  for( std::vector<vtkIdType>::iterator V=vertices.begin(); V!=vertices.end(); ++V ) {
    Body body;
    body.X[0] = 2 * drand48() - 1;
    body.X[1] = 2 * drand48() - 1;
    body.X[2] = 2 * drand48() - 1;
    body.V = 0;
    body.F = 0;
    body.I = 0;
    bodies.push_back(body);
  }
}

void readEdges(std::string fileid) {
  std::string line0, line1, line;
  std::string file="../fdgl/"+fileid+"_oEdges.csv";
  std::ifstream fid(file.c_str(),std::ios::in);
  while( !fid.eof() ) {
    fid >> line >> line0 >> line1 >> line >> line;
    edges.push_back(atoi(line1.c_str()));
    bodies[atoi(line0.c_str())].I[1] = edges.size();
  }
  fid.close();
}

void edgeRange() {
  bodies.begin()->I[0] = 0;
  for( B_iter B=bodies.begin()+1; B!=bodies.end(); ++B ) {
    if( (B-1)->I[1] == 0 ) {
      if( B-bodies.begin() == 1 ) {
        (B-1)->I[1] = 0;
      } else {
        (B-1)->I[1] = (B-2)->I[1];
      }
    }
    B->I[0] = (B-1)->I[1];
  }
}

void setEdges() {
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    for( int i=B->I[0]; i<B->I[1]; ++i ) {
      graph->AddEdge(vertices[B-bodies.begin()], vertices[edges[i]]);
    }
  }
}

void repulsion() {
  for( B_iter BI=bodies.begin(); BI!=bodies.end(); ++BI ) {
    vec<3,float> F = 0;
    for( B_iter BJ=bodies.begin(); BJ!=bodies.end(); ++BJ ) {
      vec<3,float> dist = BI->X - BJ->X;
      float R2 = norm(dist) + 1e-4;
      float R3 = std::sqrt(R2) * R2;
      F += dist / R3;
    }
    BI->F = F * 0.001;
  }
}

void spring() {
  for( B_iter BI=bodies.begin(); BI!=bodies.end(); ++BI ) {
    vec<3,float> F = 0;
    for( int i=BI->I[0]; i<BI->I[1]; ++i ) {
      B_iter BJ = bodies.begin()+edges[i];
      vec<3,float> dist = BI->X - BJ->X;
      F -= dist;
    }
    BI->F += F * 0.005;
  }
}

void central() {
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    vec<3,float> dist = B->X;
    float R2 = norm(dist) + 1e-4;
    float R3 = std::sqrt(R2) * R2;
    B->F -= dist / R3 * 0.001;
  }
}

void moveBodies() {
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    B->X += B->F;
  } 
}

void setVertices() {
  vtkPoints *points = vtkPoints::New();
  for( B_iter B=bodies.begin(); B!=bodies.end()-8; ++B ) {
    points->InsertNextPoint(B->X[0], B->X[1], B->X[2]);
  }
  points->InsertNextPoint(-1,-1,-1);
  points->InsertNextPoint( 1,-1,-1);
  points->InsertNextPoint( 1, 1,-1);
  points->InsertNextPoint(-1, 1,-1);
  points->InsertNextPoint(-1,-1, 1);
  points->InsertNextPoint( 1,-1, 1);
  points->InsertNextPoint( 1, 1, 1);
  points->InsertNextPoint(-1, 1, 1);
  graph->SetPoints(points);
}

void drawGraph(int step, int nstep) {
  view->AddRepresentationFromInput(graph);
  view->SetLayoutStrategy("Pass Through");
  view->GetInteractor()->GetRenderWindow()->SetSize(700,700);
  view->ResetCamera();
  view->Render();
  if( step == nstep ) {
    vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
    view->GetInteractor()->SetInteractorStyle(style);
    view->GetInteractor()->Initialize();
    view->GetInteractor()->Start();
  }
}

int main() {
  int nstep = 1;
  double t0, t[5] = {0,0,0,0,0};
  t0 = get_time();
  readNodes("test3");
  initBodies();
  readEdges("test3");
  edgeRange();
  setEdges();
  t[0] += get_time() - t0;
  for( int step=0; step<=nstep; ++step ) {
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
    moveBodies();
    t[4] += get_time() - t0;
    setVertices();
    if( step%1 == 0 ) drawGraph(step,nstep);
  }
  std::cout << "initialize : " << t[0] << std::endl;
  std::cout << "repulsion  : " << t[1] << std::endl;
  std::cout << "spring     : " << t[2] << std::endl;
  std::cout << "move       : " << t[3] << std::endl;
  std::cout << "central    : " << t[4] << std::endl;
}
