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
#ifndef vtk_h
#define vtk_h
#define VTK_EXCLUDE_STRSTREAM_HEADERS
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkSliderRepresentation2D.h>
#include <vtkSliderWidget.h>
#include <vtkWidgetEvent.h>
#include <vtkWidgetEventTranslator.h>
#include <vtkCallbackCommand.h>
#include <vtkCommand.h>
#include "types.h"
const int maxGroups = 1000000;

//! Interactive VTK class
class vtkSliderCallback : public vtkCommand {
public:
  vtkPoints *points[maxGroups];
  vtkPolyData *polydata;
  vtkVertexGlyphFilter *filter;
  vtkSliderCallback() {}
  static vtkSliderCallback *New() {
    return new vtkSliderCallback();
  }
  virtual void Execute(vtkObject *caller, unsigned long, void*) {
    vtkSliderWidget *widget = reinterpret_cast<vtkSliderWidget*>(caller);
    int value = static_cast<int>(static_cast<vtkSliderRepresentation *>(widget->GetRepresentation())->GetValue());
    polydata->SetPoints(points[value]);
    filter->SetInputConnection(polydata->GetProducerPort());
    filter->Update();
  }
};

//! Base VTK class
class vtkPlot {
  int I[maxGroups];
  vtkPoints *points[maxGroups];
  vtkPoints *hexPoints;
public:
  void setDomain(const real r0, const vect x0) {
    hexPoints = vtkPoints::New();
    hexPoints->SetNumberOfPoints(8);
    hexPoints->SetPoint(0, x0[0]-r0, x0[1]-r0, x0[2]-r0);
    hexPoints->SetPoint(1, x0[0]+r0, x0[1]-r0, x0[2]-r0);
    hexPoints->SetPoint(2, x0[0]+r0, x0[1]+r0, x0[2]-r0);
    hexPoints->SetPoint(3, x0[0]-r0, x0[1]+r0, x0[2]-r0);
    hexPoints->SetPoint(4, x0[0]-r0, x0[1]-r0, x0[2]+r0);
    hexPoints->SetPoint(5, x0[0]+r0, x0[1]-r0, x0[2]+r0);
    hexPoints->SetPoint(6, x0[0]+r0, x0[1]+r0, x0[2]+r0);
    hexPoints->SetPoint(7, x0[0]-r0, x0[1]+r0, x0[2]+r0);
  }

  void setGroup(const int Igroup, const int Npoints) {
    I[Igroup] = 0;
    points[Igroup] = vtkPoints::New();
    points[Igroup]->SetNumberOfPoints(Npoints);
  }

  void setPoints(const int Igroup, const vect X) {
    points[Igroup]->SetPoint(I[Igroup],X[0],X[1],X[2]);
    I[Igroup]++;
  }

  void setGroupOfPoints(Bodies &bodies, int &Ncell) {
    int begin=0, size=0;
    bigint index = bodies[0].ICELL;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      if( B->ICELL != index ) {
        setGroup(Ncell,size);
        for( int i=begin; i!=begin+size; ++i ) {
          setPoints(Ncell,bodies[i].X);
        }
        begin = B-bodies.begin();
        size = 0;
        index = B->ICELL;
        Ncell++;
        assert(Ncell < maxGroups);
      }
      size++;
    }
    setGroup(Ncell,size);
    for( int i=begin; i!=begin+size; ++i ) {
      setPoints(Ncell,bodies[i].X);
    }
    Ncell++;
    assert(Ncell < maxGroups);
  }

  void plot(const int Ngroup) {
    //Create a polygon object for points
    vtkPolyData *polydata = vtkPolyData::New();
    polydata->SetPoints(points[0]);

    //Create a filter object for points
    vtkVertexGlyphFilter *filter = vtkVertexGlyphFilter::New();
    filter->SetInputConnection(polydata->GetProducerPort());
    filter->Update();

    //Create a mapper object for points
    vtkPolyDataMapper *pointMapper = vtkPolyDataMapper::New();
    pointMapper->SetInputConnection(filter->GetOutputPort());

    //Associate the mapper to an actor object for points
    vtkActor *pointActor = vtkActor::New();
    pointActor->SetMapper(pointMapper);
    pointActor->GetProperty()->SetColor(1,0,0);

    //Create a hexahedron for cells
    vtkHexahedron *hex = vtkHexahedron::New();
    hex->GetPointIds()->SetId(0,0);
    hex->GetPointIds()->SetId(1,1);
    hex->GetPointIds()->SetId(2,2);
    hex->GetPointIds()->SetId(3,3);
    hex->GetPointIds()->SetId(4,4);
    hex->GetPointIds()->SetId(5,5);
    hex->GetPointIds()->SetId(6,6);
    hex->GetPointIds()->SetId(7,7);

    //Create a grid for cells
    vtkUnstructuredGrid *grid = vtkUnstructuredGrid::New();
    grid->Allocate(1,1);
    grid->InsertNextCell(hex->GetCellType(),hex->GetPointIds());
    grid->SetPoints(hexPoints);

    //Create a mapper object for cells
    vtkDataSetMapper *hexMapper = vtkDataSetMapper::New();
    hexMapper->SetInput(grid);

    //Associate the mapper to an actor object for cells
    vtkActor *hexActor = vtkActor::New();
    hexActor->SetMapper(hexMapper);
    hexActor->GetProperty()->SetOpacity(.1);

    //Add that actor to the renderer
    vtkRenderer *renderer = vtkRenderer::New();
    renderer->AddActor(pointActor);
    renderer->AddActor(hexActor);
    renderer->SetBackground(0,0,0);

    //Create a render window
    vtkRenderWindow *window = vtkRenderWindow::New();
    window->AddRenderer(renderer);
    window->SetSize(700,700);

    //Create an interactor and associate it to the render window
    vtkRenderWindowInteractor *interactor = vtkRenderWindowInteractor::New();
    interactor->SetRenderWindow(window);

    //Create a slider representation
    vtkSliderRepresentation2D *representation = vtkSliderRepresentation2D::New();
    representation->SetMinimumValue(0);
    representation->SetMaximumValue(Ngroup-1);
    representation->GetPoint1Coordinate()->SetCoordinateSystemToDisplay();
    representation->GetPoint1Coordinate()->SetValue(50,50);
    representation->GetPoint2Coordinate()->SetCoordinateSystemToDisplay();
    representation->GetPoint2Coordinate()->SetValue(650,50);

    //Create a slider widget
    vtkSliderWidget *widget = vtkSliderWidget::New();
    widget->SetInteractor(interactor);
    widget->SetRepresentation(representation);
    widget->SetAnimationModeToAnimate();
    widget->EnabledOn();

    //Create a slider callback
    vtkSliderCallback *callback = vtkSliderCallback::New();
    for( int i=0; i!=Ngroup; ++i ) {
      callback->points[i] = points[i];
    }
    callback->polydata = polydata;
    callback->filter = filter;
    widget->AddObserver(vtkCommand::InteractionEvent,callback);

    //Define the interacting style
    vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
    interactor->SetInteractorStyle(style);

    //Start to interact
    interactor->Initialize();
    interactor->Start();
  }
};

#endif
