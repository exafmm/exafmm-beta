#ifndef vtk_h
#define vtk_h
#define VTK_EXCLUDE_STRSTREAM_HEADERS
#include <string>
#include "types.h"
#include <vector>
#include <vtkAxis.h>
#include <vtkChartXY.h>
#include <vtkCommand.h>
#include <vtkContextScene.h>
#include <vtkContextView.h>
#include <vtkDataSetMapper.h>
#include <vtkFloatArray.h>
#include <vtkHexahedron.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkPlot.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSliderRepresentation2D.h>
#include <vtkSliderWidget.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertexGlyphFilter.h>

//! 2D plot VTK class
class vtk2DPlot {
private:
  int numColumns;

public:
  vtkTable * table;

public:
  vtk2DPlot() : numColumns(0) {
    table = vtkTable::New();
  }

  void setNumRows(int numRows) {
    table->SetNumberOfRows(numRows);
  }

  void setName(std::string name) {
    vtkFloatArray * array = vtkFloatArray::New();
    array->SetName(name.c_str());
    table->AddColumn(array);
  }

  template<typename T>
  void setData(int begin, int end, T * data) {
    for (int i=begin; i<end; i++) {
      table->SetValue(i-begin,numColumns,data[i]);
    }
    numColumns++;
  }

  void plot() {
    //Set xlabel, ylabel, title
    vtkChartXY * chart = vtkChartXY::New();
    chart->GetAxis(vtkAxis::LEFT)->SetTitle("y");
    chart->GetAxis(vtkAxis::BOTTOM)->SetTitle("x");
    chart->SetShowLegend(true);

    // Add multiple line plots, setting the colors etc
    for (int i=1; i<numColumns; i++) {
      vtkPlot * line = chart->AddPlot(vtkChart::LINE);
      line->SetInput(table, 0, i);
    }

    // Set up the view
    vtkContextView * view = vtkContextView::New();
    view->GetScene()->AddItem(chart);
    view->GetInteractor()->Initialize();
    view->GetInteractor()->Start();
  }
};

//! 3D plot VTK class
class vtk3DPlot {
private:
  std::vector<vtkPoints*> groups;
  vtkPoints * hexPoints;

  struct vtkSliderCallback : public vtkCommand {
    std::vector<vtkPoints*> groups;
    vtkPolyData * polydata;
    vtkVertexGlyphFilter * filter;
    vtkSliderCallback() {}
    static vtkSliderCallback * New() {
      return new vtkSliderCallback();
    }
    virtual void Execute(vtkObject * caller, unsigned long, void*) {
      vtkSliderWidget * widget = reinterpret_cast<vtkSliderWidget*>(caller);
      int igroup = static_cast<int>(static_cast<vtkSliderRepresentation*>(widget->GetRepresentation())->GetValue());
      polydata->SetPoints(groups[igroup]);
      filter->SetInputConnection(polydata->GetProducerPort());
      filter->Update();
    }
  };

public:
  void setBounds(const real_t r0, const vec3 x0) {
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

  void setPoints(B_iter B0, B_iter BN) {
    vtkPoints * group = vtkPoints::New();
    group->SetNumberOfPoints(BN-B0);
    for (B_iter B=B0; B!=BN; B++) {
      group->SetPoint(B-B0, B->X[0], B->X[1], B->X[2]);
    }
    groups.push_back(group);
  }

  void setGroupOfPoints(Bodies & bodies) {
    B_iter B0 = bodies.begin(); 
    B_iter BN = B0; 
    int index = B0->IBODY;
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++, BN++) {
      if (B->IBODY != index) {
        setPoints(B0, BN);
        B0 = BN = B;
        index = B->IBODY;
      }
    }
    setPoints(B0, BN);
  }

  void plot() {
    //Create a polygon object for points
    vtkPolyData * polydata = vtkPolyData::New();
    polydata->SetPoints(groups[0]);

    //Create a filter object for points
    vtkVertexGlyphFilter * filter = vtkVertexGlyphFilter::New();
    filter->SetInputConnection(polydata->GetProducerPort());
    filter->Update();

    //Create a mapper object for points
    vtkPolyDataMapper * pointMapper = vtkPolyDataMapper::New();
    pointMapper->SetInputConnection(filter->GetOutputPort());

    //Associate the mapper to an actor object for points
    vtkActor * pointActor = vtkActor::New();
    pointActor->SetMapper(pointMapper);
    pointActor->GetProperty()->SetColor(1, 0, 0);

    //Create a hexahedron for cells
    vtkHexahedron * hex = vtkHexahedron::New();
    hex->GetPointIds()->SetId(0, 0);
    hex->GetPointIds()->SetId(1, 1);
    hex->GetPointIds()->SetId(2, 2);
    hex->GetPointIds()->SetId(3, 3);
    hex->GetPointIds()->SetId(4, 4);
    hex->GetPointIds()->SetId(5, 5);
    hex->GetPointIds()->SetId(6, 6);
    hex->GetPointIds()->SetId(7, 7);

    //Create a grid for cells
    vtkUnstructuredGrid * grid = vtkUnstructuredGrid::New();
    grid->Allocate(1, 1);
    grid->InsertNextCell(hex->GetCellType(), hex->GetPointIds());
    grid->SetPoints(hexPoints);

    //Create a mapper object for cells
    vtkDataSetMapper * hexMapper = vtkDataSetMapper::New();
    hexMapper->SetInput(grid);

    //Associate the mapper to an actor object for cells
    vtkActor * hexActor = vtkActor::New();
    hexActor->SetMapper(hexMapper);
    hexActor->GetProperty()->SetOpacity(.1);

    //Add that actor to the renderer
    vtkRenderer * renderer = vtkRenderer::New();
    renderer->AddActor(pointActor);
    renderer->AddActor(hexActor);
    renderer->SetBackground(0, 0, 0);

    //Create a render window
    vtkRenderWindow * window = vtkRenderWindow::New();
    window->AddRenderer(renderer);
    window->SetSize(700, 700);

    //Create an interactor and associate it to the render window
    vtkRenderWindowInteractor * interactor = vtkRenderWindowInteractor::New();
    interactor->SetRenderWindow(window);

    //Create a slider representation
    vtkSliderRepresentation2D * representation = vtkSliderRepresentation2D::New();
    representation->SetMinimumValue(0);
    representation->SetMaximumValue(groups.size()-1);
    representation->GetPoint1Coordinate()->SetCoordinateSystemToDisplay();
    representation->GetPoint1Coordinate()->SetValue(50, 50);
    representation->GetPoint2Coordinate()->SetCoordinateSystemToDisplay();
    representation->GetPoint2Coordinate()->SetValue(650, 50);

    //Create a slider widget
    vtkSliderWidget * widget = vtkSliderWidget::New();
    widget->SetInteractor(interactor);
    widget->SetRepresentation(representation);
    widget->SetAnimationModeToAnimate();
    widget->EnabledOn();

    //Create a slider callback
    vtkSliderCallback * callback = vtkSliderCallback::New();
    callback->groups = groups;
    callback->polydata = polydata;
    callback->filter = filter;
    widget->AddObserver(vtkCommand::InteractionEvent, callback);

    //Define the interacting style
    vtkInteractorStyleTrackballCamera * style = vtkInteractorStyleTrackballCamera::New();
    interactor->SetInteractorStyle(style);

    //Start to interact
    interactor->Initialize();
    interactor->Start();
  }
};
#endif
