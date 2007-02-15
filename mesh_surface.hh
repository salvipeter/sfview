// SFView - Surface File Viewer
//
// Copyright (C) 2007 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#ifndef MESH_SURFACE_HH
#define MESH_SURFACE_HH

#include "surface.hh"

class MeshSurface : public Surface {
public:
  MeshSurface(std::string filename);
  bool showControlNet() const { return false; }
  void toggleShowControlNet() { }
  void setVisualization(Visualization const v) { vis = v; }

  void increaseDensity();
  void decreaseDensity();
  void display() const;
private:
  void approximateNormalsAndCurvatures();
  size_t decimation;
  PointVector points;
  VectorVector normals;
  ValueVector Gauss_curvature;
  double Gauss_min;
  double Gauss_max;
  ValueVector mean_curvature;
  double mean_min;
  double mean_max;
  int isophote_density;
  double slicing_density;
};

#endif // MESH_SURFACE_HH
