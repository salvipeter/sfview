// SFView - Surface File Viewer
//
// Copyright (C) 2007-2008 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#ifndef MESH_SURFACE_HH
#define MESH_SURFACE_HH

#include "surface.hh"

class MeshSurface : public Surface {
public:
  MeshSurface(std::string filename, size_t max_n_of_quads);
  bool showControlNet() const { return false; }
  void toggleShowControlNet() { }
  void setVisualization(Visualization const v) { vis = v; }

  void increaseDensity();
  void decreaseDensity();
  void display(Point const &eye_pos, bool high_density);
private:
  void approximateNormalsAndCurvatures();
  void isophoteColor(Point const &p, Vector const &n, int d,
		     Point const &eye_pos);
  void slicingColor(Point const &p, double d, Point const &eye_pos);
  void rainbowColor(double value, double min, double max);

  size_t decimation;
  size_t resx, resy;
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
