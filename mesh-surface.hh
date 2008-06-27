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
  MeshSurface(std::ifstream &in, size_t max_n_of_quads);
  static bool load(std::string const &filename, SurfacePVector &sv,
		   size_t max_n_of_quads);
  bool showControlNet() const { return false; }
  void toggleShowControlNet() { }
  void setVisualization(Visualization const v) { vis = v; }

  void increaseDensity();
  void decreaseDensity();
  void display(Point const &eye_pos, bool high_density);
private:
  void approximateNormalsAndCurvatures();
  void isophoteColor(Point const &p, Vector const &n, Point const &eye_pos);
  void slicingColor(Point const &p, Point const &eye_pos);
  static void rainbowColor(double value, double min, double max);

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
  double isophote_width, slicing_density;
};

#endif // MESH_SURFACE_HH
