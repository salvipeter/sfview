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
  ControlNetType showControlNet() const { return CN_NONE; }
  void toggleControlNet() { }
  void setVisualization(Visualization const v) { vis = v; }

  void increaseDensity();
  void decreaseDensity();
  void display(Point const &eye_pos, Vector const &eye_dir, bool high_density);
private:
  void approximateNormalsAndCurvatures();
  void isophoteColor(Point const &p, Vector const &n, Point const &eye_pos);
  void slicingColor(Point const &p, Point const &eye_pos, Vector const &eye_dir,
		    int n, unsigned char const data[][3]);
  void contourColor(Point const &p, Point const &eye_pos,
		    Vector const &eye_dir);
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
};

#endif // MESH_SURFACE_HH
