// SFView - Surface File Viewer
//
// Copyright (C) 2007 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#ifndef NURBS_SURFACE_HH
#define NURBS_SURFACE_HH

#include "surface.hh"

class NurbsSurface : public Surface {
public:
  NurbsSurface(std::string filename);
  bool showControlNet() const { return show_control_net; }
  void toggleShowControlNet() { show_control_net = !show_control_net; }
  void setVisualization(Visualization const v) {
    if(v != POINTS && v != WIREFRAME)
      vis = v;
  }

  void increaseDensity();
  void decreaseDensity();
  void display() const;
private:
  int degree_u, degree_v;
  ValueVector knots_u, knots_v;
  PointVector control_net;
  bool show_control_net;
};

#endif // NURBS_SURFACE_HH
