// SFView - Surface File Viewer
//
// Copyright (C) 2007-2008 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#ifndef NURBS_SURFACE_HH
#define NURBS_SURFACE_HH

#include <GL/glut.h>

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
  void display(Point const &eye_pos, bool); // parameter high_density ignored
private:
  void findOpenParen(std::ifstream &in) const;
  bool isCloseParen(std::ifstream &in) const;

  int degree_u, degree_v;
  std::vector<float> knots_u, knots_v, linear_cpts;
  int nu, nv;
  PointVector control_net;
  bool show_control_net;
  GLUnurbsObj *globj;
};

#endif // NURBS_SURFACE_HH