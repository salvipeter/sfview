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
  NurbsSurface(std::ifstream &in);
  ~NurbsSurface();
  static bool load(std::string const &filename, SurfacePVector &sv,
		   int texwidth, int texheight);
  void GLInit();
  bool showControlNet() const { return show_control_net; }
  void toggleShowControlNet() { show_control_net = !show_control_net; }
  void setVisualization(Visualization const v) {
    if(v != POINTS && v != WIREFRAME)
      vis = v;
  }
  void increaseDensity();
  void decreaseDensity();
  void display(Point const &eye_pos, bool high_density);
private:
  static void findOpenParen(std::ifstream &in);
  static bool isCloseParen(std::ifstream &in);
  static double readLispFloat(std::ifstream &in);
  static void ignoreWhitespaces(std::ifstream &in);
  void generateTexture(GLuint &name,
		       void (NurbsSurface::*fn)(unsigned char *) const);
  void generateIsophoteTexture(unsigned char *data) const;

  int degree_u, degree_v, texture_width, texture_height;
  double isophote_width;
  std::vector<float> knots_u, knots_v, linear_cpts;
  int nu, nv;
  PointVector control_net;
  bool show_control_net;
  GLUnurbsObj *globj, *gltextobj;
  GLuint mean_texture, gauss_texture, isophote_texture, slicing_texture;
  static GLuint default_isophote_texture;
  GLfloat texknots_u[4], texknots_v[4];
  static GLfloat texcpts[2][2][2];
};

#endif // NURBS_SURFACE_HH
