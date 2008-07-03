// SFView - Surface File Viewer
//
// Copyright (C) 2007-2008 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#ifndef GL_WINDOW_HH
#define GL_WINDOW_HH

#include "common.hh"
#include "surface.hh"

class GLWindow {
public:
  GLWindow();
  void init(int argc, char *argv[]);
  void show();
  void display();
  void keyboard(unsigned char key, int x, int y);
  void mouseButton(int button, int state, int x, int y);
  void mouseMotion(int x, int y);
  void reshape(int w, int h);
private:
  StringVector parseCommandLine(int argc, char *argv[]);
  bool loadFile(std::string filename);
  bool reloadActive();
  void zoomToBoundingBox();
  Point boundingBoxPoint(int i);
  void setClippingPlanes();
  std::string activeName(bool capital);
  void changeVisualization(Visualization v);
  void saveScreenShot(std::string filename);

  // Parameters
  int width, height, texture_width, texture_height;
  size_t max_n_of_quads;
  // Variables
  SurfacePVector surfaces;
  Box bounding_box;
  Point center, eye;
  Vector up;
  double object_width;
  bool high_density;
  size_t active;
  int mouse_start[2], next_id;
  enum MouseMode { NOTHING, ROTATION, ZOOM, PAN } mouse_mode;
  // Constants
  static std::string const help_string, copyright_string;
  static double const view_angle, znear_coefficient, zfar_coefficient;
};

#endif	// GL_WINDOW_HH
