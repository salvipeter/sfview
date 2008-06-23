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
  void zoomToBoundingBox(Box const &b);
  void updateMatrices();
  std::string activeName(bool capital);
  void changeVisualization(Visualization v);
  void saveScreenShot(std::string filename);
  Point getObjectCoordinates(int x, int y);

  // Parameters
  int width, height;
  size_t max_n_of_quads;
  // Variables
  SurfacePVector surfaces;
  Point center, eye_pos;
  double modelview[16], inverse[16], object_width;
  bool high_density;
  size_t active;
  int mouse_start[2], next_id;
  enum MouseMode { NOTHING, ROTATION, ZOOM, PAN } mouse_mode;
  // Constants
  static std::string const help_string, copyright_string;
};

#endif	// GL_WINDOW_HH
