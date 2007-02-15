// SFView - Surface File Viewer
//
// Copyright (C) 2007 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#ifndef RENDERER_HH
#define RENDERER_HH

enum Visualization { SHADED, GAUSS, MEAN, ISOPHOTE,
		     SLICING, WIREFRAME, POINTS };

class Renderer {
public:
  void saveScreenShot(std::string filename);
  void updateMatrices();
private:
  static void display();
  static void reshape(int w, int h);
  void init();
  void zoomToBoundingBox(Box const &b);

  int width, height;
  SurfacePVector surfaces;
  Point center, eye_pos;
  double modelview[16], inverse[16], object_width;
  bool high_density = false;
};

#endif // RENDERER_HH
