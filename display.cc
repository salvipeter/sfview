// SFView - Surface File Viewer
//
// Copyright (C) 2007 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#include <algorithm>
#include <fstream>

#include <GL/glut.h>

#include "globals.hh"
#include "utilities.hh"

void saveScreenShot(std::string filename)
{
  std::vector<char> buffer;
  buffer.resize(width * height * 3);

  glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, &buffer[0]);

  std::ofstream f(filename.c_str());
  f << "P6" << std::endl;
  f << width << " " << height << " 255" << std::endl;
  for(int j = height - 1; j >= 0; --j)
    for(int i = 0; i < width; ++i)
      f << buffer[3 * (j * width + i)]
	<< buffer[3 * (j * width + i) + 1]
	<< buffer[3 * (j * width + i) + 2];
  f << std::endl;
  f.close();
}

void updateMatrices()
{
  glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
  invert4x4Matrix(modelview, inverse);
  eye_pos = Point(-inverse[8], -inverse[9], -inverse[10]);
}

void isophoteColor(Point p, Vector n, int d)
{
  if((int)(std::acos((p - eye_pos).normalized() * n) * (double)d) % 2 == 0)
    glColor3d(1.0, 0.0, 0.0);
  else
    glColor3d(1.0, 1.0, 1.0);
}

void slicingColor(Point p, double d)
{
  Vector posvec(p[0], p[1], p[2]);
  Vector direction = Vector(eye_pos[0], eye_pos[1], eye_pos[2]).normalized();

  if((int)(posvec * direction * d) % 2 == 0)
    glColor3d(1.0, 0.0, 0.0);
  else
    glColor3d(0.0, 0.0, 1.0);
}

void rainbowColor(double value, double min, double max)
{
  double const d = (value - min) / (max - min);
  glColor3d(d, 1.0 - d, 0.0);
}

void display()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  for(SurfacePIterator i = surfaces.begin(); i != surfaces.end(); ++i)
    if(!(*i)->isHidden())
      (*i)->display();
  glutSwapBuffers();
}

void reshape(int w, int h)
{
  if(height == 0)
    height = 1;

  width = w;
  height = h;
  object_width = (double)w / (double)h;
  
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-object_width / 2.0, object_width / 2.0, -0.5, 0.5, -10.0, 10.0);

  glMatrixMode(GL_MODELVIEW);
}

void zoomToBoundingBox(Box const &b)
{
  center = affineCombine(b.first, 0.5, b.second);
  Vector const tmp = b.second - b.first;
  double const scaling = 0.9 / std::max(std::max(tmp[0], tmp[1]), tmp[2]);

  glLoadIdentity();
  glScaled(scaling, scaling, scaling);
  glTranslated(-center[0], -center[1], -center[2]);
}

void init()
{
  object_width = (double)width / (double)height;

  GLfloat light0_position[] = { 0.0, 10.0, 0.0, 0.0 };
  GLfloat light1_position[] = { 0.0, 0.0, -10.0, 0.0 };
  GLfloat light1_ambient[] = { 0.4, 0.4, 0.4, 1.0 };

  glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
  glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
  glLightfv(GL_LIGHT0, GL_AMBIENT, light1_ambient);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glDepthFunc(GL_LEQUAL);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_COLOR_MATERIAL);

  Box b = surfaces.front()->boundingBox();
  for(SurfacePIterator i = surfaces.begin(); i != surfaces.end(); ++i)
    b = boxUnion(b, (*i)->boundingBox());
  zoomToBoundingBox(b);

  updateMatrices();
}
