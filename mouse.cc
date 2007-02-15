// SFView - Surface File Viewer
//
// Copyright (C) 2007 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#include <GL/glut.h>

#include "display.hh"
#include "globals.hh"

Point getObjectCoordinates(int x, int y)
{
  Vector v;
  Point result;

  v[0] = ((double)x / (double)width - 0.5) * object_width;
  v[1] = (double)(height - y) / (double)height - 0.5;
  v[2] = -10.0;

  result[0] = v * Vector(inverse[0], inverse[4], inverse[8]);
  result[1] = v * Vector(inverse[1], inverse[5], inverse[9]);
  result[2] = v * Vector(inverse[2], inverse[6], inverse[10]);

  return result;
}

void mouseButton(int button, int state, int x, int y)
{
  if(state == GLUT_UP) {
    mouse_mode = NOTHING;
    display();
  }
  else if(state == GLUT_DOWN) {
    mouse_start[0] = x; mouse_start[1] = y;
    switch(button) {
    case GLUT_LEFT_BUTTON : mouse_mode = ROTATION; break;
    case GLUT_RIGHT_BUTTON : mouse_mode = ZOOM; break;
    case GLUT_MIDDLE_BUTTON : mouse_mode = PAN; break;
    }
  }
}

void mouseMotion(int x, int y)
{
  if(mouse_mode == NOTHING)
    return;

  Point oldp, p;
  Vector diff, axis, rotmp;
  double scaling, theta;

  switch(mouse_mode) {
  case ROTATION :
    rotmp = Vector(y - mouse_start[1], x - mouse_start[0], 0.0);
    theta = rotmp.length() / (double)height * 180.0;

    axis[0] = rotmp * Vector(inverse[0], inverse[4], inverse[8]);
    axis[1] = rotmp * Vector(inverse[1], inverse[5], inverse[9]);
    axis[2] = rotmp * Vector(inverse[2], inverse[6], inverse[10]);

    glTranslated(center[0], center[1], center[2]);
    glRotated(theta, axis[0], axis[1], axis[2]);
    glTranslated(-center[0], -center[1], -center[2]);
    break;
  case ZOOM :
    scaling = 1.0 / std::exp((y - mouse_start[1]) * 0.01);
    glTranslated(center[0], center[1], center[2]);
    glScaled(scaling, scaling, scaling);
    glTranslated(-center[0], -center[1], -center[2]);
    break;
  case PAN :
    oldp = getObjectCoordinates(mouse_start[0], mouse_start[1]);
    p = getObjectCoordinates(x, y);
    diff = oldp - p;
    center = center + diff;
    glTranslated(-diff[0], -diff[1], -diff[2]);
    break;
  default: break;
  }
  mouse_start[0] = x; mouse_start[1] = y;

  updateMatrices();

  display();
}
