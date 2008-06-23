// SFView - Surface File Viewer
//
// Copyright (C) 2007-2008 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#ifndef SFVIEW_HH
#define SFVIEW_HH

#include "glwindow.hh"

extern double version;
extern GLWindow window;

void display();
void keyboard(unsigned char key, int x, int y);
void mouseButton(int button, int state, int x, int y);
void mouseMotion(int x, int y);
void reshape(int w, int h);

#endif	// SFVIEW_HH
