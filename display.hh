// SFView - Surface File Viewer
//
// Copyright (C) 2007 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#ifndef DISPLAY_HH
#define DISPLAY_HH

#include <string>

#include "common.hh"

void init();
void isophoteColor(Point p, Vector n, int d);
void slicingColor(Point p, double d);
void rainbowColor(double value, double min, double max);
void display();
void reshape(int w, int h);
void updateMatrices();
void saveScreenShot(std::string filename);

#endif // DISPLAY_HH
