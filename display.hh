// SFView - Surface File Viewer
//
// Copyright (C) 2007 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#ifndef DISPLAY_HH
#define DISPLAY_HH

#include <string>

void init();
void display();
void reshape(int w, int h);
void updateMatrices();
void saveScreenShot(std::string filename);

#endif // DISPLAY_HH
