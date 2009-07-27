// SFView - Surface File Viewer
//
// Copyright (C) 2007-2008 Peter Salvi <vukung@yahoo.com>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#include "sfview.hh"

double version = 0.35;

GLWindow window;

// Dispatch functions (GLUT cannot handle member functions)
void display() { window.display(); }
void keyboard(unsigned char key, int x, int y) { window.keyboard(key, x, y); }
void mouseButton(int button, int state, int x, int y) {
  window.mouseButton(button, state, x, y);
}
void mouseMotion(int x, int y) { window.mouseMotion(x, y); }
void reshape(int w, int h) { window.reshape(w, h); }

int main(int argc, char *argv[])
{
  window.init(argc, argv);

  window.show();

  return 0;
}
