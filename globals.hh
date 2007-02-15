// SFView - Surface File Viewer
//
// Copyright (C) 2007 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#ifndef GLOBALS_HH
#define GLOBALS_HH

#include "surface.hh"

extern int width, height;
extern int nurbs_density;
extern size_t max_n_of_quads;

extern SurfacePVector surfaces;
extern Point center, eye_pos;
extern double modelview[16], inverse[16], object_width;
extern bool high_density;

extern size_t active;
extern int mouse_start[2], next_id;
extern MouseMode mouse_mode;

#endif // GLOBALS_HH
