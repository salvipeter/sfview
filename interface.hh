// SFView - Surface File Viewer
//
// Copyright (C) 2007 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#ifndef INTERFACE_HH
#define INTERFACE_HH

#include "renderer.hh"

enum MouseMode { NOTHING, ROTATION, ZOOM, PAN };

class Interface {
public:
private:
  size_t active = 0;
  int mouse_start[2], next_id = -1;
  MouseMode mouse_mode = NOTHING;
};

#endif // INTERFACE_HH
