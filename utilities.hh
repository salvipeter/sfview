// SFView - Surface File Viewer
//
// Copyright (C) 2007-2008 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#ifndef UTILITIES_HH
#define UTILITIES_HH

#include "common.hh"

inline Point minPoint(Point const &p1, Point const &p2);
inline Point maxPoint(Point const &p1, Point const &p2);
Box boxUnion(Box const &b1, Box const &b2);
inline double interpolate(double d1, double t, double d2);
Point affineCombine(Point const &p1, double t, Point const &p2);
Point rotatePoint(Point const &p, Point const &center, Vector const &axis,
		  double theta);

#endif // UTILITIES_HH
