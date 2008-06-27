// SFView - Surface File Viewer
//
// Copyright (C) 2007-2008 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#include "utilities.hh"

inline Point minPoint(Point const &p1, Point const &p2)
{
  return Point(p1[0] < p2[0] ? p1[0] : p2[0],
	       p1[1] < p2[1] ? p1[1] : p2[1],
	       p1[2] < p2[2] ? p1[2] : p2[2]);
}

inline Point maxPoint(Point const &p1, Point const &p2)
{
  return Point(p1[0] > p2[0] ? p1[0] : p2[0],
	       p1[1] > p2[1] ? p1[1] : p2[1],
	       p1[2] > p2[2] ? p1[2] : p2[2]);
}

Box boxUnion(Box const &b1, Box const &b2)
{
  Box result;
  result.first = minPoint(b1.first, b2.first);
  result.second = maxPoint(b1.second, b2.second);
  return result;
}

inline double interpolate(double d1, double t, double d2)
{
  return d1 + (d2 - d1) * t;
}

Point affineCombine(Point const &p1, double t, Point const &p2)
{
  return Point(interpolate(p1[0], t, p2[0]),
	       interpolate(p1[1], t, p2[1]),
	       interpolate(p1[2], t, p2[2]));
}

Point rotatePoint(Point const &p, Point const &center, Vector const &axis,
		  double theta)
{
  double const c = cos(theta), s = sin(theta);
  double const C = 1.0 - c;
  Vector tmp = p - center;
  return center +
    Vector(tmp[0] * (axis[0] * axis[0] * C + c) +
	   tmp[1] * (axis[0] * axis[1] * C - axis[2] * s) +
	   tmp[2] * (axis[0] * axis[2] * C + axis[1] * s),
	   tmp[0] * (axis[1] * axis[0] * C + axis[2] * s) +
	   tmp[1] * (axis[1] * axis[1] * C + c) +
	   tmp[2] * (axis[1] * axis[2] * C - axis[0] * s),
	   tmp[0] * (axis[2] * axis[0] * C - axis[1] * s) +
	   tmp[1] * (axis[2] * axis[1] * C + axis[0] * s) +
	   tmp[2] * (axis[2] * axis[2] * C + c));
}
