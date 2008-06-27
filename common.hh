// SFView - Surface File Viewer
//
// Copyright (C) 2007-2008 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#ifndef COMMON_HH
#define COMMON_HH

#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

typedef std::vector<std::string> StringVector;

class Vector {
public:
  Vector(double a = 0.0, double b = 0.0, double c = 0.0) {
    v[0] = a; v[1] = b; v[2] = c;
  }
  double &operator[](int i) { return v[i]; }
  double operator[](int i) const { return v[i]; }
  Vector operator+(Vector const &w) const {
    return Vector(v[0] + w[0], v[1] + w[1], v[2] + w[2]);
  }
  Vector operator-(Vector const &w) const {
    return Vector(v[0] - w[0], v[1] - w[1], v[2] - w[2]);
  }
  Vector operator*(double const t) const {
    return Vector(v[0] * t, v[1] * t, v[2] * t);
  }
  Vector operator/(double const t) const {
    return *this * (1.0 / t);
  }
  double operator*(Vector const &w) const {
    return v[0] * w[0] + v[1] * w[1] + v[2] * w[2];
  }
  Vector operator^(Vector const &w) const {
    return Vector(v[1] * w[2] - v[2] * w[1],
		  v[2] * w[0] - v[0] * w[2],
		  v[0] * w[1] - v[1] * w[0]);
  }
  double length() const {
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  }
  Vector normalized() const { return *this / length(); }
  friend std::ostream &operator<<(std::ostream &os, const Vector &v) {
    os << "<" << v[0] << ", " << v[1] << ", " << v[2] << ">";
    return os;
  }
private:
  double v[3];
};

typedef std::vector<Vector> VectorVector;
typedef VectorVector::iterator VectorIterator;

class Point {
public:
  Point(double a = 0.0, double b = 0.0, double c = 0.0) {
    p[0] = a; p[1] = b; p[2] = c;
  }
  double &operator[](int i) { return p[i]; }
  double operator[](int i) const { return p[i]; }
  Point operator+(Vector const &v) const {
    return Point(p[0] + v[0], p[1] + v[1], p[2] + v[2]);
  }
  Vector operator-(Point const &q) const {
    return Vector(p[0] - q[0], p[1] - q[1], p[2] - q[2]);
  }
  friend std::ostream &operator<<(std::ostream &os, const Point &p) {
    os << "<" << p[0] << ", " << p[1] << ", " << p[2] << ">";
    return os;
  }
private:
  double p[3];
};
  
typedef std::vector<Point> PointVector;
typedef PointVector::iterator PointIterator;

typedef double Value;
typedef std::vector<Value> ValueVector;
typedef ValueVector::iterator ValueIterator;

typedef std::pair<Point, Point> Box;

Vector const slicing_direction(0.0, 0.0, 1.0);

#endif // COMMON_HH
