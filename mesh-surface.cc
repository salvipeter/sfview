// SFView - Surface File Viewer
//
// Copyright (C) 2007-2008 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#include <algorithm>
#include <fstream>
#include <iostream>

#include <GL/glut.h>

#include "mesh-surface.hh"

void MeshSurface::approximateNormalsAndCurvatures()
{
  VectorVector u_derivatives, v_derivatives;

  normals.clear();
  mean_curvature.clear();
  Gauss_curvature.clear();

  for(size_t i = 0; i < resx; ++i)
    for(size_t j = 0; j < resy; ++j) {
      Vector du, dv;
      if(i == 0)
	du = points[resy + j] - points[j];
      else
	du = points[i * resy + j] - points[(i - 1) * resy + j];
      if(j == 0)
	dv = points[i * resy + 1] - points[i * resy];
      else
	dv = points[i * resy + j] - points[i * resy + j - 1];
      u_derivatives.push_back(du);
      v_derivatives.push_back(dv);
    }
  for(size_t i = 0; i < resx; ++i)
    for(size_t j = 0; j < resy; ++j) {
      Vector duu, duv, dvv;
      Vector const &du = u_derivatives[i * resy + j];
      Vector const &dv = v_derivatives[i * resy + j];
      if(i == 0)
	duu = u_derivatives[2 * resy + j] - u_derivatives[resy + j];
      else
	duu = du - u_derivatives[(i - 1) * resy + j];
      if(j == 0) {
	duv = u_derivatives[i * resy + 2] - u_derivatives[i * resy + 1];
	dvv = v_derivatives[i * resy + 2] - v_derivatives[i * resy + 1];
      } else {
	duv = du - u_derivatives[i * resy + j - 1];
	dvv = dv - v_derivatives[i * resy + j - 1];
      }
      Vector const normal = (du ^ dv).normalized();
      double const E = du * du;
      double const F = du * dv;
      double const G = dv * dv;
      double const L = normal * duu;
      double const M = normal * duv;
      double const N = normal * dvv;
      double const m = N * E - 2 * F * M + L * G;
      double const g = 2 * (L * N - M * M);
      double const a = 2 * (E * G - F * F);
      double const d = std::sqrt(m * m - a * g);
      double const k1 = (m + d) / a;
      double const k2 = (m - d) / a;

      // OpenGL seems to want the opposite direction... strange
      normals.push_back(normal * -1.0);
      Gauss_curvature.push_back(k1 * k2);
      mean_curvature.push_back((k1 + k2) * 0.5);
    }
  Gauss_min =
    *std::min_element(Gauss_curvature.begin(), Gauss_curvature.end());
  Gauss_max =
    *std::max_element(Gauss_curvature.begin(), Gauss_curvature.end());
  mean_min =
    *std::min_element(mean_curvature.begin(), mean_curvature.end());
  mean_max =
    *std::max_element(mean_curvature.begin(), mean_curvature.end());
}

bool MeshSurface::load(std::string const &filename, SurfacePVector &sv,
		       size_t max_n_of_quads)
{
  std::ifstream in(filename.c_str());
  std::cout << "Loading file `" << filename << "'... " << std::flush;
  MeshSurface *sf = new MeshSurface(in, max_n_of_quads);
  in.close();
  if(!sf->error) {
    sf->filename = filename;
    sv.push_back(sf);
  }
  return !sf->error;
}

MeshSurface::MeshSurface(std::ifstream &in, size_t max_n_of_quads) :
  isophote_width(2.0),
  slicing_density(0.5)
{
  Point p;

  in >> resx >> resy;
  for(size_t j = 0; j < resy; ++j)
    for(size_t i = 0; i < resx; ++i) {
      if(in.eof()) {
	error = true;
	return;
      }
      in >> p[0] >> p[1] >> p[2];
      points.push_back(p);
    }

  in.close();
  std::cout << "ok (" << resx << "x" << resy << " points)" << std::endl;

  std::cout << "Calculating bounding box... " << std::flush;
  bounding_box.first = bounding_box.second = points.front();
  for(PointIterator i = points.begin(); i != points.end(); ++i) {
    if((*i)[0] < bounding_box.first[0]) bounding_box.first[0] = (*i)[0];
    if((*i)[1] < bounding_box.first[1]) bounding_box.first[1] = (*i)[1];
    if((*i)[2] < bounding_box.first[2]) bounding_box.first[2] = (*i)[2];
    if((*i)[0] > bounding_box.second[0]) bounding_box.second[0] = (*i)[0];
    if((*i)[1] > bounding_box.second[1]) bounding_box.second[1] = (*i)[1];
    if((*i)[2] > bounding_box.second[2]) bounding_box.second[2] = (*i)[2];
  }
  std::cout << "ok (" << bounding_box.first << " - "
	    << bounding_box.second << ")" << std::endl;

  decimation = (size_t)std::ceil(std::sqrt((double)(resx * resy) /
					   (double)max_n_of_quads));

  std::cout << "Using a decimation factor of " << decimation
	    << " while moving." << std::endl;

  std::cout << "Approximating normals and curvatures... " << std::flush;
  approximateNormalsAndCurvatures();
  std::cout << "ok" << std::endl;

  error = false;
}

void MeshSurface::increaseDensity()
{
  switch(vis) {
  case SLICING : slicing_density *= 2.0; break;
  case ISOPHOTE : isophote_width /= 2.0; break;
  default: ;
  }
}

void MeshSurface::decreaseDensity()
{
  switch(vis) {
  case SLICING : slicing_density /= 2.0; break;
  case ISOPHOTE : isophote_width *= 2.0; break;
  default: ;
  }
}

void MeshSurface::isophoteColor(Point const &p, Vector const &n, double d,
				Point const &eye_pos)
{
  double const angle = std::acos((eye_pos - p).normalized() * n) * 180.0 / M_PI;
  if(static_cast<int>(std::floor(angle / isophote_width)) % 2 == 0)
    glColor3d(1.0, 0.0, 0.0);
  else
    glColor3d(1.0, 1.0, 1.0);
}

void MeshSurface::slicingColor(Point const &p, double d, Point const &eye_pos)
{
  Vector posvec(p[0], p[1], p[2]);
  Vector direction = Vector(eye_pos[0], eye_pos[1], eye_pos[2]).normalized();

  if((int)(posvec * direction * d) % 2 == 0)
    glColor3d(1.0, 0.0, 0.0);
  else
    glColor3d(0.0, 0.0, 1.0);
}

void MeshSurface::rainbowColor(double value, double min, double max)
{
  double const d = std::min(std::max((value - min) / (max - min), min), max);
  glColor3d(d, 1.0 - d, 0.0);
}

void MeshSurface::display(Point const &eye_pos, bool high_density)
{
  if(hidden)
    return;

  size_t dec_incr = high_density ? 1 : decimation;

  PointVector const &p = points;
  VectorVector const &n = normals;
  ValueVector const &g = Gauss_curvature;
  ValueVector const &m = mean_curvature;

  if(vis == GAUSS || vis == MEAN)
    glDisable(GL_LIGHTING);
  else
    glEnable(GL_LIGHTING);

  for(size_t j = 0; j < resy-1; j += dec_incr)
    for(size_t i = 0; i < resx-1; i += dec_incr) {
      size_t const next_i = std::min(i + dec_incr, resx - 1);
      size_t const next_j = std::min(j + dec_incr, resy - 1);
      size_t const i1 = i * resy + j;
      size_t const i2 = i * resy + next_j;
      size_t const i3 = next_i * resy + next_j;
      size_t const i4 = next_i * resy + j;

      switch(vis) {
      case SHADED :
	glBegin(GL_QUADS);
	glColor3d(1.0, 1.0, 1.0);
	glVertex3d(p[i1][0], p[i1][1], p[i1][2]);
	glNormal3d(n[i1][0], n[i1][1], n[i1][2]);
	glVertex3d(p[i2][0], p[i2][1], p[i2][2]);
	glNormal3d(n[i2][0], n[i2][1], n[i2][2]);
	glVertex3d(p[i3][0], p[i3][1], p[i3][2]);
	glNormal3d(n[i3][0], n[i3][1], n[i3][2]);
	glVertex3d(p[i4][0], p[i4][1], p[i4][2]);
	glNormal3d(n[i4][0], n[i4][1], n[i4][2]);
	glEnd();
	break;
      case GAUSS :
	glBegin(GL_QUADS);
	rainbowColor(g[i1], Gauss_min, Gauss_max);
	glVertex3d(p[i1][0], p[i1][1], p[i1][2]);
	glNormal3d(n[i1][0], n[i1][1], n[i1][2]);
	rainbowColor(g[i2], Gauss_min, Gauss_max);
	glVertex3d(p[i2][0], p[i2][1], p[i2][2]);
	glNormal3d(n[i2][0], n[i2][1], n[i2][2]);
	rainbowColor(g[i3], Gauss_min, Gauss_max);
	glVertex3d(p[i3][0], p[i3][1], p[i3][2]);
	glNormal3d(n[i3][0], n[i3][1], n[i3][2]);
	rainbowColor(g[i4], Gauss_min, Gauss_max);
	glVertex3d(p[i4][0], p[i4][1], p[i4][2]);
	glNormal3d(n[i4][0], n[i4][1], n[i4][2]);
	glEnd();
	break;
      case MEAN :
	glBegin(GL_QUADS);
	rainbowColor(m[i1], mean_min, mean_max);
	glVertex3d(p[i1][0], p[i1][1], p[i1][2]);
	glNormal3d(n[i1][0], n[i1][1], n[i1][2]);
	rainbowColor(m[i2], mean_min, mean_max);
	glVertex3d(p[i2][0], p[i2][1], p[i2][2]);
	glNormal3d(n[i2][0], n[i2][1], n[i2][2]);
	rainbowColor(m[i3], mean_min, mean_max);
	glVertex3d(p[i3][0], p[i3][1], p[i3][2]);
	glNormal3d(n[i3][0], n[i3][1], n[i3][2]);
	rainbowColor(m[i4], mean_min, mean_max);
	glVertex3d(p[i4][0], p[i4][1], p[i4][2]);
	glNormal3d(n[i4][0], n[i4][1], n[i4][2]);
	glEnd();
	break;
      case ISOPHOTE :
	glBegin(GL_QUADS);
	isophoteColor(p[i1], n[i1], isophote_width, eye_pos);
	glVertex3d(p[i1][0], p[i1][1], p[i1][2]);
	glNormal3d(n[i1][0], n[i1][1], n[i1][2]);
	isophoteColor(p[i2], n[i2], isophote_width, eye_pos);
	glVertex3d(p[i2][0], p[i2][1], p[i2][2]);
	glNormal3d(n[i2][0], n[i2][1], n[i2][2]);
	isophoteColor(p[i3], n[i3], isophote_width, eye_pos);
	glVertex3d(p[i3][0], p[i3][1], p[i3][2]);
	glNormal3d(n[i3][0], n[i3][1], n[i3][2]);
	isophoteColor(p[i4], n[i4], isophote_width, eye_pos);
	glVertex3d(p[i4][0], p[i4][1], p[i4][2]);
	glNormal3d(n[i4][0], n[i4][1], n[i4][2]);
	glEnd();
	break;
      case SLICING :
	glBegin(GL_QUADS);
	slicingColor(p[i1], slicing_density, eye_pos);
	glVertex3d(p[i1][0], p[i1][1], p[i1][2]);
	glNormal3d(n[i1][0], n[i1][1], n[i1][2]);
	slicingColor(p[i2], slicing_density, eye_pos);
	glVertex3d(p[i2][0], p[i2][1], p[i2][2]);
	glNormal3d(n[i2][0], n[i2][1], n[i2][2]);
	slicingColor(p[i3], slicing_density, eye_pos);
	glVertex3d(p[i3][0], p[i3][1], p[i3][2]);
	glNormal3d(n[i3][0], n[i3][1], n[i3][2]);
	slicingColor(p[i4], slicing_density, eye_pos);
	glVertex3d(p[i4][0], p[i4][1], p[i4][2]);
	glNormal3d(n[i4][0], n[i4][1], n[i4][2]);
	glEnd();
	break;
      case WIREFRAME :
	glBegin(GL_LINE_STRIP);
	glColor3d(1.0, 1.0, 1.0);
	glVertex3d(p[i1][0], p[i1][1], p[i1][2]);
	glVertex3d(p[i2][0], p[i2][1], p[i2][2]);
	glVertex3d(p[i3][0], p[i3][1], p[i3][2]);
	glVertex3d(p[i4][0], p[i4][1], p[i4][2]);
	glEnd();
	break;
      case POINTS :
	glBegin(GL_POINTS);
	glColor3d(1.0, 1.0, 1.0);
	glVertex3d(p[i1][0], p[i1][1], p[i1][2]);
	glVertex3d(p[i2][0], p[i2][1], p[i2][2]);
	glVertex3d(p[i3][0], p[i3][1], p[i3][2]);
	glVertex3d(p[i4][0], p[i4][1], p[i4][2]);
	glEnd();
	break;
      }
    }
}
