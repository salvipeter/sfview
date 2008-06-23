// SFView - Surface File Viewer
//
// Copyright (C) 2007-2008 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#include <fstream>

#include "nurbs-surface.hh"

void NurbsSurface::findOpenParen(std::ifstream &in)
{
  char c;
  do {
    in >> c;
  } while(c != '(');
}

bool NurbsSurface::isCloseParen(std::ifstream &in)
{
  ignoreWhitespaces(in);
  char c;
  bool result = (in.peek() == ')');
  if(result)
    in >> c;
  return result;
}

double NurbsSurface::readLispFloat(std::ifstream &in)
{
  double result;
  char c;
  in >> result;
  c = in.peek();
  if(c == 'd' || c == 'f') {
    int exponent;
    in >> c;
    in >> exponent;
    result *= std::pow(10.0, exponent);
  }
  return result;
}

void NurbsSurface::ignoreWhitespaces(std::ifstream &in)
{
  char c = in.peek();
  while(!in.eof() && (c == ' ' || c == '\t' || c == '\n' || c == '\r')) {
    in.get(c);
    c = in.peek();
  }
}

NurbsSurface::NurbsSurface(std::ifstream &in) :
  Surface(),
  show_control_net(false)
{
  std::string s;
  float x, y, z;

  nu = 0;

  findOpenParen(in);
  in >> s;
  if(s != ":bspline-surface") {
    error = true;
    return;
  }
  while(!isCloseParen(in)) {
    in >> s;
    if(s == ":degrees") {
      findOpenParen(in);
      in >> degree_u >> degree_v;
      isCloseParen(in);
    } else if(s == ":knot-vectors") {
      findOpenParen(in);
      findOpenParen(in);
      do {
	knots_u.push_back(readLispFloat(in));
      } while(!isCloseParen(in));
      findOpenParen(in);
      do {
	knots_v.push_back(readLispFloat(in));
      } while(!isCloseParen(in));
      isCloseParen(in);
    } else if(s == ":control-net") {
      findOpenParen(in);
      do {
	++nu;
	findOpenParen(in);
	do {
	  findOpenParen(in);
	  x = readLispFloat(in);
	  y = readLispFloat(in);
	  z = readLispFloat(in);
	  control_net.push_back(Point(x, y, z));
	  linear_cpts.push_back(x);
	  linear_cpts.push_back(y);
	  linear_cpts.push_back(z);
	  isCloseParen(in);
	} while(!isCloseParen(in));
      } while(!isCloseParen(in));
    } else {
      error = true;
      return;
    }
  }
  nv = control_net.size() / nu;

  std::cout << "ok (" << nu << "x" << nv << " control points)" << std::endl;
  std::cout << "Calculating bounding box... " << std::flush;
  bounding_box.first = bounding_box.second = control_net.front();
  for(PointIterator i = control_net.begin(); i != control_net.end(); ++i) {
    if((*i)[0] < bounding_box.first[0]) bounding_box.first[0] = (*i)[0];
    if((*i)[1] < bounding_box.first[1]) bounding_box.first[1] = (*i)[1];
    if((*i)[2] < bounding_box.first[2]) bounding_box.first[2] = (*i)[2];
    if((*i)[0] > bounding_box.second[0]) bounding_box.second[0] = (*i)[0];
    if((*i)[1] > bounding_box.second[1]) bounding_box.second[1] = (*i)[1];
    if((*i)[2] > bounding_box.second[2]) bounding_box.second[2] = (*i)[2];
  }
  std::cout << "ok (" << bounding_box.first << " - "
	    << bounding_box.second << ")" << std::endl;

  globj = gluNewNurbsRenderer();

  error = false;
}

bool NurbsSurface::load(std::string const &filename, SurfacePVector &sv)
{
  bool error = false;

  std::ifstream in(filename.c_str());
  std::cout << "Loading file `" << filename << "'... " << std::flush;

  while(!in.eof()) {
    NurbsSurface *sf = new NurbsSurface(in);
    error |= sf->error;
    if(!sf->error) {
      sf->filename = filename;
      sv.push_back(sf);
    }
    ignoreWhitespaces(in);
  }
  in.close();

  return !error;
}

void NurbsSurface::increaseDensity()
{
}

void NurbsSurface::decreaseDensity()
{
}

void NurbsSurface::display(Point const &eye_pos, bool)
{
  if(show_control_net) {
    glDisable(GL_LIGHTING);
    glColor3d(0.0, 0.0, 0.0);
    for(int i = 0; i < nu; ++i)
      for(int j = 0; j < nv; ++j) {
	glBegin(GL_LINE_STRIP);
	if(i != 0) {
	  glVertex3d(control_net[(i - 1) * nv + j][0],
		     control_net[(i - 1) * nv + j][1],
		     control_net[(i - 1) * nv + j][2]);
	  glVertex3d(control_net[i * nv + j][0],
		     control_net[i * nv + j][1],
		     control_net[i * nv + j][2]);
	}
	if(j != 0) {
	  glVertex3d(control_net[i * nv + j - 1][0],
		     control_net[i * nv + j - 1][1],
		     control_net[i * nv + j - 1][2]);
	  glVertex3d(control_net[i * nv + j][0],
		     control_net[i * nv + j][1],
		     control_net[i * nv + j][2]);
	}
	glEnd();
      }
  }

  if(hidden)
    return;

  glEnable(GL_LIGHTING);
  glColor3d(1.0, 1.0, 1.0);

  glEnable(GL_AUTO_NORMAL);
  glEnable(GL_NORMALIZE);
  gluBeginSurface(globj);
  gluNurbsSurface(globj,
		  knots_u.size(), &knots_u[0],
		  knots_v.size(), &knots_v[0],
		  nv * 3, 3, &linear_cpts[0],
		  degree_u + 1, degree_v + 1,
		  GL_MAP2_VERTEX_3);
  gluEndSurface(globj);
}
