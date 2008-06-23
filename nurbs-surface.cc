// SFView - Surface File Viewer
//
// Copyright (C) 2007-2008 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#include <fstream>

#include "nurbs-surface.hh"

void NurbsSurface::findOpenParen(std::ifstream &in) const
{
  char c;
  do {
    in >> c;
  } while(c != '(');
}

bool NurbsSurface::isCloseParen(std::ifstream &in) const
{
  char c;
  bool result = (in.peek() == ')');
  if(result)
    in >> c;
  return result;
}

NurbsSurface::NurbsSurface(std::string filename) :
  Surface(filename),
  show_control_net(false)
{
  std::ifstream in(filename.c_str());

  std::cout << "Loading file `" << filename << "'... " << std::flush;

  std::string s;
  float d, x, y, z;

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
	in >> d;
	knots_u.push_back(d);
      } while(!isCloseParen(in));
      findOpenParen(in);
      do {
	in >> d;
	knots_v.push_back(d);
      } while(!isCloseParen(in));
      isCloseParen(in);
    } else if(s == ":control-net") {
      findOpenParen(in);
      do {
	++nu;
	findOpenParen(in);
	do {
	  findOpenParen(in);
	  in >> x >> y >> z;
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
  in.close();
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

  double const axis = (bounding_box.second - bounding_box.first).length();
  globj = gluNewNurbsRenderer();
  gluNurbsProperty(globj, GLU_SAMPLING_TOLERANCE, axis / 20.0);

  error = false;
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

  glEnable(GL_AUTO_NORMAL);
  glEnable(GL_LIGHTING);
  glColor3d(1.0, 1.0, 1.0);

  gluBeginSurface(globj);
  gluNurbsSurface(globj,
		  knots_u.size(), &knots_u[0],
		  knots_v.size(), &knots_v[0],
		  nv * 3, 3, &linear_cpts[0],
		  degree_u + 1, degree_v + 1,
		  GL_MAP2_VERTEX_3);
  gluEndSurface(globj);
}
