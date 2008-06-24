// SFView - Surface File Viewer
//
// Copyright (C) 2007-2008 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#include <fstream>

#include "nurbs-surface.hh"

GLuint NurbsSurface::default_isophote_texture = 0;

GLfloat NurbsSurface::texcpts[2][2][2] =
  { {{0.0, 0.0}, {0.0, 1.0}}, {{1.0, 0.0}, {1.0, 1.0}} };

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
  Surface(), show_control_net(false)
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

  error = false;
}

NurbsSurface::~NurbsSurface()
{
//   glDeleteTextures(1, &mean_texture);
//   glDeleteTextures(1, &gauss_texture);
//   glDeleteTextures(1, &default_isophote_texture);
  glDeleteTextures(1, &isophote_texture);
//   glDeleteTextures(1, &slicing_texture);
}

bool NurbsSurface::load(std::string const &filename, SurfacePVector &sv,
			int texwidth, int texheight)
{
  bool error = false;

  std::ifstream in(filename.c_str());
  std::cout << "Loading file `" << filename << "'... " << std::flush;

  while(!in.eof()) {
    NurbsSurface *sf = new NurbsSurface(in);
    error |= sf->error;
    if(!sf->error) {
      sf->filename = filename;
      sf->texture_width = texwidth;
      sf->texture_height = texheight;
      sv.push_back(sf);
    }
    ignoreWhitespaces(in);
  }
  in.close();

  return !error;
}

void NurbsSurface::GLInit()
{
  globj = gluNewNurbsRenderer();

  std::cout << "Generating textures... " << std::flush;

  texknots_u[0] = texknots_u[1] = knots_u[degree_u];
  texknots_u[2] = texknots_u[3] = knots_u[knots_u.size() - degree_u - 1];
  texknots_v[0] = texknots_v[1] = knots_v[degree_v];
  texknots_v[2] = texknots_v[3] = knots_v[knots_v.size() - degree_v - 1];

//   void generateTexture(GLuint &name, void (NurbsSurface::*)(void) const fn);

//   if(default_isophote_texture == 0)
//     generateTexture(default_isophote_texture, generateIsophoteTexture);
//   isophote_texture = default_isophote_texture;

  // update:
//   if(isophote_texture != default_isophote_texture)
//     glDeleteTextures(1, &isophote_texture);
//   generateTexture(isophote_texture, generateIsophoteTexture);

  glGenTextures(1, &isophote_texture);
  glBindTexture(GL_TEXTURE_2D, isophote_texture);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  unsigned char *data = new unsigned char[texture_width * texture_height * 3];
  double const isophote_width = 5.0;
  // Sphere mapping gives the 3-dimensional reflection (unit) vector
  // in 2-dimensional form: let m = 2sqrt(x^2+y^2+(1+z)^2), then
  // the coordinates are (x/m+0.5, y/m+0.5).
  // We want to get the reflection angle in the line of sight, whose
  // cosine is the z coordinate of the reflection vector, since the
  // reflection vector is originally given in eye coordinates.
  // Using the equation x^2+y^2+z^2=1, it turns out that
  //   m = 8 * sqrt(1/4 - (x/m)^2 - (y/m)^2).
  for(int i = 0, index = 0; i < texture_width; ++i) {
    double const xm = (double)i / (double)(texture_width - 1) - 0.5;
    for(int j = 0; j < texture_height; ++j) {
      double const ym = (double)j / (double)(texture_height - 1) - 0.5;
      double const length2 = xm * xm + ym * ym;
      if(length2 <= 0.25) {
	double const m = 8.0 * std::sqrt(0.25 - length2);
	double const x = xm * m;
	double const y = ym * m;
	double const z = std::sqrt(1.0 - x * x - y * y);
	double const angle = std::acos(z) * 180.0 / M_PI;
	bool color =
	  static_cast<int>(std::floor(angle / isophote_width)) % 2 == 0;
	data[index++] = 255;
	data[index++] = color ? 0 : 255;
	data[index++] = color ? 0 : 255;
      } else { // outside the interesting region, just fill with zeroes
	data[index++] = 0;
	data[index++] = 0;
	data[index++] = 0;
      }
    }
  }
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texture_width, texture_height,
	       0, GL_RGB, GL_UNSIGNED_BYTE, data);
  delete[] data;

  std::cout << "ok" << std::endl;
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

  if(vis != SHADED)
    glEnable(GL_TEXTURE_2D);

  switch(vis) {
//   case MEAN:     glBindTexture(GL_TEXTURE_2D, mean_texture); break;
//   case GAUSS:    glBindTexture(GL_TEXTURE_2D, gauss_texture); break;
  case ISOPHOTE:
    glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
    glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
    glEnable(GL_TEXTURE_GEN_S);
    glEnable(GL_TEXTURE_GEN_T);	
    glBindTexture(GL_TEXTURE_2D, isophote_texture);
    break;
//   case SLICING:  glBindTexture(GL_TEXTURE_2D, slicing_texture); break;
  default: ;
  }

  glEnable(GL_LIGHTING);
  glEnable(GL_AUTO_NORMAL);
  glEnable(GL_NORMALIZE);
  gluBeginSurface(globj);
  if(vis != SHADED)
    gluNurbsSurface(globj, 4, &texknots_u[0], 4, &texknots_v[0], 2 * 2, 2,
		    &texcpts[0][0][0], 2, 2, GL_MAP2_TEXTURE_COORD_2);
  glColor3d(1.0, 1.0, 1.0);
  gluNurbsSurface(globj,
		  knots_u.size(), &knots_u[0],
		  knots_v.size(), &knots_v[0],
		  nv * 3, 3, &linear_cpts[0],
		  degree_u + 1, degree_v + 1,
		  GL_MAP2_VERTEX_3);
  gluEndSurface(globj);
  glDisable(GL_NORMALIZE);
  glDisable(GL_AUTO_NORMAL);

  if(vis == ISOPHOTE) {
    glDisable(GL_TEXTURE_GEN_S);
    glDisable(GL_TEXTURE_GEN_T);
  }

  if(vis != SHADED)
    glDisable(GL_TEXTURE_2D);
}
