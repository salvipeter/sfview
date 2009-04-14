// SFView - Surface File Viewer
//
// Copyright (C) 2007-2008 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#include <algorithm>
#include <fstream>

#include <boost/algorithm/string.hpp>

#include "utilities.hh"

#include "nurbs-surface.hh"

// Fix parameters
int const NurbsSurface::texture_width_low = 64;
int const NurbsSurface::texture_height_low = 64;
int const NurbsSurface::isophote_map_size = 1024;
int const NurbsSurface::slicing_map_size = 16;
GLuint NurbsSurface::default_isophote_texture = 0;
GLuint NurbsSurface::slicing_texture = 0;
size_t NurbsSurface::isophote_users = 0;
size_t NurbsSurface::slicing_users = 0;

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
  Surface(), high_quality_textures(false), cn_vis(CN_NONE)
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
    if(boost::to_lower_copy(s) == ":degrees") {
      findOpenParen(in);
      in >> degree_u >> degree_v;
      isCloseParen(in);
    } else if(boost::to_lower_copy(s) == ":knot-vectors") {
      findOpenParen(in);
      findOpenParen(in);
      do {
	knots_u.push_back(readLispFloat(in));
	fknots_u.push_back(knots_u.back());
      } while(!isCloseParen(in));
      findOpenParen(in);
      do {
	knots_v.push_back(readLispFloat(in));
	fknots_v.push_back(knots_v.back());
      } while(!isCloseParen(in));
      isCloseParen(in);
    } else if(boost::to_lower_copy(s) == ":control-net") {
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
  --slicing_users;

  if(isophote_texture == default_isophote_texture)
    --isophote_users;
  else
    glDeleteTextures(1, &isophote_texture);

  if(isophote_users == 0) {
    glDeleteTextures(1, &default_isophote_texture);
    default_isophote_texture = 0;
  }
  if(slicing_users == 0) {
    slicing_texture = 0;
  }

  glDeleteTextures(1, &gauss_texture);
  glDeleteTextures(1, &mean_texture);
  glDeleteTextures(1, &slicing_texture);
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

void NurbsSurface::toggleControlNet()
{
  switch(cn_vis) {
  case CN_NONE:    cn_vis = CN_LINES;   break;
  case CN_LINES:   cn_vis = CN_FULL;    break;
  case CN_FULL:    cn_vis = CN_BORDERS; break;
  case CN_BORDERS: cn_vis = CN_NONE;    break;
  }
}

void NurbsSurface::texturePrologue(GLuint &name)
{
  glGenTextures(1, &name);
  glBindTexture(GL_TEXTURE_2D, name);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
}

void NurbsSurface::generateIsophoteTexture(GLuint &name) const
{
  texturePrologue(name);
  int const &w = isophote_map_size, &h = isophote_map_size;
  unsigned char *data = new unsigned char[w * h * 3];

  // Sphere mapping gives the 3-dimensional reflection (unit) vector
  // in 2-dimensional form: let m = 2sqrt(x^2+y^2+(1+z)^2), then
  // the coordinates are (x/m+0.5, y/m+0.5).
  // We want to get the reflection angle in the line of sight, whose
  // (doubled) cosine can be approximated by the z coordinate of the
  // reflection vector, since the reflection vector is given in eye
  // coordinates.
  // Using the equation x^2+y^2+z^2=1, it turns out that
  //   m^2 = 64 * (1/4 - (x/m)^2 - (y/m)^2), from which
  //   z = sqrt(m^2 * m^2/64) - 1 = m^2 / 8 - 1.
  for(int i = 0, index = 0; i < w; ++i) {
    double const xm = (double)i / (double)(w - 1) - 0.5;
    for(int j = 0; j < h; ++j) {
      double const ym = (double)j / (double)(h - 1) - 0.5;
      double const length2 = xm * xm + ym * ym;
      if(length2 <= 0.25) {
	double const z = 8.0 * (0.25 - length2) - 1.0;
	double const angle = (std::acos(z) * 180.0 / M_PI) / 2.0;
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

  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE,
	       data);
  delete[] data;  
}

void NurbsSurface::generateSlicingTexture()
{
  // Create a 1-dimensional map with its left half red, right half blue.
  // This will be repeatedly used with the GL_OBJECT_LINEAR map, that
  // selects a position in the map according to the point's distance
  // from a reference plane (= the eye plane in our case).
  int const &w = slicing_map_size;
  unsigned char *data = new unsigned char[w * 3];
  glGenTextures(1, &slicing_texture);
  glBindTexture(GL_TEXTURE_1D, slicing_texture);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  int index = 0;
  for(int i = 0; i < w; ++i) {
    bool color = i < w/2;
    data[index++] = color ? 255 : 0;
    data[index++] = 0;
    data[index++] = color ? 0 : 255;
  }
  glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, w, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
}

void NurbsSurface::fillRainbow(DoubleMatrix const &m, int w, int h,
			       unsigned char *output)
{
  double min = m[0][0], max = m[0][0];
  for(int i = 0; i < w; ++i)
    for(int j = 0; j < h; ++j)
      if(m[i][j] < min)
	min = m[i][j];
      else if(m[i][j] > max)
	max = m[i][j];
  double const len = max - min;

  int index = 0;
  for(int j = 0; j < h; ++j)
    for(int i = 0; i < w; ++i) {
      int const color =
	static_cast<int>(interpolate(0, (m[i][j] - min) / len, 510));
      output[index++] = static_cast<unsigned char>(std::max(color - 255, 0));
      output[index++] = static_cast<unsigned char>(255 - std::abs(color - 255));
      output[index++] = static_cast<unsigned char>(std::max(255 - color, 0));
    }
}

void NurbsSurface::generateEvaluatedTextures()
{
  int const w = high_quality_textures ? texture_width : texture_width_low;
  int const h = high_quality_textures ? texture_height : texture_height_low;
  unsigned char *data = new unsigned char[w * h * 3];
  DoubleMatrix gauss(w, DoubleVector(h));
  DoubleMatrix mean(w, DoubleVector(h));

  // Compute points, derivatives and curvatures
  double const lu = lowerBoundU(), uu = upperBoundU();
  double const lv = lowerBoundV(), uv = upperBoundV();
  for(int i = 0; i < w; ++i) {
    if(high_quality_textures && (i+1) % (w/10) == 0)
      std::cout << '.' << std::flush;
    double const u = interpolate(lu, (double)i / (double)(w - 1), uu);
    for(int j = 0; j < h; ++j) {
      double const v = interpolate(lv, (double)j / (double)(h - 1), uv);
      VectorMatrix const deriv = derivatives(u, v, 2);
      Vector const normal = (deriv[1][0] ^ deriv[0][1]).normalized();
      double const E = deriv[1][0] * deriv[1][0];
      double const F = deriv[1][0] * deriv[0][1];
      double const G = deriv[0][1] * deriv[0][1];
      double const L = normal * deriv[2][0];
      double const M = normal * deriv[1][1];
      double const N = normal * deriv[0][2];
      double const divisor = E * G - F * F;
      gauss[i][j] = divisor == 0.0 ? 0.0 : (L * N - M * M) / divisor;
      mean[i][j] = divisor == 0.0 ? 0.0 : (N*E - 2*M*F + L*G) / (2 * divisor);
    }
  }

  texturePrologue(gauss_texture);
  fillRainbow(gauss, w, h, data);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE,
	       data);

  texturePrologue(mean_texture);
  fillRainbow(mean, w, h, data);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE,
	       data);

  delete[] data;
}

void NurbsSurface::GLInit()
{
  globj = gluNewNurbsRenderer();
  gluNurbsProperty(globj, GLU_SAMPLING_TOLERANCE, 50.0);  // low quality
  gluNurbsProperty(globj, GLU_PARAMETRIC_TOLERANCE, 0.05); // high quality
  gluNurbsProperty(globj, GLU_CULLING, GL_TRUE);

  std::cout << "Generating textures... " << std::flush;

  texknots_u[0] = texknots_u[1] = lowerBoundU();
  texknots_u[2] = texknots_u[3] = upperBoundU();
  texknots_v[0] = texknots_v[1] = lowerBoundV();
  texknots_v[2] = texknots_v[3] = upperBoundV();

  if(default_isophote_texture == 0)
    generateIsophoteTexture(default_isophote_texture);
  isophote_texture = default_isophote_texture;
  ++isophote_users;

  if(slicing_texture == 0)
    generateSlicingTexture();
  ++slicing_users;

  generateEvaluatedTextures();

  std::cout << "ok" << std::endl;
}

void NurbsSurface::calculateLargeMaps()
{
  if(!high_quality_textures) {
    high_quality_textures = true;

    glDeleteTextures(1, &gauss_texture);
    glDeleteTextures(1, &mean_texture);

    std::cout << "Generating high quality textures" << std::flush;
    generateEvaluatedTextures();
    std::cout << " ok" << std::endl;
  }
}

void NurbsSurface::increaseDensity()
{
  switch(vis) {
  case SLICING: slicing_density /= 2.0; break;
  case ISOPHOTE:
    isophote_width /= 2.0;
    if(isophote_texture != default_isophote_texture)
      glDeleteTextures(1, &isophote_texture);
    else
      --isophote_users;
    generateIsophoteTexture(isophote_texture);
    break;
  default: ;
  }
}

void NurbsSurface::decreaseDensity()
{
  switch(vis) {
  case SLICING: slicing_density *= 2.0; break;
  case ISOPHOTE:
    isophote_width *= 2.0;
    if(isophote_texture != default_isophote_texture)
      glDeleteTextures(1, &isophote_texture);
    else
      --isophote_users;
    generateIsophoteTexture(isophote_texture);
    break;
  default: ;
  }
}

void NurbsSurface::display(Point const &eye_pos, Vector const &eye_dir,
			   bool high_density)
{
  if(cn_vis == CN_FULL || cn_vis == CN_LINES) {
    glDisable(GL_LIGHTING);
    for(int i = 0; i < nu; ++i)
      for(int j = 0; j < nv; ++j) {
	if(cn_vis == CN_FULL) {
	  // Points
	  glColor3d(0.8, 0.0, 0.5);
	  glPointSize(6.0);
	  glBegin(GL_POINTS);
	  glVertex3d(control_net[i * nv + j][0],
		     control_net[i * nv + j][1],
		     control_net[i * nv + j][2]);
	  glEnd();
	}
	// Lines
	glColor3d(0.0, 0.0, 0.0);
	glLineWidth(2.0);
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
  } else if(cn_vis == CN_BORDERS) {
    glDisable(GL_LIGHTING);
    glColor3d(0.0, 0.0, 0.0);
    glLineWidth(2.0);
    glBegin(GL_LINE_STRIP);
    for(int i = 0; i < nu; ++i) {
      Point &p = control_net[i * nv];
      glVertex3d(p[0], p[1], p[2]);
    }
    glEnd();
    glBegin(GL_LINE_STRIP);
    for(int i = 0; i < nu; ++i) {
      Point &p = control_net[i * nv + nv - 1];
      glVertex3d(p[0], p[1], p[2]);
    }
    glEnd();
    glBegin(GL_LINE_STRIP);
    for(int j = 0; j < nv; ++j) {
      Point &p = control_net[j];
      glVertex3d(p[0], p[1], p[2]);
    }
    glEnd();
    glBegin(GL_LINE_STRIP);
    for(int j = 0; j < nv; ++j) {
      Point &p = control_net[(nu-1) * nv + j];
      glVertex3d(p[0], p[1], p[2]);
    }
    glEnd();
  }

  if(hidden)
    return;

  if(high_density)
    gluNurbsProperty(globj, GLU_SAMPLING_METHOD, GLU_PARAMETRIC_ERROR);
  else
    gluNurbsProperty(globj, GLU_SAMPLING_METHOD, GLU_PATH_LENGTH);

  GLfloat plane[4];
  switch(vis) {
  case MEAN:
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, mean_texture);
    break;
  case GAUSS:
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, gauss_texture);
    break;
  case ISOPHOTE:
    glEnable(GL_TEXTURE_2D);
    glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
    glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
    glEnable(GL_TEXTURE_GEN_S);
    glEnable(GL_TEXTURE_GEN_T);
    glBindTexture(GL_TEXTURE_2D, isophote_texture);
    break;
  case SLICING:
    glEnable(GL_TEXTURE_1D);
    glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);
    plane[0] = eye_dir[0] / slicing_density;
    plane[1] = eye_dir[1] / slicing_density;
    plane[2] = eye_dir[2] / slicing_density;
    plane[3] = -(eye_dir * (eye_pos - Point(0.0, 0.0, 0.0)));
    glTexGenfv(GL_S, GL_OBJECT_PLANE, plane);
    glEnable(GL_TEXTURE_GEN_S);
    glBindTexture(GL_TEXTURE_1D, slicing_texture);
    break;
  default: ;
  }

  glEnable(GL_LIGHTING);
  glEnable(GL_AUTO_NORMAL);
  glEnable(GL_NORMALIZE);
  gluBeginSurface(globj);
  if(vis == MEAN || vis == GAUSS)
    gluNurbsSurface(globj, 4, &texknots_u[0], 4, &texknots_v[0], 2 * 2, 2,
		    &texcpts[0][0][0], 2, 2, GL_MAP2_TEXTURE_COORD_2);
  glColor3d(1.0, 1.0, 1.0);
  gluNurbsSurface(globj,
		  fknots_u.size(), &fknots_u[0],
		  fknots_v.size(), &fknots_v[0],
		  nv * 3, 3, &linear_cpts[0],
		  degree_u + 1, degree_v + 1,
		  GL_MAP2_VERTEX_3);
  gluEndSurface(globj);
  glDisable(GL_NORMALIZE);
  glDisable(GL_AUTO_NORMAL);

  switch(vis) {
  case MEAN:
  case GAUSS:
    glDisable(GL_TEXTURE_2D);
    break;
  case ISOPHOTE:
    glDisable(GL_TEXTURE_GEN_S);
    glDisable(GL_TEXTURE_GEN_T);
    glDisable(GL_TEXTURE_2D);
    break;
  case SLICING:
    glDisable(GL_TEXTURE_GEN_S);
    glDisable(GL_TEXTURE_1D);
  default: ;
  }
}

// Evaluation as in The NURBS Book, algorithms A2.1, A2.3, A3.6.

int NurbsSurface::findSpan(DoubleVector const &knots, double t, int n)
{
  if(t == knots[n])
    return n - 1;
  return std::upper_bound(knots.begin(), knots.end(), t) - knots.begin() - 1;
}

DoubleMatrix NurbsSurface::basisDerivatives(DoubleVector const &knots,
					    int i, int p, double u, int n)
{
  DoubleMatrix ndu(p + 1, DoubleVector(p + 1));
  DoubleMatrix a(2, DoubleVector(p + 1));
  DoubleMatrix ders(n + 1, DoubleVector(p + 1));
  DoubleVector left(p + 1), right(p + 1);

  ndu[0][0] = 1.0;
  for(int j = 1; j <= p; ++j) {
    left[j] = u - knots[i + 1 - j];
    right[j] = knots[i + j] - u;
    double saved = 0.0;
    for(int r = 0; r < j; ++r) {
      ndu[j][r] = right[r+1] + left[j-r];
      double temp = ndu[r][j-1] / ndu[j][r];

      ndu[r][j] = saved + right[r+1] * temp;
      saved = left[j-r] * temp;
    }
    ndu[j][j] = saved;
  }
  
  for(int j = 0; j <= p; ++j)
    ders[0][j] = ndu[j][p];

  for(int r = 0; r <= p; ++r) {
    int s1 = 0, s2 = 1;
    a[0][0] = 1.0;
    for(int k = 1; k <= n; ++k) {
      double d = 0.0;
      int rk = r - k, pk = p - k;
      if(r >= k) {
	a[s2][0] = a[s1][0] / ndu[pk+1][rk];
	d = a[s2][0] * ndu[rk][pk];
      }
      for(int j = (rk >= -1 ? 1 : -rk); j <= (r-1 <= pk ? k-1 : p-r); ++j) {
	a[s2][j] = (a[s1][j] - a[s1][j-1]) / ndu[pk+1][rk+j];
	d += a[s2][j] * ndu[rk+j][pk];
      }
      if(r <= pk) {
	a[s2][k] = -a[s1][k-1] / ndu[pk+1][r];
	d += a[s2][k] * ndu[r][pk];
      }
      ders[k][r] = d;
      std::swap(s1, s2);
    }
  }
  
  int r = p;
  for(int k = 1; k <= n; ++k) {
    for(int j = 0; j <= p; ++j)
      ders[k][j] *= r;
    r *= p - k;
  }

  return ders;
}

VectorMatrix NurbsSurface::derivatives(double u, double v, int d) const
{
  int const p = degree_u, q = degree_v;
  int const du = std::min(d, p), dv = std::min(d, q);
  VectorMatrix ders(d + 1, VectorVector(d + 1, Vector(0, 0, 0)));
  VectorVector temp(q + 1, Vector(0.0, 0.0, 0.0));

  int uspan = findSpan(knots_u, u, nu), vspan = findSpan(knots_v, v, nv);
  DoubleMatrix nders_u = basisDerivatives(knots_u, uspan, degree_u, u, du);
  DoubleMatrix nders_v = basisDerivatives(knots_v, vspan, degree_v, v, dv);

  for(int k = 0; k <= du; ++k) {
    for(int s = 0; s <= q; ++s)
      for(int r = 0; r <= p; ++r) {
	Vector const &pvector =
	  control_net[(uspan-p+r) * nv + (vspan-q+s)] - Point(0.0, 0.0, 0.0);
	temp[s] = temp[s] + pvector * nders_u[k][r];
      }
    int const dd = std::min(d - k, dv);
    for(int l = 0; l <= dd; ++l) {
      for(int s = 0; s <= q; ++s)
	ders[k][l] = ders[k][l] + temp[s] * nders_v[l][s];
    }
  }

  return ders;
}
