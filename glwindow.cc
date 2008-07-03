// SFView - Surface File Viewer
//
// Copyright (C) 2007-2008 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#include <getopt.h>
#include <sys/types.h>
#include <unistd.h>

#include <algorithm>
#include <fstream>
#include <sstream>

#include <boost/lambda/bind.hpp>
#include <boost/lambda/construct.hpp>
#include <boost/lambda/lambda.hpp>

#include <GL/glut.h>

#include "sfview.hh"
#include "mesh-surface.hh"
#include "nurbs-surface.hh"
#include "utilities.hh"

using namespace boost::lambda;
using std::for_each;

std::string const GLWindow::help_string =
  "Usage: sfview [OPTION]... FILE [FILE]...\n"
  "Loads the file(s) designated in the command line, displaying their\n"
  "contents in one common window. It recognizes the PTS an RBN file formats."
  "\n\n"
  "Mandatory arguments to long options are mandatory for short options too.\n"
  "  -g, --geometry=WIDTHxHEIGHT\n"
  "                             Geometry of the window at startup.\n"
  "                               800x600 by default.\n"
  "  -h, --help                 Display this help and exit.\n"
  "  -q, --quads=NUMBER         Maximum number of quads in low-density mode.\n"
  "                               2500 by default.\n"
  "  -t, --texture-size=WIDTHxHEIGHT\n"
  "                             Texture size of the curvature maps for\n"
  "                               B-Spline surfaces. 256x256 by default.\n"
  "  -v, --version              Print version information and exit."
  "\n\n"
  "See the man page for more information."
  "\n\n"
  "Report bugs to <vukung@yahoo.com>.";

std::string const GLWindow::copyright_string =
  "Copyright (C) 2007-2008 Peter Salvi <vukung@yahoo.com>\n"
  "sfview comes with NO WARRANTY, to the extent permitted by law.\n"
  "This is free software.  You may redistribute copies of it under the terms\n"
  "of the GNU General Public License <http://www.gnu.org/licenses/gpl.html>.\n"
  "For more information about these matters, see the files named COPYING.";

double const GLWindow::view_angle = 45.0; // field of view in the y direction
double const GLWindow::znear_coefficient = 0.4;
double const GLWindow::zfar_coefficient = 1.7;

GLWindow::GLWindow() :
  width(800), height(600), texture_width(256), texture_height(256),
  max_n_of_quads(2500), high_density(false),
  active(0), next_id(-1), mouse_mode(NOTHING)
{
}

void GLWindow::init(int argc, char *argv[])
{
  StringVector files = parseCommandLine(argc, argv);
  if(files.empty()) {
    std::cerr << argv[0] << ": no files" << std::endl;
    std::cerr << "Try " << argv[0]
	      << " --help for more information." << std::endl;
    exit(1);
  }

  for(StringVector::const_iterator i = files.begin(); i != files.end(); ++i)
    if(!loadFile(*i)) {
      std::cerr << "Cannot load file " << *i << std::endl;
      exit(2);
    }
}

void GLWindow::show()
{
  int argc = 1;			// dummy arguments
  char *argv[] = { "sfview" };	//  for glutInit
  glutInit(&argc, argv);
  glutInitWindowSize(width, height);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
  glutCreateWindow("SFView");

  object_width = (double)width / (double)height;

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);

  GLfloat ambientLight[] = { 0.2, 0.2, 0.2, 1.0 };
  GLfloat diffuseLight[] = { 0.7, 0.7, 0.7, 1.0 };
  GLfloat specularLight[] = { 0.7, 0.7, 0.7, 1.0 };
  GLfloat position[] = { -4.5, 4.0, 4.0, 1.0 };

  glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
  glLightfv(GL_LIGHT0, GL_POSITION, position);

  glDepthFunc(GL_LEQUAL);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_COLOR_MATERIAL);
  glShadeModel(GL_SMOOTH);

//   GLfloat ambientMaterialFront[] = { ... };
//   GLfloat diffuseMaterialFront[] = { ... };
//   GLfloat specularMaterialFront[] = { ... };
//   GLfloat shininessMaterialFront = ...;

//   glMaterialfv(GL_FRONT, GL_AMBIENT, ambientMaterialFront);
//   glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuseMaterialFront);
//   glMaterialfv(GL_FRONT, GL_SPECULAR, specularMaterialFront);
//   glMaterialfv(GL_FRONT, GL_SHININESS, &shininessMaterialFront);

//   GLfloat ambientMaterialBack[] = { ... };
//   GLfloat diffuseMaterialBack[] = { ... };
//   GLfloat specularMaterialBack[] = { ... };
//   GLfloat shininessMaterialBack = ...;

//   glMaterialfv(GL_BACK, GL_AMBIENT, ambientMaterialBack);
//   glMaterialfv(GL_BACK, GL_DIFFUSE, diffuseMaterialBack);
//   glMaterialfv(GL_BACK, GL_SPECULAR, specularMaterialBack);
//   glMaterialfv(GL_BACK, GL_SHININESS, &shininessMaterialBack);

  glClearColor(0.7, 0.7, 0.7, 1.0);

  Box b = surfaces.front()->boundingBox();
  for(SurfacePIterator i = surfaces.begin(); i != surfaces.end(); ++i)
    b = boxUnion(b, (*i)->boundingBox());
  for(SurfacePIterator i = surfaces.begin(); i != surfaces.end(); ++i)
    (*i)->setSlicingDensity((b.second - b.first).length() / 50);
  bounding_box = b;
  zoomToBoundingBox();

  glutDisplayFunc(::display);
  glutKeyboardFunc(::keyboard);
  glutMouseFunc(::mouseButton);
  glutMotionFunc(::mouseMotion);
  glutReshapeFunc(::reshape);

  for_each(surfaces.begin(), surfaces.end(), bind(&Surface::GLInit, _1));

  glutMainLoop();
}

void GLWindow::display()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  gluLookAt(eye[0], eye[1], eye[2],
	    center[0], center[1], center[2],
	    up[0], up[1], up[2]);

  Vector const eye_dir = (center - eye).normalized();
  for(SurfacePIterator i = surfaces.begin(); i != surfaces.end(); ++i)
    (*i)->display(eye, eye_dir, mouse_mode == NOTHING && high_density);

  glFlush();
  glutSwapBuffers();
}

bool GLWindow::reloadActive()
{
  SurfacePVector old = surfaces;
  surfaces.clear();
  bool ok = true;
  if(active != 0) {
    std::string fname = surfaces[active-1]->fileName();
    std::cout << "Reloading " << fname << "..." << std::endl;
    if(loadFile(fname)) {
      for_each(surfaces.begin(), surfaces.end(), bind(&Surface::GLInit, _1));
      for(SurfacePIterator i = old.begin(); i != old.end();) {
	if((*i)->fileName() == fname) {
	  delete *i;
	  i = old.erase(i); 
	} else
	  ++i;
      }
      surfaces.insert(surfaces.begin(), old.begin(), old.end());
      active = old.size() + 1;
    } else
      ok = false;
  } else {
    std::cout << "Reloading all surfaces..." << std::endl;
    for(SurfacePIterator i = old.begin(); ok && i != old.end(); ++i) {
      if(!loadFile((*i)->fileName()))
	ok = false;
    }
    if(ok) {
      for_each(old.begin(), old.end(), bind(delete_ptr(), _1));
      for_each(surfaces.begin(), surfaces.end(), bind(&Surface::GLInit, _1));
    }
  }
  if(!ok) {
    for_each(surfaces.begin(), surfaces.end(), bind(delete_ptr(), _1));
    surfaces = old;
    return false;
  }
  return true;
}

void GLWindow::keyboard(unsigned char key, int x, int y)
{
  std::string filename;
  std::stringstream str;

  switch(tolower(key)) {
  case '0' : // T
  case '1' : // |
  case '2' : // |
  case '3' : // |
  case '4' : // |
  case '5' : // |
  case '6' : // |
  case '7' : // |
  case '8' : // |
  case '9' : // V
    active = key - '0';
    if(active > surfaces.size())
      active = surfaces.size();
    std::cout << "Changing surface to: " << activeName(false) << std::endl;
    if(active != 0) {
      surfaces[active-1]->toggleHidden();
      display();
      usleep(10000);
      surfaces[active-1]->toggleHidden();
      display();
    }
    break;
  case 'c' :
    if(active != 0)
      surfaces[active - 1]->toggleShowControlNet();
    else
      for(SurfacePIterator i = surfaces.begin(); i != surfaces.end(); ++i)
	(*i)->toggleShowControlNet();
    std::cout << "Toggling control net for " << activeName(false) << std::endl;
    display();
    break;
  case 'g' :
    changeVisualization(GAUSS);
    std::cout << activeName(true) << ": Gauss map" << std::endl;
    display();
    break;
  case 'h' :
    if(active != 0)
      surfaces[active - 1]->toggleHidden();
    else
      for(SurfacePIterator i = surfaces.begin(); i != surfaces.end(); ++i)
	(*i)->toggleHidden();
    std::cout << "Toggling visibility for " << activeName(false) << std::endl;
    display();
    break;
  case 'i' :
    changeVisualization(ISOPHOTE);
    std::cout << activeName(true) << ": isophote map" << std::endl;
    display();
    break;
  case 'l' :
    changeVisualization(SLICING);
    std::cout << activeName(true) << ": slicing map" << std::endl;
    display();
    break;
  case 'm' :
    changeVisualization(MEAN);
    std::cout << activeName(true) << ": mean map" << std::endl;
    display();
    break;
  case 'n' :
    changeVisualization(SHADED);
    std::cout << activeName(true) << ": no map" << std::endl;
    display();
    break;
  case 'p' :
    changeVisualization(POINTS);
    std::cout << activeName(true) << ": points" << std::endl;
    display();
    break;
  case 'r' :
    if(reloadActive()) {
      std::cout << "Reloading done." << std::endl;  
      display();
    } else
      std::cout << "Reloading failed." << std::endl;
    break;
  case 'w' :
    changeVisualization(WIREFRAME);
    std::cout << activeName(true) << ": wireframe" << std::endl;
    display();
    break;
  case 'y' :
    if(active != 0)
      surfaces[active - 1]->calculateLargeMaps();
    else
      for(SurfacePIterator i = surfaces.begin(); i != surfaces.end(); ++i)
	(*i)->calculateLargeMaps();
    display();
    break;
  case '+' :
    if(active != 0)
      surfaces[active - 1]->increaseDensity();
    else
      for(SurfacePIterator i = surfaces.begin(); i != surfaces.end(); ++i)
	(*i)->increaseDensity();
    std::cout << "Increasing density for " << activeName(false) << std::endl;
    display();
    break;
  case '-' :
    if(active != 0)
      surfaces[active - 1]->decreaseDensity();
    else
      for(SurfacePIterator i = surfaces.begin(); i != surfaces.end(); ++i)
	(*i)->decreaseDensity();
    std::cout << "Decreasing density for " << activeName(false) << std::endl;
    display();
    break;
  case 'd' :
    high_density = !high_density;
    std::cout << "High density: ";
    if(high_density)
      std::cout << "on";
    else
      std::cout << "off";
    std::cout << std::endl;
    display();
    break;
  case 's' :
    str << "sfview-" << getpid() << next_id-- << ".ppm";
    filename = str.str();
    std::cout << "Saving screenshot in " << filename << "..." << std::flush;
    saveScreenShot(filename);
    std::cout << "ok" << std::endl;
    break;
  case 'q' :
    for_each(surfaces.begin(), surfaces.end(), bind(delete_ptr(), _1));
    std::cout << "Bye!" << std::endl;
    exit(0);
  }
}

void GLWindow::mouseButton(int button, int state, int x, int y)
{
  if(state == GLUT_UP) {
    mouse_mode = NOTHING;
    display();
  }
  else if(state == GLUT_DOWN) {
    mouse_start[0] = x; mouse_start[1] = y;
    switch(button) {
    case GLUT_LEFT_BUTTON : mouse_mode = ROTATION; break;
    case GLUT_RIGHT_BUTTON : mouse_mode = ZOOM; break;
    case GLUT_MIDDLE_BUTTON : mouse_mode = PAN; break;
    }
  }
}

void GLWindow::mouseMotion(int x, int y)
{
  if(mouse_mode == NOTHING)
    return;

  double dx = (double)(x - mouse_start[0]) / (double)width;
  double dy = (double)(y - mouse_start[1]) / (double)height;
  double length;
  Vector axis;

  switch(mouse_mode) {
  case ROTATION :
    axis = (up ^ (center - eye)).normalized();
    eye = rotatePoint(eye, center, up, -dx * M_PI);
    eye = rotatePoint(eye, center, axis, dy * M_PI);
    up = (rotatePoint((center + up), center, axis, dy * M_PI) - center)
      .normalized();
    break;
  case ZOOM :
    eye = center + (eye - center) * (dy + 1.0);
    break;
  case PAN :
    length = (eye - center).length() * tan(view_angle * M_PI / 360.0) * 2.0;
    axis = (up ^ (center - eye)).normalized();
    center = center + axis * dx * length * object_width;
    center = center + up * dy * length;
    break;
  default: break;
  }
  mouse_start[0] = x; mouse_start[1] = y;

  setClippingPlanes();

  display();
}

void GLWindow::reshape(int w, int h)
{
  if(height == 0)
    height = 1;

  width = w;
  height = h;
  object_width = (double)w / (double)h;
  
  glViewport(0, 0, w, h);

  setClippingPlanes();
}

StringVector GLWindow::parseCommandLine(int argc, char *argv[])
{
  StringVector sv;
  
  static struct option long_options[] = {
    {"geometry", required_argument, NULL, 'g'},
    {"help", no_argument, NULL, 'h'},
    {"quads", required_argument, NULL, 'q'},
    {"texture-size", required_argument, NULL, 't'},
    {"version", no_argument, NULL, 'v'},
    {0, 0, 0, 0}
  };

  std::string geom;
  int c, index;

  while((c = getopt_long(argc, argv, "g:hq:t:v", long_options, NULL)) != -1) {
    switch(c) {
    case 'g' :
      geom = optarg;
      index = geom.find('x');
      width = std::atoi(geom.substr(0, index).c_str());
      height = std::atoi(geom.substr(index + 1).c_str());
      if(width <= 0 || height <= 0) {
	std::cerr << argv[0] << ": geometry should be of "
	  "the format WIDTHxHEIGHT." << std::endl;
	exit(2);
      } else
	std::cout << "Geometry: " << width << "x" << height << std::endl;
      break;
    case 'h' :
      std::cout << help_string << std::endl;
      exit(0);
    case 'q' :
      max_n_of_quads = std::atoi(optarg);
      if(max_n_of_quads != 0)
	std::cout << "# of quads in low-density-mode: "
		  << max_n_of_quads << std::endl;
      else {
	std::cerr << argv[0] << ": " << optarg
		  << " is not a valid number." << std::endl;
	exit(2);
      }
      break;
    case 't' :
      geom = optarg;
      index = geom.find('x');
      texture_width = std::atoi(geom.substr(0, index).c_str());
      texture_height = std::atoi(geom.substr(index + 1).c_str());
      if(texture_width <= 0 || texture_height <= 0) {
	std::cerr << argv[0] << ": texture size should be of "
	  "the format WIDTHxHEIGHT." << std::endl;
	exit(2);
      } else
	std::cout << "Texture size: " << width << "x" << height << std::endl;
      break;
    case 'v' :
      std::cout << "SFView " << version << std::endl;
      std::cout << copyright_string << std::endl;
      exit(0);
    case '?' :
      std::cerr << "Try " << argv[0]
		<< " --help for more information." << std::endl;
      exit(1);
    }
  }

  while(optind < argc)
    sv.push_back(argv[optind++]);

  return sv;
}

bool GLWindow::loadFile(std::string filename)
{
  std::string str;
  bool result;

  std::ifstream in(filename.c_str());
  if(in.fail())
    return false;
  in >> str;
  in.close();

  // Simplified test: if it begins with a number, it is PTS, otherwise RBN.
  if(std::atoi(str.c_str()) > 0)
    result = MeshSurface::load(filename, surfaces, max_n_of_quads);
  else
    result = NurbsSurface::load(filename, surfaces,
				texture_width, texture_height);

  return result;
}

void GLWindow::zoomToBoundingBox()
{
  center = affineCombine(bounding_box.first, 0.5, bounding_box.second);
  up = Vector(0.0, 1.0, 0.0);

  // Eye position is calculated from the bounding box and the view angle.
  // The eye is located behind the center of the bounding box.
  // Also give some extra space on the sides (the factor is now 1.2).
  double const bwidth = bounding_box.second[0] - bounding_box.first[0];
  double const bheight = bounding_box.second[1] - bounding_box.first[1];
  double const length = 1.2 * std::max(bwidth, bheight);
  double const eye_z =
    bounding_box.second[2] + length / (2.0 * tan(view_angle * M_PI / 360.0));
  eye = Point(center[0], center[1], eye_z);

  setClippingPlanes();
}

Point GLWindow::boundingBoxPoint(int i)
{
  Point result;
  Point const &b1 = bounding_box.first, &b2 = bounding_box.second;
  switch(i) {
  case 0: result = b1; break;
  case 1: result = Point(b2[0], b1[1], b1[2]); break;
  case 2: result = Point(b1[0], b2[1], b1[2]); break;
  case 3: result = Point(b2[0], b2[1], b1[2]); break;
  case 4: result = Point(b1[0], b1[1], b2[2]); break;
  case 5: result = Point(b2[0], b1[1], b2[2]); break;
  case 6: result = Point(b1[0], b2[1], b2[2]); break;
  case 7: result = b2; break;
  default: ;
  }
  return result;
}

void GLWindow::setClippingPlanes()
{
  double znear = 1, zfar = 10;

  // Get the eye distance of the nearest and farest bounding box vertex
  Vector const eye_dir = (center - eye).normalized();
  for(int i = 0; i < 8; ++i) {
    double const d = (boundingBoxPoint(i) - eye) * eye_dir;
    znear = std::min(znear, d);
    zfar = std::max(zfar, d);
  }
  znear *= znear_coefficient;
  zfar *= zfar_coefficient;

  if(znear < 0.0)		// It should be positive, but we have no
    znear = 0.0;		// better alternative...

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(view_angle, object_width, znear, zfar);
  glMatrixMode(GL_MODELVIEW);
}

std::string GLWindow::activeName(bool capital)
{
  std::stringstream str;

  if(active == 0)
    str << "all surfaces";
  else
    str << "surface " << active;

  std::string result(str.str());

  if(capital)
    result[0] = toupper(result[0]);
  return result;
}

void GLWindow::changeVisualization(Visualization v)
{
  if(active != 0)
    surfaces[active - 1]->setVisualization(v);
  else
    for(SurfacePIterator i = surfaces.begin(); i != surfaces.end(); ++i)
      (*i)->setVisualization(v);
}

void GLWindow::saveScreenShot(std::string filename)
{
  std::vector<char> buffer;
  buffer.resize(width * height * 3);

  glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, &buffer[0]);

  std::ofstream f(filename.c_str());
  f << "P6" << std::endl;
  f << width << " " << height << " 255" << std::endl;
  for(int j = height - 1; j >= 0; --j)
    for(int i = 0; i < width; ++i)
      f << buffer[3 * (j * width + i)]
	<< buffer[3 * (j * width + i) + 1]
	<< buffer[3 * (j * width + i) + 2];
  f << std::endl;
  f.close();
}
