// SFView - Surface File Viewer
//
// Copyright (C) 2007 Peter Salvi <vukung@yahoo.com>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#include <getopt.h>
#include <sys/types.h>
#include <unistd.h>

#include <GL/glut.h>

#include "display.hh"
#include "globals.hh"
#include "keyboard.hh"
#include "mouse.hh"

// **************
// * Parameters *
// **************

int nurbs_density = 30;
int width = 800;
int height = 600;
size_t max_n_of_quads = 2500;

// ********************
// * Global variables *
// ********************

SurfacePVector surfaces;
Point center, eye_pos;
double modelview[16], inverse[16], object_width;
bool high_density = false;

size_t active = 0;
int mouse_start[2], next_id = -1;
MouseMode mouse_mode = NOTHING;

// *************
// * Constants *
// *************

std::string const help_string =
  "Usage: sfview [OPTION]... FILE [FILE]...\n"
  "Loads the file(s) designated in the command line, displaying their\n"
  "contents in one common window. It recognizes the PTS an RBN file formats."
  "\n\n"
  "Mandatory arguments to long options are mandatory for short options too.\n"
  "  -d, --density=NUMBER       Sampling density of RBN files. It represents\n"
  "                               the average number of sampled points\n"
  "                               between two control points. 30 by default.\n"
  "  -g, --geometry=WIDTHxHEIGHT\n"
  "                             Geometry of the window at startup.\n"
  "                               800x600 by default.\n"
  "  -h, --help                 Display this help and exit.\n"
  "  -q, --quads=NUMBER         Maximum number of quads in low-density mode.\n"
  "                               2500 by default.\n"
  "  -v, --version              Print version information and exit."
  "\n\n"
  "See the man page for more information."
  "\n\n"
  "Report bugs to <vukung@yahoo.com>.";

std::string const version_string =
  "SFView 0.22\n"
  "Copyright (C) 2007 Peter Salvi <vukung@yahoo.com>\n"
  "sfview comes with NO WARRANTY, to the extent permitted by law.\n"
  "This is free software.  You may redistribute copies of it under the terms\n"
  "of the GNU General Public License <http://www.gnu.org/licenses/gpl.html>.\n"
  "For more information about these matters, see the files named COPYING.";

// ****************
// * Command line *
// ****************

StringVector parseCommandLine(int argc, char *argv[])
{
  StringVector sv;
  
  static struct option long_options[] = {
    {"density", required_argument, NULL, 'd'},
    {"geometry", required_argument, NULL, 'g'},
    {"help", no_argument, NULL, 'h'},
    {"quads", required_argument, NULL, 'q'},
    {"version", no_argument, NULL, 'v'},
    {0, 0, 0, 0}
  };

  std::string geom;
  int c, index;

  while((c = getopt_long(argc, argv, "d:g:hq:v", long_options, NULL)) != -1) {
    switch(c) {
    case 'd' :
      nurbs_density = std::atoi(optarg);
      if(nurbs_density != 0)
	std::cout << "Density: " << nurbs_density << std::endl;
      else {
	std::cerr << argv[0] << ": " << optarg
		  << " is not a valid number." << std::endl;
	exit(2);
      }
      break;
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
    case 'v' :
      std::cout << version_string << std::endl;
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

// ********
// * Main *
// ********

int main(int argc, char *argv[])
{
  StringVector files = parseCommandLine(argc, argv);
  if(files.empty()) {
    std::cerr << argv[0] << ": no files" << std::endl;
    std::cerr << "Try " << argv[0]
	      << " --help for more information." << std::endl;
    return 1;
  }

  for(StringVector::const_iterator i = files.begin(); i != files.end(); ++i)
    if(!loadFile(*i)) {
      std::cerr << "Cannot load file " << *i << std::endl;
      return 2;
    }

  glutInit(&argc, argv);
  glutInitWindowSize(width, height);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
  glutCreateWindow("SFView");

  init();

  glutDisplayFunc(display);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouseButton);
  glutMotionFunc(mouseMotion);
  glutReshapeFunc(reshape);

  glutMainLoop();

  return 0;
}
