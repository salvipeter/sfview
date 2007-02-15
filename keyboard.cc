// SFView - Surface File Viewer
//
// Copyright (C) 2007 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#include <cstdlib>
#include <sstream>

#include <unistd.h>

#include "display.hh"
#include "globals.hh"
#include "mesh_surface.hh"
#include "nurbs_surface.hh"

std::string activeName(bool capital)
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

void changeVisualization(Visualization v)
{
  if(active != 0)
    surfaces[active - 1]->setVisualization(v);
  else
    for(SurfacePIterator i = surfaces.begin(); i != surfaces.end(); ++i)
      (*i)->setVisualization(v);
}

bool loadFile(std::string filename)
{
  Surface *s = new MeshSurface(filename);
  if(s->noError())
    surfaces.push_back(s);
  else
    return false;
  return true;
}

void keyboard(unsigned char key, int x, int y)
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
  case 'w' :
    changeVisualization(WIREFRAME);
    std::cout << activeName(true) << ": wireframe" << std::endl;
    display();
    break;
  case 'r' :
    std::cout << "Reloading " << activeName(false) << "..." << std::endl;
    if(active != 0) {
      if(loadFile(surfaces[active-1]->fileName())) {
	surfaces.back()->setVisualization(surfaces[active-1]->visualization());
	delete surfaces[active-1];
	surfaces[active-1] = surfaces.back();
	surfaces.pop_back();
      }
    } else {
      for(int i = 0, ie = surfaces.size(); i != ie; ++i) {
	if(loadFile(surfaces[i]->fileName())) {
	  surfaces.back()->setVisualization(surfaces[i]->visualization());
	  delete surfaces[i];
	  surfaces[i] = surfaces.back();
	  surfaces.pop_back();
	}
      }
    }
    std::cout << "Reloading done." << std::endl;
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
    for(SurfacePIterator i = surfaces.begin(); i != surfaces.end(); ++i)
      delete *i;
    std::cout << "Bye!" << std::endl;
    exit(0);
  }
}
