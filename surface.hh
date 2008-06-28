// SFView - Surface File Viewer
//
// Copyright (C) 2007-2008 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#ifndef SURFACE_HH
#define SURFACE_HH

#include "common.hh"

enum Visualization { SHADED, GAUSS, MEAN, ISOPHOTE,
		     SLICING, WIREFRAME, POINTS };

class Surface {
public:
  Surface() : vis(SHADED), isophote_width(5.0), hidden(false) { }
  virtual ~Surface() { }
  virtual void GLInit() { }

  virtual bool showControlNet() const = 0;
  virtual void toggleShowControlNet() = 0;
  std::string fileName() const { return filename; }
  Box boundingBox() const { return bounding_box; }
  Visualization visualization() const { return vis; }
  virtual void setVisualization(Visualization const v) = 0;
  void setSlicingDensity(double d) { slicing_density = d; }
  virtual void calculateLargeMaps() { }
  bool isHidden() const { return hidden; }
  void toggleHidden() { hidden = !hidden; }
  bool noError() const { return !error; }

  virtual void increaseDensity() = 0;
  virtual void decreaseDensity() = 0;
  virtual void display(Point const &eye_pos, Vector const &eye_dir,
		       bool high_density) = 0;
protected:
  std::string filename;
  Box bounding_box;
  Visualization vis;
  double isophote_width, slicing_density;
  bool hidden, error;
};

typedef std::vector<Surface *> SurfacePVector;
typedef SurfacePVector::iterator SurfacePIterator;

#endif // SURFACE_HH
