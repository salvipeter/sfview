// SFView - Surface File Viewer
//
// Copyright (C) 2007-2008 Peter Salvi <vukung@yahoo.com>
//
// See the file `sfview.cc' for copyright details.

#include "utilities.hh"

inline Point minPoint(Point const &p1, Point const &p2)
{
  return Point(p1[0] < p2[0] ? p1[0] : p2[0],
	       p1[1] < p2[1] ? p1[1] : p2[1],
	       p1[2] < p2[2] ? p1[2] : p2[2]);
}

inline Point maxPoint(Point const &p1, Point const &p2)
{
  return Point(p1[0] > p2[0] ? p1[0] : p2[0],
	       p1[1] > p2[1] ? p1[1] : p2[1],
	       p1[2] > p2[2] ? p1[2] : p2[2]);
}

Box boxUnion(Box const &b1, Box const &b2)
{
  Box result;
  result.first = minPoint(b1.first, b2.first);
  result.second = maxPoint(b1.second, b2.second);
  return result;
}

inline double interpolate(double d1, double t, double d2)
{
  return d1 + (d2 - d1) * t;
}

Point affineCombine(Point const &p1, double t, Point const &p2)
{
  return Point(interpolate(p1[0], t, p2[0]),
	       interpolate(p1[1], t, p2[1]),
	       interpolate(p1[2], t, p2[2]));
}

void invert4x4Matrix(double *m, double *result)
{
  // Taken from GLT ZPR (http://www.nigels.com/glt/gltzpr/)
  
  double d12, d13, d23, d24, d34, d41;

  d12 = (m[2] * m[7] - m[3] * m[6]);
  d13 = (m[2] * m[11] - m[3] * m[10]);
  d23 = (m[6] * m[11] - m[7] * m[10]);
  d24 = (m[6] * m[15] - m[7] * m[14]);
  d34 = (m[10] * m[15] - m[11] * m[14]);
  d41 = (m[14] * m[3] - m[15] * m[2]);

  result[0] =  (m[5] * d34 - m[9] * d24 + m[13] * d23);
  result[1] = -(m[1] * d34 + m[9] * d41 + m[13] * d13);
  result[2] =  (m[1] * d24 + m[5] * d41 + m[13] * d12);
  result[3] = -(m[1] * d23 - m[5] * d13 + m[9] * d12);

  double const determinant = m[0] * result[0] + m[4] * result[1] +
    m[8] * result[2] + m[12] * result[3];

  if(determinant == 0.0)
    std::cerr << "Inversion of a singular matrix!" << std::endl;

  double const det_inv = 1.0 / determinant;

  result[0] *= det_inv;
  result[1] *= det_inv;
  result[2] *= det_inv;
  result[3] *= det_inv;

  result[4] = -(m[4] * d34 - m[8] * d24 + m[12] * d23) * det_inv;
  result[5] =  (m[0] * d34 + m[8] * d41 + m[12] * d13) * det_inv;
  result[6] = -(m[0] * d24 + m[4] * d41 + m[12] * d12) * det_inv;
  result[7] =  (m[0] * d23 - m[4] * d13 + m[8] * d12) * det_inv;

  d12 = m[0] * m[5] - m[1] * m[4];
  d13 = m[0] * m[9] - m[1] * m[8];
  d23 = m[4] * m[9] - m[5] * m[8];
  d24 = m[4] * m[13] - m[5] * m[12];
  d34 = m[8] * m[13] - m[9] * m[12];
  d41 = m[12] * m[1] - m[13] * m[0];

  result[8] =  (m[7] * d34 - m[11] * d24 + m[15] * d23) * det_inv;
  result[9] = -(m[3] * d34 + m[11] * d41 + m[15] * d13) * det_inv;
  result[10] =  (m[3] * d24 + m[7] * d41 + m[15] * d12) * det_inv;
  result[11] = -(m[3] * d23 - m[7] * d13 + m[11] * d12) * det_inv;
  result[12] = -(m[6] * d34 - m[10] * d24 + m[14] * d23) * det_inv;
  result[13] =  (m[2] * d34 + m[10] * d41 + m[14] * d13) * det_inv;
  result[14] = -(m[2] * d24 + m[6] * d41 + m[14] * d12) * det_inv;
  result[15] =  (m[2] * d23 - m[6] * d13 + m[10] * d12) * det_inv;
}

void invert4x4MatrixSlow(double *m, double *result)
{
  double const determinant =
    m[0] * m[5] * m[10] * m[15] - m[0] * m[5] * m[14] * m[11] -
    m[0] * m[9] * m[6] * m[15] + m[0] * m[9] * m[14] * m[7] +
    m[0] * m[13] * m[6] * m[11] - m[0] * m[13] * m[10] * m[7] -
    m[4] * m[1] * m[10] * m[15] + m[4] * m[1] * m[14] * m[11] +
    m[4] * m[9] * m[2] * m[15] - m[4] * m[9] * m[14] * m[3] -
    m[4] * m[13] * m[2] * m[11] + m[4] * m[13] * m[10] * m[3] +
    m[8] * m[1] * m[6] * m[15] - m[8] * m[1] * m[14] * m[7] -
    m[8] * m[5] * m[2] * m[15] + m[8] * m[5] * m[14] * m[3] +
    m[8] * m[13] * m[2] * m[7] - m[8] * m[13] * m[6] * m[3] -
    m[12] * m[1] * m[6] * m[11] + m[12] * m[1] * m[10] * m[7] +
    m[12] * m[5] * m[2] * m[11] - m[12] * m[5] * m[10] * m[3] -
    m[12] * m[9] * m[2] * m[7] + m[12] * m[9] * m[6] * m[3];

  if(determinant == 0.0)
    std::cerr << "Inversion of a singular matrix!" << std::endl;

  result[0] = (m[5] * m[10] * m[15] - m[5] * m[14] * m[11] -
	       m[9] * m[6] * m[15] + m[9] * m[14] * m[7] +
	       m[13] * m[6] * m[11] - m[13] * m[10] * m[7]) / determinant;
  result[1] = -(m[1] * m[10] * m[15] - m[1] * m[14] * m[11] -
		m[9] * m[2] * m[15] + m[9] * m[14] * m[3] +
		m[13] * m[2] * m[11] - m[13] * m[10] * m[3]) / determinant;
  result[2] = (m[1] * m[6] * m[15] - m[1] * m[14] * m[7] -
	       m[5] * m[2] * m[15] + m[5] * m[14] * m[3] +
	       m[13] * m[2] * m[7] - m[13] * m[6] * m[3]) / determinant;
  result[3] = -(m[1] * m[6] * m[11] - m[1] * m[10] * m[7] -
		m[5] * m[2] * m[11] + m[5] * m[10] * m[3] +
		m[9] * m[2] * m[7] - m[9] * m[6] * m[3]) / determinant;
  result[4] = -(m[4] * m[10] * m[15] - m[4] * m[14] * m[11] -
		m[8] * m[6] * m[15] + m[8] * m[14] * m[7] +
		m[12] * m[6] * m[11] - m[12] * m[10] * m[7]) / determinant;
  result[5] = (m[0] * m[10] * m[15] - m[0] * m[14] * m[11] -
	       m[8] * m[2] * m[15] + m[8] * m[14] * m[3] +
	       m[12] * m[2] * m[11] - m[12] * m[10] * m[3]) / determinant;
  result[6] = -(m[0] * m[6] * m[15] - m[0] * m[14] * m[7] -
		m[4] * m[2] * m[15] + m[4] * m[14] * m[3] +
		m[12] * m[2] * m[7] - m[12] * m[6] * m[3]) / determinant;
  result[7] = (m[0] * m[6] * m[11] - m[0] * m[10] * m[7] -
	       m[4] * m[2] * m[11] + m[4] * m[10] * m[3] +
	       m[8] * m[2] * m[7] - m[8] * m[6] * m[3]) / determinant;
  result[8] = (m[4] * m[9] * m[15] - m[4] * m[13] * m[11] -
	       m[8] * m[5] * m[15] + m[8] * m[13] * m[7] +
	       m[12] * m[5] * m[11] - m[12] * m[9] * m[7]) / determinant;
  result[9] = -(m[0] * m[9] * m[15] - m[0] * m[13] * m[11] -
		m[8] * m[1] * m[15] + m[8] * m[13] * m[3] +
		m[12] * m[1] * m[11] - m[12] * m[9] * m[3]) / determinant;
  result[10] = (m[0] * m[5] * m[15] - m[0] * m[13] * m[7] -
		m[4] * m[1] * m[15] + m[4] * m[13] * m[3] +
		m[12] * m[1] * m[7] - m[12] * m[5] * m[3]) / determinant;
  result[11] = -(m[0] * m[5] * m[11] - m[0] * m[9] * m[7] -
		 m[4] * m[1] * m[11] + m[4] * m[9] * m[3] +
		 m[8] * m[1] * m[7] - m[8] * m[5] * m[3]) / determinant;
  result[12] = -(m[4] * m[9] * m[14] - m[4] * m[13] * m[10] -
		 m[8] * m[5] * m[14] + m[8] * m[13] * m[6] +
		 m[12] * m[5] * m[10] - m[12] * m[9] * m[6]) / determinant;
  result[13] = (m[0] * m[9] * m[14] - m[0] * m[13] * m[10] -
		m[8] * m[1] * m[14] + m[8] * m[13] * m[2] +
		m[12] * m[1] * m[10] - m[12] * m[9] * m[2]) / determinant;
  result[14] = -(m[0] * m[5] * m[14] - m[0] * m[13] * m[6] -
		 m[4] * m[1] * m[14] + m[4] * m[13] * m[2] +
		 m[12] * m[1] * m[6] - m[12] * m[5] * m[2]) / determinant;
  result[15] = (m[0] * m[5] * m[10] - m[0] * m[9] * m[6] -
		m[4] * m[1] * m[10] + m[4] * m[9] * m[2] +
		m[8] * m[1] * m[6] - m[8] * m[5] * m[2]) / determinant;
}
