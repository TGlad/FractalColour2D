// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>
#include <iostream>
#include <vector>
#include <windows.h>
#include <strstream>

#include "/Code/Eigen/Eigen"
#include "/Code/Eigen/StdVector"
using namespace std;
using namespace Eigen;
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Vector2d);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Vector3d);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Vector4d);
const double pi = 3.14159265;
inline double sqr(double x)
{
  return x*x;
}
inline double random(double from, double to)
{
  return from + (to - from)*(rand() % 10000) / 10000.0;
}