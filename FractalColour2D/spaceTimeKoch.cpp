// Folding.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <vector>
#include "/Code/Eigen/Eigen"
#include "/Code/Eigen/StdVector"
using namespace std;
using namespace Eigen;
//EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Vector2d);
//EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Vector3d);


#include <fstream>
#include <windows.h>
#include <stdio.h>       // for memset
static ofstream svg;
double scale = 800.0;
void openSVG(const string &fileName)
{
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)scale << "\" height = \"" << (int)(scale*1.05) << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;
}
static int icount = 0;
void saveSVG(const Vector2d &offset, const vector<Vector2d> &points)
{
  svg << "<path d = \"M " << scale*(points[0][1] + offset[1]) << " " << scale*(points[0][0] + offset[0]);
  for (int i = 1; i < (int)points.size(); i++)
    svg << " L " << scale*(points[i][1] + offset[1]) << " " << scale*(1.0 - points[i][0] - offset[0]);
  string cols[] = { "blue", "red", "grey", "orange", "purple" };
  svg << "\" fill=\"transparent\" stroke=\"" << cols[icount++] << "\" />\n";

  vector<Vector2d> ps;
  ps.push_back(Vector2d(0, 0));   ps.push_back(Vector2d(1, 0));
  ps.push_back(Vector2d(0, -0.5));   ps.push_back(Vector2d(0, 0.5));
  ps.push_back(Vector2d(0, 0));   ps.push_back(Vector2d(0.5, -0.5));
  ps.push_back(Vector2d(0, 0));   ps.push_back(Vector2d(0.5, 0.5));

  for (int i = 0; i < ps.size(); i+=2)
  {
    // axes
    svg << "<path d = \"M " << scale*(ps[i][1] + offset[1]) << " " << scale*(1.0 - ps[i][0] - offset[0]);
    svg << " L " << scale*(ps[i+1][1] + offset[1]) << " " << scale*(1.0 - ps[i+1][0] - offset[0]);
    if (i<=2)
      svg << "\" fill=\"transparent\" stroke=\"black\" />\n";
    else
      svg << "\" fill=\"transparent\" stroke=\"green\" />\n";
  }
}

void closeSVG()
{
  svg << "</svg>" << endl;
  svg.close();
}

vector<Vector2d> curve;
static double roothalf = 1.0 / sqrt(2.0);
static Vector2d u(roothalf, -roothalf);
static Vector2d v(roothalf, roothalf);
static double downScale = 0.0;
static double angle = 1.2;
void addKochChild(int order, const Vector2d &p0, const Vector2d &p1, bool flip = false)
{
#if 0 // koch
  const double lift = 0.3;
  Vector2d dir(p1[1] - p0[1], p0[0] - p1[0]);
  if (flip)
    dir = -dir;
  Vector2d mid = (p0 + p1)*0.5 + dir * lift;
#else
  Vector2d dir = (p1 - p0)/2.0;
  if (flip)
    dir = u*u.dot(dir) / angle + v*v.dot(dir)*angle;
  else
    dir = u*u.dot(dir) * angle + v*v.dot(dir)/angle;
  dir *= downScale;
  Vector2d mid = p0 + dir;
#endif

  if (order > 0)
    addKochChild(order - 1, p0, mid, !flip);
  curve.push_back(mid);
  if (order > 0)
    addKochChild(order - 1, mid, p1, !flip);
}
int _tmain(int argc, _TCHAR* argv[])
{
  openSVG("kochy.svg");
  double dims[] = { 0.75, 0.5, 0, -1 };
  for (int d = 0; d < 4; d++)
  {
    double dimension = dims[d]; // 0.75 - 0.25*(double)d;
    double time0 = pow(2.0, 1.0 - dimension);
    // time0 = half/angle + half*angle
    // half*angle^2 - time0*angle + half = 0
    double a = 0.5;
    double b = -time0;
    double c = 0.5;
    double b24ac = b*b - 4.0*a*c;
    angle = (-b + sqrt(b24ac)) / (2.0*a);

    Vector2d time(1, 0);
    time = u*u.dot(time) / angle + v*v.dot(time)*angle;
    downScale = 1.0 / time[0];

    curve.push_back(Vector2d(0, 0));
    Vector2d curveEnd(1, 0);
    addKochChild(10, curve[0], curveEnd, true);
    curve.push_back(curveEnd);
    saveSVG(Vector2d(0, 0.5), curve);
  }
  closeSVG();
  return 1;
}