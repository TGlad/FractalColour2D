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
void saveSVG(const Vector2d &offset, const vector<Vector3d> &points)
{
  double s = 2.0;
  svg << "<path d = \"M " << scale*(s*points[0][2] + offset[1]) << " " << scale*(1.0 - s*points[0][1] - offset[0]);
  for (int i = 1; i < (int)points.size(); i++)
    svg << " L " << scale*(s*points[i][2] + offset[1]) << " " << scale*(1.0 - s*points[i][1] - offset[0]);
  string cols[] = { "blue", "red", "grey", "orange", "purple" };
  svg << "\" fill=\"transparent\" stroke=\"" << cols[icount++] << "\" />\n";
}

void closeSVG()
{
  svg << "</svg>" << endl;
  svg.close();
}

vector<Vector3d> curve;
static double roothalf = 1.0 / sqrt(2.0);
static Vector3d u(roothalf, -roothalf, 0);
static Vector3d v(roothalf, roothalf, 0);
static double downScale = 0.0;
static double angle = 1.2;
static double dimension = 0.0;
void addKochChild(int order, const Vector3d &p0, const Vector3d &p1, bool flip = false)
{
  double time0 = pow(2.0, 1.0 - dimension);
  // time0 = half/angle + half*angle
  // time0 = roothalf*u.dot(dt,dist)/angle + roothalf*v.dot(dt,dist)*angle
  //       = 
  // roothalf*v.dot(dt,dist)*angle^2 - time0*angle + roothalf*u.dot(dt,dist) = 0
  Vector3d dir = p1 - p0;
  double dt = dir[0];
  Vector3d dirn = dir;
  dirn[0] = 0;
  double dist = dirn.norm();
  Vector3d proj(dt, dist, 0);
  
  
  double a = roothalf*v.dot(proj);
  double b = -time0*dt;
  double c = roothalf*u.dot(proj);
  double b24ac = b*b - 4.0*a*c;
  angle = (-b + sqrt(b24ac)) / (2.0*a);
  
  Vector3d proja = u*u.dot(proj) / angle + v*v.dot(proj)*angle;
  double downScale = dt / proja[0];
  double len = proja[1] * downScale;
  if (len < dist)
    cout << "bad" << endl;
  double lift = 0.5*sqrt(len*len - dist*dist);




  Vector3d orth(0.0, -dir[2], dir[1]);
  if (dist == 0)
    orth = Vector3d(0, 1, 0);
  else
    orth.normalize();
  if (flip)
    orth = -orth;
  Vector3d mid = p0 + dirn*0.5 + orth*lift + Vector3d(0.5*proja[0]*downScale, 0, 0);

  if (order > 0)
    addKochChild(order - 1, p0, mid, flip);
  curve.push_back(mid);
  if (order > 0)
    addKochChild(order - 1, mid, p1, flip);
}
int _tmain(int argc, _TCHAR* argv[])
{
  openSVG("levy2D2.svg");
  double dims[] = { 0.99, 0.95, 0.75, 0.5, 0 };
//  double dims[] = { 0.99, 0.95, 0.9, 0.8 };
  for (int d = 0; d < 5; d++)
  {
    dimension = dims[d]; //  0.75 - 0.25*(double)d;
    double time0 = pow(2.0, 1.0 - dimension);
    // time0 = half/angle + half*angle
    // half*angle^2 - time0*angle + half = 0
    double a = 0.5;
    double b = -time0;
    double c = 0.5;
    double b24ac = b*b - 4.0*a*c;
    angle = (-b + sqrt(b24ac)) / (2.0*a);

    Vector3d time(1, 0, 0);
    time = u*u.dot(time) / angle + v*v.dot(time)*angle;
    downScale = 1.0 / time[0];

    curve.push_back(Vector3d(0, 0, 0));
    Vector3d curveEnd(1, 0, 0);
    addKochChild(8, curve[0], curveEnd, true);
    curve.push_back(curveEnd);
    saveSVG(Vector2d(0.95, 0.4), curve);
  }
  closeSVG();
  return 1;
}