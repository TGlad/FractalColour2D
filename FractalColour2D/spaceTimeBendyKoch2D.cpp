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
double scale = 1500.0;
void openSVG(const string &fileName)
{
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)scale << "\" height = \"" << (int)(scale*1.05) << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;
}
static int icount = 0;
void saveSVG(const vector<Vector3d> &points, const vector<Vector3d> &circs)
{
  Vector2d mins(1e10,1e10), maxs(-1e10,-1e10);
 // for (int i = 0; i < points.size(); i++)
  {
    mins[0] = -1;// min(mins[0], points[i][0]);
    mins[1] = -1.5;// min(mins[1], points[i][1]);
    maxs[0] = 2.0;// max(maxs[0], points[i][0]);
    maxs[1] = 1.5;// max(maxs[1], points[i][1]);
  }
  double mult = 1.0 / max(maxs[0] - mins[0], maxs[1] - mins[1]);
  svg << "<path d = \"M " << scale*(points[0][0] - mins[0]) *mult << " " << scale*(points[0][1] - mins[1]) *mult;
  for (int i = 1; i < (int)points.size(); i++)
    svg << " L " << scale*(points[i][0] - mins[0]) *mult << " " << scale*(points[i][1] - mins[1]) *mult;
  string cols[] = { "blue", "red", "grey", "orange", "purple" };
  svg << "\" fill=\"transparent\" stroke=\"" << cols[icount++] << "\" />\n";
  /*
  for (int i = 0; i < circs.size(); i++)
  {
    svg << "<circle cx = \"" << scale*(circs[i][0] - mins[0])*mult <<
      "\" cy = \"" << scale*(circs[i][1] - mins[1])*mult <<
      "\" r = \"" << scale*circs[i][2] * mult << "\" stroke = \"black\" stroke-width = \"1\" fill=\"none\" />\n";
  }*/
}

void closeSVG()
{
  svg << "</svg>" << endl;
  svg.close();
}
static double scales[20] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

#if 0
double findLength(double R, double angle, double maxAngle, Vector3d &sideVec, double &height)
{
  Vector2d m = R * Vector2d(sin(maxAngle - angle), cos(maxAngle - angle));
  double endAngle = maxAngle - angle / 2.0;
  Vector2d vec = R*Vector2d(sin(endAngle), cos(endAngle));
  vec /= vec[0];
  vec *= 2 * m[0];
  double h = 2 * m[1] - vec[1];
  Vector2d mid = m - Vector2d(0, h);
  double rMid = mid.norm();
  double angMid = atan2(mid[0], mid[1]);
  Vector2d end = vec - m;
  double angEnd = atan2(end[0], end[1]) - endAngle;
  double rEnd = end.norm();
  height = h;
  sideVec = Vector3d(vec[0], vec[1], 0);
  return 2.0 * rMid*angMid + 4.0 * rEnd*angEnd;
}
#else
static double bendAngle = 0;
double findLength(double R, double angle, double maxAngle, Vector3d &sideVec, Vector3d &vecMid, double &height, double &rr1, double &rr2)
{
  double phi = maxAngle - angle;
  double xi = phi + bendAngle;
  double r1 = R*(sin(phi) / sin(xi));
  double X = 2.0*bendAngle - angle;
  double r2 = R*(sin(angle / 2.0) / sin(X / 2.0));
  double length = 2.0*(abs(r1*xi) + abs(r2*X));
  
  double endAngle = angle / 2.0;
  sideVec = (R*cos(endAngle) + r2*cos(X/2.0))*Vector3d(sin(phi + endAngle), cos(phi + endAngle), 0);
  vecMid = R*Vector3d(sin(phi), cos(phi), 0);
  height = R*cos(phi) - r1*cos(xi);

  rr1 = abs(r2);
  rr2 = abs(r1);
  return abs(length);
}
#endif
// simple binary search...
// returns the three nodes, the two vertices are midway between
void findNodes(const Vector3d &a, const Vector3d &b, Vector3d *nodes, Vector3d *verts, bool flip, int index)
{
  double r1, r2;
  double R = a.norm();
  double dot = a.dot(b);
  double crss = a.dot(Vector3d(-b[1], b[0], 0));
  double maxAngle = atan2(crss, dot) / 2.0;
  if (index == 1 && maxAngle < 0.0 && !flip)
    maxAngle += pi;
  if (crss == 0.0)
    maxAngle = pi;
  if (flip)
    maxAngle = -maxAngle;
  double length = abs(R * 2.0*maxAngle);
  bendAngle = 0.5;
  double scaleUp = scales[index]; //  bendAngle / sin(bendAngle);
  for (int i = 0; i < 100; i++)
  {
    double ang = 3.14*((double)i) / 100.0;
    if ((ang / sin(ang)) > scales[index])
    {
      bendAngle = ang;
      break;
    }
  }
  double desiredLength = length*scaleUp;

  double height;
  Vector3d vecSide, vecMid;
#if 1
  const int count = 50;
  double startLengths[count];
  double startAngles[count];
  double minLength = 1e10;
  int minI = 0;
  double minAng = 0;
  for (int i = 0; i < count; i++)
  {
    startAngles[i] = maxAngle * (0.5 + (double)i) / (double)count;
  //  startAngles[i] = (0.5*startAngles[i]) + 0.25*maxAngle;
    startLengths[i] = abs(findLength(R, startAngles[i], maxAngle, vecSide, vecMid, height, r1, r2) - desiredLength);
    if (startLengths[i] < minLength)
    {
      minI = i;
      minLength = startLengths[i];
    }
  }
  findLength(R, startAngles[minI], maxAngle, vecSide, vecMid, height, r1, r2);
  /*
  double angles[3];
  double lengths[3];
  for (int i = -1; i < 2; i++)
  {
    angles[i + 1] = startAngles[minI + i];
    lengths[i + 1] = startLengths[minI + i];
  }
  for (int i = 0; i < 10; i++)
  {
    angles[1] = (angles[0] + angles[2]) / 2.0;
    double diff = (angles[2] - angles[0]) / 2.0;
    //    angles[2] += diff*0.1;
    //    angles[0] -= diff*0.1;
    //    angles[0] = max(angles[0], 1e-6);
    //    angles[2] = min(angles[2], maxAngle - 1e-6);
    lengths[1] = findLength(R, angles[1], maxAngle, vecSide, vecMid, height, r1, r2);
    lengths[1] = abs(lengths[1] - desiredLength);

    if (lengths[0] < lengths[2])
    {
      lengths[2] = lengths[1];
      angles[2] = angles[1];
    }
    else
    {
      lengths[0] = lengths[1];
      angles[0] = angles[1];
    }
  }*/
#else
  findLength(R, maxAngle * 2.0/3.0, maxAngle, vecSide, vecMid, height, r1, r2);
#endif
  Vector3d mid = (a + b).normalized();
  double dir = 1.0;
  if (index == 1.0 && crss < 0.0 && !flip)
    mid = -mid;
  if (crss == 0.0)
    mid = -mid;
  if (flip)
    dir = -1;
  Vector3d midside(mid[1], -mid[0], 0);
  nodes[0] = mid*vecSide[1] - dir*midside*vecSide[0];
  nodes[1] = mid*height;
  nodes[2] = mid*vecSide[1] + dir*midside*vecSide[0];
  nodes[0][2] = r1;
  nodes[2][2] = r1;
  nodes[1][2] = r2;
  verts[0] = mid*vecMid[1] - dir*midside*vecMid[0];
  verts[1] = mid*vecMid[1] + dir*midside*vecMid[0];
}
static vector<Vector3d> circs;
void arcRecurse(vector<Vector3d> &points, const Vector3d &a, const Vector3d &o, const Vector3d &b, int level, bool flip, int index)
{
  const double minLength = 0.005;
  Vector3d ns[3];
  Vector3d vs[2];
  Vector3d o2 = o;
  o2[2] = 0.0;
  findNodes(a - o2, b - o2, ns, vs, flip, index);
  for (int i = 0; i<3; i++)
    ns[i] += o2;
  for (int i = 0; i<2; i++)
    vs[i] += o2;
  if (level>0)//(vs[0] - a).norm() > minLength)
    arcRecurse(points, a, ns[0], vs[0], level - 1, !flip, index+1);
  else
    circs.push_back(ns[0]);
  points.push_back(vs[0]);
  if (level>0)//(vs[1] - vs[0]).norm() > minLength)
    arcRecurse(points, vs[0], ns[1], vs[1], level-1, flip, index+1);
  else
    circs.push_back(ns[1]);
  points.push_back(vs[1]);
  if (level>0)//(b - vs[1]).norm() > minLength)
    arcRecurse(points, vs[1], ns[2], b, level-1, !flip, index+1);
  else
    circs.push_back(ns[2]);
}


vector<Vector3d> curve;
static double roothalf = 1.0 / sqrt(2.0);
static Vector3d u(roothalf, -roothalf, 0);
static Vector3d v(roothalf, roothalf, 0);
static double downScale = 0.0;
static double boostAngle = 1.2;
static double dimension = 0.0;
/*
void addKochChild(int order, const Vector3d &p0, const Vector3d &p1, const Vector3d &centre, bool flip = false)
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
    addKochChild(order - 1, p0, mid, !flip);
  curve.push_back(mid);
  if (order > 0)
    addKochChild(order - 1, mid, p1, !flip);
}*/
int _tmain(int argc, _TCHAR* argv[])
{
  openSVG("bendy2.svg");
  double dims[] = { 0.75, 0.5, 0, -1 };
//  double dims[] = { 0.99, 0.95, 0.9, 0.8 };
  for (int d = 1; d < 4; d++)
  {
    dimension = dims[d];
    double time0 = pow(2.0, 1.0 - dimension);

    Vector3d proj(1.0, 0, 0);
    double psize = 0;
    for (int i = 0; i <= 20; i++)
    {
      double a = roothalf*v.dot(proj);
      double b = -time0*proj[0];
      double c = roothalf*u.dot(proj);
      double b24ac = b*b - 4.0*a*c;
      double angle = (-b + sqrt(b24ac)) / (2.0*a);

      Vector3d proja = u*u.dot(proj) / angle + v*v.dot(proj)*angle;
      proja /= proja[0];
      if (i)
        scales[i - 1] = proja[1] / proj[1];
      else
        psize = proja[1];
      proj = proja;
    }

    curve.clear();
    curve.push_back(Vector3d(-1, 0, 0));
    Vector3d curveEnd(-1, 0, 0);
    Vector3d centre(0, 0, 0);
    //addKochChild(8, curve[0], curveEnd, centre, true);
    arcRecurse(curve, curve[0], centre, curveEnd, 8, false, 0);
    curve.push_back(curveEnd);
 
    for (int i = 0; i < curve.size(); i++)
      curve[i] = Vector3d(-1, 0, 0) + (curve[i] - Vector3d(-1, 0, 0))*psize;
    saveSVG(curve, circs);
  }
  closeSVG();
  return 1;
}