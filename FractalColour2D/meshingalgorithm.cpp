#include "stdafx.h"
#include "bmp.h"
#include <fstream>

double scale = 760.0;
Vector2d offset(0, 0);
static vector<Vector2d> tri;

void saveSVG(const string &fileName, const vector<Vector2d> &points, double k = -1, int drawType = 0)
{
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)(2.0*scale) << "\" height = \"" << (int)scale+2 << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;

  double q = scale * 0.5;
  if (drawType == 1) // parabola
  {
    int n = 100;
    for (int j = -1; j < 7; j++)
    {
      for (int i = 0; i < n; i++)
      {
        double x = (double)i / (double)n;
        Vector2d pos(x, (double)j * 0.2 + k * sqr(x - 0.5));
        Vector2d p = scale*pos + offset;
        if (!i)
          svg << "<path d = \"M " << p[0] << " " << scale - p[1];
        else
          svg << " L " << p[0] << " " << scale - p[1];
      }
      svg << "\" stroke=\"grey\" fill = \"transparent\" stroke-width=\"1\" />\n";
    }
  }
  else if (drawType==2) // circles
  {
    svg << "<circle cx = \"" << q + offset[0] << "\" cy = \"" << scale - q + offset[1] << "\" r = \"" << scale*0.5 << "\" stroke = \"grey\" fill = \"transparent\" stroke-width = \"1\" />\n";
    double end = k == 0 ? 1.4 : 0.8;
    for (double r = 0.2; r <= end; r += 0.2)
    {
      double radius = 1.0 - (1.0 - r) / (1.0 + k*r*r);
      svg << "<circle cx = \"" << q + offset[0] << "\" cy = \"" << scale - q + offset[1] << "\" r = \"" << scale * radius*0.5 << "\" stroke = \"grey\" fill = \"transparent\" stroke-width = \"1\" />\n";
    }
  }
  for (auto &point : points)
  {
    Vector2d p = scale*point + offset;
    svg << "<circle cx = \"" << p[0] << "\" cy = \"" << scale - p[1] << "\" r = \"1\" stroke = \"black\" fill = \"solid\" stroke-width = \"1\" />\n";
  }
  svg << "</svg>" << endl;
  svg.close();
}
void saveSVG(const string &fileName, const vector<Vector2d> &points, const vector<Vector2d> &linestrip, const vector<Vector2d> &linelist)
{
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)(2.0*scale) << "\" height = \"" << (int)scale+2 << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;
  svg << "<defs>" << endl;
  svg << "  <marker id = \"arrowhead\" markerWidth = \"8\" markerHeight = \"5\" refX = \"0\" refY = \"2.5\" orient = \"auto\" stroke = \"green\" fill = \"green\" >" << endl;
  svg << "    <polygon points = \"0 0, 8 2.5, 0 5\" />" << endl;
  svg << "  </marker>" << endl;
  svg << "</defs>" << endl;

  for (auto &point : points)
  {
    Vector2d p = scale*point + offset;
    svg << "<circle cx = \"" << p[0] << "\" cy = \"" << scale - p[1] << "\" r = \"1\" stroke = \"black\" fill = \"solid\" stroke-width = \"1\" />\n";
  }

  for (int i = 0; i < (int)linestrip.size(); i++)
  {
    Vector2d p = scale*linestrip[i] + offset;
    if (!i)
      svg << "<path d = \"M " << p[0] << " " << scale - p[1];
    else
      svg << " L " << p[0] << " " << scale - p[1];
  }
  svg << "\" stroke=\"black\" fill = \"transparent\" stroke-width=\"2\" />\n"; 

  for (int i = 0; i < (int)linelist.size(); i+=2)
  {
    Vector2d p = scale*linelist[i] + offset;
    Vector2d p2 = scale*linelist[i+1] + offset;
    svg << "<path d = \"M " << p[0] << " " << scale - p[1] << " L " << p2[0] << " " << scale - p2[1] << "\" stroke=\"green\" stroke-width=\"2\" marker-end=\"url(#arrowhead)\" />\n";
  }

  svg << "</svg>" << endl;
  svg.close();
}

void wobblyLine(vector<Vector2d> &floor, const Vector2d &a, const Vector2d &b, int level)
{
  if (level == 0)
    return;
  Vector2d toB = b - a;
  Vector2d normal(toB[1], -toB[0]);
  Vector2d mid = a + 0.5*toB;
  Vector2d point = mid + normal*random(-0.2, 0.2);

  wobblyLine(floor, a, point, level - 1);
  floor.push_back(point);
  wobblyLine(floor, point, b, level - 1);
}

int _tmain(int argc, _TCHAR* argv[])
{
  srand(0);
/*  vector<Vector2d> poins;
  double k = 1.0;
  double r0 = 0.05;
  double d = 0.04;
  double sr = 0.1;
  double r = r0;
  for (int i = 0; i < 500; i++)
  {
    poins.push_back(Vector2d(r, sr));

    double l = sqrt(sqr(r + k*d*d) + sqr(d));
    double w = l / (r + k*d*d);
    sr *= w;
    r = l;
  }
  d = 0.01;
  sr = 0.1;
  r = 0.1;
  for (int i = 0; i < 500; i++)
  {
    poins.push_back(Vector2d(r, sr));

    double l = sqrt(sqr(r + k*d*d) + sqr(d));
    double w = l / (r + k*d*d);
    sr *= w;
    r = l;
  }
  for (double x = 0.0; x< 1.5; x += 0.03)
  {
    poins.push_back(Vector2d(x, 1.0 - 1.0/(pow(x, 2.0) + 1.0)));
  }
  saveSVG("testSphere.svg", poins);
  return 1;*/
#define TEST_SPHERE
#if defined TEST_SPHERE
  vector<Vector2d> poins;
  double k = 8.0;
  for (int j = 0; j < 9; j++)
  {
    for (double x = -0.25; x < 0.25; x += 0.005)
    {
      double r2 = (double)j * 0.1;
      double y = k*x*x + r2;
      Vector2d p(x, y);
      double r = p.norm();
      p.normalize();
      double R = 1.0;
      // r2 = r;
      double sr = r*(1.0+k*R)/(1.0 + k*r);
      Vector2d pos = p*sr;
      poins.push_back(Vector2d(0.5, 0.5) + 0.5*pos);
    }
  }
  saveSVG("testSphere.svg", poins);
  return 1;
#endif
  //#define SPHERICAL
#if defined SPHERICAL
  Vector3d rombuses1[3] = { Vector3d(0.1, 0, 0.9), Vector3d(0.0, 0.35, 1.0), Vector3d(0.3, 0.6, 0.4) };
  Vector3d rombuses2[3] = { Vector3d(0.1, 0.35, 0.9), Vector3d(0.25, 0.6, 0.75), Vector3d(0.3, 0.7, 0.4) };
  Vector2d mid(0.5, 0.3);
  vector<Vector2d> points;
  int count = 3000;
  for (int c = 0; c < count; c++)
  {
    Vector2d pos(random(0.0, 1.0), random(0.0, 0.7));
    bool inside = false;
    for (int i = 0; i < 3; i++)
    {
      if (pos[1] > rombuses2[i][1] || pos[1] < rombuses1[i][1])
        continue;
      double blend = (pos[1] - rombuses1[i][1]) / (rombuses2[i][1] - rombuses1[i][1]);
      double xMin = rombuses1[i][0] * (1.0 - blend) + rombuses2[i][0] * blend;
      double xMax = rombuses1[i][2] * (1.0 - blend) + rombuses2[i][2] * blend;
      if (pos[0] > xMin && pos[0] < xMax)
        inside = true;
    }
    if (inside)
      points.push_back(pos);
  }
  vector<Vector2d> scaled1;
  for (auto &p : points)
    scaled1.push_back(Vector2d(0.5, 0.5) + 0.75*(p - mid));
  saveSVG("house.svg", scaled1, 0, 2);
  double k = 20.0;

  for (auto &p : points)
  {
    p -= mid;
    double r = p.norm();
    Vector2d surface = p / r;
    Vector2d distorted = surface + (p - surface) / (1 + k*r*r);
    p = Vector2d(0.5, 0.5) + 0.5*distorted;
  }
  saveSVG("houseBent.svg", points, k, 2);

  // Next we need the internal case, where a room is seen from the inside
  Vector2d boxMin[3] = { Vector2d(0.15, 0.15), Vector2d(0.8, 0.2), Vector2d( 0.2, 0.4) };
  Vector2d boxMax[3] = { Vector2d(0.8, 0.55), Vector2d(0.9, 0.45), Vector2d(0.3, 0.5) };
  points.clear();
  for (int c = 0; c < count; c++)
  {
    Vector2d pos(random(0.0, 1.0), random(0.0, 0.7));
    double s = max(max(0.0, 0.1 - pos[0]), max(0.0, 0.1 - pos[1]));
    double t = max(max(0.0, pos[0] - 0.9), max(0.0, pos[1] - 0.6));
    if (random(0.0, 0.1) < max(s,t))
      continue;
    bool inside = false;
    for (int i = 0; i < 2; i++)
    {
      if (pos[0] > boxMin[i][0] && pos[0] < boxMax[i][0] && pos[1] > boxMin[i][1] && pos[1] < boxMax[i][1])
        inside = true;
    }
    if (pos[0] > boxMin[2][0] && pos[0] < boxMax[2][0] && pos[1] > boxMin[2][1] && pos[1] < boxMax[2][1])
      inside = false;
    if (!inside)
      points.push_back(pos);
  }
  vector<Vector2d> scaled;
  for (auto &p : points)
    scaled.push_back(Vector2d(0.5, 0.5) + 0.75*(p - mid));
  saveSVG("room.svg", scaled, 0, 2);
  mid = Vector2d(0.45, 0.35);
  for (auto &p : points)
  {
    p -= mid;
    p *= 6.0;
    double r = p.norm();
    p /= (r*r);
    r = 1.0 / r;
    Vector2d surface = p / r;
    Vector2d distorted = surface + (p - surface) / (1 + k*r*r);
    p = Vector2d(0.5, 0.5) + 0.5*distorted;
  }
  saveSVG("roomBent.svg", points, k, 2);
#else

  // First, make the scene:
  vector<Vector2d> floor;
  floor.push_back(Vector2d(0, 0));
  Vector2d end(1, 0);
  wobblyLine(floor, floor[0], end, 4);
  floor.push_back(end);
  vector<Vector2d> points;
  int numPoints = 600;
  double groundThickness = 0.1;
  for (int i = 0; i < numPoints; i++)
  {
    double k = random(0.0, 0.999)*(floor.size() - 1);
    int id = (int)k;
    double blend = k - (double)id;
    Vector2d pos = floor[id] * (1.0 - blend) + floor[id + 1] * blend;
    double y = random(0.0, 1.0);
    pos[1] += (0.5*y + y*y)*groundThickness * 0.75;
    points.push_back(pos);
  }
  // now add tree shapes
  Vector2d trunkMin[3], trunkMax[3];
  double p[3] = { 0.3, 0.5, 0.7 };
  double trunkThickness = 0.03;
  double trunkHeight[3] = { 0.2, 0.25, 0.2 };
  double topRadius[3] = { 0.2, 0.25, 0.22 };
  for (int i = 0; i < 3; i++)
  {
    double k = p[i] * (floor.size() - 1);
    int id = (int)k;
    double blend = k - (double)id;
    Vector2d pos = floor[id] * (1.0 - blend) + floor[id + 1] * blend;
    pos[1] += groundThickness;
    trunkMin[i][0] = pos[0] - trunkThickness;
    trunkMax[i][0] = pos[0] + trunkThickness;
    trunkMin[i][1] = pos[1];
    trunkMax[i][1] = pos[1] + trunkHeight[i];
  }
  int numSamples = 5000;
  for (int i = 0; i < numSamples; i++)
  {
    Vector2d pos(random(0, 1), random(0, 1));
    bool found = false;
    for (int j = 0; j < 3; j++)
    {
      if (pos[0] > trunkMin[j][0] && pos[0] < trunkMax[j][0] && pos[1] > trunkMin[j][1] && pos[1] < trunkMax[j][1])
        found = true; 
      Vector2d dif = Vector2d((trunkMin[j][0] + trunkMax[j][0])*0.5, trunkMax[j][1] + topRadius[j]) - pos;
      if (dif.norm() < topRadius[j])
        found = true;
    }
    if (found)
      points.push_back(pos);
  }

  // Then save the easy ones
  saveSVG("terrain.svg", points, 0, 1);
  vector<Vector2d> points2 = points;
  double curvature = 1.5;
  for (auto &p : points2)
    p[1] += curvature * sqr(p[0] - 0.5);
  saveSVG("terrainBentUp.svg", points2, curvature, 1);
  curvature = -1.0;
  points2 = points;
  for (auto &p : points2)
    p[1] += curvature * sqr(p[0] - 0.5);
  saveSVG("terrainBentDown.svg", points2, curvature, 1);

  for (int step = 4; step > 1; step /= 2)
  {
    // Next I want to draw the terrain with low level vertex normals
    vector<Vector2d> linelist;
    vector<Vector2d> line;
    vector<Vector2d> faceNormals;
    if (step == 4)
    {
      for (int i = 0; i < (int)floor.size(); i += step)
        line.push_back(floor[i]);
    }
    else
    {
      line.push_back(Vector2d(0, 0));
      line.push_back(Vector2d(0.1, 0.042));
      line.push_back(Vector2d(0.2, 0.003));
      line.push_back(Vector2d(0.31, 0));

      line.push_back(Vector2d(0.485, 0.03));
      line.push_back(Vector2d(0.55, 0));

      line.push_back(Vector2d(0.67, 0.03));
      line.push_back(Vector2d(0.77, 0.0));
      line.push_back(Vector2d(0.836, 0));

      line.push_back(Vector2d(1.0, 0.048));
      line.push_back(Vector2d(1.123, 0));
    }
    for (int i = 1; i < line.size(); i++)
    {
      Vector2d up(line[i - 1][1] - line[i][1], line[i][0] - line[i - 1][0]);
      faceNormals.push_back(up.normalized());
    }
    vector<Vector2d> vertexNormals;
    vertexNormals.push_back(faceNormals[0]);
    for (int i = 1; i < faceNormals.size(); i++)
      vertexNormals.push_back((faceNormals[i] + faceNormals[i - 1]).normalized());
    vertexNormals.push_back(faceNormals.back());
    double normalLength = 0.1;
    for (int i = 0; i < vertexNormals.size(); i++)
    {
      linelist.push_back(line[i]);
      linelist.push_back(line[i] + vertexNormals[i] * normalLength);
    }
    stringstream st;
    st << "lineAndVertices" << step << ".svg";
    saveSVG(st.str(), points, line, linelist);

    // Now I want to bend down the points
    vector<Vector2d> line2;
    line2.push_back(line[0]);
    for (int i = 1; i < line.size(); i++)
    {
      double dist = (line[i] - line[i - 1]).norm();
      line2.push_back(Vector2d(line2.back()[0] + dist, 0));
    }
    vector<Vector2d> linelist2;
    for (auto &p : line2)
    {
      linelist2.push_back(p);
      linelist2.push_back(p + Vector2d(0, normalLength));
    }
    points2.clear();
    for (auto &p : points)
    {
      for (int i = 0; i < faceNormals.size(); i++)
      {
        double height = (p - line[i]).dot(faceNormals[i]);
        if (height > 0.0 && height < (double)step * 0.1)
        {
          Vector2d side1(vertexNormals[i][1], -vertexNormals[i][0]);
          Vector2d side2(vertexNormals[i + 1][1], -vertexNormals[i + 1][0]);
          if ((p - line[i]).dot(side1) > 0.0 && (p - line[i + 1]).dot(side2) < 0.0)
          {
            // transform:
            Vector2d p0 = line[i] + height*vertexNormals[i] / vertexNormals[i].dot(faceNormals[i]);
            Vector2d p1 = line[i + 1] + height*vertexNormals[i + 1] / vertexNormals[i + 1].dot(faceNormals[i]);

            double blend = (p - p0).dot(p1 - p0) / (p1 - p0).squaredNorm();
            points2.push_back(line2[i] * (1.0 - blend) + line2[i + 1] * blend + Vector2d(0, height));
          }
        }
      }
    }
    stringstream st2;
    st2 << "lineAndVerticesBent" << step << ".svg";
    saveSVG(st2.str(), points2, line2, linelist2);
    points = points2;
  }
#endif
}
