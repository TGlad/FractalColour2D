#include "stdafx.h"
#include "bmp.h"
#include <set>
#include <sstream>
#include <fstream>
#include "Core\EigenBase.h"
static int width = 2000; // 3200
static int height = width; // keep the same

void setPixel(vector<BYTE> &out, const Vector2i &pos, double col)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= height)
    return;
  int ind = 3 * (pos[0] + width * (height - 1 - pos[1]));
  out[ind + 0] = (int)(255.0*col);
  out[ind + 1] = (int)(255.0*col);
  out[ind + 2] = (int)(255.0*col);
}
void setPixel(vector<BYTE> &out, const Vector2i &pos, const Vector3d &col)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= height)
    return;
  int ind = 3 * (pos[0] + width * (height - 1 - pos[1]));
  out[ind + 0] = (int)(255.0*col[0]);
  out[ind + 1] = (int)(255.0*col[1]);
  out[ind + 2] = (int)(255.0*col[2]);
}

void drawCircle(vector<BYTE> &out, const Vector2d &pos, double radius, double col)
{
  for (int x = (int)(pos[0] - radius); x <= (int)(pos[0] + radius); x++)
  {
    for (int y = (int)(pos[1] - radius); y <= (int)(pos[1] + radius); y++)
    {
      double dist = sqrt(sqr((double)x - pos[0]) + sqr((double)y - pos[1]));
      if (dist >= radius - 2.0 && dist <= radius)
        setPixel(out, Vector2i(x, y), Vector3d(col, 0.2, 1.0 - col));
    }
  }
}
void drawDisk(vector<BYTE> &out, const Vector2d &pos, double radius, double col)
{
  for (int x = (int)(pos[0] - radius); x <= (int)(pos[0] + radius); x++)
  {
    for (int y = (int)(pos[1] - radius); y <= (int)(pos[1] + radius); y++)
    {
      double dist = sqrt(sqr((double)x - pos[0]) + sqr((double)y - pos[1]));
      if (dist <= radius)
        setPixel(out, Vector2i(x, y), Vector3d(col, 0.2, 1.0 - col));
    }
  }
}

double mod(double x)
{
  if (x > 0.0)
    return (double)((int) x);
  return (double)((int)(1.0 - x));
}

double distanceToDisks(Vector3d p)
{
  const double rad = 0.5;
  Vector3d s0(0.0, 1.0, rad);
  Vector3d s1(sqrt(3.0) / 2.0, -0.5, rad);
  Vector3d s2(-sqrt(3.0) / 2.0, -0.5, rad);
  Vector3d t0(0.0, 1.0, 0.0);
  Vector3d t1(sqrt(3.0) / 2.0, -0.5, 0.0);
  Vector3d t2(-sqrt(3.0) / 2.0, -0.5, 0.0);
  Vector3d n0(1.0, 0.0, 0.0);
  Vector3d n1(-0.5, -sqrt(3.0) / 2.0, 0.0);
  Vector3d n2(-0.5, sqrt(3.0) / 2.0, 0.0);

  double scale = 1.0;

  int numIterations = 8;
  const double innerScale = sqrt(3.0) / (1.0 + sqrt(3.0));
  const double threshold = sqr(innerScale * 0.5);
  const double root3 = sqrt(3.0);
  double k = sqrt(0.5);

  // now modolu the space so we move to being in just the central hexagon, inner radius 0.5
  double h = p[2];
  double x = p.dot(-n2) * 2.0 / root3;
  double y = p.dot(-n1) * 2.0 / root3;
  x -= floor(x);
  y -= floor(y);
  if (x + y > 1.0)
  {
    x = 1.0 - x;
    y = 1.0 - y;
  }
  p = x*t1 - y*t2;

  // fold the space to be in a kite
  double l0 = p.dot(p);
  double l1 = (p - t1).dot(p - t1);
  double l2 = (p + t2).dot(p + t2);
  if (l1 < l0 && l1 < l2)
    p -= t1 * (2.0*t1.dot(p) - 1.0);
  else if (l2 < l0 && l2 < l1)
    p -= t2 * (2.0 * p.dot(t2) + 1.0);
  p[2] = h;

  for (int i = 0; i < numIterations; i++)
  {
    if ((p - Vector3d(0, 0, 0.5)).norm() < 0.5 || p[2] > 0.5)
      break;
    double p2 = p.squaredNorm();
    scale *= p2;
    p /= p2; // no, has to be further back 

    // now modolu the space so we move to being in just the central hexagon, inner radius 0.5
    double h = p[2];
    double x = p.dot(-n2) * 2.0 / root3;
    double y = p.dot(-n1) * 2.0 / root3;
    x -= floor(x);
    y -= floor(y);
    if (x + y > 1.0)
    {
      x = 1.0 - x;
      y = 1.0 - y;
    }
    p = x*t1 - y*t2;

    // fold the space to be in a kite
    double l0 = p.dot(p);
    double l1 = (p - t1).dot(p - t1);
    double l2 = (p + t2).dot(p + t2);
    if (l1 < l0 && l1 < l2)
      p -= t1 * (2.0*t1.dot(p) - 1.0);
    else if (l2 < l0 && l2 < l1)
      p -= t2 * (2.0 * p.dot(t2) + 1.0);
    p[2] = h;
  }
  double d = ((p - Vector3d(0, 0, 0.5*k)).norm() - 0.5*k) * scale; // the 0.4 is slightly more averaging than 0.5
  return d;
}

int _tmain(int argc, _TCHAR* argv[])
{
  long s2;
  vector<BYTE> out(width*width * 3); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is white

  double count = 0.0;
  double xyMin = -1.0;
  double xyMax = 1.0;
  double zMin = -0.25;
  double zMax = 1.75;
  double r = 0.0;// 0.025;// / 16.0;
  if (1)
  {
    for (int i = 0; i < width; i++)
    {
      double x = xyMin + (xyMax - xyMin)*(((double)i) + 0.5) / (double)width;
      for (int k = 0; k < width; k++)
      {
        double z = zMin + (zMax - zMin)*(((double)k) + 0.5) / (double)width;
        Vector3d p(0, x, z);
        double d = distanceToDisks(p);
        if (d < r)
          setPixel(out, Vector2i(i, k), 0);
      }
    }
  }
  if (0)
  {
    Vector2d rayStart(0.41, 1.0);
    Vector2d rayDir = Vector2d(0.225, -0.9).normalized();
    double depth = 0.0;
    for (int i = 0; i < 9; i++)
    {
      Vector2d pos = rayStart + rayDir*depth;
      double d = distanceToDisks(Vector3d(0, pos[0], pos[1]));
      Vector2d p((pos[0] - xyMin)*(double)width / (xyMax - xyMin), (pos[1] - zMin)*(double)width / (zMax - zMin));
      drawDisk(out, p, 6.0, (double)i / 8.0);
      drawCircle(out, p, d * (double)width / (xyMax - xyMin), (double)i / 8.0);
      depth += d*0.9;
    }
  }
  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"fordSpheres.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}
