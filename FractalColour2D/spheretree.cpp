#include "stdafx.h"
#include "bmp.h"
#include <set>
#include <sstream>
#include <fstream>
#include "Core\EigenBase.h"
static int width = 100;
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

double distanceToTree(Vector3d p)
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
  for (int i = 0; i < numIterations; i++)
  {
    if ((p - Vector3d(0, 0, innerScale*0.5)).squaredNorm() < threshold)
      break; // definitely inside
    double maxH = 0.4;
    if (i == 0)
      maxH = -100;
    if (p[2] > maxH && (p - Vector3d(0, 0, 0.5*1.1)).squaredNorm() > sqr(0.5*1.1))
      break; // definitely outside
    double p2 = p.squaredNorm();
    if (p[2] < maxH && (p - Vector3d(0, 0, 0.5)).squaredNorm() > 0.5*0.5)
    {
      // needs a sphere inverse
      scale *= p2;
      p /= p2; // no, has to be further back 
    }
    else
    {
      // stretch onto a plane at zero 
      scale *= p2;
      p /= p2;
      p[2] -= 1.0;
      p[2] *= -1.0;
      p *= root3;
      scale /= root3;
      p[2] += 1.0;

      // and rotate it a twelfth
      double a = 3.1415 / 6.0;
      double sina = sin(a);
      double cosa = cos(a);
      double xx = p[0]*cosa + p[1]*sina;
      double yy = -p[0]*sina + p[1]*cosa;
      p[0] = xx;
      p[1] = yy;
    }
    // now modolu the space so we move to being in just the central hexagon, inner radius 0.5
    double h = p[2];
    //   p.x -= 1.0/sqrt(3.0);
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
//#define SURFACE
#define POINTS
#if defined SURFACE
  double d = ((p - Vector3d(0, 0, 0.4)).norm() - 0.4) * scale; // the 0.4 is slightly more averaging than 0.5
#elif defined POINTS
  double d = (p - Vector3d(0, 0, 0.4)).norm() * scale; // the 0.4 is slightly more averaging than 0.5
#elif 
  double d = ((p - Vector3d(0, 0, 0.4)).norm() - 0.4) * scale; // the 0.4 is slightly more averaging than 0.5
#endif

  return d;
}

int _tmain(int argc, _TCHAR* argv[])
{
//  double rs[] = { 0.025, 0.0125, 0.00625, 0.003125, 0.0015625, 0.00104167, 0.0005208 };
  double rs[] = { 0.025, 0.0125, 0.00625, 0.003125, 0.0015625, 0.00078125 };
//  double vs[] = { 0.329, 0.229, 0.160, 0.117, 0.0863, 0.0723, 0.0534 };
  double vs[] = { 0.217, 0.158, 0.11685, 0.0861, 0.0637, 0.047 };
  const int sz = 6;
  double ms[sz];

  long s2;
  vector<BYTE> out(width*width * 3); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is white

  if (1)
  {
    int end = sz - 1;
  double minx = min(log(rs[0]), log(rs[end])), maxx = max(log(rs[0]), log(rs[end]));
  double miny = min(log(vs[0]), log(vs[end])), maxy = max(log(vs[0]), log(vs[end]));

  for (int i = 0; i < sz; i++)
  {
    double x = log(rs[i]);
    ms[i] = 1.0 / rs[i];
    double y = log(vs[i]);
    double X = 5.0 + (double)(width - 10)*(x - minx) / (maxx - minx);
    double Y = 5.0 + (double)(width - 10)*(y - miny) / (maxy - miny);
    setPixel(out, Vector2i(X, Y), 0);
  }
  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"spheretree_graph_points.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;

  int a = 2; // 0
  int b = sz-1; // 3
  double gradient = (log((ms[a] * ms[a] * ms[a]) * vs[a]) - log((ms[b] * ms[b] * ms[b]) * vs[b])) / (log(ms[a]) - log(ms[b]));
  cout << "dimension: " << gradient << endl;

  for (int i = 0; i < sz; i++)
  {
    double c = ((ms[i] * ms[i] * ms[i]) * vs[i]) / pow(ms[i], gradient);
    cout << "c: " << c << endl;
  }
  }

  double count = 0.0;
  double xyMin = -0.525;
  double xyMax = 0.525;
  double zMin = -0.025;
  double zMax = 1.025;
  double r = 0.00625 / 8.0;
  for (int i = 0; i < width; i++)
  {
    double x = xyMin + (xyMax - xyMin)*(((double)i) + 0.5) / (double)width;
    for (int j = 0; j < width; j++)
    {
      double y = xyMin + (xyMax - xyMin)*(((double)j) + 0.5) / (double)width;
      for (int k = 0; k < width; k++)
      {
        double z = zMin + (zMax - zMin)*(((double)k) + 0.5) / (double)width;
        Vector3d p(x, y, z);
        double d = distanceToTree(p);
#if defined SURFACE
        if (abs(d) < r)
#elif defined POINTS
        if (d < r)
#elif 
        if (d < 0.0)
#endif
        {
          count++;
          if (j == width / 2)
            setPixel(out, Vector2i(i, k), 0);
        }
      }
    }
  }
  double step = (xyMax - xyMin) / (double)width;
  double volume = count * step*step*step;
  cout << "r: " << r << ", volume: " << volume << endl;

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"spheretreepoints.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}
