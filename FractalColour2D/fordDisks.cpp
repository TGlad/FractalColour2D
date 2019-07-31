#include "stdafx.h"
#include "bmp.h"
#include <set>
#include <sstream>
#include <fstream>
#include "Core\EigenBase.h"
static int width = 3200; // 3200
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

double distanceToDisks(Vector2d p)
{
  double scale = 1.0; 

  int numIterations = 8;
  int i = 0;
  for (i = 0; i < numIterations; i++)
  {
    p[0] -= floor(p[0]);
    if (p[0] > 0.5)
      p[0] = 1.0 - p[0];
    if ((p - Vector2d(0, 0.5)).squaredNorm() < 0.5*0.5)
      break; // definitely inside
    if (p[1] < 0.0 || p[1] > 0.5)
      break;
    double p2 = p.squaredNorm();
    // needs a sphere inverse
    scale *= p2;
    p /= p2; // no, has to be further back 
  }
//#define SURFACE
#define POINTS
#if defined POINTS
  double d = (p - Vector2d(0, 0.5)).norm() * scale; // the 0.4 is slightly more averaging than 0.5
#else
  double d2 = abs(p[1]) * scale;
  double d3 = abs(p[1] - 1.0) * scale;
  double d = ((p - Vector2d(0, 0.5)).norm() - 0.5) * scale; // the 0.4 is slightly more averaging than 0.5
  d = min(d, d2);
  if (i != 0)
    d = min(d, d3);
#endif

  return d;
}

int _tmain(int argc, _TCHAR* argv[])
{
  double rs[] = { 0.025, 0.0125, 0.00625, 0.003125, 0.0015625, 0.00078125, 0.0003609};
  double vs[] = { 0.204, 0.1158, 0.0655, 0.0366, 0.020175, 0.01113, 0.006095};
  const int sz = 7;
  double ms[sz];

  long s2;
  vector<BYTE> out(width*width * 3); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is white

  if (0)
  {
    int end = sz-1;
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
  LPCTSTR file = L"fordDisksPoints_graph.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;

  int a = 1; // 0
  int b = sz-1; // 3
  double gradient = (log((ms[a] * ms[a]) * vs[a]) - log((ms[b] * ms[b]) * vs[b])) / (log(ms[a]) - log(ms[b]));
  cout << "dimension: " << gradient << endl;

  for (int i = 0; i < sz; i++)
  {
    double c = ((ms[i] * ms[i]) * vs[i]) / pow(ms[i], gradient);
    cout << "c: " << c << endl;
  }
  }

  double count = 0.0;
  double xyMin = 0.0;
  double xyMax = 1.0;
  double zMin = -0.025;
  double zMax = 1.025;
  double r = 0.025;// / 16.0;
  for (int i = 0; i < width; i++)
  {
    double x = xyMin + (xyMax - xyMin)*(((double)i) + 0.5) / (double)width;
    for (int k = 0; k < width; k++)
    {
      double z = zMin + (zMax - zMin)*(((double)k) + 0.5) / (double)width;
      Vector2d p(x, z);
      double d = distanceToDisks(p);
#if defined SURFACE
      if (abs(d) < r)
#elif defined POINTS
      if (d < r)
#else 
      if (d < 0.0)
#endif
      {
        count++;
        setPixel(out, Vector2i(i, k), 0);
      }
    }
  }
  double stepX = (xyMax - xyMin) / (double)width;
  double stepY = (zMax - zMin) / (double)width;
  double area = count * stepX*stepY;
  cout << "r: " << r << ", area: " << area << endl;

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"fordDisksPoints.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}
