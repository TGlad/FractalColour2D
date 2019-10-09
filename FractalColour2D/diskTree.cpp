#include "stdafx.h"
#include "bmp.h"
#include <set>
#include <sstream>
#include <fstream>
#include "Core\EigenBase.h"
static int width = 4000;
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

double distanceToTree(Vector2d p)
{
  double scale = 1.0;
  double ratio = 2.0 / 3.0;
  for (int i = 0; i < 8; i++)
  {
    if ((p - Vector2d(0, 0.5)).norm() > 0.5)
      if (i == 0)
        break; // not in set
    if ((p - Vector2d(0, 1.0 / 3.0)).norm() < 1.0 / 3.0)
      break; // is in set
    double ss = p.squaredNorm();
    p /= ss;
    scale /= ss;

    if (p[1] > 1.0)
    {
      p = Vector2d(0, 3.0) - 2.0 * p;
      scale *= 2.0;
    }
    p[0] = fmod(100 + p[0] + 0.5, 1.0) - 0.5;
  }
//#define SURFACE
#define POINTS
#if defined SURFACE
  double d = ((p - Vector2d(0, 1.0 / 3.0)).norm() - 1.0/3.0) / scale; 
#elif defined POINTS
  double d = (p - Vector2d(0, 1.0 / 3.0)).norm() / scale; 
#else
  double d = ((p - Vector2d(0, 1.0 / 3.0)).norm() - 1.0/3.0) / scale; // the 0.4 is slightly more averaging than 0.5
#endif
  return d;
}

int _tmain(int argc, _TCHAR* argv[])
{
  double rs[] = { 0.025, 0.0125, 0.00625, 0.003125, 0.0015625, 0.00078125 };
//  double vs[] = { 0.34, 0.255, 0.192, 0.144, 0.109, 0.0818 };
  double vs[] = { 0.213, 0.160, 0.121, 0.0908, 0.0684, 0.0515 };
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
  LPCTSTR file = L"disktree_graph_points.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;

  int a = 0; // 0
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
  double xyMin = -0.6;
  double xyMax = 0.6;
  double zMin = -0.1;
  double zMax = 1.1;
  double r = 0.025 / 32.0;
  for (int i = 0; i < width; i++)
  {
    double x = xyMin + (xyMax - xyMin)*(((double)i) + 0.5) / (double)width;
    for (int j = 0; j < width; j++)
    {
      double y = zMin + (zMax - zMin)*(((double)j) + 0.5) / (double)width;
      Vector2d p(x, y);
      double d = distanceToTree(p);
#if defined SURFACE
      if (abs(d) < r)
#elif defined POINTS
      if (d < r)
#else 
      if (d < 0.0)
#endif
      {
        count++;
        setPixel(out, Vector2i(i, j), 0);
      }
    }
  }
  double step = (xyMax - xyMin) / (double)width;
  double area = count * step*step;
  cout << "r: " << r << ", area: " << area << endl;

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"disktreepoints.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}
