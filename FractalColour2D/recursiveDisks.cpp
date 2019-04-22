#include "stdafx.h"
#include "bmp.h"
#include <fstream>


static int width = 1024;
static int height = 1024;

void putpixel(vector<BYTE> &out, const Vector2i &pos, int shade)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= height)
    return;
  int ind = 3 * (pos[0] + width*(height-pos[1]));
  out[ind + 0] = out[ind + 1] = out[ind + 2] = shade;
}
static int w = 120;
static int h = 120;
void putpixel(vector<BYTE> &out, const Vector2i &pos, int red, int green, int blue)
{
  if (pos[0] < 0 || pos[0] >= w || pos[1] < 0 || pos[1] >= h)
    return;
  int ind = 3 * (pos[0] + w*(h - 1 - pos[1]));
  out[ind + 0] = red;
  out[ind + 1] = green;
  out[ind + 2] = blue;
}
void drawPoly(const vector<double> &poly, LPCTSTR file)
{
  long s2;
  vector<BYTE> out(w*h * 3); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey

  for (int i = 0; i < h; i++)
  {
    putpixel(out, Vector2i(w / 2 - 1, i), 0, 0, 255);
    putpixel(out, Vector2i(w / 2, i), 0, 0, 255);
    putpixel(out, Vector2i(w / 2 + 1, i), 0, 0, 255);
  }
  for (int i = 0; i < w; i++)
  {
    double x = -1.0 + 2.0*(double)i / (double)(w - 1);
    double y = 0.0;
    double xn = x;
    for (int j = 0; j < poly.size(); j++)
    {
      y += poly[j] * xn;
      xn *= x;
    }
    putpixel(out, Vector2i(i, h / 2 + 1), 0, 0, 255);
    putpixel(out, Vector2i(i, h / 2), 0, 0, 255);
    putpixel(out, Vector2i(i, h / 2 - 1), 0, 0, 255);
    int Y = (double)h / 2 + y*(double)(h - 1) / 2.0;
    putpixel(out, Vector2i(i, Y + 1), 0, 255, 0);
    putpixel(out, Vector2i(i, Y), 0, 255, 0);
    putpixel(out, Vector2i(i, Y - 1), 0, 255, 0);
  }
  BYTE* c = ConvertRGBToBMPBuffer(&out[0], w, h, &s2);
  SaveBMP(c, w, h, s2, file);
  delete[] c;
}

int _tmain(int argc, _TCHAR* argv[])
{
#define DISKS
#if defined DISKS
  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey
   
  for (int x = 0; x < width; x++)
  {
    for (int y = 0; y < height; y++)
    {
      Vector2d p;
      p[0] = ((double)(x - width / 2)) / ((double)width);
      p[1] = ((double)y) / (double)height;
      double dScale = 1.0;
      double ratio = 2.0/3.0;
      for (int i = 0; i < 4; i++)
      {
  //      if ((p - Vector2d(0, 0.5)).norm() > 0.5)
  //        if (i == 0)
  //          break; // not in set
        if ((p - Vector2d(0, 1.0/3.0)).norm() < 1.0/3.0)
          break; // is in set
        p /= p.squaredNorm();
        if (p[1] > 1.0)
          p = Vector2d(0, 3.0) - 2.0 * p;
        p[0] = fmod(100 + p[0] + 0.5, 1.0) - 0.5;
      }
      if ((p - Vector2d(0, 0.5)).norm() < 0.5)
        putpixel(out, Vector2i(x, y), 0);
    }
  }

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"recursiveFordDisks2.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
#else // make little graphs
  vector<double> polynomial;
  {
    polynomial.clear();
    polynomial.push_back(0.0);
    LPCTSTR file = L"poly0.bmp";
    drawPoly(polynomial, file);
  }
  {
    polynomial.clear();
    polynomial.push_back(1.0);
    LPCTSTR file = L"poly1.bmp";
    drawPoly(polynomial, file);
  }
  {
    polynomial.clear();
    polynomial.push_back(1.0);
    polynomial.push_back(-1.0);
    LPCTSTR file = L"poly2.bmp";
    drawPoly(polynomial, file);
  }
  {
    polynomial.clear();
    polynomial.push_back(-1.0/2.0);
    polynomial.push_back(2.0);
    LPCTSTR file = L"poly3.bmp";
    drawPoly(polynomial, file);
  }
  {
    polynomial.clear();
    polynomial.push_back(1.0);
    polynomial.push_back(-1.0);
    polynomial.push_back(2.0);
    LPCTSTR file = L"poly4.bmp";
    drawPoly(polynomial, file);
  }
#endif
}
