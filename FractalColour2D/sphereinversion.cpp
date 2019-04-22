#include "stdafx.h"
#include "bmp.h"
#include <sstream>
#include <fstream>
#include <complex>
typedef complex<double> Complex;
static const Complex I(0, 1);
static int width = 1000;
static int height = 1000;

void putpixel(vector<BYTE> &out, const Vector2i &pos, Vector3d &col)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= height)
    return;
  int ind = 3 * (pos[0] + width * pos[1]);
  out[ind + 0] = (int)(255.0*min(col[0], 1.0));
  out[ind + 1] = (int)(255.0*min(col[1], 1.0));
  out[ind + 2] = (int)(255.0*min(col[2], 1.0));
}

void scalepixel(vector<BYTE> &out, const Vector2i &pos, double scale)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= height)
    return;
  int ind = 3 * (pos[0] + width * pos[1]);
  out[ind + 0] = (int)((double)out[ind + 0] * scale);
  out[ind + 1] = (int)((double)out[ind + 1] * scale);
  out[ind + 2] = (int)((double)out[ind + 2] * scale);
}

int _tmain(int argc, _TCHAR* argv[])
{
  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  double res = 10.0;
#if 1
  memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey
  for (int x = 0; x < width; x++)
  {
    for (int y = 0; y < height; y++)
    {
      double X = res*(-1.0 + 2.0*(double)x / (double)width);
      double Y = res*(-1.0 + 2.0*(double)y / (double)height);
      Vector3d col;
      col[0] = 0.5;
      Complex z(X*2.0/res, Y*2.0/res);
      col[1] = max(0.0, 0.5 + 0.5*z.real() / 2.0);
      col[2] = max(0.0, 0.5 + 0.5*z.imag() / 2.0);
      col *= 0.85;
      double xx = fmod(X + 10.0*res + 0.5, 1.0);
      double xx2 = fmod(X + 10.0*res + 0.5, 5.0);        
      if (abs(xx2 - 0.5) < 0.05 || abs(xx - 0.5) < 0.02)
        putpixel(out, Vector2i(x, y), col);
      double yy = fmod(Y + 10.0*res + 0.5, 1.0);
      double yy2 = fmod(Y + 10.0*res + 0.5, 5.0);
      if (abs(yy2 - 0.5) < 0.05 || abs(yy - 0.5) < 0.02)
        putpixel(out, Vector2i(x, y), col);
    }
  }
  {
    BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    LPCTSTR file = L"sphereinversion1.bmp";
    SaveBMP(c, width, height, s2, file);
    delete[] c;
  }
#endif
#if 1
  memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey
  for (int x = 0; x < width; x++)
  {
    for (int y = 0; y < height; y++)
    {
      Vector3d col;
      col[0] = 0.5;
      double pale = 0.8;
      double X = res*(-1.0 + 2.0*(double)x / (double)width);
      double Y = res*(-1.0 + 2.0*(double)y / (double)height);
      Complex z(X*2.0 / res, Y*2.0 / res);
      z = conj(1.0 / z);
      {
        X = z.real() * res / 2.0;
        Y = z.imag() * res / 2.0;
        col[1] = max(0.0, 0.5 + 0.5*z.real() / 2.0);
        col[2] = max(0.0, 0.5 + 0.5*z.imag() / 2.0);
        double xx = fmod(X + 10.0*res + 0.5, 1.0);
        double xx2 = fmod(X + 10.0*res + 0.5, 5.0);
        if (abs(xx2 - 0.5) < 0.05 || abs(xx - 0.5) < 0.02)
          putpixel(out, Vector2i(x, y), col);
        double yy = fmod(Y + 10.0*res + 0.5, 1.0);
        double yy2 = fmod(Y + 10.0*res + 0.5, 5.0);
        if (abs(yy2 - 0.5) < 0.05 || abs(yy - 0.5) < 0.02)
          putpixel(out, Vector2i(x, y), col);
      }
    }
  }
  {
    BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    LPCTSTR file = L"sphereinversion2.bmp";
    SaveBMP(c, width, height, s2, file);
    delete[] c;
  }
#endif
}
