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
      if (abs(Y)<0.06 && X < 0.0 && fmod(-X, 0.2)<0.1)
        putpixel(out, Vector2i(x, y), Vector3d(0, 0, 0));
      if (X < 0.0 && Y < -X*0.08 && Y>0.0)
        scalepixel(out, Vector2i(x, y), 0.9);
    }
  }
  {
    BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    LPCTSTR file = L"square1.bmp";
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
      double mag = abs(z);
      double angle = arg(z);
      bool swappy = angle < 0.0;
      double power = 1.0 / 1.5;
      angle *= power;
      mag = pow(mag, power);
      Complex zs[2];
      zs[0] = mag*cos(angle) + I*mag*sin(angle);
      int top = 1;
      angle = arg(z);
      if (abs(angle) > pi / 2.0)
      {
        if (angle < 0.0)
          angle += 2.0*pi;
        else
          angle -= 2.0*pi;
        angle *= power;
        zs[1] = mag*cos(angle) + I*mag*sin(angle);
        top = 2;
        if (swappy)
          swap(zs[0], zs[1]);
      }
      if (X < 0.0)
        putpixel(out, Vector2i(x, y), Vector3d(0.99, 0.99, 0.99));
      else if (Y>0.0 && X < 0.14*Y)
        putpixel(out, Vector2i(x, y), Vector3d(0.9, 0.9, 0.9));
      //  z = pow(z, 1.0/1.5);
      for (int i = 0; i < top; i++)
      {
        X = zs[i].real() * res / 2.0;
        Y = zs[i].imag() * res / 2.0;
        col[1] = max(0.0, 0.5 + 0.5*zs[i].real() / 2.0);
        col[2] = max(0.0, 0.5 + 0.5*zs[i].imag() / 2.0);
        col *= 0.9;
        Vector3d cols[2] = { col*(1.0 - pale) + Vector3d(1, 1, 1)*pale, col };

        double xx = fmod(X + 10.0*res + 0.5, 1.0);
        double xx2 = fmod(X + 10.0*res + 0.5, 5.0);
        if (abs(xx2 - 0.5) < 0.05 || abs(xx - 0.5) < 0.02)
          putpixel(out, Vector2i(x, y), top==2 ? cols[i] : cols[1]);
        double yy = fmod(Y + 10.0*res + 0.5, 1.0);
        double yy2 = fmod(Y + 10.0*res + 0.5, 5.0);
        if (abs(yy2 - 0.5) < 0.05 || abs(yy - 0.5) < 0.02)
          putpixel(out, Vector2i(x, y), top == 2 ? cols[i] : cols[1]);
        if (abs(Y)<0.06 && X < 0.0 && fmod(-X, 0.2)<0.1)
          putpixel(out, Vector2i(x, y), top == 2 && i == 0 ? Vector3d(pale, pale, pale) : Vector3d(0.4, 0.4, 0.4));
      }
    }
  }
  {
    BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    LPCTSTR file = L"square2.bmp";
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
      double mag = abs(z);
      double angle = arg(z);
      double off = sin(angle/2.0) * mag * 0.01;
      double power = 0.5;
      angle *= power;
      mag = pow(mag, power);
      Complex zs[2];
      zs[0] = mag*cos(angle) + I*mag*sin(angle);
      zs[0] += I*off;
      int top = 1;
      int k = 0;
      angle = arg(z);
      if (angle < 0.0)
        k = 1;
      {
        if (angle < 0.0)
          angle += 2.0*pi;
        else
          angle -= 2.0*pi;
        off = sin(angle / 2.0) * mag * 0.01;
        angle *= power;
        zs[1] = mag*cos(angle) + I*mag*sin(angle);
        zs[1] += off;
        top = 2;
      }
      for (int i = 0; i < top; i++)
      {
        X = zs[i].real() * res / 2.0;
        Y = zs[i].imag() * res / 2.0;
        col[1] = max(0.0, 0.5 + 0.5*zs[(i+k)%2].real() / 2.0);
        col[2] = max(0.0, 0.5 + 0.5*zs[(i + k) % 2].imag() / 2.0);
        Vector3d cols[2] = { col*(1.0 - pale) + Vector3d(1, 1, 1)*pale, col };
        double xx = fmod(X + 10.0*res + 0.5, 1.0);
        double xx2 = fmod(X + 10.0*res + 0.5, 5.0);
        if (abs(xx2 - 0.5) < 0.05 || abs(xx - 0.5) < 0.02)
          putpixel(out, Vector2i(x, y), cols[i]);
        double yy = fmod(Y + 10.0*res + 0.5, 1.0);
        double yy2 = fmod(Y + 10.0*res + 0.5, 5.0);
        if (abs(yy2 - 0.5) < 0.05 || abs(yy - 0.5) < 0.02)
          putpixel(out, Vector2i(x, y), cols[i]);
        if (abs(Y)<0.03 && X < 0.0 && fmod(-X, 0.2)<0.1)
          putpixel(out, Vector2i(x, y), top == 2 && i == 0 ? Vector3d(0.5, 0.5, 0.5) : Vector3d(0, 0, 0));
      }
    }
  }
  {
    BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    LPCTSTR file = L"square3.bmp";
    SaveBMP(c, width, height, s2, file);
    delete[] c;
  }
#endif
}
