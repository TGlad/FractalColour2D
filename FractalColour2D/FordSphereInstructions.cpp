#include "stdafx.h"
#include "bmp.h"
#include <fstream>


static int width = 1024;
static int height = 1024;

void putpixel(vector<BYTE> &out, const Vector2i &pos, int shade)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= height)
    return;
  int ind = 3 * (pos[0] + width*(height - 1 - pos[1]));
  out[ind + 0] = out[ind + 1] = out[ind + 2] = shade;
}

void drawCircle(vector<BYTE> &out, const Vector2d &pos, double radius, int shade)
{
  for (int i = (int)(pos[0] - radius); i < (int)(pos[0] + radius); i++)
    for (int j = (int)(pos[1] - radius); j < (int)(pos[1] + radius); j++)
      if ((Vector2d(i, j) - pos).squaredNorm() < sqr(radius))
        putpixel(out, Vector2i(i, j), shade);
}
int _tmain(int argc, _TCHAR* argv[])
{
  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey
  Vector2d vs[3] = { Vector2d(0, 1), Vector2d(sqrt(3) / 2.0, -0.5), Vector2d(-sqrt(3) / 2.0, -0.5) };
#if 1 // inverse
  for (int xx = 0; xx < width; xx++)
  {
    for (int yy = 0; yy < width; yy++)
    {
      for (int i = -3; i <= 3; i++)
      {
        for (int j = -3; j <= 3; j++)
        {
          if (i == 0 && j == 0)
            continue;
          double x = (double)i*sqrt(3.0) / 2.0;
          double y = (double)j;
          if ((i + 100) % 2)
            y -= 0.5;
          if (Vector2d(x, y).squaredNorm() > sqr(3.1))
            continue;

          double s = 2.25;
          Vector2d p(xx, yy);
          p /= (double)width * 0.5;
          p -= Vector2d(1, 1);
          p *= s;
          p *= 0.75 / p.squaredNorm();
          if ((p - Vector2d(x, y)).squaredNorm() < sqr(0.5))
            putpixel(out, Vector2i(xx, yy), 196);
        }
      }
    }
  }
#else // forwards
  for (int i = -3; i <= 3; i++)
  {
    for (int j = -3; j <= 3; j++)
    {
      if (i == 0 && j == 0)
        continue;
      double x = (double)i*sqrt(3.0)/2.0;
      double y = (double)j;
      if ((i + 100) % 2)
        y -= 0.5;
      double s = 2.25;
      drawCircle(out, (double)width * 0.5 * ((Vector2d(x, y) / s) + Vector2d(1.0, 1.0)), (double)width * 0.5*0.5 / s, 196);
    }
  }
#endif

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"disksInv.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}
