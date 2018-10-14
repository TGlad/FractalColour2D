#include "stdafx.h"
#include "bmp.h"
#include <fstream>


static int width = 1024;
static int height = 1024;

void putpixel(vector<BYTE> &out, const Vector2i &pos, int shade)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= width)
    return;
  int ind = 3 * (pos[0] + width*pos[1]);
  out[ind + 0] = out[ind + 1] = out[ind + 2] = shade;
}

int _tmain(int argc, _TCHAR* argv[])
{
  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey

  for (int x = 0; x < width; x++)
  {
    for (int y = 0; y < height; y++)
    {
      bool set = false;
      Vector2d p;
      p[0] = ((double)(x - width / 2)) / ((double)width / 3.0);
      p[1] = ((double)(y - height / 2)) / ((double)height / 3.0);
      Vector2d n1(-1, -sqrt(3.0));
      Vector2d n2(1, -sqrt(3.0));
      for (int i = 0; i < 6; i++)
      {
        if (p.dot(n1) < 0.0 || p.dot(n2) < 0.0)
        {
          double rot = p[0] > 0.0 ? pi*2.0 / 3.0 : -pi*2.0 / 3.0;
          p = Vector2d(p[0]*cos(rot) + p[1]*sin(rot), -p[0]*sin(rot) + p[1]*cos(rot));
        }
        double w = 0.58;
        if (p[0] > w || p[0] < -w)
        {
          Vector2d dir(sqrt(3.0) / 2.0, -0.5);
          if (p[0] < -w)
            dir[0] = -dir[0];
          p += dir * (2 * sqrt(3.0) - 1.0);
          p /= 2.0*sqrt(3.0); // so it is at distance of 1
          p /= p.squaredNorm();
          double dot = p.dot(dir);
          p += dir * (1.0 - 2.0*dot);
          double scale = 4.0;
          p *= scale * 2.0*sqrt(3.0);
        }
        else
        {
          p[1] -= 1.0;
          double r = sqrt(3.0);
          double thickness = /*5.0 / (double)width; */  0.1 * (1.0 - sqr(p[0] / w));
          if (p.norm() > r - thickness && p.norm() < r + thickness)
            set = true;
          break;
        }
      }
      if (set)
        putpixel(out, Vector2i(x, y), 0);
    }
  }

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"bubbles_small.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}
