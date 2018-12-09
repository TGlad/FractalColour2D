#include "stdafx.h"
#include "bmp.h"
#include <fstream>


static int width = 2048;
static int height = 2048;

void putpixel(vector<BYTE> &out, const Vector2i &pos, int shade)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= height)
    return;
  int ind = 3 * (pos[0] + width*(height-pos[1]));
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
      Vector2d p;
      p[0] = ((double)(x - width / 2)) / ((double)width);
      p[1] = ((double)y) / (double)height;
      double dScale = 1.0;
      double ratio = 2.0/3.0;
      for (int i = 0; i < 2; i++)
      {
        if ((p - Vector2d(0, 0.5)).norm() > 0.5)
          if (i == 0 || p[1] > 0.5)
            break; // not in set
        if ((p - Vector2d(0, 0.5*ratio)).norm() < 0.5*ratio)
          break; // is in set
        if ((p - Vector2d(0, 0.5)).norm() < 0.5)
        {
          // otherwise we need to transform the thing
          p /= p.squaredNorm();
          p[1] -= 1.0;
          double s = 1.0 / ((1.0 / ratio) - 1.0);
          p *= s;
          p[1] = 1.0 - p[1];
        }
        else 
        {
          p /= p.squaredNorm();
//          p[0] = abs(p[0]);
//          p[0] -= 1.0;
//          p /= p.squaredNorm();
//          p[0] += 1.0;
        }
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
}
