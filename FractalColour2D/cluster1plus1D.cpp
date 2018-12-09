#include "stdafx.h"
#include "bmp.h"
#include <fstream>


static int width = 1024;
static int height = 1024;// 512;

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
      Vector2d p((double)x / (double)width * 3.0 - 1.5, (double)y / (double)height * 3.0 - 1.5);
      for (int it = 0; it < 4; it++)
      {
        if (p[0]>-0.5 && p[0]<0.5 && p[1]>-0.5 && p[1] < 0.5) // found it
        {
       //   if (abs(p[1]) < (0.4 - 2.0*sqr(p[0])))
          if (p.norm() < 0.4)
          {
            putpixel(out, Vector2i(x, y), 0);
            break;
          }
        }
        else
        {
          p[0] = abs(p[0]);
          p[1] = abs(p[1]);
          if (p[0] > 0.5 && p[1] > 0.5)
            p -= Vector2d(1.0, 1.0);
          else if (p[0] > 0.5)
            p -= Vector2d(1.0, 0.0);
          else
            p -= Vector2d(0.0, 1.0);
          p *= 3.0;
        }
      }
    }
  }

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"cluster2D.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}
