#include "stdafx.h"
#include "bmp.h"
#include <fstream>


static int width = 2048;
static int height = 2048;

void putpixel(vector<BYTE> &out, const Vector2i &pos, int shade)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= height)
    return;
  int ind = 3 * (pos[0] + width*(height - 1 - pos[1]));
  out[ind + 0] = out[ind + 1] = out[ind + 2] = shade;
}

int _tmain(int argc, _TCHAR* argv[])
{
  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey
  Vector2d vs[3] = { Vector2d(0, 1), Vector2d(sqrt(3) / 2.0, -0.5), Vector2d(-sqrt(3) / 2.0, -0.5) };
 
  /*
  Vector2d x = vs[1];// -vs[2] * (sqrt(3.0) - 1.0);
  double l = 2.0;// 2.0*sqrt(3.0);
  Vector2d o = -vs[1] * l;
  Vector2d pos = x - o;
  pos *= sqr(l + 1.0) / pos.squaredNorm();
  pos += o;
  double d = pos.dot(vs[1]) - 1.0;
  pos -= vs[1] * 2.0*d;
  pos[0] -= vs[1][0];
  */
  
  
  double scale = 4.0;
  for (int x = 0; x < width; x++)
  {
    for (int y = 0; y < height; y++)
    {
      Vector2d p;
      p[0] = ((double)(x - width / 2)) / ((double)width / 4.0);
      p[1] = ((double)(y - height / 2)) / ((double)height / 4.0);
      double dScale = 1.0;
      for (int i = 0; i < 8; i++)
      {
        double bendDist = 3.0;
        double l0 = (p - vs[0]).squaredNorm();
        double l1 = (p - vs[1]).squaredNorm();
        double l2 = (p - vs[2]).squaredNorm();
        int vi = 2;
        if (l0 < l1 && l0 < l2)
          vi = 0;
        else if (l1 < l0 && l1 < l2)
          vi = 1;
        p += vs[vi] * (bendDist - 1.0);
        double s = sqr(bendDist) / p.squaredNorm();
        p *= s;
        dScale *= s;
        p -= vs[vi] * (p.dot(vs[vi]) * 2.0 - bendDist);
        p *= scale;
        dScale *= scale;
      }
      double l0 = (p - vs[0]).squaredNorm();
      double l1 = (p - vs[1]).squaredNorm();
      double l2 = (p - vs[2]).squaredNorm();
      int vi = 2;
      if (l0 < l1 && l0 < l2)
        vi = 0;
      else if (l1 < l0 && l1 < l2)
        vi = 1;
      double t = p.dot(vs[vi]);
      double dist = (p - t*vs[vi]).norm() / dScale;
      if (dist < 0.004 )
        putpixel(out, Vector2i(x, y), 0);
    }
  }

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"bubble1.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}
