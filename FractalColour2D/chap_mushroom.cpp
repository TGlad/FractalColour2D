// experiments with structural coloration of fractal curves, in particular Koch, Levy, Dragon and random curves.
#include "stdafx.h"
#include "bmp.h"
#include "spectrumToRGB.h"
#include <set>

static int width = 4096;
static int height = 1664;

void setpixel(vector<BYTE> &out, const Vector2i &pos, const Vector3d &colour)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= height)
    return;
  int ind = 3 * (pos[0] + width*(height - 1 - pos[1]));
  out[ind + 0] = max(0, min((int)(255.0*colour[0]), 255));
  out[ind + 1] = max(0, min((int)(255.0*colour[1]), 255));
  out[ind + 2] = max(0, min((int)(255.0*colour[2]), 255));
}

void drawBox(vector<BYTE> &out, const Vector2d &start, const Vector2d &end, const Vector3d &colour)
{
  for (int x = start[0]; x <= end[0]; x++)
  {
    for (int y = start[1]; y <= end[1]; y++)
    {
      setpixel(out, Vector2i(x, y), colour);
    }
  }
}
 
void recurse(vector<BYTE> &out, const Vector2d &start, const Vector2d &end)
{
  double mag = (end - start).norm();
  if (mag <= 1.0)
  {
    Vector2d bmin(min(start[0], end[0]), min(start[1], end[1]));
    Vector2d bmax(max(start[0], end[0]), max(start[1], end[1]));
    double w = 0.5;
    drawBox(out, bmax - Vector2d(w, w), bmax + Vector2d(w, w), Vector3d(0.4, 0.0, 0.0));
    return;
  }
  double length = mag / 4.0;
  Vector2d along = (end - start) / 4.0;
  Vector2d up(-along[1], along[0]);
  recurse(out, start, start + along);
  recurse(out, start + along + up, start + along);
  recurse(out, start + along + up, start + along + 2.0*up);
  recurse(out, start + along + 2.0*up, start + 2.0*along + 2.0*up);
  recurse(out, start + 2.0*along + 2.0*up, start + 3.0*along + 2.0*up);
  recurse(out, start + 3.0*along + 2.0*up, start + 3.0*along + up);
  recurse(out, start + 3.0*along, start + 3.0*along + up);
  recurse(out, start + 3.0*along, start + 4.0*along);
}
int _tmain(int argc, _TCHAR* argv[])
{
  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey

  for (int x = 0; x < width; x++)
    for (int y = 0; y < height; y++)
      setpixel(out, Vector2i(x, y), Vector3d(0.7, 0.8, 0.9));

  double length = (double)height*0.97;
  recurse(out, Vector2d(length*0.75, (double)height - length), Vector2d(length * 0.75, (double)height));
  recurse(out, Vector2d(length*0.75, (double)height - length), Vector2d(length*1.75, (double)height - length));
  recurse(out, Vector2d(length*1.75, (double)height), Vector2d(length*1.75, (double)height - length));

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"chap_mushroom.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}

