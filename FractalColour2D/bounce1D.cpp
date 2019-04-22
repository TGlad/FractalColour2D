#include "stdafx.h"
#include "bmp.h"
#include <sstream>
#include <fstream>
static int width = 1024;
static int height = 512;

void setPixel(vector<BYTE> &out, int x, int y, const Vector3i &col)
{
  if (x < 0 || x >= width || y < 0 || y >= height)
    return;
  int ind = 3 * (x + width * (height - 1 - y));
  out[ind + 0] = col[0];
  out[ind + 1] = col[1];
  out[ind + 2] = col[2];
}

void drawThing(vector<BYTE> &out, double scale, double h, Vector3i &col)
{
  for (int y = (int)(h*(double)height); y < (int)((h + scale)*(double)height); y++)
  {
    double Y = 1.0 - abs((double)y - (h + 0.5*scale)*(double)height) / (0.5*scale*(double)height);
    Y = max(3.0 / 6.0, (1.0 / 6.0) + 4.0*Y / 6.0);
    double startX = (double)width * scale * Y;
    for (int x = startX; x < startX + (double)width*scale*1.0 / 6.0; x++)
      setPixel(out, x, y, col);
  }
}

void main()
{
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey
  long s2;
  int i = 0;
  for (double scale = 1.0; scale > 0.001; scale /= 2.0)
  {
    Vector3i col = i == 0 ? Vector3i(20, 60, 20) : Vector3i(20, 20, 60);
    i = 1 - i;
    for (double y = -0.5; y < 1.0; y+=scale)
      drawThing(out, scale, y, col);
  }

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"bounce1D.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}
