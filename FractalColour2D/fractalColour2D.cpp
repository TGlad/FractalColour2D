// experiments with structural coloration of fractal curves, in particular Koch, Levy, Dragon and random curves.
#include "stdafx.h"
#include "bmp.h"
#include "spectrumToRGB.h"
#include <set>

static int width = 729;
static int height = 729;
static double drawHeight = (double)height * 0.8;

void putpixel(vector<BYTE> &out, const Vector2i &pos, int shade)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= width)
    return;
  int ind = 3 * (pos[0] + width*pos[1]);
  out[ind + 0] = out[ind + 1] = out[ind + 2] = shade;
}
int getpixel(vector<BYTE> &out, const Vector2i &pos)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= width)
    return 0;
  int ind = 3 * (pos[0] + width*pos[1]);
  return out[ind + 0];
}
static double ccount = 0;
//#define BASIC_CANTOR2D
//#define DENSE
#if defined BASIC_CANTOR2D
void addPoint(vector<BYTE> &out, const Vector2i &bMin, int length)
{
  if (length == 1)
  {
    putpixel(out, bMin, 0);
    return;
  }
  int step = length / 3;
  for (int x = 0; x < 2; x++)
    for (int y = 0; y < 2; y++)
      addPoint(out, bMin + step*Vector2i(2 * x, 2 * y), step);
}
#elif defined DENSE
void addPoint(vector<BYTE> &out, const Vector2i &bMin, int length)
{
  if (length == 1)
  {
    putpixel(out, bMin, 0);
    return;
  }
  int step = length / 9;
  for (int xx = 0; xx < 2; xx++)
    for (int yy = 0; yy < 2; yy++)
      for (int x = 0; x < 2; x++)
        for (int y = 0; y < 2; y++)
          addPoint(out, bMin + step*Vector2i(7 * xx + x, 7 * yy + y), step);
}
static const int sequenceLength = 5;
static double scales[] = { 1.0, 2.0 / 9.0, 1.0 / 9.0, 4.0 / 81.0, 2.0 / 81.0 };
static double boxCounts[] = { 4.0, 16.0, 36.0, 64.0, 37.0*4.0 };
#else
void addPoint(vector<BYTE> &out, const Vector2i &bMin, int length)
{
  if (length == 1)
  {
    putpixel(out, bMin, 0);
    return;
  }
  int step = length / 27;
  for (int xx = 0; xx < 2; xx++)
    for (int yy = 0; yy < 2; yy++)
      for (int x = 0; x < 4; x++)
        for (int y = 0; y < 4; y++)
          addPoint(out, bMin + step*Vector2i(20 * xx + 2*x, 20 * yy + 2*y), step);
}
static const int sequenceLength = 5;
static double scales[] = { 1.0, 1.0 / 3.0, 1.0 / 9.0, 1.0/27.0, 1.0/81.0};
static double boxCounts[] = { 4.0, 16.0, 64.0, 64.0, 256.0 };
#endif

int _tmain(int argc, _TCHAR* argv[])
{
  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey

  addPoint(out, Vector2i(0, 0), width);
  int i = 0;
  for (int x = 0; x < width; x++)
  {
    double scale = pow(3.0, -3.0*(double)x / (double)width);
    if (scale < scales[i] && i < sequenceLength-1)
      i++;
    putpixel(out, Vector2i(x, height - 1 - boxCounts[i]), 127);
  }

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);

#if defined BASIC_CANTOR2D
  LPCTSTR file = L"cantor2Dgrowth_basic.bmp";
#elif defined DENSE
  LPCTSTR file = L"cantor2Dgrowth_dense.bmp";
#else
  LPCTSTR file = L"cantor2Dgrowth.bmp";
#endif
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}

