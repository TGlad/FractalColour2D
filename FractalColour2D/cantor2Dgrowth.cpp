// experiments with structural coloration of fractal curves, in particular Koch, Levy, Dragon and random curves.
#include "stdafx.h"
#include "bmp.h"
#include "spectrumToRGB.h"
#include <set>

static int width = 729;
static int height = 729 +160;
void putpixel(vector<BYTE> &out, const Vector2i &pos, int shade)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= height)
    return;
  int ind = 3 * (pos[0] + width*pos[1]);
  out[ind + 0] = out[ind + 1] = out[ind + 2] = shade;
}
void putpixel(vector<BYTE> &out, const Vector2d &pos, int shade)
{
  putpixel(out, Vector2i(pos[0], pos[1]), shade);
}
int getpixel(vector<BYTE> &out, const Vector2i &pos)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= height)
    return 0;
  int ind = 3 * (pos[0] + width*pos[1]);
  return out[ind + 0];
}
static double ccount = 0;
#define BASIC_CANTOR2D
//#define TRIPLE
//#define FUNNY
#if defined BASIC_CANTOR2D
void addPoint(vector<BYTE> &out, const Vector2d &bMin, double length)
{
  if (length <= 1 + 1e-10)
  {
    putpixel(out, bMin, 0);
    return;
  }
  double step = length / 3.0;
  for (int x = 0; x < 2; x++)
    for (int y = 0; y < 2; y++)
      addPoint(out, bMin + step*Vector2d(2 * x, 2 * y), step);
}
static const int sequenceLength = 5;
static double scales[] = { 1.0, 3.0, 9.0, 27.0, 81.0 };
static double boxCounts[] = { 4.0, 16.0, 64.0, 256.0, 4.0*256.0 };
#elif defined TRIPLE
void addPoint(vector<BYTE> &out, const Vector2d &bMin, double length)
{
  if (length <= 1 + 1e-10)
  {
    putpixel(out, bMin, 0);
    return;
  }
  double step = length * 0.5* 4.7 / 5.7;

  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++)
      addPoint(out, bMin + step*Vector2d(x, y), length / 5.7);
}
static const int sequenceLength = 6;
static double scales[] = { 1.0, 2.0, 5.7, 5.7*2.0, 5.7*5.7, 5.7*5.7 * 2.0 };
static double boxCounts[] = { 4.0, 9.0, 36.0, 81.0, 81.0*4.0, 81.0*9.0 };
#elif defined FUNNY
void addPoint(vector<BYTE> &out, const Vector2d &bMin, double length)
{
  if (length <= 1 + 1e-10)
  {
    putpixel(out, bMin, 0);
    return;
  }
  double step = length * 1.0/9.0;

  for (int x = 0; x < 2; x++)
  {
    for (int y = 0; y < 2; y++)
    {
      addPoint(out, bMin + step*Vector2d(5 * x, 5 * y), length / 9.0);
      addPoint(out, bMin + step*Vector2d(5 * x + 3, 5 * y), length / 9.0);
      addPoint(out, bMin + step*Vector2d(5 * x, 5 * y + 3), length / 9.0);
      addPoint(out, bMin + step*Vector2d(5 * x + 3, 5 * y + 3), length / 9.0);
    }
  }
}
static const int sequenceLength = 7;
static double scales[] = { 1.0, 9.0 / 4.0, 3.0, 9.0, 9.0*9.0 / 4.0, 27.0, 81.0 };
static double boxCounts[] = { 4.0, 9.0, 16.0, 4.0*16.0, 9.0*16.0, 16.0*16.0, 4.0*16.0*16.0 };

#else
void addPoint(vector<BYTE> &out, const Vector2d &bMin, double length)
{
  if (length <= 1 + 1e-10)
  {
    putpixel(out, bMin, 0);
    return;
  }
  double step = length * 8 / 27;

  for (int x = 0; x < 4; x++)
    for (int y = 0; y < 4; y++)
      addPoint(out, bMin + step*Vector2d(x, y), length / 9.0);
}
static const int sequenceLength = 7;
static double scales[] = { 1.0, 27.0 / 11.0, 3.0, 9.0, 9.0*27.0 / 11.0, 27.0, 81.0 };
static double boxCounts[] = { 4.0, 9.0, 16.0, 4.0*16.0, 9.0*16.0, 16.0*16.0, 4.0*16.0*16.0 };
#endif

int _tmain(int argc, _TCHAR* argv[])
{
  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey

  addPoint(out, Vector2d(0, 0), width);
/*  int i = -1;
  double lastY = -1;
  for (int x = 0; x < width; x++)
  {
    double scale = pow(3.0, 4.0*(double)(x-20) / (double)width);
    double y = 1.0;
    if (scale > scales[i+1] && i+1 < sequenceLength-1)
      i++;
    if (i >= 0)
      y = boxCounts[i];
    if (y != lastY)
      cout << "x: " << x << ", y: " << y << endl;
    lastY = y;
    putpixel(out, Vector2i(x, height - 1 - 0.5*y), 127);
  }*/

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);

#if defined BASIC_CANTOR2D
  LPCTSTR file = L"cantor2D_basic2.bmp";
#elif defined TRIPLE
  LPCTSTR file = L"cantor2Dgrowth_triple.bmp";
#elif defined FUNNY
  LPCTSTR file = L"cantor2Dgrowth_funny.bmp";
#else
  LPCTSTR file = L"cantor2Dgrowth_dense.bmp";
#endif
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}

