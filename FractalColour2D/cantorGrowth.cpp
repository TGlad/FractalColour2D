// experiments with structural coloration of fractal curves, in particular Koch, Levy, Dragon and random curves.
#include "stdafx.h"
#include "bmp.h"
#include "spectrumToRGB.h"
#include <set>

static int width = 1024;
static int height = 1024;
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
//#define POSITIVE_DIMENSIONAL
void addPoint(vector<BYTE> &out, int x, double yMin, double yMax)
{
#if !defined POSITIVE_DIMENSIONAL
  if (yMax < -drawHeight / 2.0)
    return;
  if (yMin > drawHeight / 2.0)
    return;
#endif
  if (yMax - yMin < 1.0)
  {
    if (yMax > -drawHeight / 2.0)
    {
      putpixel(out, Vector2i(x, (double)drawHeight * 0.5 - yMin), 0);
      putpixel(out, Vector2i(x, (double)drawHeight * 0.5 - yMax), 0);
    }
    ccount++;
  }
  else
  {
    addPoint(out, x, yMin, yMin + (yMax - yMin)*1.0 / 3.0);
    addPoint(out, x, yMin + (yMax - yMin)*2.0 / 3.0, yMax);
  }
}
int _tmain(int argc, _TCHAR* argv[])
{
  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey

  for (int x = 0; x < width; x++)
  {
#if defined POSITIVE_DIMENSIONAL
    double scale = pow(2.0, (double)x / 220.0);
#else
    double scale = pow(2.0, (double)x / 220.0);
#endif
    ccount = 0;
#if defined POSITIVE_DIMENSIONAL
    addPoint(out, x, -20.0*scale, 20.0*scale);
#else
    addPoint(out, x, -scale*drawHeight/4.0, scale*3.0*drawHeight/4.0);
#endif

    // now count up the pixels and plot the line:
    putpixel(out, Vector2i(x, height - 1 - 1.5 * ccount), 127);
    putpixel(out, Vector2i(x, height - 2 - 1.5 * ccount), 127);
  }


  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
#if defined POSITIVE_DIMENSIONAL
  LPCTSTR file = L"cantorgrowthPositive.bmp";
#else
  LPCTSTR file = L"cantorgrowthSigned.bmp";
#endif

  SaveBMP(c, width, height, s2, file);
  delete[] c;
}

