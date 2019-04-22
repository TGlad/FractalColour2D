// experiments with structural coloration of fractal curves, in particular Koch, Levy, Dragon and random curves.
#include "stdafx.h"
#include "bmp.h"
#include "spectrumToRGB.h"
#include <set>

static int width = 1000;
static int height = 270;

void putpixel(vector<BYTE> &image, const Vector2i &pos, const Vector3d &col)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= height)
    return;
  int ind = 3 * (pos[0] + width*pos[1]);
  for (int i = 0; i < 3; i++)
    image[ind + i] = (BYTE)(255.0*col[i]);
}
void drawLine(vector<BYTE> &out, double x, double y, double x2, double y2, const Vector3d &col)
{
  if (abs(x2 - x) > abs(y2 - y))
  {
    if (x > x2)
    {
      swap(x, x2);
      swap(y, y2);
    }
    for (int i = x; i < x2; i++)
    {
      double j = y + (y2 - y)*((double)(i - x) / (x2 - x));
      putpixel(out, Vector2i(i, j), col);
      putpixel(out, Vector2i(i, j + 1), col);
    }
  }
  else
  {
    if (y > y2)
    {
      swap(x, x2);
      swap(y, y2);
    }
    for (int j = y; j < y2; j++)
    {
      double i = x + (x2 - x)*((double)(j - y) / (y2 - y));
      putpixel(out, Vector2i(i, j), col);
      putpixel(out, Vector2i(i+1, j), col);
    }
  }
}
int _tmain(int argc, _TCHAR* argv[])
{
  long s2;
  vector<BYTE> out(3*width*height); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); 
  // right, not keep an array of radii
  vector<double> rs(width);
  vector<double> ins(width);
  for (auto &r : rs)
    r = 0.0;
  for (auto &i : ins)
    i = 0.0;
  double tones[] = { 1.0, 16.0 / 15.0, 9.0 / 8.0, 6.0 / 5.0, 5.0 / 4.0, 4.0 / 3.0, 7.0 / 5.0, 3.0 / 2.0, 8.0 / 5.0, 5.0 / 3.0, /*9.0 / 5.0*/16.0/9.0, 15.0 / 8.0, 2.0 };
  for (int i = 0; i < 13; i++)
  {
    double x = (double)(width - 1) * (tones[i]/2.0);
    ins[x] = 1.0;
  }

  for (int numerator = 0; numerator <= 16; numerator++)
  {
    for (int denominator = numerator/2; denominator <= 15; denominator++)
    {
      if (denominator == 0)
        continue;
      double radius = 1.0 / sqr(denominator);
      double location = (double)(width-1) * 0.5 * (double)numerator / (double)denominator;
      if (location >= 0 && location < width)
        rs[location] = max(rs[location], radius);
    }
  }

  for (int i = 0; i < width; i++)
  {
    if (rs[i] <= 0.0)
      continue;
    double v = rs[i] * 64.0;
    double s = 1.0 - min(v, 1.0);
    double s2 = ins[i] ? 1.0 : s;

    for (int j = 0; j < 6; j++)
      for (int xx = i - 1; xx <= i + 1; xx++)
        putpixel(out, Vector2i(xx, 0 + j), Vector3d(s2, s, s));
  }


  // now try a stronger one with widths:
  vector<double> avoid(width);
  for (auto &v : avoid)
    v = 0.0;
  for (int i = 0; i < width; i++)
  {
    if (rs[i] < 0.002)
      continue;
    double w = rs[i] * 30.0;
    for (int x = (int)max(0.0, (double)i - w); x <= (int)min((double)width - 1, (double)i + w); x++)
      avoid[x] = 1.0;
  }
  // now try the opposite, staying away from the strong ones
  for (int y = -20; y < 20; y++)
  {
    for (int i = width/2; i < width; i++)
    {
      double rad = width/2;
      double x = (sqrt(sqr(rad) - sqr(8.0 * y)) - rad);
      if (avoid[i] == 0.0)
        putpixel(out, Vector2i(1 + i + (int)x, height - 21 + y), Vector3d(168, 139, 53)/255.0);
    }
  }

  // now draw some diagonal lines
  double vs[6] = { 1.0, 1 + 1.0/3.0, 3.0 / 2.0, 1.0 + 2.0/3.0, 2.0 - 2.0/(double)width, (1 + sqrt(5.0)) / 2.0 };
  for (int i = 0; i < 6; i++)
  {
    double w = (double)(width / 2);
    if (i == 5)
      drawLine(out, -1 + w*vs[i], height - 51, -1 + w*vs[i], height - 21, Vector3d(0.0, 0.0, 1.0));
    else
      drawLine(out, -1 + w*vs[i], 6, -1 + w*vs[i], height - 21, Vector3d(0.85, 0.85, 0.85));
  }
  // now try some music:
  double semitone = pow(2.0, 1.0 / 12.0);
  for (double d = 1.0; d <= 2.0001; d *= semitone)
  {
    double x = d * (double)(width - 1) / 2.0;
    double s = 0;
    for (int j = 0; j < 6; j++)
    {
      putpixel(out, Vector2i(x, 13 + j), Vector3d(s, s, s));
      putpixel(out, Vector2i(x - 1, 13 + j), Vector3d(s, s, s));
      putpixel(out, Vector2i(x + 1, 13 + j), Vector3d(s, s, s));
    }
  }

  double vals[8] = { 1.0 / 1.0, 2.0 / 1.0, 3.0 / 2.0, 5.0 / 3.0, 8.0 / 5.0, 13.0 / 8.0, 21.0 / 13.0, (1 + sqrt(5.0)) / 2.0 };
  double y = 70.0;
  for (int i = 0; i < 6; i++)
  {
    double w = (double)(width / 2);
    double newY = y + 20;
    drawLine(out, -1 + w*vals[i], y, -1 + w*vals[i + 1], newY, Vector3d(0,0,1));
    y = newY;
  }

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"conclusion.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}

