// experiments with structural coloration of fractal curves, in particular Koch, Levy, Dragon and random curves.
#include "stdafx.h"
#include "bmp.h"
#include "spectrumToRGB.h"
#include <set>

static int width = 4000;
static int height = 1600;
static int gap = 50;

void setpixel(vector<BYTE> &out, const Vector2i &pos, const Vector3d &colour)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= height)
    return;
  int ind = 3 * (pos[0] + width*(height - 1 - pos[1]));
  out[ind + 0] = max(0, min((int)(255.0*colour[0]), 255));
  out[ind + 1] = max(0, min((int)(255.0*colour[1]), 255));
  out[ind + 2] = max(0, min((int)(255.0*colour[2]), 255));
}

void drawLineX(vector<BYTE> &out, double y, const Vector3d &colour, bool thick)
{
  for (int x = 0; x < gap; x++)
  {
    setpixel(out, Vector2i(x, y-1), Vector3d(0, 0, 0));
    setpixel(out, Vector2i(x, y), Vector3d(0, 0, 0));
    setpixel(out, Vector2i(x, y+1), Vector3d(0, 0, 0));
  }
  for (int x = gap; x < width; x++)
  {
    setpixel(out, Vector2i(x, y), colour);
    if (thick)
    {
      setpixel(out, Vector2i(x, y-1), colour);
      setpixel(out, Vector2i(x, y+1), colour);
    }
  }
}
void drawLineY(vector<BYTE> &out, double x, const Vector3d &colour, bool thick)
{
  for (int y = 0; y < gap; y++)
  {
    setpixel(out, Vector2i(x-1, y), Vector3d(0, 0, 0));
    setpixel(out, Vector2i(x, y), Vector3d(0, 0, 0));
    setpixel(out, Vector2i(x+1, y), Vector3d(0, 0, 0));
  }
  for (int y = gap; y < height; y++)
  {
    setpixel(out, Vector2i(x, y), colour);
    if (thick)
    {
      setpixel(out, Vector2i(x-1, y), colour);
      setpixel(out, Vector2i(x+1, y), colour);
    }
  }
}
void drawLine(vector<BYTE> &out, double y0, double y1, const Vector3d &colour)
{
  for (int x = gap; x < width; x++)
  {
    double y = y0 + (y1 - y0) * (double)(x - gap) / (double)(width - gap);
    setpixel(out, Vector2i(x, y - 2.0), colour);
    setpixel(out, Vector2i(x, y - 1.0), colour);
    setpixel(out, Vector2i(x, y), colour);
    setpixel(out, Vector2i(x, y+1.0), colour);
    setpixel(out, Vector2i(x, y+2.0), colour);
  }
}
void drawTri(vector<BYTE> &out, double y0, double y1, double x, double y, const Vector3d &colour)
{
  double top = y0 + (y1 - y0) * (double)(x - gap) / (double)(width - gap);
  double left = gap + (y - y0)*(double)(width - gap) / (y1 - y0);
  for (int i = left; i < x; i++)
  {
    setpixel(out, Vector2i(i, y-1), colour);
    setpixel(out, Vector2i(i, y), colour);
    setpixel(out, Vector2i(i, y+1), colour);
  }
  for (int j = y; j < top; j++)
  {
    setpixel(out, Vector2i(x-1, j), colour);
    setpixel(out, Vector2i(x, j), colour);
    setpixel(out, Vector2i(x+1, j), colour);
  }
}

void drawPoint(vector<BYTE> &out, double x, double y, double radius, const Vector3d &colour)
{
  for (int i = x - radius; i <= x + radius; i++)
  {
    for (int j = y - radius; j <= y + radius; j++)
    {
      double r2 = sqr((double)i - x) + sqr((double)j - y);
      if (r2 <= sqr(radius))
        setpixel(out, Vector2i(i, j), colour);
    }
  }
}
 
int _tmain(int argc, _TCHAR* argv[])
{
  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey

  for (int x = gap; x < width; x++)
    for (int y = gap; y < height; y++)
      setpixel(out, Vector2i(x, y), Vector3d(0.6, 0.8, 1.0));
  int s = 1;
  int cols = 5;
  for (int j = 0; j < cols; j++)
  {
    for (int i = 1; i <= 10; i++)
    {
      double pos = log10((double)(i*s)) * (double)(width-gap) / (double)cols;
      double shade = 0.5;
      drawLineY(out, pos + (double)gap, Vector3d(shade, shade, shade), i==1);
    }
    s *= 10;
  }
  int rows = 2;
  s = 1;
  for (int j = 0; j < rows; j++)
  {
    for (int i = 1; i <= 10; i++)
    {
      double pos = log10((double)(i*s)) * (double)(height-gap) / (double)rows;
      double shade = 0.5;
      drawLineX(out, pos + (double)gap, Vector3d(shade, shade, shade), i==1);
    }
    s *= 10;
  }
  double start = height / 10;
  double end = 1.1*(double)height;
  drawLine(out, start, end, Vector3d(0, 0, 0.6));
  drawTri(out, start, end, width / 2, height / 2, Vector3d(0, 0, 0));
  double alongs[] = { 0.15, 0.3, 0.4, 0.54, 0.6, 0.75 };
  for (int i = 0; i < 6; i++)
  {
    double x = alongs[i] * (double)width;
    double y = start + (end - start) * (double)(x - gap) / (double)(width - gap);
    y += random(-0.05, 0.05) * (double)height;
    drawPoint(out, x, y, 15.0, Vector3d(0.9, 0, 0));
  }

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"graph.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}

