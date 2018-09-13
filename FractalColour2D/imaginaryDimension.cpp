// experiments with structural coloration of fractal curves, in particular Koch, Levy, Dragon and random curves.
#include "stdafx.h"
#include "bmp.h"
#include "spectrumToRGB.h"

static const int width = 256*3;
static const int height = width;
static const double factor = 3;
static vector<Vector2d> ps;
static const double widthThreshold = 1.0 / factor*factor; // 1 pixel after zooming by factor^2
static double minWidth = 0.0;

void recurse(const Vector2d &pos, double scale)
{
  if (scale <= widthThreshold)
  {
    ps.push_back(pos);
    minWidth = scale;
    return;
  }
  scale /= factor*factor;
  for (int i = 0; i < factor; i++)
  {
    for (int j = 0; j < factor; j++)
    {
      recurse(pos + scale*(Vector2d(1,1) + 3.0*Vector2d(i, j)), scale);
    }
  }
}

int _tmain(int argc, _TCHAR* argv[])
{
  long s2;
  recurse(Vector2d(0, 0), width);

  const int extraHeight = 40*3;
  const int fullHeight = height + extraHeight;
  const int maxFrames = 6;// 0;
  vector<double> pixelCount(maxFrames+1);
  for (int run = 0; run < 2; run++)
  {
    for (int i = 0; i < maxFrames+1; i++)
    {
      if (run == 0)
        pixelCount[i] = 0;
      vector<BYTE> out(width*fullHeight * 3); // .bmp pixel buffer
      memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is white

      double scale = pow(factor*factor, (double)i / (double)maxFrames);

      // now just fill in colour
      for (int p = 0; p < ps.size(); p++)
      {
        Vector2d offset = Vector2d(1.0, 1.0)*(double)width / 2.0;
        Vector2d pos = Vector2d(width/2 + 0.38, height/2+0.61) + (ps[p] - offset) * scale;
        double minScale = minWidth * scale;
        int x = (int)floor(pos[0] + minScale / (factor*factor));
        int y = (int)floor(pos[1] + minScale / (factor*factor));
        int x2 = (int)floor(pos[0] + minScale*(1.0 - 1.0 / (factor*factor)));
        int y2 = (int)floor(pos[1] + minScale*(1.0 - 1.0 / (factor*factor)));

        for (int X = x; X <= x2; X++)
        {
          for (int Y = y; Y <= y2; Y++)
          {
            if (X >= 0 && X < width && Y >= 0 && Y < height)
            {
              out[3 * (X + Y*width) + 0] = 0;
              out[3 * (X + Y*width) + 1] = 0;
              out[3 * (X + Y*width) + 2] = 0;
            }
          }
        }
      }
      if (run == 0)
      {
        for (int x = 0; x < width; x++)
          for (int y = 0; y < height; y++)
            if (out[3 * (x + y*width)] == 0)
              pixelCount[i]++;
      }
      if (run == 1)
      {
        for (int xx = 0; xx < width; xx++)
        {
          for (int yy = height; yy < fullHeight; yy++)
          {
            out[3 * (xx + yy*width) + 0] = 220;
            out[3 * (xx + yy*width) + 1] = 220;
            out[3 * (xx + yy*width) + 2] = 220;
          }
        }
        double maxCount = 0;
        for (int p = 0; p < maxFrames+1; p++)
          maxCount = max(maxCount, pixelCount[p]);
        pixelCount[maxFrames] = pixelCount[0]; // set exactly equal (should be very close)
        double midY = 0;
        double startX = -i;

        for (int p = 0; p < maxFrames+1; p++)
        {
          double y = (pixelCount[p] / maxCount) * (double)(extraHeight-10);
          int Y = 3*((fullHeight-5 - (int)y)/3);
          double x = (double)width * (1.0 + (double)(p-i) / (double)maxFrames);
          int X = 3*((((int)x + width / 2) % width)/3);

          for (int xx = 0; xx < 3; xx++)
          {
            for (int yy = 0; yy < 3; yy++)
            {
              int ind = 3 * (X + xx + (Y + yy)*width);
              out[ind + 0] = 0;
              out[ind + 1] = 0;
              out[ind + 2] = 50;
            }
          }
          if (p == i)
            midY = Y;
        }
        for (int xx = -4; xx < 5; xx++)
        {
          for (int yy = -4; yy < 5; yy++)
          {
            int ind = 3 * (width/2 + xx + (midY + yy)*width);
            out[ind + 0] = 128;
            out[ind + 1] = 0;
            out[ind + 2] = 0;
          }
        }

        BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, fullHeight, &s2);
        wstringstream str;
        str << L"imaginary/triple" << i << L".bmp";
        wstring strng = str.str();
        const TCHAR * bll = strng.c_str();
        LPCTSTR file = bll;
        SaveBMP(c, width, fullHeight, s2, file);
        delete[] c;
      }
    }
  }
}
