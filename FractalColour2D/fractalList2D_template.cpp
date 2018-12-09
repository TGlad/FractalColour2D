#include "stdafx.h"
#include "bmp.h"
#include <fstream>

int blockWidth = 1000;
int gap = 50;
int startGap = 250;
int endGap = 5;
int aboveGap = 5;
int belowGap = 100;
int thickness = 5;

static int width = startGap + 1*blockWidth + 4*gap + endGap;
static int height = aboveGap + 5*blockWidth + 4*gap + belowGap;

void putpixel(vector<BYTE> &out, const Vector2i &pos, const Vector3i &colour)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= height)
    return;
  int ind = 3 * (pos[0] + width*pos[1]);
  out[ind + 0] = colour[0];
  out[ind + 1] = colour[1]; 
  out[ind + 2] = colour[2];
}

void drawSquare(vector<BYTE> &out, int i, int j, int w, const Vector3i &colour)
{
  for (int x = i; x < i + w; x++)
  {
    for (int y = j; y < j + w; y++)
      putpixel(out, Vector2i(x, y), colour);
  }
}

int _tmain(int argc, _TCHAR* argv[])
{
  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); 
  
  for (int i = 0; i < 1; i++)
  {
    for (int j = 0; j < 5; j++)
    {
      drawSquare(out, startGap + (blockWidth + gap)*i - thickness, aboveGap+(blockWidth+gap)*j - thickness, blockWidth + 2 * thickness, Vector3i(127,127,255));
      drawSquare(out, startGap + (blockWidth + gap)*i, aboveGap + (blockWidth + gap)*j, blockWidth, Vector3i(0, 255, 0));
    }
  }

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"list_template.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}
