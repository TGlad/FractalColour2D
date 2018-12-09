#include "stdafx.h"
#include "bmp.h"
#include <fstream>

static int width = 5*5*5*5*5;
static int height = width;

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

void drawSquareRec(vector<BYTE> &out, int i, int j, int w, int shade)
{
  drawSquare(out, i, j, 5 * w, Vector3i(shade, shade, shade));
  if (w == 1)
    return;
  for (int x = 0; x < 2; x++)
  {
    for (int y = 0; y < 2; y++)
    {
      drawSquareRec(out, i+w+2*w*x, j+w+2*w*y, w/5, 255 - shade);
    }
  }
}
int _tmain(int argc, _TCHAR* argv[])
{
  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); 
  
  drawSquareRec(out, 0, 0, 5 * 5 * 5 * 5, 255);

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"tableSpongeCluster.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}
