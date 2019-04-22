// experiments with structural coloration of fractal curves, in particular Koch, Levy, Dragon and random curves.
#include "stdafx.h"
#include "bmp.h"
#include "spectrumToRGB.h"
#include <set>

static vector<Vector2i> pointset, pointset2;
static int gap = 4;
static int vgap = 14;
static int eachWidth = 2700;
static const int cols = 1;
static const int rows = 1;
static int width = (eachWidth + gap) * rows + 1 - gap;
static int height = (eachWidth + vgap) * cols;
static int X = 0;
static int Y = 0;

Vector3i colsx[5] = { Vector3i(0, 0, 0), Vector3i(255, 64, 44), Vector3i(64, 255, 44), Vector3i(64, 44, 255), Vector3i(230, 210, 44) };
void putpix(vector<BYTE> &out, const Vector2i &pos, int shade)
{
  if (pos[0] < -eachWidth / 2 || pos[0] >= eachWidth / 2 || pos[1] < -eachWidth / 2 || pos[1] > eachWidth/2)
    return;
  int ind = 3 * (X*(eachWidth + gap) + pos[0] + eachWidth/2 + width*(Y*(eachWidth + vgap) + pos[1] + eachWidth/2));
  for (int i = 0; i < 3; i++)
    out[ind + i] = colsx[shade][i];
}
int getpix(vector<BYTE> &out, const Vector2i &pos)
{
  if (pos[0] < -eachWidth / 2 || pos[0] >= eachWidth / 2 || pos[1] < -eachWidth / 2 || pos[1] > eachWidth / 2)
    return 1;
  int ind = 3 * (X*(eachWidth + gap) + pos[0] + eachWidth / 2 + width*(Y*(eachWidth + vgap) + pos[1] + eachWidth / 2));
  for (int i = 0; i<5; i++)
    if (out[ind + 0] == colsx[i][0] && out[ind + 1] == colsx[i][1] && out[ind + 2] == colsx[i][2])
      return i;
  return 0;
}
static double level = 2;

Vector2i toVector2i(const Vector2d &pos)
{
  Vector2i vec;
  vec[0] = (int)pos[0];
  if ((double)vec[0] > pos[0])
    vec[0]--;
  vec[1] = (int)pos[1];
  if ((double)vec[1] > pos[1])
    vec[1]--;
  return vec;
}

bool addMid(vector<BYTE> &out, const Vector2d &pos0, const Vector2d &pos1, const Vector2d &pos2)
{
  Vector2i mid = toVector2i((pos0 + pos1 + pos2) / 3.0);
  if (mid[0] < -eachWidth / 2 || mid[0] >= eachWidth / 2 || mid[1] < -eachWidth / 2 || mid[1] >= eachWidth/2)
    return false;
  if (getpix(out, mid) > 0)
    return false;
  int col0 = getpix(out, toVector2i(pos0));
  int col1 = getpix(out, toVector2i(pos1));
  int col2 = getpix(out, toVector2i(pos2));

  if (col0 == 0 || col1 == 0 || col2 == 0)
  {
    return false;
  }
  int colour;
  if (col0 == col1 && col0 == col2) // positivity
    colour = col0;
  // ambiguous case
  else if (col0 == col1)
    colour = col0;
  else if (col0 == col2)
    colour = col0;
  else if (col1 == col2)
    colour = col1;
  // if there's not 3 the same, or 2 the same, then there must be 1 of each
  else // 012=3, 013=4, 023=5 123=6 -->  3, 2, 1, 0, i.e. 6 - sum. .: 10-sum
    colour = 10 - (col0 + col1 + col2);
  putpix(out, mid, colour);
  return true;
}

int _tmain(int argc, _TCHAR* argv[])
{
  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  int startWidth = eachWidth;
  memset(&out[0], 0, out.size() * sizeof(BYTE)); 
  Vector2d xaxis[2] = { Vector2d(1, 0), Vector2d(cos(30.0*pi / 180.0), sin(30.0*pi / 180.0)) };
  Vector2d yaxis[2] = { Vector2d(cos(60.0*pi / 180.0), sin(60.0*pi / 180.0)), Vector2d(cos(30.0*pi / 180.0), -sin(30.0*pi / 180.0)) };
  srand(120);
  bool found = true;
  double scale = 81.0 * 3.0 *3.0;
  Vector2d offset(0.6, 0.6);
  // initialise
  int count = (int)(2.0*(double)eachWidth / scale);
  for (int i = -count; i < count; i++)
    for (int j = -count; j < count; j++)
    {
      int col = 0;
    /*  if ((i == 0 && j == 0) || (i == 1 && j == 0) || (i == 0 && j == 1))
        col = 0;
      else */
      /*if (i>0 && (i + j) > 0)
        col = 1;
      else if (j > 0)
        col = 2;
      else 
        col = 3;*/
      if (i == 0 && abs(j) == 1)
        col = 2;
      else if (abs(i) == 1 && j == 0)
        col = 3;
      else if (i == 1 && j == -1)
        col = 4;
      else if (i == -1 && j == 1)
        col = 4;
      else 
        col = 1;
   /*   if ((i == 0 && j == 1) || (i == 1 && j == 0))
        col = 2;
      else if ((i == 1 && j == -1) || (i == 0 && j == -1))
        col = 3;
      else if ((i == -1 && j == 0) || (i == -1 && j == 1))
        col = 4;
      else
        col = 1;*/
      putpix(out, toVector2i(offset + (xaxis[0] * (double)i + yaxis[0] * (double)j)*scale), col);// rand() % 4 + 1);
    }
  do
  {
    found = false;
    int count = (int)(2.0*(double)eachWidth / scale);
    for (int i = -count; i < count; i++)
    {
      for (int j = -count; j < count; j++)
      {
        Vector2d pos0 = offset + (xaxis[0] * (double)i + yaxis[0] * (double)j)*scale;
        Vector2d pos1 = offset + (xaxis[0] * (double)(i + 1) + yaxis[0] * (double)j)*scale;
        Vector2d pos2 = offset + (xaxis[0] * (double)i + yaxis[0] * (double)(j + 1))*scale;
        Vector2d pos3 = offset + (xaxis[0] * (double)(i+1) + yaxis[0] * (double)(j+1))*scale;
        found |= addMid(out, pos0, pos1, pos2);
        found |= addMid(out, pos1, pos2, pos3);
      }
    }
    if (!found)
      break;
    scale /= sqrt(3.0);
    count = (int)(2.0*(double)eachWidth / scale);
    for (int i = -count; i < count; i++)
    {
      for (int j = -count; j < count; j++)
      {
        Vector2d pos0 = offset + (xaxis[1] * (double)i + yaxis[1] * (double)j)*scale;
        Vector2d pos1 = offset + (xaxis[1] * (double)(i + 1) + yaxis[1] * (double)j)*scale;
        Vector2d pos2 = offset + (xaxis[1] * (double)i + yaxis[1] * (double)(j + 1))*scale;
        Vector2d pos3 = offset + (xaxis[1] * (double)(i+1) + yaxis[1] * (double)(j + 1))*scale;
        found |= addMid(out, pos0, pos1, pos2);
        found |= addMid(out, pos1, pos2, pos3);
      }
    }
    scale /= sqrt(3.0);
  } while (found); 
  
  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"thingo.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}

