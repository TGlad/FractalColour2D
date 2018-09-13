// experiments with structural coloration of fractal curves, in particular Koch, Levy, Dragon and random curves.
#include "stdafx.h"
#include "bmp.h"
#include "spectrumToRGB.h"
#include <set>
//#define DIAGONAL
//#define ADJACENT
//#define FIVE_SIDE
//#define THREE_SIDE
#define THREE_SIDE_VOID
//#define THREE_SIDE2
//#define SPONGE
//#define SPONGE_CLUSTER
//#define TREE
//#define TREE2
//#define AUTOMATA
//#define AUTOMATA3COLOUR
//#define AUTOMATA_TRIANGLE 

static vector<Vector2i> pointset, pointset2;
static int gap = 4;
static int vgap = 14;
#if defined(FIVE_SIDE)
static int depth = 3;
static int eachWidth = (int)pow(5, depth);
static const int rows = 27;
static const int cols = 27;
#elif defined(THREE_SIDE_VOID)
static int depth = 5;
static int eachWidth = (int)pow(3, depth);
static const int rows = 4;
static const int cols = 2;
#elif defined(THREE_SIDE)
static int depth = 5;
static int eachWidth = (int)pow(3, depth);
static const int rows = 9;
static const int cols = 3;
#elif defined(THREE_SIDE2)
static int depth = 6;
static int eachWidth = (int)pow(3, depth);
static const int rows = 8;
static const int cols = 8;
#elif defined(SPONGE)
static int depth = 4;
static int eachWidth = 2 * (int)pow(5, depth);
static const int rows = 1;
static const int cols = 1;
#elif defined(SPONGE_CLUSTER)
static int depth = 4;
static int eachWidth = 4*(int)pow(6, depth);
static const int rows = 1;
static const int cols = 1;
#elif defined(TREE)
static int depth = 5;
static int eachWidth = 4 * (int)pow(4, depth);
static const int rows = 1;
static const int cols = 1;
#elif defined(TREE2)
static int depth = 5;
static int eachWidth =  4*(int)pow(4, depth);
static const int rows = 1;
static const int cols = 1;
#elif defined(AUTOMATA)
static int depth = 7;
static int eachWidth = 3 * (1 << depth);
static const int rows = 8;
static const int cols = 8;
#elif defined(AUTOMATA3COLOUR)
static int depth = 10;
static int eachWidth = 3 * (1 << depth);
static const int cols = 1;
static const int rows = 1;
#elif defined(AUTOMATA_TRIANGLE)
static int eachWidth = 2700;
static const int cols = 1;
static const int rows = 1;
#else
static int depth = 8;
static int eachWidth = 1 << depth;
static const int rows = 8;
static const int cols = 8;
#endif
static int width = (eachWidth + gap) * rows + 1 - gap;
static int height = (eachWidth + vgap) * cols;
static int X = 0;
static int Y = 0;

void putpixel(vector<BYTE> &out, const Vector2i &pos, int shade)
{
  if (pos[0] < -eachWidth / 2 || pos[0] >= eachWidth / 2 || pos[1] < -eachWidth / 2 || pos[1] >= eachWidth / 2)
    return;
  int ind = 3 * (X*(eachWidth + gap) + pos[0] + eachWidth / 2 + width*(Y*(eachWidth + vgap) + pos[1] + eachWidth / 2));
  out[ind + 0] = out[ind + 1] = out[ind + 2] = shade;
}
int getpixel(vector<BYTE> &out, const Vector2i &pos)
{
  if (pos[0] < -eachWidth / 2 || pos[0] >= eachWidth / 2 || pos[1] < -eachWidth / 2 || pos[1] >= eachWidth / 2)
    return 0;
  int ind = 3 * (X*(eachWidth + gap) + pos[0] + eachWidth / 2 + width*(Y*(eachWidth + vgap) + pos[1] + eachWidth / 2));
  return out[ind + 0];
}

#if defined(AUTOMATA3COLOUR)
void putpix(vector<BYTE> &out, const Vector2i &pos, int shade)
{
  if (pos[0] < 0 || pos[0] > eachWidth || pos[1] < 0 || pos[1] > eachWidth)
    return;
  int ind = 3 * (X*(eachWidth + gap) + pos[0] + width*(Y*(eachWidth + vgap) + pos[1]));
  out[ind + 0] = out[ind + 1] = out[ind + 2] = 64;
  out[ind + shade] = 255;
}
int getpix(vector<BYTE> &out, const Vector2i &pos)
{
  if (pos[0] < 0 || pos[0] > eachWidth || pos[1] < 0 || pos[1] > eachWidth)
    return 0;
  int ind = 3 * (X*(eachWidth + gap) + pos[0] + width*(Y*(eachWidth + vgap) + pos[1]));
  if (out[ind + 0] > out[ind + 1])
    return out[ind + 0] > out[ind + 2] ? 0 : 2;
  return out[ind + 1] > out[ind + 2] ? 1 : 2;
}
#elif defined(AUTOMATA_TRIANGLE)
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
#else
void putpix(vector<BYTE> &out, const Vector2i &pos, int shade)
{
  if (pos[0] < 0 || pos[0] > eachWidth || pos[1] < 0 || pos[1] > eachWidth)
    return;
  int ind = 3 * (X*(eachWidth + gap) + pos[0] + width*(Y*(eachWidth + vgap) + pos[1]));
  out[ind + 0] = out[ind + 1] = out[ind + 2] = shade;
}
int getpix(vector<BYTE> &out, const Vector2i &pos)
{
  if (pos[0] < 0 || pos[0] > eachWidth || pos[1] < 0 || pos[1] > eachWidth)
    return 0;
  int ind = 3 * (X*(eachWidth + gap) + pos[0] + width*(Y*(eachWidth + vgap) + pos[1]));
  return out[ind + 0];
}
#endif
static double level = 2;

void colour(vector<BYTE> &out, const Vector2i &pos, const Vector2i &xAxis, const Vector2i &yAxis, int length, bool flip = false)
{
  // draw diagonals
#if defined(DIAGONAL)
  length /= 2;
  int minx = min(xAxis[0], yAxis[0])*length;
  int maxx = max(xAxis[0], yAxis[0])*length;
  int miny = min(xAxis[1], yAxis[1])*length;
  int maxy = max(xAxis[1], yAxis[1])*length;
  for (int lx = minx; lx < maxx; lx++)
    for (int ly = miny; ly < maxy; ly++)
      putpixel(out, pos + Vector2i(lx, ly), 0);
  for (int lx = -maxx; lx < -minx; lx++)
    for (int ly = -maxy; ly < -miny; ly++)
      putpixel(out, pos + Vector2i(lx,ly), 255);
  if (length > 1)
  {
    Vector2i xAxis2 = X < 4 ? xAxis : yAxis;
    Vector2i yAxis2 = X < 4 ? yAxis : xAxis;
    if (X % 2)
      xAxis2 = -xAxis2;
    if ((X / 2) % 2)
      yAxis2 = -yAxis2;
    colour(out, pos - xAxis*length / 2 + yAxis*length / 2, xAxis2, yAxis2, length);

    Vector2i xAxis3 = Y < 4 ? xAxis : yAxis;
    Vector2i yAxis3 = Y < 4 ? yAxis : xAxis;
    if (Y % 2)
      yAxis3 = -yAxis3;
    if ((Y / 2) % 2)
      xAxis3 = -xAxis3;
    colour(out, pos + xAxis*length / 2 - yAxis*length / 2, xAxis3, yAxis3, length);
  }
#elif defined(ADJACENT)
  length /= 2;
  int minx = min(-xAxis[0], yAxis[0])*length;
  int maxx = max(-xAxis[0], yAxis[0])*length;
  int miny = min(-xAxis[1], yAxis[1])*length;
  int maxy = max(-xAxis[1], yAxis[1])*length;
  int minx2 = min(-xAxis[0], -yAxis[0])*length;
  int maxx2 = max(-xAxis[0], -yAxis[0])*length;
  int miny2 = min(-xAxis[1], -yAxis[1])*length;
  int maxy2 = max(-xAxis[1], -yAxis[1])*length;
  for (int lx = minx; lx < maxx; lx++)
    for (int ly = miny; ly < maxy; ly++)
      putpixel(out, pos + Vector2i(lx, ly), 0);
  for (int lx = minx2; lx < maxx2; lx++)
    for (int ly = miny2; ly < maxy2; ly++)
      putpixel(out, pos + Vector2i(lx, ly), 255);
  if (length > 1)
  {
    Vector2i xAxis2 = X < 4 ? xAxis : yAxis;
    Vector2i yAxis2 = X < 4 ? yAxis : xAxis;
    if (X % 2)
      xAxis2 = -xAxis2;
    if ((X / 2) % 2)
      yAxis2 = -yAxis2;
    colour(out, pos + xAxis*length / 2 + yAxis*length / 2, xAxis2, yAxis2, length);

    Vector2i xAxis3 = Y < 4 ? xAxis : yAxis;
    Vector2i yAxis3 = Y < 4 ? yAxis : xAxis;
    if (Y % 2)
      yAxis3 = -yAxis3;
    if ((Y / 2) % 2)
      xAxis3 = -xAxis3;
    colour(out, pos + xAxis*length / 2 - yAxis*length / 2, xAxis3, yAxis3, length);
  }
#elif defined(FIVE_SIDE)
  length /= 5;
  int types[3][3] = { { X%3, (X/3)%3, (X/9)%3 }, { 1, Y%3, (Y/3)%3 }, { 1, 1, (Y/9)%3 } };
  for (int i = -2; i <= 2; i++)
  {
    for (int j = -2; j <=2; j++)
    {
      int x = pos[0] + length*i;
      int y = pos[1] + length*j;

      int I = abs(i);
      int J = abs(j);
      if (I > J)
        swap(I, J);
      if (types[I][J] < 2)
      {
        for (int lx = -length/2; lx <= length/2; lx++)
          for (int ly = -length/2; ly <= length/2; ly++)
            putpixel(out, Vector2i(x + lx, y + ly), types[I][J]*255);
      }
      else if (length > 1)
      {
        colour(out, Vector2i(x, y), xAxis, yAxis, length);
      }
    }
  }
#elif defined(THREE_SIDE_VOID)
  length /= 3;
  int types[2][2] = { { X % 2, (X / 2) % 2 }, { 1, Y % 2 } }; 
  for (int i = -1; i <= 1; i++)
  {
    for (int j = -1; j <= 1; j++)
    {
      int x = pos[0] + length*i;
      int y = pos[1] + length*j;

      int I = abs(i);
      int J = abs(j);
      if (I > J)
        swap(I, J);
      if (types[I][J] < 1)
      {
        for (int lx = -length / 2; lx <= length / 2; lx++)
          for (int ly = -length / 2; ly <= length / 2; ly++)
            putpixel(out, Vector2i(x + lx, y + ly), types[I][J] * 255);
      }
      else if (length > 1)
      {
        colour(out, Vector2i(x, y), xAxis, yAxis, length);
      }
    }
  }
#elif defined(THREE_SIDE)
length /= 3;
int types[2][2] = { { X % 3, (X / 3) % 3 }, { 1, Y % 3 } }; 
for (int i = -1; i <= 1; i++)
{
  for (int j = -1; j <= 1; j++)
  {
    int x = pos[0] + length*i;
    int y = pos[1] + length*j;

    int I = abs(i);
    int J = abs(j);
    if (I > J)
      swap(I, J);
    if (types[I][J] < 2)
    {
      for (int lx = -length / 2; lx <= length / 2; lx++)
        for (int ly = -length / 2; ly <= length / 2; ly++)
          putpixel(out, Vector2i(x + lx, y + ly), types[I][J] * 255);
    }
    else if (length > 1)
    {
      colour(out, Vector2i(x, y), xAxis, yAxis, length);
    }
  }
}
#elif defined(THREE_SIDE2)
  length /= 3;
  int ind = X + 8 * Y;
  int types[2][2] = { { ind%4, (ind/4)%4 }, { 1, (ind/16)%4 } };
  for (int i = -1; i <= 1; i++)
  {
    for (int j = -1; j <= 1; j++)
    {
      int x = pos[0] + length*i;
      int y = pos[1] + length*j;

      int I = abs(i);
      int J = abs(j);
      if (I > J)
        swap(I, J);
      if (types[I][J] < 2)
      {
        int z = types[I][J] * 255;
        if (flip)
          z = 255 - z;
        for (int lx = -length / 2; lx <= length / 2; lx++)
          for (int ly = -length / 2; ly <= length / 2; ly++)
            putpixel(out, Vector2i(x + lx, y + ly), z);
      }
      else if (length > 1)
      {
        bool newflip = flip;
        if (types[I][J] == 3)
          newflip = !flip;
        colour(out, Vector2i(x, y), xAxis, yAxis, length, newflip);
      }
    }
  }
#elif defined(SPONGE)
int radius = length / 10;
// for each vector2i, draw a 2*radius ring...
vector<Vector2i> newPoints;
for (int i = 0; i < (int)pointset.size(); i++)
{
  Vector2i pos = pointset[i];
  // draw 4 rectangles
  for (int x = pos[0] - 5 * radius; x <= pos[0] + 5 * radius; x+=2*radius)
  {
    for (int y = pos[1] - 5 * radius; y <= pos[1] + 5 * radius; y+=2*radius)
    {
      bool inside = false;
      if (abs(x - pos[0]) < 5 * radius && abs(y - pos[1]) < 5 * radius)
        inside = true;
//      if (inside)
//        continue;
      for (int x2 = -radius; x2 < radius; x2++)
        for (int y2 = -radius; y2 < radius; y2++)
        {
          if (!inside)
            putpixel(out, Vector2i(x + x2, y + y2), 255);
          else if (getpixel(out, Vector2i(x + x2, y + y2)) != 255)
            putpixel(out, Vector2i(x + x2, y + y2), (int)(255.0 - 255.0*0.5*level));
        }
      if (radius > 1 && !inside)
      {
        for (int x2 = -2 * radius; x2 <= 2 * radius; x2 += 2 * radius)
        {
          for (int y2 = -2 * radius; y2 <= 2 * radius; y2 += 2 * radius)
          {
            if (abs(y2) != 2 * radius && abs(x2) != 2 * radius)
              continue;
            Vector2i newpos(x + x2, y + y2);
            if (find(newPoints.begin(), newPoints.end(), newpos) == std::end(newPoints))
              newPoints.push_back(newpos);
          }
        }
      }
    }
  }
}
int oldCol = (int)(255.0 - 255.0*0.5*level);
level *= 0.6;
int newCol = (int)(255.0 - 255.0*0.5*level * 0.5);
pointset = newPoints;
if (radius > 1)
  colour(out, Vector2i(0, 0), Vector2i(0, 0), Vector2i(0, 0), length / 5);
else 
{
  // finish off in grey
  for (int x = 0; x < eachWidth; x += 2)
  {
    for (int y = 0; y < eachWidth; y += 2)
    {
      bool whiteFound = false;
      for (int i = -1; i <= 1; i++)
      {
        for (int j = -1; j <= 1; j++)
          if (getpix(out, Vector2i(x + 2 * i, y + 2 * j)) == 255)
            whiteFound = true;
      }
      if (whiteFound && getpix(out, Vector2i(x, y)) !=255)
      {
        int add = (depth % 2) ? -1 : 1;
        putpix(out, Vector2i(x, y), newCol);
        putpix(out, Vector2i(x + add, y), newCol);
        putpix(out, Vector2i(x, y + add), newCol);
        putpix(out, Vector2i(x + add, y + add), newCol);
      }
    }
  }
}
#elif defined(SPONGE_CLUSTER)
int radius = length / 12;
// for each vector2i, draw a 2*radius ring...
vector<Vector2i> newPoints;
for (int i = 0; i < (int)pointset.size(); i++)
{
  Vector2i pos = pointset[i];
  // draw 4 rectangles
  if (radius > 1)
    newPoints.push_back(pos);
  for (int x = pos[0] - 6 * radius; x <= pos[0] + 6 * radius; x += 2 * radius)
  {
    for (int y = pos[1] - 6 * radius; y <= pos[1] + 6 * radius; y += 2 * radius)
    {
      int dx = abs(x - pos[0]);
      int dy = abs(y - pos[1]);
      if (dx < 6 * radius && dy < 6 * radius)
        continue;
      for (int x2 = -radius; x2 < radius; x2++)
        for (int y2 = -radius; y2 < radius; y2++)
          putpixel(out, Vector2i(x + x2, y + y2), 255);
      if (radius > 1)
      {
        for (int x2 = -2 * radius; x2 <= 2 * radius; x2 += 2 * radius)
        {
          for (int y2 = -2 * radius; y2 <= 2 * radius; y2 += 2 * radius)
          {
            if (abs(y2) != 2 * radius && abs(x2) != 2 * radius)
              continue;
            Vector2i newpos(x + x2, y + y2);
            if (find(newPoints.begin(), newPoints.end(), newpos) == std::end(newPoints))
              newPoints.push_back(newpos);
          }
        }
      }
    }
  }
}
pointset = newPoints;
if (radius > 1)
colour(out, Vector2i(0, 0), Vector2i(0, 0), Vector2i(0, 0), length / 6);
else
{
  // finish off in grey
  for (int x = 0; x < eachWidth; x += 2)
  {
    for (int y = 0; y < eachWidth; y += 2)
    {
      bool whiteFound = false;
      for (int i = -1; i <= 1; i++)
      {
        for (int j = -1; j <= 1; j++)
          if (getpix(out, Vector2i(x + 2 * i, y + 2 * j)) == 255)
            whiteFound = true;
      }
      if (whiteFound && getpix(out, Vector2i(x, y)) == 0)
      {
        int add = -1;
        putpix(out, Vector2i(x, y), 127);
        putpix(out, Vector2i(x + add, y), 127);
        putpix(out, Vector2i(x, y + add), 127);
        putpix(out, Vector2i(x + add, y + add), 127);
      }
    }
  }
  // now fill in little centre dots
  for (int x = 0; x < eachWidth; x += 2)
  {
    for (int y = 0; y < eachWidth; y += 2)
    {
      int greyFound = 0;
      for (int i = -2; i <= 2; i+=2)
      {
        for (int j = -2; j <= 2; j+=2)
          if (getpix(out, Vector2i(x + 2 * i, y + 2 * j)) == 127)
            greyFound++;
      }
      if (greyFound==8 && getpix(out, Vector2i(x, y)) == 0)
      {
        int add = -1;
        putpix(out, Vector2i(x, y), 127);
        putpix(out, Vector2i(x + add, y), 127);
        putpix(out, Vector2i(x, y + add), 127);
        putpix(out, Vector2i(x + add, y + add), 127);
      }
    }
  }
}
#elif defined(TREE) || defined (TREE2)
int radius = length / 4;
// for each vector2i, draw a 2*radius ring...
const int numChildren = 16;
Vector2d children[numChildren] = { Vector2d(-2.5, -1.5), Vector2d(-1.5, -2.5), Vector2d(1.5, -2.5), Vector2d(2.5, -1.5),
Vector2d(-2.5, 1.5), Vector2d(-1.5, 2.5), Vector2d(1.5, 2.5), Vector2d(2.5, 1.5),
Vector2d(-1.5, -0.5), Vector2d(-1.5, 0.5), Vector2d(1.5, -0.5), Vector2d(1.5, 0.5),
Vector2d(-0.5, -1.5), Vector2d(0.5, -1.5), Vector2d(-0.5, 1.5), Vector2d(0.5, 1.5) };
vector<Vector2i> newPoints;
for (int i = 0; i < (int)pointset.size(); i++)
{
  Vector2i pos = pointset[i];
  for (int x2 = -length/2; x2 < length/2; x2++)
    for (int y2 = -length / 2; y2 < length / 2; y2++)
    {
//#define TREE_CLUSTER
#if defined(TREE_CLUSTER)
      if ((x2 >= -length/4) && (x2 < length/4) && (y2 >= -length/4) && (y2 < length/4))
        continue;
#endif
#if defined(TREEVOID)
      if (getpixel(out, Vector2i(pos[0] + x2, pos[1] + y2))==0)
        putpixel(out, Vector2i(pos[0] + x2, pos[1] + y2), (int)(255.0 * 0.5 * level ));
#else
      putpixel(out, Vector2i(pos[0] + x2, pos[1] + y2), 255);
#endif
    }
  
  if (radius >= 1)
  {
    for (int c = 0; c < numChildren; c++)
    {
      Vector2d p = children[c] * radius;
      newPoints.push_back(pos + Vector2i(p[0], p[1]));
    }
  }
}
level*=0.6;
#if defined(TREE2)
const int numClusters = 4;
Vector2d clusters[numClusters] = {
//  Vector2d(-5, -3), Vector2d(-3, -5), Vector2d(3, -5), Vector2d(5, -3),
//  Vector2d(-5, 3), Vector2d(-3, 5), Vector2d(3, 5), Vector2d(5, 3),
//  Vector2d(-5, 0), Vector2d(0, -5), Vector2d(5, 0), Vector2d(0, 5),
  Vector2d(-5, -5), Vector2d(5, -5), Vector2d(5, 5), Vector2d(-5, 5),
};
vector<Vector2i> newPoints2;
for (int i = 0; i < (int)pointset2.size(); i++)
{
  Vector2i pos = pointset2[i];
  for (int x2 = -length / 2; x2 < length / 2; x2++)
    for (int y2 = -length / 2; y2 < length / 2; y2++)
      putpixel(out, Vector2i(pos[0] + x2, pos[1] + y2), 255);

  if (radius >= 1)
  {
    for (int c = 0; c < numChildren; c++)
    {
      Vector2d p = children[c] * radius;
      newPoints.push_back(pos + Vector2i(p[0], p[1]));
    }
    for (int c = 0; c < numClusters; c++)
    {
      Vector2d p = clusters[c] * radius;
      newPoints2.push_back(pos + Vector2i(p[0], p[1]));
    }
  }
}
pointset2 = newPoints2;
#endif
pointset = newPoints;
if (radius >= 1)
  colour(out, Vector2i(0, 0), Vector2i(0, 0), Vector2i(0, 0), radius); 
#endif
}
#if defined(AUTOMATA)
void square(vector<BYTE> &out, int x0, int y0, int x2, int y2)
{
  int on[] = { X % 2, (X / 2) % 2, (X / 4) % 2, Y % 2, (Y / 2) % 2, (Y / 4) % 2 };
  int x1 = (x0 + x2) / 2;
  int y1 = (y0 + y2) / 2;
  int count = 0; // will be 0 to 5
  count += getpix(out, Vector2i(x0, y0)) / 255;
  count += getpix(out, Vector2i(x0, y2)) / 255;
  count += getpix(out, Vector2i(x2, y0)) / 255;
  count += getpix(out, Vector2i(x2, y2)) / 255;
  if (count == 2)
  {
    if ((getpix(out, Vector2i(x0, y0)) > 0) == (getpix(out, Vector2i(x2, y2)) > 0)) // extra case
      count = 5;
  }
  putpix(out, Vector2i(x1, y1), on[count] ? 255 : 0);
}

void diamond(vector<BYTE> &out, int x0, int y0, int x2, int y2)
{
  int on[] = { X % 2, (X / 2) % 2, (X / 4) % 2, Y % 2, (Y / 2) % 2, (Y / 4) % 2 };
  int x1 = (x0 + x2) / 2;
  int y1 = (y0 + y2) / 2;
  int count = 0; // will be 0 to 5
  count += getpix(out, Vector2i(x0, y1)) / 255;
  count += getpix(out, Vector2i(x2, y1)) / 255;
  count += getpix(out, Vector2i(x1, y0)) / 255;
  count += getpix(out, Vector2i(x1, y2)) / 255;
  if (count == 2)
  {
    if ((getpix(out, Vector2i(x0, y1)) > 0) == (getpix(out, Vector2i(x2, y1)) > 0)) // extra case
      count = 5;
  }
  putpix(out, Vector2i(x1, y1), on[count] ? 255 : 0);
}
#elif defined(AUTOMATA3COLOUR)
int getcol(int count1, int count2)
{
  int count0 = 4 - count1 - count2;
  // positivity
  int col;
  if (count0 == 4)
    col = 0;
  else if (count1 == 4)
    col = 1;
  else if (count2 == 4)
    col = 2;
  // unambiguous
  else if (count0 == count1)
    col = 2;
  else if (count0 == count2)
    col = 1;
  else if (count1 == count2)
    col = 0;
  // decision:
  else
  {
    int zeroCol = (count0 == 0 ? 0 : (count1 == 0 ? 1 : 2));
    int oneCol = (count0 == 1 ? 0 : (count1 == 1 ? 1 : 2));
    int threeCol = (count0 == 3 ? 0 : (count1 == 3 ? 1 : 2));
    col = X == 0 ? zeroCol : (X == 1 ? oneCol : threeCol);
  }
  return col;
}
void square(vector<BYTE> &out, int x0, int y0, int x2, int y2)
{ 
  int x1 = (x0 + x2) / 2;
  int y1 = (y0 + y2) / 2;
  int count1 = 0, count2 = 0; // will be 0 to 5
  count1 += getpix(out, Vector2i(x0, y0)) == 1;
  count1 += getpix(out, Vector2i(x0, y2)) == 1;
  count1 += getpix(out, Vector2i(x2, y0)) == 1;
  count1 += getpix(out, Vector2i(x2, y2)) == 1;
  count2 += getpix(out, Vector2i(x0, y0)) == 2;
  count2 += getpix(out, Vector2i(x0, y2)) == 2;
  count2 += getpix(out, Vector2i(x2, y0)) == 2;
  count2 += getpix(out, Vector2i(x2, y2)) == 2;
  putpix(out, Vector2i(x1, y1), getcol(count1, count2));
}

void diamond(vector<BYTE> &out, int x0, int y0, int x2, int y2)
{ 
  int x1 = (x0 + x2) / 2;
  int y1 = (y0 + y2) / 2;
  int count1 = 0, count2 = 0; // will be 0 to 5
  count1 += getpix(out, Vector2i(x0, y1)) == 1;
  count1 += getpix(out, Vector2i(x2, y1)) == 1;
  count1 += getpix(out, Vector2i(x1, y0)) == 1;
  count1 += getpix(out, Vector2i(x1, y2)) == 1;
  count2 += getpix(out, Vector2i(x0, y1)) == 2;
  count2 += getpix(out, Vector2i(x2, y1)) == 2;
  count2 += getpix(out, Vector2i(x1, y0)) == 2;
  count2 += getpix(out, Vector2i(x1, y2)) == 2;
  putpix(out, Vector2i(x1, y1), getcol(count1, count2));
}
#endif
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
#if defined(SPONGE) || defined(SPONGE_CLUSTER) || defined(TREE) || defined(TREE2)
  memset(&out[0], 0, out.size() * sizeof(BYTE)); // background is grey
  pointset.push_back(Vector2i(0, 0));
  pointset2.push_back(Vector2i(0, 0));
 // startWidth /= 2; // to see whoe thing (optional)
#if defined(TREE2)
  startWidth /= 2;
#endif
#elif defined(AUTOMATA) || defined(AUTOMATA3COLOUR) || defined(AUTOMATA_TRIANGLE)
  memset(&out[0], 0, out.size() * sizeof(BYTE)); 
#else
  memset(&out[0], 127, out.size() * sizeof(BYTE)); // background is grey
#endif
#if !defined(AUTOMATA_TRIANGLE)
  for (X = 0; X < rows; X++)
  {
    for (Y = 0; Y < cols; Y++)
    {
#if defined(AUTOMATA) || defined(AUTOMATA3COLOUR)
//      int start[4][4] = { { 0, 0, 0, 1 },
//      { 0, 0, 1, 0 },
//      { 0, 0, 1, 1 },
//      { 1, 1, 1, 1 } };
      int start[4][4] = { { 0, 2, 0, 1 },
      { 1, 1, 2, 1 },
      { 1, 2, 0, 2 },
      { 0, 1, 2, 2 } };
      int w = eachWidth / 3;
      for (int x = 0; x < 4; x++)
        for (int y = 0; y < 4; y++)
#if defined(AUTOMATA3COLOUR)
          putpix(out, Vector2i(x, y)*w, start[y][x]);
#else
          putpix(out, Vector2i(x, y)*w, start[y][x]*255);
#endif
      while (w > 1)
      {
        for (int x = 0; x < eachWidth; x += w)
          for (int y = 0; y < eachWidth; y += w)
            square(out, x, y, x + w, y + w);
        // diamonds are more tricky...
        for (int x = 0; x < eachWidth+1; x += w)
        {
          for (int y = 0; y < eachWidth+1; y += w)
          {
            if (y<eachWidth)
              diamond(out, x - w / 2, y, x + w - w / 2, y + w);
            if (x<eachWidth)
              diamond(out, x, y - w / 2, x + w, y + w - w / 2);
          }
        }
        w = w / 2;
      }
#else
      colour(out, Vector2i(0, 0), Vector2i(1, 0), Vector2i(0, 1), startWidth, false);
#endif
    }
  }
#else
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
  /*    if (i == 0 && abs(j) == 1)
        col = 2;
      else if (abs(i) == 1 && j == 0)
        col = 3;
      else if (i == 1 && j == -1)
        col = 4;
      else if (i == -1 && j == 1)
        col = 4;
      else 
        col = 1;*/
  /*    if ((i == 0 && j == 1) || (i == 1 && j == 0))
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
  
  #endif


  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"basictypes_void.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}

