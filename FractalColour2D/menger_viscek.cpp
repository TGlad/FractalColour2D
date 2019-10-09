#include "stdafx.h"
#include "bmp.h"
#include <sstream>
#include <fstream>
struct Node
{
  Vector2d pos;
  Vector2d up;
  double length;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#define MENGER
//#define VISCEK

static int width = 3*3*3*3*3*3;
static int height = width;
static double minLength = 0.0015;

void recurse(vector<Node> &list, Node &node)
{
  if (node.length <= 1.001)
  {
    list.push_back(node);
    return;
  }
  Node child;
  child.length = node.length / 3.0;
  for (int x = -1; x <= 1; x++)
  {
    for (int y = -1; y <= 1; y++)
    {
      if (x == 0 && y == 0)
        continue;
      child.pos = node.pos + child.length * Vector2d(x,y);
      recurse(list, child);
    }
  }
}
void recursev(vector<Node> &list, Node &node)
{
  if (node.length <= 3*3*3 + 1e-5)
  {
    list.push_back(node);
//    recurse(list, node);
    return;
  }
  Node child;
  child.length = node.length / 3.0;
  child.pos = node.pos + Vector2d(-child.length, 0);
  recursev(list, child);
  child.pos = node.pos + Vector2d(0, 0);
  recursev(list, child);
  child.pos = node.pos + Vector2d(child.length, 0);
  recursev(list, child);
  child.pos = node.pos + Vector2d(0, -child.length);
  recursev(list, child);
  child.pos = node.pos + Vector2d(0, child.length);
  recursev(list, child);
}

void putpixel(vector<BYTE> &out, const Vector2i &pos, int shade)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= height)
    return;
  int ind = 3 * (pos[0] + width*pos[1]);
  out[ind + 0] = out[ind + 1] = out[ind + 2] = shade;
}

int _tmain(int argc, _TCHAR* argv[])
{
  Node base;
  base.pos = Vector2d(0, 0);
  base.up = Vector2d(0, 1);
  base.length = width;
  vector<Node> list;
  recursev(list, base);

  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  vector<double> graph;
  double c = 0;
  {
    memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey
    for (auto &node : list)
    {
      for (int x = -2; x <= 2; x++)
        for (int y = -2; y <= 2; y++)
          if (abs(x) <= 1 || abs(y)<=1)
            putpixel(out, Vector2i(x + node.pos[0]*0.6 + width / 2, y + node.pos[1]*0.6 + width / 2), 0);
    }
    double d = log(8.0) / log(3.0); 
    double d2 = log(5.0) / log(3.0);
    cout << "d: " << d << ", " << d2 << endl;

    BYTE* buf = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    LPCTSTR file = L"viscekb.bmp";
    SaveBMP(buf, width, height, s2, file);
    delete[] buf;
  }
}
