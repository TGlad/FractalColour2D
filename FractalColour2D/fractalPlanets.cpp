#include "stdafx.h"
#include "bmp.h"
#include <fstream>

static int width = 1024;
static int height = 1024;

void putpixel(vector<BYTE> &out, const Vector2i &pos, int shade)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= width)
    return;
  int ind = 3 * (pos[0] + width*pos[1]);
  out[ind + 0] = out[ind + 1] = out[ind + 2] = shade;
}
void drawDisk(const Vector2d &pos, vector<BYTE> &out, double rad, int shade)
{
  for (int x = (int)(pos[0] - rad); x <= (int)(pos[0] + rad); x++)
    for (int y = (int)(pos[1] - rad); y <= (int)(pos[1] + rad); y++)
      if (sqr(x - pos[0]) + sqr(y - pos[1]) <= sqr(rad))
        putpixel(out, Vector2i(x, y), 255 - shade);
}
struct Node
{
  Vector2d pos;
  double angle;
  Vector2d xAxis() const { return Vector2d(cos(angle), sin(angle)); }
  Vector2d yAxis() const { return Vector2d(-sin(angle), cos(angle)); }
  double radius;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

double time = 250.0;


void buildCluster(vector<Node> &cluster, const Node &node)
{
  cluster.push_back(node);
  if (node.radius < 1.0)
    return;

  double scale1 = 1.0 / 3.41828;
  double orbitalDist1 = 3.0;
  double scale2 = 1.0 / 4.317;
  double orbitalDist2 = 6.315;

  Node child1;
  child1.angle = time / (node.radius * orbitalDist1);
  child1.pos = node.pos + child1.yAxis() * node.radius * orbitalDist1;
  child1.radius = node.radius * scale1;
  buildCluster(cluster, child1);

  Node child2;
  child2.angle = time / (node.radius * orbitalDist2);
  child2.pos = node.pos + child2.yAxis() * node.radius * orbitalDist2;
  child2.radius = node.radius * scale2;
  buildCluster(cluster, child2);
}

int _tmain(int argc, _TCHAR* argv[])
{

  Node base;
  base.angle = 0;
  base.radius = 150;
  base.pos = Vector2d(700, -50);

  vector<Node> cluster;

  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey

  for (int shade = 0; shade < 256; shade += 1)
  {
//    time = 240.0 + 75.0 * (double)shade / 255.0;
    time = 190.0 + 125.0 * (double)shade / 255.0;
    cluster.clear();
    buildCluster(cluster, base);
    for (auto &planet : cluster)
    {
      if (shade == 255)
        drawDisk(planet.pos, out, planet.radius, 255);
      else
        drawDisk(planet.pos, out, planet.radius / 2.0, min(shade, 127));
    }
  }

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"planets5.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}
