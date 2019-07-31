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
  if (node.radius < 0.5)
    return;

  Vector2d offset(512, 512);
  double scale1 = 0.8;
  Node child1;
  Vector2d vec = node.pos - offset;
  double ang = 0.1; // 0.0 for unrotated version.
  Vector2d rotated(cos(ang)*vec[0] + sin(ang)*vec[1], -sin(ang)*vec[0] + cos(ang)*vec[1]);
  child1.angle = 0;
  child1.pos = offset + rotated * scale1;
  child1.radius = node.radius * scale1;
  buildCluster(cluster, child1);

 // Node child2;
 // child2.angle = time / (node.radius * orbitalDist2);
 // child2.pos = node.pos + child2.yAxis() * node.radius * orbitalDist2;
 // child2.radius = node.radius * scale2;
 // buildCluster(cluster, child2);
}

void buildCluster2(vector<Node> &cluster, const Node &node)
{
  cluster.push_back(node);
  if (node.radius < 0.5)
    return;

  Vector2d offset(512, 512);
  double scale1 = 0.5;
  Node child1;
  Vector2d vec = node.pos - offset;
  double ang = 0.5; // 0.0 for unrotated version.

  child1.angle = node.angle -0.6;
  child1.pos = node.pos + child1.yAxis()*node.radius*2.5;
  child1.radius = node.radius * scale1;
  buildCluster2(cluster, child1);

  Node child2;
  child2.angle = node.angle + 0.35;
  child2.pos = node.pos + child2.yAxis()*node.radius*3.0;
  child2.radius = node.radius * scale1;
  buildCluster2(cluster, child2);
}

int _tmain(int argc, _TCHAR* argv[])
{
  Node base;
  base.angle = 3.1415;
  base.radius = 120;
//  base.pos = Vector2d(0, 95);
  base.pos = Vector2d(512, 900);

  vector<Node> cluster;

  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey

  buildCluster2(cluster, base);
  for (auto &planet : cluster)
    drawDisk(planet.pos, out, planet.radius, 255);
/*
  base.pos = Vector2d(1024, 500);
  vector<Node> cluster2;
  buildCluster(cluster2, base);
  base.pos = Vector2d(354, 1024);
  vector<Node> cluster3;
  buildCluster(cluster3, base);
  for (auto &planet : cluster2)
    drawDisk(planet.pos, out, planet.radius, 255);
  for (auto &planet : cluster3)
    drawDisk(planet.pos, out, planet.radius, 255);
  */
  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"discreteFractalSS2.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}
