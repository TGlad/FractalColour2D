#include "stdafx.h"
#include "bmp.h"
#include <fstream>
struct Node
{
  Vector2d pos;
  double radius;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

double scale = 750.0;
static double minLength = 2.0;

void saveSVG(const string &fileName, const vector<Node> &list)
{
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)scale << "\" height = \"" << (int)scale << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;

  for (int i = 0; i < (int)list.size(); i++)
  {
    double r2 = list[i].radius*1.1; 
    double r1 = list[i].radius*0.9;
//    double r2 = list[i].radius + 1.0;
//    double r1 = list[i].radius - 1.0;
    svg << "<circle cx = \"" << list[i].pos[0] << "\" cy = \"" << list[i].pos[1] << "\" r = \"" << r2 << "\" fill=\"black\" />" << endl;
    svg << "<circle cx = \"" << list[i].pos[0] << "\" cy = \"" << list[i].pos[1] << "\" r = \"" << r1 << "\" fill=\"white\" />" << endl;
  }
  svg << "</svg>" << endl;
  svg.close();
}

void split(vector<Node> &list, const Node &node, bool flip)
{
  list.push_back(node);
  if (node.radius <= minLength)
  {
    return;
  }
  // flip = !flip;
  double dir = flip ? -1.0 : 1.0;
  Node child1;
  double s = 1.151;
  child1.pos = node.pos + s*Vector2d(sqrt(3.0) / 2.0, dir*0.5)*node.radius;
  child1.radius = node.radius / 3.0;
  split(list, child1, flip);
  Node child2;
  child2.pos = node.pos + s*Vector2d(-sqrt(3.0) / 2.0, dir*0.5)*node.radius;
  child2.radius = node.radius / 3.0;
  split(list, child2, flip);
  Node child3;
  child3.pos = node.pos + s*Vector2d(0, -1.0*dir)*node.radius;
  child3.radius = node.radius / 3.0;
  split(list, child3, flip);
}

int _tmain(int argc, _TCHAR* argv[])
{
  Node base;
  base.pos = Vector2d(scale/2.0, 50 + scale/2.0);
  base.radius = 200.0;
  vector<Node> list;
  split(list, base, false);

  saveSVG("spongeo.svg", list);
}
