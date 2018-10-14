#include "stdafx.h"
#include "bmp.h"
#include <fstream>
static double gradient = 0.15;
struct Node
{
  Vector2d pos;
  double angle;
  Vector2d xAxis() const { return Vector2d(cos(angle), sin(angle)); }
  Vector2d yAxis() const { return Vector2d(-sin(angle), cos(angle)); }
  double length;
  bool flip;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

double scale = 750.0;
Vector2d offset(0.6, 0);
static double widthPerLength = 0.2;

void saveSVG(const string &fileName, const vector<Node> &tree)
{
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)scale << "\" height = \"" << (int)scale << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;

  for (auto &node : tree)
  {
    Vector2d start = node.pos;
    Vector2d x = node.xAxis() * node.length * widthPerLength * 0.5;
    Vector2d y = node.yAxis() * node.length * 1.1;

    Vector2d corners[] = { start - x, start + x, start + 0.8*x + y, start - 0.8*x + y };
    double gscale = 0.9*scale;
    svg << "<path d = \"M " << gscale*(corners[3][0] + offset[0]) << " " << scale - gscale*(corners[3][1] + offset[1]);
    for (int i = 0; i < 4; i++)
      svg << " L " << gscale*(corners[i][0] + offset[0]) << " " << scale - gscale*(corners[i][1] + offset[1]);
    svg << "\" fill=\"black\" stroke-width = \"1\" stroke=\"black\" />\n";
  }

  svg << "</svg>" << endl;
  svg.close();
}
double time = 124.0;
double omega = 2.0;

void buildTree(vector<Node> &tree, const Node &node)
{
  tree.push_back(node);
  if (node.length < 0.005)
    return;

//  double scale1 = 0.75;  double angle1 = -0.4;  double scale2 = 0.53828;  double angle2 = 0.7;
  double scale1 = 0.75;  double angle1 = -0.4;  double scale2 = 0.53828;  double angle2 = 0.6;
  double sway = 0.25;
  if (node.flip)
  {
    angle1 *= -1;
    angle2 *= -1;
    swap(angle1, angle2);
    swap(scale1, scale2);
    sway = -sway;
  }

  Node child1;
  child1.pos = node.pos + node.yAxis() * node.length;
  child1.length = node.length * scale1;
  child1.angle = node.angle + angle1 + sway * sin(time*omega / child1.length);
  child1.flip = !node.flip;
  buildTree(tree, child1);

  Node child2;
  child2.pos = node.pos + node.yAxis() * node.length;
  child2.length = node.length * scale2;
  child2.angle = node.angle + angle2 + sway * sin(time*omega / child1.length);
  child2.flip = !node.flip;
  buildTree(tree, child2);
}

int _tmain(int argc, _TCHAR* argv[])
{
  Node base;
  base.angle = 0.2;
  base.length = 0.28;
  base.pos = Vector2d(0, 0);
  base.flip = false;

  vector<Node> tree;
  buildTree(tree, base);

  saveSVG("asymmetric_tree.svg", tree);
}
