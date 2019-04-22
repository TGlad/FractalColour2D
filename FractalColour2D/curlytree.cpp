#include "stdafx.h"
#include "bmp.h"
#include <fstream>
static double gradient = 0.07;
struct Node
{
  Vector2d pos;
  double angle;
  double width;
  double dir;
  double fwd;
  int flip;
  double length;
  vector<Node> children;
  void split();
  void draw(ofstream &svg);
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

double scale = 700.0;
Vector2d offset(0.5, 0);
static vector<Vector2d> leaves;
static double minLength = 0.0006;

void Node::draw(ofstream &svg)
{
  Vector2d x(cos(angle), -sin(angle));
  Vector2d y = Vector2d(sin(angle), cos(angle)) * length;
  Vector2d start = pos;

  Vector2d corners[] = { start - x*width, start + x*width, start + x*(width-gradient*length) + y, start - x*(width-gradient*length) + y};
  double gscale = scale;
  svg << "<path d = \"M " << gscale*(corners[3][0] + offset[0]) << " " << scale-gscale*(corners[3][1] + offset[1]);
  for (int i = 0; i < 4; i++)
    svg << " L " << gscale*(corners[i][0] + offset[0]) << " " << scale-gscale*(corners[i][1] + offset[1]);
  svg << "\" fill=\"black\" stroke-width = \"0\" stroke=\"black\" />\n";

  for (auto &c : children)
    c.draw(svg);
}

void saveSVG(const string &fileName, Node &tree)
{
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)scale << "\" height = \"" << (int)scale << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;

  tree.draw(svg);
  svg << "</svg>" << endl;
  svg.close();
}
static double tscale = 0.7;
static double ang = 0.6;

void Node::split()
{
  length = (1.0 - tscale) * width / gradient;
  if (width <= minLength)
    return;
  Node child1, child2;

  child1.width = width * tscale;
  child1.dir = -dir; 
  child1.angle = angle - dir*ang;
  child1.pos = pos + length * Vector2d(sin(angle), cos(angle));
  child1.fwd = fwd;
  child1.flip = flip+1;

  child2.width = width * tscale;
  child2.angle = angle + dir*ang;
  child2.pos = pos + length * Vector2d(sin(angle), cos(angle));
  child2.flip = flip + 1;

  child2.fwd = fwd;
  child2.dir = -dir;

  children.push_back(child1);
  children.push_back(child2);

  for (auto &c : children)
    c.split();
}

int _tmain(int argc, _TCHAR* argv[])
{
  {
    Node base;
    base.angle = ang*3.0;
    base.width = 0.04 / (tscale*tscale*tscale);
    base.pos = Vector2d(-0.9, 0.1);
    base.dir = 1.0;
    base.fwd = 1.0;
    base.flip = 0;

    base.split();

    saveSVG("curlytree1.svg", base);
  }
  {
    Node base;
    base.angle = 0.0;
    base.width = 0.04;
    base.pos = Vector2d(0, 0.0);
    base.dir = 1.0;
    base.fwd = 1.0;
    base.flip = 0;

    base.split();

    saveSVG("curlytree2.svg", base);
  }
  {
    minLength = 0.006;
    Node base;
    base.angle = 0.0;
    base.width = 0.04;
    base.pos = Vector2d(0, 0.0);
    base.dir = 1.0;
    base.fwd = 1.0;
    base.flip = 0;

    base.split();

    saveSVG("curlytree3.svg", base);
  }
}
