#include "stdafx.h"
#include "bmp.h"
#include <fstream>
static double gradient = 0.0001;
struct Node
{
  Vector2d pos;
  Vector2d xAxis, yAxis;
  double width;
  vector<Node> children;
  void split(int index);
  void draw(ofstream &svg, const Vector2d &origin, const Vector2d &xAx, const Vector2d &yAx);
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

double scale = 750.0;
Vector2d offset(0.55, 0);
static vector<Vector2d> leaves;
static double minLength = 0.001;
static double area = 0.0;
void Node::draw(ofstream &svg, const Vector2d &origin, const Vector2d &xAx, const Vector2d &yAx)
{
  Vector2d start = origin + xAx*pos[0] + yAx*pos[1];
  Vector2d x = xAx*xAxis[0] + yAx*xAxis[1];
  Vector2d y = xAx*yAxis[0] + yAx*yAxis[1];
  double length = width * 0.4;

  if (children.empty())
  {
    Vector2d corners[] = { start - x*width*0.5, start + x*width*0.5, start + x*0.5*width + y*length, start + y*(1.5*length), start - x*0.5*width + y*length };
    double gscale = 0.9*scale;
    svg << "<path d = \"M " << gscale*(corners[4][0] + offset[0]) << " " << scale - gscale*(corners[4][1] + offset[1]);
    for (int i = 0; i < 5; i++)
      svg << " L " << gscale*(corners[i][0] + offset[0]) << " " << scale - gscale*(corners[i][1] + offset[1]);
//    svg << "\" fill=\"grey\" stroke-width = \"0\" stroke=\"grey\" />\n";
    svg << "\" fill=\"black\" stroke-width = \"1\" stroke=\"black\" />\n";
  }
  Vector2d p = start + y*length;
  for (auto &c : children)
    c.draw(svg, p, x, y);
 // if (children.empty())
 //   leaves.push_back(corners[2]); // pos is in local space!
}

void saveSVGLeaves(const string &fileName)
{
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)scale << "\" height = \"" << (int)scale << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;
  double gscale = 0.9*scale;

  for (auto &p : leaves)
    svg << "<circle cx = \"" << gscale*(p[0] + offset[0]) << "\" cy = \"" << scale - gscale*(p[1] + offset[1]) << "\" r = \"2\" stroke = \"black\" stroke-width = \"1\" fill = \"black\" />" << endl;
  svg << "</svg>" << endl;
  svg.close();
}
void saveSVG(const string &fileName, Node &tree)
{
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)scale << "\" height = \"" << (int)scale << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;

  tree.draw(svg, Vector2d(0, 0), Vector2d(1, 0), Vector2d(0, 1));
  svg << "</svg>" << endl;
  svg.close();
}

void Node::split(int index)
{
  if (index == 0)
  {
    return;
  }
  Node child1, child2;

  double scale = 0.5;
  double ang = -0.5;
  child1.pos = width*Vector2d(-0.4, 0.5);
  child1.yAxis = Vector2d(sin(ang), cos(ang));
  child1.xAxis = Vector2d(cos(ang), -sin(ang));
  child1.width = width * scale;

  ang = 0.7;
  child2.pos = width*Vector2d(0.5, 0.35);
  child2.yAxis = Vector2d(sin(ang), cos(ang));
  child2.xAxis = Vector2d(cos(ang), -sin(ang));
  child2.width = width * scale;

  children.push_back(child1);
  children.push_back(child2);
  for (auto &c : children)
    c.split(index - 1);
}

int _tmain(int argc, _TCHAR* argv[])
{
  Node base;
  base.xAxis = Vector2d(1, 0);
  base.yAxis = Vector2d(0, 1);
  base.width = 0.35; 
  base.pos = Vector2d(0, 0.05);

//  base.split(8);
//  saveSVG("multiplicationrule_example2.svg", base);
  base.split(7);
  saveSVG("substitutionrule_example3.svg", base);
}
