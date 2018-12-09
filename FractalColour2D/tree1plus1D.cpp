#include "stdafx.h"
#include "bmp.h"
#include <fstream>
static double gradient = 0.07;
struct Node
{
  Vector2d pos;
  double velocity;
  double width;
  double dir;
  double fwd;
  int flip;
  double time;
  vector<Node> children;
  void split();
  void draw(ofstream &svg);
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

double scale = 750.0;
Vector2d offset(0.55, 0);
static vector<Vector2d> leaves;
static double minLength = 0.0005;

void Node::draw(ofstream &svg)
{
  Vector2d x(-1, 0);
  Vector2d y = Vector2d(velocity, fwd) * time;
  Vector2d start = pos;

  Vector2d corners[] = { start - x*width, start + x*width, start + x*(width-gradient*time) + y, start - x*(width-gradient*time) + y};
  double gscale = 0.9*scale;
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

void Node::split()
{
  double scale1 = 0.7;
  time = (1.0 - scale1) * width / gradient;
  if (width <= minLength)
    return;
  Node child1, child2;

  child1.pos = pos + Vector2d(velocity*time, time*fwd);
  child1.width = width * scale1;
  child1.dir = -dir; 
  child1.velocity = velocity - dir*0.3;
  child1.fwd = fwd;
  child1.flip = flip+1;

  double scale = 0.4;
  child2.pos = pos + Vector2d(velocity*time + dir*width*(scale1 - scale), fwd*time);
  child2.width = width * scale;
  child2.velocity = velocity + dir*0.5;
  child2.flip = flip+1;
#define BIDIRECTIONAL
#if defined BIDIRECTIONAL
  child2.dir = -dir;
  child2.fwd = fwd;
  if (velocity > 0.0 && dir > 0.0)
    child2.fwd = -fwd, child2.dir = dir;
  if (velocity < 0.0 && dir < 0.0)
    child2.fwd = -fwd, child2.dir = dir;
//  child2.fwd = dir > 0 ? fwd : -fwd;
//  child2.dir = dir > 0 ? dir : dir;
#else
  child2.fwd = fwd;
  child2.dir = -dir;
#endif

  children.push_back(child1);
  children.push_back(child2);

  for (auto &c : children)
    c.split();
}

int _tmain(int argc, _TCHAR* argv[])
{
  Node base;
  base.velocity =0.15;
  base.width = 0.14;// 0.07;
  base.pos = Vector2d(0, -1.0);
  base.dir = 1.0;
  base.fwd = 1.0;
  base.flip = 0;

  base.split();

#if defined BIDIRECTIONAL
  saveSVG("tree1plus1Dbidirectional.svg", base);
#else
  saveSVG("tree1plus1D.svg", base);
#endif
}
