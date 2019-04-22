#include <vector>
#include <fstream>
#include "stdafx.h"
static double gradient = 0.1;

struct Section
{
  Vector2d pos, peak;
  Vector2d xAxis, yAxis;
  double width, width2, length;
  double dir;
  vector<Section> children;
  void split(int level);
  void draw(ofstream &svg, const Vector2d &origin, const Vector2d &xAx, const Vector2d &yAx);
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

double scale = 750.0;
Vector2d offset(0.0, 0);
static double xVal = 0;
static double weight = 0;
void Section::draw(ofstream &svg, const Vector2d &origin, const Vector2d &xAx, const Vector2d &yAx)
{
  Vector2d start = origin + xAx*pos[0] + yAx*pos[1];
  Vector2d x = xAx*xAxis[0] + yAx*xAxis[1];
  Vector2d y = xAx*yAxis[0] + yAx*yAxis[1];

  Vector2d corners[] = { start - x*width*0.5, start + x*width*0.5, start + x*0.5*width2 + y*length, start + x*peak[0] + y*(peak[1] + length), start - x*0.5*width2 + y*length };
  double gscale = 0.9*scale;
  svg << "<path d = \"M " << gscale*(corners[4][0] + offset[0]) << " " << scale-gscale*(corners[4][1] + offset[1]);
  for (int i = 0; i < 5; i++)
    svg << " L " << gscale*(corners[i][0] + offset[0]) << " " << scale-gscale*(corners[i][1] + offset[1]);
  svg << "\" fill=\"black\" stroke=\"black\" stroke-width=\"1\" />\n"; \

  Vector2d p = start + y*length;
  for (auto &c : children)
    c.draw(svg, p, x, y);
}

void saveSVG(const string &fileName, vector<Section> &forest)
{
  xVal = 0;
  weight = 0;
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)(2.0*scale) << "\" height = \"" << (int)scale << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;

  for (auto &tree: forest)
    tree.draw(svg, Vector2d(0,0), Vector2d(1,0), Vector2d(0,1));
  svg << "</svg>" << endl;
  svg.close();
}

static double cantorScale = 0.2;// 1.0 / 3.0;// 3.0 / 5.0 - 0.08;
static double p = 1.25;// 1.55;// 2.0;
void Section::split(int level)
{
  double minLength = 0.005;
  if (length <= minLength)
    return;
  Section child1, child2;

  double w2 = pow(width, p);
  double w2b = pow(width2, p);
  double areaLoss = (w2 - w2b) * cantorScale;

  double a = pow(areaLoss, 1.0/p);
  double b = pow((w2 + w2b - areaLoss) / 2.0, 1.0 / p);
  double c = pow((w2 + w2b + areaLoss) / 2.0, 1.0 / p);
  double theta = acos((sqr(b) + sqr(c) - sqr(a)) / (2.0*b*c));
  Vector2d apex(c*0.5 - b*cos(theta), b*sin(theta));
  Vector2d B(-c/2.0, 0);
  Vector2d A(c/2.0, 0);
  Vector2d fromB = apex - B;
  Vector2d fromA = apex - A;
//  length *= pow((1.0 - cantorScale) / 2.0, 1.0/p);
  length = (length - a / 2.0) / 2.0;

  child1.pos = B + 0.5*fromB;
  child1.yAxis = Vector2d(-fromB[1], fromB[0]).normalized();
  child1.xAxis = fromB.normalized();
  child1.pos[0] *= dir; child1.yAxis[0] *= dir; child1.xAxis[1] *= dir;
  child1.length = a / gradient;
  child1.width = a;
  child1.width2 = 0;
  child1.dir = -dir; // remove negation for curly
  child1.peak = Vector2d(0, 0);

  child2.pos = A + 0.5*fromA;
  child2.yAxis = Vector2d(fromA[1], -fromA[0]).normalized();
  child2.xAxis = -fromA.normalized();
  child2.pos[0] *= dir; child2.yAxis[0] *= dir; child2.xAxis[1] *= dir;
  child2.length = length;
  child2.width = b;
  child2.width2 = width2;
  child2.children = children;
  child2.dir = -dir;
  child2.peak = peak;

  peak = apex; 
  peak[0] *= dir;
  width2 = c;

  children.clear();
  if (level > 3)
  {
    children.push_back(child1);
    children.push_back(child2);
  }

  for (auto &c : children)
    c.split(level+1);
  dir = -dir;
  split(level+1);
}

int _tmain(int argc, _TCHAR* argv[])
{
  gradient = (1.0 - 2.0*pow((1.0 - cantorScale) / 2.0, 1.0 / p)) / (0.5*pow(cantorScale, 1.0 / p));
  if (gradient < 0.0)
  {
    cout << "gradient less than zero." << endl;
    exit(1);
  }
  vector<Section> forest;
  Section base;
  base.xAxis = Vector2d(0, 1);
  base.yAxis = Vector2d(-1, 0);
  base.length = 50.0;
  base.width = base.length * gradient;
  base.width2 = 0;
  base.pos = Vector2d(2.3, -6.6); //Vector2d(2.0, -0.7);
  base.dir = 1.0;
  base.split(0);
  forest.push_back(base);

  saveSVG("forest3.svg", forest);
}
