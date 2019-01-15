#include "stdafx.h"
#include "bmp.h"
#include <fstream>
struct Node
{
  Vector2d left, right;
  Vector2d getTop() const
  {
    Vector2d mid = (left + right) / 2.0;
    Vector2d dif = (right - left) * sqrt(3.0) / 2.0;
    return mid + Vector2d(-dif[1], dif[0]);
  }
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

double scale = 750.0;
Vector2d offset(0.5, 0);
static double minLength = 0.001;
static vector<Vector2d> leaves;

void saveTree(const string &fileName, vector<Node> &tree)
{
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)scale << "\" height = \"" << (int)scale << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;

  for (auto &node : tree)
  {
    Vector2d top = node.getTop();
    svg << "<path d = \"M " << scale*(node.left[0] + offset[0]) << " " << scale - scale*(node.left[1] + offset[1]);
    svg << " L " << scale*(top[0] + offset[0]) << " " << scale - scale*(top[1] + offset[1]);
    svg << " L " << scale*(node.right[0] + offset[0]) << " " << scale - scale*(node.right[1] + offset[1]);
    svg << "\" fill=\"black\" stroke-width = \"1\" stroke=\"black\" />\n";
  }
  svg << "</svg>" << endl;
  svg.close();
}

void saveOutline(const string &fileName, vector<Node> &tree)
{
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)scale << "\" height = \"" << (int)scale << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;

  for (auto &node : tree)
  {
    Vector2d top = node.getTop();
    Vector2d midLeft = (node.left + top) / 2.0;
    Vector2d midRight = (node.right + top) / 2.0;
    svg << "<path d = \"M " << scale*(node.left[0] + offset[0]) << " " << scale - scale*(node.left[1] + offset[1]);
    svg << " L " << scale*(midLeft[0] + offset[0]) << " " << scale - scale*(midLeft[1] + offset[1]);
    svg << "\" fill=\"none\" stroke-width = \"2\" stroke=\"black\" />\n";
    svg << "<path d = \"M " << scale*(node.right[0] + offset[0]) << " " << scale - scale*(node.right[1] + offset[1]);
    svg << " L " << scale*(midRight[0] + offset[0]) << " " << scale - scale*(midRight[1] + offset[1]);
    svg << "\" fill=\"none\" stroke-width = \"2\" stroke=\"black\" />\n";
  }
  svg << "</svg>" << endl;
  svg.close();
}

void saveLeaves(const string &fileName, vector<Node> &tree)
{
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)scale << "\" height = \"" << (int)scale << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;

  for (auto &p : leaves)
  {
    svg << "<circle cx = \"" << scale*(p[0] + offset[0]) << "\" cy = \"" << scale - scale*(p[1] + offset[1]) << "\" r = \"2\" stroke = \"black\" stroke-width = \"0\" fill = \"black\" />" << endl;
  }
  svg << "</svg>" << endl;
  svg.close();
}

void saveSkeleton(const string &fileName, vector<Node> &tree)
{
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)scale << "\" height = \"" << (int)scale << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;

  for (auto &node : tree)
  {
    Vector2d top = node.getTop();
    Vector2d base = (node.left + node.right) / 2.0;
    Vector2d peak = base + (top - base)*2.0 / 3.0;
    base = base - (top - base)*1.0 / 3.0;
    svg << "<path d = \"M " << scale*(base[0] + offset[0]) << " " << scale - scale*(base[1] + offset[1]);
    svg << " L " << scale*(peak[0] + offset[0]) << " " << scale - scale*(peak[1] + offset[1]);
    svg << "\" fill=\"none\" stroke-width = \"2\" stroke=\"black\" />\n";
  }
  svg << "</svg>" << endl;
  svg.close();
}

void saveBranchPoints(const string &fileName, vector<Node> &tree)
{
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)scale << "\" height = \"" << (int)scale << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;

  for (auto &node : tree)
  {
    Vector2d top = node.getTop();
    Vector2d base = (node.left + node.right) / 2.0;
    Vector2d p = base + (top - base)*2.0 / 3.0;
    svg << "<circle cx = \"" << scale*(p[0] + offset[0]) << "\" cy = \"" << scale - scale*(p[1] + offset[1]) << "\" r = \"2\" stroke = \"black\" stroke-width = \"0\" fill = \"black\" />" << endl;
  }
  svg << "</svg>" << endl;
  svg.close();
}

void recurse(vector<Node> &list, const Node &node)
{
  list.push_back(node);
  double length = (node.left - node.right).norm();
  if (length <= minLength)
  {
    leaves.push_back(node.getTop());
    return;
  }

  Vector2d top = node.getTop();
  Node leftNode;
  leftNode.left = (node.left + top) / 2.0;
  leftNode.right = top;
  recurse(list, leftNode);
  Node rightNode;
  rightNode.left = top; 
  rightNode.right = (node.right + top) / 2.0;
  recurse(list, rightNode);
}

int _tmain(int argc, _TCHAR* argv[])
{
  Node base;
  base.left = Vector2d(-0.32, 0.0);
  base.right = Vector2d(0.32, 0.0);
  vector<Node> list;
  recurse(list, base);

  saveTree("triangletree1.svg", list);
  saveOutline("triangletree2.svg", list);
  saveLeaves("triangletree3.svg", list); \
  saveSkeleton("triangletree4.svg", list); \
  saveBranchPoints("triangletree5.svg", list); \
}
