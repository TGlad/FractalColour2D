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
//#define MENGER
//#define VISCEK
//#define MUSHROOM
//#define TRITREE_SOLID
//#define TRITREE_OUTLINE
//#define TRITREE_LEAVES
//#define TRITREE_SKELETON
//#define TRITREE_BRANCHES
//#define BINARY_TREE1
#define BINARY_TREE2
//#define BINARY_TREE3
//#define BINARY_TREE_LEAVES

static int width = 1024;
static int height = 750;
Vector2d offset(500.0, 600.0);
static vector<Vector2d> leaves;
#if defined TRITREE_SOLID || defined TRITREE_OUTLINE || defined TRITREE_SKELETON
static double minLength = 0.01;
#elif defined BINARY_TREE1 || defined BINARY_TREE2 || defined BINARY_TREE3 
static double minLength = 0.0025;
#else
static double minLength = 0.0015;
#endif
static double area = 0.0;

void drawGraph(const string &fileName, vector<double> &graph, double yscale, double minY)
{
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << width << "\" height = \"" << height << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;

  double xscale = 0.9*(double)width;
  svg << "<path d = \"M " << 0 << " " << (double)height - yscale*(graph[0] - minY);
  for (int i = 1; i < (int)graph.size(); i++)
    svg << " L " << xscale * (double)i / (double)graph.size() << " " << (double)height - yscale*(graph[i] - minY);
  svg << "\" fill = \"none\" stroke-width = \"2\" stroke=\"black\" />\n";
  svg << "</svg>" << endl;
  svg.close();
}
static Vector2d minVec(1e10, 1e10);
static Vector2d maxVec(-1e10, -1e10);

#if defined MENGER
void recurse(vector<Node> &list, Node &node)
{
  if (node.length < minLength)
  {
    minVec = Vector2d(min(minVec[0], node.pos[0]), min(minVec[1], node.pos[1]));
    maxVec = Vector2d(max(maxVec[0], node.pos[0]), max(maxVec[1], node.pos[1]));
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
#elif defined VISCEK
void recurse(vector<Node> &list, Node &node)
{
  if (node.length < minLength)
  {
    minVec = Vector2d(min(minVec[0], node.pos[0]), min(minVec[1], node.pos[1]));
    maxVec = Vector2d(max(maxVec[0], node.pos[0]), max(maxVec[1], node.pos[1]));
    list.push_back(node);
    return;
  }
  Node child;
  child.length = node.length / 3.0;
  child.pos = node.pos + Vector2d(-child.length, 0);
  recurse(list, child);
  child.pos = node.pos + Vector2d(0, 0);
  recurse(list, child);
  child.pos = node.pos + Vector2d(child.length, 0);
  recurse(list, child);
  child.pos = node.pos + Vector2d(0, -child.length);
  recurse(list, child);
  child.pos = node.pos + Vector2d(0, child.length);
  recurse(list, child);
}
#elif defined MUSHROOM
void recurse(vector<Node> &list, Node &node)
{
  if (node.length < minLength)
  {
    minVec = Vector2d(min(minVec[0], node.pos[0]), min(minVec[1], node.pos[1]));
    maxVec = Vector2d(max(maxVec[0], node.pos[0]), max(maxVec[1], node.pos[1]));
    list.push_back(node);
    return;
  }
  Vector2d side(-node.up[1], node.up[0]);
  Node child;
  child.length = node.length / 4.0;
  child.pos = node.pos - side * child.length*1.5;
  child.up = node.up;
  recurse(list, child);
  child.pos = node.pos - side * child.length + node.up*child.length*0.5;
  child.up = side;
  recurse(list, child);
  child.pos = node.pos - side * child.length + node.up*child.length*1.5;
  child.up = -side;
  recurse(list, child);

  child.pos = node.pos - side * child.length*0.5 + node.up*child.length*2.0;
  child.up = node.up;
  recurse(list, child);
  child.pos = node.pos + side * child.length*0.5 + node.up*child.length*2.0;
  child.up = node.up;
  recurse(list, child);

  child.pos = node.pos + side * child.length*1.5;
  child.up = node.up;
  recurse(list, child);
  child.pos = node.pos + side * child.length + node.up*child.length*0.5;
  child.up = -side;
  recurse(list, child);
  child.pos = node.pos + side * child.length + node.up*child.length*1.5;
  child.up = side;
  recurse(list, child);
}
#elif defined TRITREE_SOLID || defined TRITREE_OUTLINE || defined TRITREE_LEAVES || defined TRITREE_SKELETON || defined TRITREE_BRANCHES
static const double cos30 = cos(pi / 6.0);
static const double rt75 = sqrt(0.75);
void recurse(vector<Node> &list, Node &node)
{
  minVec = Vector2d(min(minVec[0], node.pos[0]), min(minVec[1], node.pos[1]));
  maxVec = Vector2d(max(maxVec[0], node.pos[0]), max(maxVec[1], node.pos[1]));
#if defined TRITREE_LEAVES
  if (node.length < minLength)
  {
    list.push_back(node);
    return;
  }
#elif defined TRITREE_BRANCHES
  Node node2 = node;
  node2.pos += node.up*node.length*rt75 / 6.0;
  list.push_back(node2);
  if (node.length < minLength)
    return;
#else
  list.push_back(node);
  if (node.length < minLength)
    return;
#endif
  Vector2d side(-node.up[1], node.up[0]);


  Node child;
  child.length = node.length / 2.0;
  Vector2d top = node.pos + node.up*rt75*node.length;

  child.pos = (node.pos + 0.5*side*node.length) * 0.25 + top*0.75;
  child.up = node.up * 0.5 + side*cos30;
  recurse(list, child);
  child.pos = (node.pos - 0.5*side*node.length) * 0.25 + top*0.75;
  child.up = node.up * 0.5 - side*cos30;
  recurse(list, child);
}
#elif defined BINARY_TREE1 || defined BINARY_TREE2 || defined BINARY_TREE3 || defined BINARY_TREE_LEAVES
static const double rt5 = sqrt(0.5);
void recurse(vector<Node> &list, Node &node)
{
  minVec = Vector2d(min(minVec[0], node.pos[0]), min(minVec[1], node.pos[1]));
  maxVec = Vector2d(max(maxVec[0], node.pos[0]), max(maxVec[1], node.pos[1]));
#if defined BINARY_TREE_LEAVES
  if (node.length < minLength)
  {
    list.push_back(node);
    return;
  }
#else
  list.push_back(node);
  if (node.length < minLength)
    return;
#endif
  Vector2d side(-node.up[1], node.up[0]);


  Node child;
#if defined BINARY_TREE1
  child.length = node.length * 5.0 / 12.0;
#elif defined BINARY_TREE2
  child.length = node.length * 6.0 / 12.0;
#elif defined BINARY_TREE3
  child.length = node.length * 7.0 / 12.0;
#endif
  Vector2d top = node.pos + node.up*node.length;

  child.pos = top;
  child.up = node.up * rt5 + side*rt5;
  recurse(list, child);
  child.pos = top;
  child.up = node.up * rt5 - side*rt5;
  recurse(list, child);
}
#endif

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
  base.length = 1.0;
  vector<Node> list;
  recurse(list, base);
  cout << "min: " << minVec.transpose() << ", max: " << maxVec.transpose() << endl;

  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  vector<double> graph;
  double c = 0;
  {
    memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey
#if defined TRITREE_SOLID
    double r = 0.002;
    double w = 0.002; 
#else
    double r = 0.1;
    double w = 0.0015; 
#endif
    double m = 1.0 / r;
    double minX = minVec[0] - r;
    double minY = minVec[1] - r;
    double maxX = maxVec[0] + r;
    double maxY = maxVec[1] + r;
    double radius2 = r*r;
    double Nmink = 0;
    Vector2d di(sin(22.5*pi / 180.0), cos(22.5*pi / 180.0));
    for (double x = minX; x < maxX; x += w)
    {
      for (double y = minY; y < maxY; y += w)
      {
        Vector2d p(x, y);
        for (auto &node : list)
        {
          Vector2d nodePos = node.pos;
#if defined TRITREE_SOLID // intersect with triangle
          Vector2d side(-node.up[1], node.up[0]);
          Vector2d ps[3] = { node.pos - 0.5*side*node.length, node.pos + 0.5*side*node.length, node.pos + node.up*rt75*node.length };
          bool inside = true;
          Vector2d posInTriangle = p;
          for (int i = 0; i < 3; i++)
          {
            Vector2d dir = ps[(i + 1) % 3] - ps[i];
            Vector2d side(-dir[1], dir[0]);
            if ((p - ps[i]).dot(side) > 0.0)
            {
              inside = false;
              //         double d = (p - ps[i]).dot(dir) / dir.squaredNorm();
              //         nodePos = ps[i] + dir*max(0.0, min(d, 1.0));
            }
          }
          if (inside)
#elif defined TRITREE_OUTLINE // intersect with two lines
          Vector2d side(-node.up[1], node.up[0]);
          Vector2d ps[3] = { node.pos - 0.5*side*node.length, node.pos + 0.5*side*node.length, node.pos + node.up*rt75*node.length };
          Vector2d lines[2][2] = { { ps[0], ( ps[0] + ps[2]) / 2.0 }, { ps[1], (ps[1] + ps[2]) / 2.0 } };
          double minDist = 1e10;
          for (int i = 0; i < 2; i++)
          {
            Vector2d dir = lines[i][1] - lines[i][0];
            double d = (p - lines[i][0]).dot(dir) / dir.squaredNorm();
            Vector2d nearest = lines[i][0] + dir*max(0.0, min(d, 1.0));
            minDist = min(minDist, (nearest - p).squaredNorm());
          }
          if (minDist < radius2)
#elif defined TRITREE_SKELETON // intersect with one lines
          Vector2d side(-node.up[1], node.up[0]);
          double k = rt75*1.0/3.0;
          double j = rt75*2.0/3.0;
          Vector2d ps[2] = { node.pos - k*node.up*node.length, node.pos + j*node.up*node.length };
          Vector2d dir = ps[1] - ps[0];
          double d = (p - ps[0]).dot(dir) / dir.squaredNorm();
          Vector2d nearest = ps[0] + dir*max(0.0, min(d, 1.0));
          double minDist = (nearest - p).squaredNorm();
          if (minDist < radius2)
#elif defined BINARY_TREE1 || defined BINARY_TREE2 || defined BINARY_TREE3 // intersect with one lines
          Vector2d p(x, y);
          if (y < 0.5 && (y < x || p.dot(di) < 0.0))
            continue;
          Vector2d side(-node.up[1], node.up[0]);
          Vector2d ps[2] = { node.pos, node.pos + node.up*node.length };
          Vector2d dir = ps[1] - ps[0];
          double d = (p - ps[0]).dot(dir) / dir.squaredNorm();
          Vector2d nearest = ps[0] + dir*max(0.0, min(d, 1.0));
          double minDist = (nearest - p).squaredNorm();
          if (minDist < radius2)
#else
          if ((p - nodePos).squaredNorm() < radius2)
#endif
          {
            Nmink++;
            double W = w*2;
            putpixel(out, Vector2i(0.01 + (x - minX) / W, 0.01 + (y - minY) / W), 192);
            break;
          }
        }
      }
    }
    for (auto &node : list)
    {
      double W = w;
      putpixel(out, Vector2i(0.01 + (node.pos[0] - minX) / W, 0.01 + (node.pos[1] - minY) / W), 0);
    }
    Nmink *= w*w; // now a volume
    Nmink *= m*m; // now a Minkowski number
#if defined MENGER
    double d = log(8.0) / log(3.0); 
#elif defined VISCEK
    double d = log(5.0) / log(3.0);
#elif defined MUSHROOM
    double d = 1.5;
#elif defined TRITREE_SOLID
    double d = 2.0;
#elif defined TRITREE_LEAVES || defined TRITREE_BRANCHES || defined BINARY_TREE_LEAVES
    double d = 1.0;
#endif
#if defined TRITREE_OUTLINE || defined TRITREE_SKELETON  
    c = Nmink / (m * log(m));
    cout << "radius: " << r << ", c: " << c << endl;
#elif defined BINARY_TREE1
    cout << "N" << Nmink;
    cout << "r scaling process:" << endl;
    double d = log(2.0) / log(12.0 / 5.0);
    for (int i = 0; i < 10; i++)
    {
      Nmink *= 2.0;
      m *= 12.0 / 5.0;
      Nmink += 2.0*m;
      Nmink -= 1.5*tan(22.5*pi/180.0)/(m*m);
      Nmink -= 0.5*tan(45.0*pi / 180.0) / (m*m);
      cout << "m: " << m << ", Nmink: " << Nmink << endl;
      cout << "estimated c: " << (12.0*m - Nmink) / pow(m, d) << endl;
    }

#elif defined BINARY_TREE2
    cout << "N" << Nmink;
    cout << "r scaling process:" << endl;
    double d = log(2.0) / log(12.0 / 6.0);
    for (int i = 0; i < 10; i++)
    {
      Nmink *= 2.0;
      m *= 12.0 / 6.0;
      Nmink += 2.0*m;
      Nmink -= 1.5*tan(22.5*pi / 180.0) / (m*m);
      Nmink -= 0.5*tan(45.0*pi / 180.0) / (m*m);
      cout << "m: " << m << ", Nmink: " << Nmink << endl;
      double est = m*log(m) * 2.0 / log(2.0);
      cout << "estimated c: " << (Nmink - est) / pow(m, d) << endl;
    }
#elif defined BINARY_TREE3
      cout << "N" << Nmink;
    cout << "r scaling process:" << endl;
    double d = log(2.0) / log(12.0 / 7.0);
    for (int i = 0; i < 10; i++)
    {
      Nmink *= 2.0;
      m *= 12.0 / 7.0;
      Nmink += 2.0*m;
      Nmink -= 1.5*tan(22.5*pi / 180.0) / (m*m);
      Nmink -= 0.5*tan(45.0*pi / 180.0) / (m*m);
      cout << "m: " << m << ", Nmink: " << Nmink << endl;
      cout << "estimated c: " << (12.0*m + Nmink) / pow(m, d) << endl;
    }
#else
    c = Nmink / pow(m, d);
    cout << "radius: " << r << ", d: " << d << ", c: " << c << endl;
#endif

    graph.push_back(c);

    BYTE* buf = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    LPCTSTR file = L"blah.bmp";
    SaveBMP(buf, width, height, s2, file);
    delete[] buf;
  }
#if 0
  double dif = graph.back() - graph[0];
  for (int i = 0; i < (int)graph.size(); i++)
  {
    graph[i] -= dif * (double)i / (double)(graph.size() - 1);
  }
  cout << "dif " << dif << ", " << graph[0] << ", " << graph.back() << endl;
  
  graph.insert(graph.end(), graph.begin(), graph.end());
  double yscale = 2.0*(double)height;
  double minY = 9.6;

  drawGraph("vicsekGraph.svg", graph, yscale, minY);
#endif
}
