#include "stdafx.h"
#include "bmp.h"
#include <sstream>
#include <fstream>
static double gradient = 0.0001;
struct Node
{
  Vector2d pos;
  Vector2d xAxis, yAxis;
  double nodeWidth;
  vector<Node> children;
  void split(int index);
  void draw(ofstream &svg, const Vector2d &origin, const Vector2d &xAx, const Vector2d &yAx);
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

static int width = 1024;
static int height = 750;
Vector2d offset(500.0, 600.0);
Vector2d minVec(1e10, 1e10);
Vector2d maxVec(-1e10, -1e10);
static vector<Vector2d> leaves;
static double minLength = 0.001;
static double area = 0.0;
void Node::draw(ofstream &svg, const Vector2d &origin, const Vector2d &xAx, const Vector2d &yAx)
{
  Vector2d minV(-1.26, 0.62);
  Vector2d maxV(1.36, 1.85);
  double mult = (maxV[0] - minV[0]);
  Vector2d start = origin + xAx*pos[0] + yAx*pos[1];
  Vector2d x = xAx*xAxis[0] + yAx*xAxis[1];
  Vector2d y = xAx*yAxis[0] + yAx*yAxis[1];
  double length = nodeWidth * 0.4;

  if (children.empty())
  {
    Vector2d corners[] = { start - x*nodeWidth*0.5, start + x*nodeWidth*0.5, start + x*0.5*nodeWidth + y*length, start + y*(1.5*length), start - x*0.5*nodeWidth + y*length };
    leaves.push_back(corners[2]);
    minVec[0] = min(minVec[0], corners[2][0]);
    minVec[1] = min(minVec[1], corners[2][1]);
    maxVec[0] = max(maxVec[0], corners[2][0]);
    maxVec[1] = max(maxVec[1], corners[2][1]);
    for (int i = 0; i < 5; i++)
      corners[i] = minV + (double)width * (corners[i] - minV) / mult;
    svg << "<path d = \"M " << corners[4][0] << " " << corners[4][1];
    for (int i = 0; i < 5; i++)
      svg << " L " << corners[i][0]  << " " << corners[i][1];
//    svg << "\" fill=\"grey\" stroke-width = \"0\" stroke=\"grey\" />\n";
    svg << "\" fill=\"black\" stroke-width = \"1\" stroke=\"black\" />\n";
  }
  Vector2d p = start + y*length;
  for (auto &c : children)
    c.draw(svg, p, x, y);
}

void saveSVGLeaves(const string &fileName)
{
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << width << "\" height = \"" << height << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;
  double gscale = 0.9*(double)height;

  for (auto &p : leaves)
    svg << "<circle cx = \"" << gscale*(p[0] + offset[0]) << "\" cy = \"" << (double)height - gscale*(p[1] + offset[1]) << "\" r = \"2\" stroke = \"black\" stroke-width = \"1\" fill = \"black\" />" << endl;
  svg << "</svg>" << endl;
  svg.close();
}
void saveSVG(const string &fileName, Node &tree)
{
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << width << "\" height = \"" << height << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;

  tree.draw(svg, Vector2d(0, 0), Vector2d(1, 0), Vector2d(0, 1));
  svg << "</svg>" << endl;
  svg.close();
}

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

void Node::split(int index)
{
  if (index == 0)
  {
    return;
  }
  Node child1, child2;

  double scale = 0.5;
  double ang = -0.5;
  child1.pos = nodeWidth*Vector2d(-0.4, 0.5);
  child1.yAxis = Vector2d(sin(ang), cos(ang));
  child1.xAxis = Vector2d(cos(ang), -sin(ang));
  child1.nodeWidth = nodeWidth * scale;

  ang = 0.7;
  child2.pos = nodeWidth*Vector2d(0.5, 0.35);
  child2.yAxis = Vector2d(sin(ang), cos(ang));
  child2.xAxis = Vector2d(cos(ang), -sin(ang));
  child2.nodeWidth = nodeWidth * scale;

  children.push_back(child1);
  children.push_back(child2);
  for (auto &c : children)
    c.split(index - 1);
}



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
  base.xAxis = Vector2d(1, 0);
  base.yAxis = Vector2d(0, 1);
  base.nodeWidth = 1.0;
  base.pos = Vector2d(0, 0.05);

  base.split(12);
  saveSVG("substitutionrule_thing.svg", base);
  cout << "min: " << minVec.transpose() << ", max: " << maxVec.transpose() << endl;

  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  vector<double> graph;
  vector<double> dimgraph; 
  double c = 0;
//  for (double r = 0.25; r > 1.0 / 8.0; r /= 1.03)
  for (double r = 0.191; r > 0.0191; r /= 2.0)
  {
    memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey
    double w = 2e-3; 
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
        for (auto &leaf : leaves)
        {
          if ((p - leaf).squaredNorm() < radius2)
          {
            Nmink++;
            double W = 2.0*w;
            putpixel(out, Vector2i(0.01 + (x - minX) / W, 0.01 + (y - minY) / W), 192);
            break;
          }
        }
      }
    }
    for (auto &leaf : leaves)
    {
      double W = 2.0*w;
      for (int i = -1; i < 2; i++)
        for (int j = -1; j < 2; j++)
          putpixel(out, Vector2i((double)i + 0.01 + (leaf[0] - minX) / W, (double) j + 0.01 + (leaf[1] - minY) / W), 0);
    }
    Nmink *= w*w; // now a volume
    Nmink *= m*m; // now a Minkowski number
    double d = 1.0;
    if (Nmink / pow(m, d) < c)
      cout << "blah at radius " << r << endl;
    c = Nmink / pow(m, d);
    cout << "radius: " << r << ", c: " << c << endl;
    cout << "log(N): " << log(Nmink) << ", log m: " << log(m) << ", ratio: " << log(Nmink) / log(m) << endl;

    graph.push_back(c);
    dimgraph.push_back(log(Nmink));


    BYTE* buf = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    LPCTSTR file = L"blah.bmp";
    SaveBMP(buf, width, height, s2, file);
    delete[] buf;
  }
  double dif = graph.back() - graph[0];
  for (int i = 0; i < (int)graph.size(); i++)
    graph[i] -= dif * (double)i / (double)(graph.size() - 1);
  cout << "dif " << dif << ", " << graph[0] << ", " << graph.back() << endl;
  graph.insert(graph.end(), graph.begin(), graph.end());
  double yscale = 2.0*(double)height;
  double minY = 9.6;
//  drawGraph("substitutionGraph.svg", graph, yscale, minY);
  int dimsize = dimgraph.size();
  double diff = dimgraph.back() - dimgraph[0];
  for (int i = 1; i < dimsize; i++)
    dimgraph.push_back(dimgraph[i] + diff);
  yscale = (double)height / 2.5;
  minY = 3.0;
//  drawGraph("substitutionDimGraph.svg", dimgraph, yscale, minY);
}
