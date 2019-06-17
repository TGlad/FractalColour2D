#include <vector>
#include <fstream>
#include "stdafx.h"

double scale = 750.0;
struct Square
{
  Vector3d pos;
  double width;
  int layer;
  string colour;
};

void saveSVG(const string &fileName, vector<Square> &squares)
{
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)(2.0*scale) << "\" height = \"" << (int)scale << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;

  for (auto &square : squares)
  {
    Vector3d ps[4] = { square.pos, square.pos + Vector3d(square.width, 0, 0), square.pos + Vector3d(square.width, square.width, 0), square.pos + Vector3d(0, square.width, 0) };
    Vector3d mid = (ps[0] + ps[1] + ps[2] + ps[3]) / 4.0;
    double dist = max(abs(mid[0]), abs(mid[1]));
    double type = log(dist) / log(2.0);
    double pitch = 0.25;
    Vector3d side(1, 0, 0);
    Vector3d up(0.0, -sin(pitch), cos(pitch));
    Vector3d fwd(0.0, cos(pitch), sin(pitch));
    for (int i = 0; i < 5; i++)
    {
      Vector3d w = ps[i % 4];
      Vector3d p = w[0] * side + w[2] * up + w[1] * fwd;
      p /= (2.5 + p[1]);
      if (i==0)
        svg << "<path d = \"M " << scale*(p[0]+1.0) << " " << scale * (0.5 - p[2]);
      else
        svg << " L " << scale*(p[0] + 1.0) << " " << scale * (0.5 - p[2]);
    }
    string cols[] = { "black", "dimgrey", "grey", "darkgrey", "lightgrey", "white"};
    double wid = 1.0;
    string fill = cols[square.layer + 1];
    if (square.colour != "")
      fill = square.colour;
    svg << "\" opacity=\"0.9\" fill=\"" << fill << "\" stroke=\"" << cols[max(0, square.layer-1)] << "\" stroke-width=\"" << wid << "\" />\n";
  }
  svg << "</svg>" << endl;
  svg.close();
}


int _tmain(int argc, _TCHAR* argv[])
{
  double d = 0.0;
  double scale = 2.0;
  double length = 0.0;
  double numSegments = 1.0;
  for (int i = 0; i < 20; i++)
  {
    double s = (scale + sqrt(sqr(scale) + d*d - 1.0)) / (d + 1.0);
    d = (s*s*(d + 1.0) + d - 1.0) / (s*s*(d + 1.0) + 1.0 - d);
    numSegments *= 2.0;
    if (i > 0)
    {
      double newLength = d / numSegments;
      double h = sqrt(sqr(newLength) - sqr(length / 2.0));
      double curvature = h / sqr(length / 2.0);
      cout << "i: " << i << ", d: " << d << ", curv: " << curvature << endl;

 //     double r = 1.0 / (2.0*sqrt(3.0));
 //     double h1 = r - sqrt(r*r - sqr(length/2.0));
 //     cout << "h1: " << h << ", h2: " << h1 << endl;
    }
    else
      cout << "i: " << i << ", d: " << d << endl;
    length = d / numSegments;
  }

  exit(1);
  vector<Square> squares;
  int layer = 0;
  int nums = 8;
  double w = 0.1;
  for (int i = 0; i<4; i++)
  {
    for (int x = 0; x < nums; x++)
    {
      for (int y = 0; y < nums; y++)
      {
        Square square;
        square.pos[0] = -0.06 + ((double)x - ((double)nums)/2.0) * w;
        square.pos[1] = ((double)y - ((double)nums) / 2.0) * w;
        square.pos[2] = -0.35 + 0.17*(double)layer;
        square.layer = layer;
        square.width = w;
        square.colour = "";
  /*      if (layer == 0 && x >= 2 && x <= 3 && y >= 2 && y <= 3)
          square.colour = "purple";
        if (layer == 1 && x == 1 && y == 1)
          square.colour = "red";
        else if (layer == 1 && x < 3 && y < 3)
          square.colour = "orange";
        else if (layer == 2)
          square.colour = "yellow";*/
        squares.push_back(square);
      }
    }
    layer++;
    nums /= 2;
    w *= 2.0;
  }

  saveSVG("automaton.svg", squares);
}
