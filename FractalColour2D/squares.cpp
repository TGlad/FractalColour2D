#include <vector>
#include <fstream>
#include "stdafx.h"

double scale = 750.0;
struct Square
{
  Vector3d pos;
  double width;
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
    double pitch = 1.0;
    Vector3d side(1, 0, 0);
    Vector3d up(0.0, -sin(pitch), cos(pitch));
    Vector3d fwd(0.0, cos(pitch), sin(pitch));
    for (int i = 0; i < 5; i++)
    {
      Vector3d w = ps[i % 4];
      Vector3d p = w[0] * side + w[2] * up + w[1] * fwd;
      p /= (1.5 + p[1]);
      if (i==0)
        svg << "<path d = \"M " << scale*(p[0]+1.0) << " " << scale * (0.5 - p[2]);
      else
        svg << " L " << scale*(p[0] + 1.0) << " " << scale * (0.5 - p[2]);
    }
    string cols[] = { "grey", "purple", "blue", "yellow", "red", "orange", "green", "purple", "blue", "lightgreen" };
    double wid = 8.0*square.width;
    svg << "\" opacity=\"0.75\" fill=\"" << cols[(int)(7.5 + type)] << "\" stroke=\"black\" stroke-width=\"" << wid << "\" />\n";
  }
  svg << "</svg>" << endl;
  svg.close();
}


int _tmain(int argc, _TCHAR* argv[])
{
  vector<Square> squares;
  double layer = 0.0;
  for (double w = 0.02; w < 8.0; w *= 2.0)
  {
    for (int x = 0; x < 16; x++)
    {
      for (int y = 8; y < 16; y++)
      {
        if (x>=7 && x < 9 && y>=7 && y < 9)
          continue;
        Square square;
        square.pos[0] = ((double)x - 8.0) * w / 8.0;
        square.pos[1] = ((double)y - 8.0) * w / 8.0;
        square.pos[2] = -0.35 + 0.1*w;
        square.width = w / 8.0;
        if (abs(square.pos[1]) < 1.0)
//        if (square.pos[1] > 0.0)
          squares.push_back(square);
      }
    }
    layer++;
  }

  saveSVG("squares.svg", squares);
}
