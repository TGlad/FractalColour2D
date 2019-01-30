#include "stdafx.h"
#include "bmp.h"
#include <fstream>
#include <complex>
typedef complex<double> Complex;

double width = 700;
double height = 700;

void saveSVG(const string &fileName, vector<vector<Vector2d>> &lists, vector<int> &colours)
{
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)width << "\" height = \"" << (int)height << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;

  for (int j = 0; j < lists.size(); j++)
  {
    auto &list = lists[j];
    int colour = colours[j];
    Vector2d end = list.back();
    svg << "<path d = \"M " << end[0] << " " << height - end[1];
    for (int i = 0; i < list.size(); i++)
      svg << " L " << list[i][0] << " " << height - list[i][1];
    if (colour == 1)
      svg << "\" fill=\"darkred\" stroke-width = \"2\" stroke=\"black\" />\n";
    else if (colour == 2)
      svg << "\" fill=\"green\" stroke-width = \"2\" stroke=\"black\" />\n";
    else if (colour == 3)
      svg << "\" fill=\"grey\" stroke-width = \"2\" stroke=\"black\" />\n";
    else
      svg << "\" fill=\"none\" stroke-width = \"2\" stroke=\"black\" />\n";
  }

  svg << "</svg>" << endl;
  svg.close();
}

double h = sqrt(3.0 / 4.0);
int _tmain(int argc, _TCHAR* argv[])
{
#if 0 // hexagonal
  vector<vector<Vector2d> >grids;
  vector<int> colours;
  double yOff = 0.0;
  double number = 20;
  int X = 0;
  for (double x = 0; x < number; x += 1.5)
  {
    yOff = h - yOff;
    int Y = 0;
    for (double y = yOff; y < number; y += 2.0*h)
    {
      int colour = 0;

 /*     if (!(X % 2) && !(Y % 3))
        colour = 3;
      if (!((X + 1) % 2) && !((Y + 1) % 3))
        colour = 3;
   */   
  /*    if (!(X % 6) && !(Y % 3))
        colour = 2;
      if (!((X + 3) % 6) && !((Y + 1) % 3))
        colour = 2;
        */
      if (!(X % 6) && !(Y % 9))
        colour = 1;
      if (!((X + 3) % 6) && !((Y + 4) % 9))
        colour = 1;

      vector<Vector2d> grid(6);
      Vector2d hex[] = { Vector2d(-0.5, h), Vector2d(0.5, h), Vector2d(1, 0), Vector2d(0.5, -h), Vector2d(-0.5, -h), Vector2d(-1, 0) };
      for (int i = 0; i < 6; i++)
        grid[i] = (Vector2d(x, y) + hex[i])*(double)width/number;
      grids.push_back(grid);
      colours.push_back(colour);
      Y++;
    }
    X++;
  }
  saveSVG("hexagontile0.svg", grids, colours);
#else
  vector<vector<Vector2d> >grids;
  vector<int> colours;
  double yOff = 0.0;
  double number = 10;
  int X = 0;
  for (double x = 0; x < number; x += 1)
  {
    int Y = 0;
    for (double y = 0; y < number; y += 1)
    {
      int colour = 0;

  /*    if (!(X % 2) && !(Y % 2))
        colour = 3;
      if (!((X + 1) % 2) && !((Y + 1) % 2))
        colour = 3;
    */  
 //     if (!(X % 2) && !(Y % 2))
 //       colour = 2;

      if (!(X % 4) && !(Y % 4))
        colour = 1;
      if (!((X + 2) % 4) && !((Y + 2) % 4))
        colour = 1;

      vector<Vector2d> grid(4);
      Vector2d hex[] = { Vector2d(-0.5, 0.5), Vector2d(0.5, 0.5), Vector2d(0.5, -0.5), Vector2d(-0.5, -0.5) };
      for (int i = 0; i < 4; i++)
        grid[i] = (Vector2d(x, y) + hex[i])*(double)width / number;
      grids.push_back(grid);
      colours.push_back(colour);
      Y++;
    }
    X++;
  }
  saveSVG("squaretile0.svg", grids, colours);
#endif
}
