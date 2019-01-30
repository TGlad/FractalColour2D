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
      svg << "\" fill=\"darkgrey\"";
    else if (colour == 2)
      svg << "\" fill=\"grey\"";
    else if (colour == 3)
      svg << "\" fill=\"lightgrey\"";
    else if (colour == 4)
      svg << "\" fill=\"darkgreen\"";
    else if (colour == 5)
      svg << "\" fill=\"green\"";
    else if (colour == 6)
      svg << "\" fill=\"lightgreen\"";
    else
      svg << "\" fill=\"none\"";
    if (colour < 4)
      svg << " fill-opacity = \"0.5\"";
    svg << " stroke-width = \"2\" stroke=\"black\" />\n";
  }

  svg << "</svg>" << endl;
  svg.close();
}

double h = sqrt(3.0 / 4.0);
int _tmain(int argc, _TCHAR* argv[])
{
  // I have to build a grid, cut off the long diagonal plane, and put lines in place, then render in 3dish, 
  // this is a bit tricky!
  double pitch = 0.5;
  double yaw = 1.05;
  Vector2d X(cos(yaw), -sin(yaw)*sin(pitch));
  Vector2d Y(-sin(yaw), -cos(yaw)*sin(pitch));
  Vector2d Z(0.0, cos(pitch));
  const int n = 2;
  int grid[n][n][n];
  Vector3d sides[3][4] = {{ Vector3d(1, -1, -1), Vector3d(1, -1, 1), Vector3d(1, 1, 1), Vector3d(1, 1, -1) },
                          { Vector3d(-1, 1, -1), Vector3d(-1, 1, 1), Vector3d(1, 1, 1), Vector3d(1, 1, -1) },
                          { Vector3d(-1, -1, 1), Vector3d(-1, 1, 1), Vector3d(1, 1, 1), Vector3d(1, -1, 1) } };
  vector<vector<Vector2d> > lists;
  vector<int> colours;
  for (int i = 0; i < 3; i++)
  {
    vector<Vector2d> list;
    for (int j = 0; j < 4; j++)
    {
      Vector3d pos = Vector3d(1, 1, 1) - 0.5*(double)n * sides[i][j];
      double w = (double)width;
      Vector2d pixel = 0.5*Vector2d(w, w) + (pos[0] * X + pos[1] * Y + pos[2] * Z)*0.5*w / (double)n;
      list.push_back(pixel);
    }
    lists.push_back(list);
    colours.push_back(i+1);
  }
  for (int i = 0; i < 3; i++)
  {
    vector<Vector2d> list;
    for (int j = 0; j < 4; j++)
    {
      Vector3d pos = Vector3d(0.75, 0.75, 0.75) + 0.25*sides[i][j];
      double w = (double)width;
      Vector2d pixel = 0.5*Vector2d(w, w) + (pos[0] * X + pos[1] * Y + pos[2] * Z)*0.5*w / (double)n;
      list.push_back(pixel);
    }
    lists.push_back(list);
    colours.push_back(i+4);
  }
  for (int x = 0; x < n; x++)
    for (int y = 0; y < n; y++)
      for (int z = 0; z < n; z++)
      {
        if (x==1 && y==1 && z==1)
          continue;
        for (int i = 0; i < 3; i++)
        {
          vector<Vector2d> list;
          for (int j = 0; j < 4; j++)
          {
            Vector3d pos = Vector3d(x+0.5, y+0.5, z+0.5) + 0.5*sides[i][j];
            double w = (double)width;
            Vector2d pixel = 0.5*Vector2d(w,w) + (pos[0] * X + pos[1] * Y + pos[2] * Z)*0.5*w/(double)n;
            list.push_back(pixel);
          }
          lists.push_back(list);
          colours.push_back(i+1);
        }
      }

  saveSVG("voxelparent.svg", lists, colours);
}
