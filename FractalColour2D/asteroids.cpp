#include "stdafx.h"
#include "bmp.h"
#include <fstream>
double width = 440.0;
double height = 800.0;

struct Disk
{
  double radius;
  double x, y;
};

void saveSVG(const string &fileName, vector<Disk> &disks)
{
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)width << "\" height = \"" << (int)height << "\" xmlns = \"http://www.w3.org/2000/svg\" style=\'background-color: black;\'>" << endl;

  double gMinX = 0.2;
  double gMinY = 0.1;
  double gMaxX = 0.98;
  double gMaxY = 0.98;

  double minY = 0.0;
  double maxY = 8.0;
  double minX = 3.0;
  double maxX = -1.0;

  for (double x = maxX; x < minX; x++)
  {
    for (double x2 = 1; x2 < 10; x2++)
    {
      double X = log(pow(10.0, x) * x2) / log(10.0);
      X = (X - minX) / (maxX - minX);
      X = X * (gMaxX - gMinX) + gMinX;
      svg << "<line x1 = \"" << width*X << "\" y1 = \"" << height - height*gMinY << "\" x2 = \"" << width*X << "\" y2 = \"" << height - height*gMaxY << "\"";
      if (x2 == 1)
        svg << " style = \"stroke:rgb(100, 100, 100); stroke - width:2\" />" << endl;
      else
        svg << " style = \"stroke:rgb(50, 50, 50); stroke - width:1\" />" << endl;
    }
  }
  for (double y = minY; y < maxY; y++)
  {
    for (double y2 = 1; y2 < 10; y2++)
    {
      double Y = log(pow(10.0, y) * y2) / log(10.0);
      Y = (Y - minY) / (maxY - minY);
      Y = Y * (gMaxY - gMinY) + gMinY;
      svg << "<line x1 = \"" << width*gMinX << "\" y1 = \"" << height - height*Y << "\" x2 = \"" << width*gMaxX << "\" y2 = \"" << height - height*Y << "\"";
      if (y2 == 1)
        svg << " style = \"stroke:rgb(100, 100, 100); stroke - width:2\" />" << endl;
      else
        svg << " style = \"stroke:rgb(50, 50, 50); stroke - width:1\" />" << endl;
    }
  }

  double avX = 0.0, avY = 0.0;
  double num = 14.0;
  for (auto &d : disks)
  {
    avX += d.x;
    avY += d.y;
  }
  avX /= num;
  avY /= num;
  double numer = 0, denom = 0;
  for (auto &d : disks)
  {
    numer += (d.x - avX)*(d.y - avY);
    denom += sqr(d.x - avX);
  }

  double gradient = -numer / denom;
  avX = 3.0 - (avX + 1.0);
  double a = avY - gradient*(avX - maxX);
  double b = a + gradient*(minX - maxX);
  a = (a - minY) / (maxY - minY);
  b = (b - minY) / (maxY - minY);
  svg << "<line stroke-dasharray=\"8, 8\" x1 = \"" << width*gMinX << "\" y1 = \"" << height - height*((a* (gMaxY - gMinY)) + gMinY) << "\" x2 = \"" << width*gMaxX << "\" y2 = \"" << height - height*((b*(gMaxY - gMinY)) + gMinY) << "\" style = \"stroke:rgb(20, 170, 20); stroke - width:6\" />" << endl;

  for (int i = disks.size() - 1; i >= 0; i--)
  {
    auto &d = disks[i];
    svg << "<circle cx = \"" << width * ((((d.x - minX) / (maxX - minX))* (gMaxX - gMinX)) + gMinX) << "\" cy = \"" << height - height*((((d.y - minY) / (maxY - minY))* (gMaxY - gMinY)) + gMinY) << "\" r = \"" << max(d.radius, 2.0) + 1.0 << "\" stroke=\"black\" fill=\"burlywood\" />" << endl;
  }

  double y = height - height*(gMinY);
  svg << "<rect x = \"" << 0 << "\" y = \"" << y << "\" width = \"152\" height = \"50\" style = \"fill:rgb(0,0,0);stroke-width:0\" />" << endl;
  double w = width*gMinX;
  svg << "<rect x = \"" << 0 << "\" y = \"" << y - 55 << "\" width = \"" << w << "\" height = \"55\" style = \"fill:rgb(0,0,0);stroke-width:0\" />" << endl;

  svg << "<line x1 = \"" << width*gMinX << "\" y1 = \"" << height - height*gMinY << "\" x2 = \"" << width*gMaxX << "\" y2 = \"" << height - height*gMinY << "\" style = \"stroke:rgb(255, 255, 255); stroke - width:4\" />" << endl;
  svg << "<line x1 = \"" << width*gMinX << "\" y1 = \"" << height - height*gMinY << "\" x2 = \"" << width*gMinX << "\" y2 = \"" << height - height*gMaxY << "\" style = \"stroke:rgb(255, 255, 255); stroke - width:4\" />" << endl;

  string nums[] = { "1", "10", "100", "1,000", "10,000", "100,000", "1 million", "10 million", "100 million" };
  for (int i = 0; i <= maxY; i++)
  {
    double x = gMinX * width;
    double y = ((double)i - minY) / (maxY - minY);
    svg << "<text x = \"" << x-4 << "\" y = \"" << 4 + height - height*(y*(gMaxY - gMinY) + gMinY) << "\" text-anchor=\"end\" fill = \"white\">" << nums[i] << "</text>" << endl;
  }
  string diams[] = { "1000", "100", "10", "1", "0.1"};
  for (int i = 0; i <= (minX - maxX); i++)
  {
    double x = ((double)i) / (minX - maxX);
    double y = 16 + height - gMinY * height;
    svg << "<text x = \"" << width*(x*(gMaxX - gMinX) + gMinX) << "\" y = \"" << y << "\" text-anchor=\"middle\" fill = \"white\">" << diams[i] << "</text>" << endl;
  }
  svg << "<text x = \"" << width*(gMaxX + gMinX) / 2.0 << "\" y = \"" << height - height*(0.5*gMinY) << "\" text-anchor=\"middle\" fill = \"white\">Asteroid diameter in km</text>" << endl;
  svg << "<text transform=\"rotate(-90)\" x = \"" << -height*(1.0 - (gMinY + gMaxY)/2.0) << "\" y = \"" << width*gMinX / 4.0 << "\" text-anchor=\"middle\" fill = \"white\">Number of asteroids larger than diameter</text>" << endl;

  svg << "</svg>" << endl;
  svg.close();
}


int _tmain(int argc, _TCHAR* argv[])
{
  double xs[] = { 0.1, 0.3, 0.5, 1, 3, 5, 10, 30, 50, 100, 200, 300, 500, 900};
  double ys[] = { 25000000, 4000000, 2000000, 750000, 200000, 90000, 10000, 1100, 600, 200, 30, 5, 3, 1};
  vector<Disk> disks(14);
  for (int i = 0; i < 14; i++)
  {
    disks[i].radius = sqrt(xs[13-i] * 1000.0) * 0.04;
    disks[i].x = log(xs[13 - i]) / log(10.0);
    disks[i].y = log(ys[13 - i]) / log(10.0);
  }
  saveSVG("asteroids1.svg", disks);
}
