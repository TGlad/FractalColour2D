#include "stdafx.h"
#include "bmp.h"
#include <fstream>
#include <complex>
typedef complex<double> Complex;

static double width = 1024;
static double height = 512;
struct Disk
{
  double angle;
  Vector2d pos;
  double radius;
};
void saveSVG(const string &fileName, vector<Disk> &disks)
{
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)width << "\" height = \"" << (int)height << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;

  for (auto &disk : disks)
  {
    svg << "<circle cx = \"" << disk.pos[0] << "\" cy = \"" << height - disk.pos[1] << "\" r = \"" << disk.radius << "\" stroke = \"black\" stroke-width = \"3\" fill = \"lightsteelblue\" />" << endl;
    Vector2d top = disk.pos + Vector2d(-sin(disk.angle), cos(disk.angle))*disk.radius*0.8;
    svg << "<path d = \"M " << disk.pos[0] << " " << height - disk.pos[1] << " L " << top[0] << " " << height - top[1] << "\" fill=\"none\" stroke-width = \"3\" stroke=\"black\" />\n";
  }
  Disk parentDisk;
  parentDisk.angle = 0.0;
  for (auto &disk : disks)
  {
    Vector2d pos = disk.pos + Vector2d(cos(disk.angle), sin(disk.angle))*disk.radius*0.75;
    // do an arc...
    bool found = false;
    for (double a = parentDisk.angle; a < parentDisk.angle + (disk.angle - parentDisk.angle)/3.0; a += 0.025)
    {
      Vector2d p = parentDisk.pos + Vector2d(cos(a), sin(a))*parentDisk.radius;
      if (a == parentDisk.angle)
        svg << "<path d = \"M " << p[0] << " " << height - p[1];
      else
        svg << " L " << p[0] << " " << height - p[1];
      found = true;
    }
    for (double a = parentDisk.angle; a > parentDisk.angle + (disk.angle - parentDisk.angle) / 3.0; a -= 0.025)
    {
      Vector2d p = parentDisk.pos + Vector2d(cos(a), sin(a))*parentDisk.radius;
      if (a == parentDisk.angle)
        svg << "<path d = \"M " << p[0] << " " << height - p[1];
      else
        svg << " L " << p[0] << " " << height - p[1];
      found = true;
    }
    if (found)
    {
      svg << "\" fill=\"none\"";
      svg << " stroke-width = \"3\" stroke=\"red\" />\n";
    }
    found = false;
    for (double a = parentDisk.angle; a < disk.angle; a += 0.025)
    {
      Vector2d p = disk.pos + Vector2d(-sin(a), cos(a))*disk.radius * 0.75;
      if (a == parentDisk.angle)
        svg << "<path d = \"M " << p[0] << " " << height - p[1];
      else
        svg << " L " << p[0] << " " << height - p[1];
      found = true;
    }
    for (double a = parentDisk.angle; a > disk.angle; a -= 0.025)
    {
      Vector2d p = disk.pos + Vector2d(-sin(a), cos(a))*disk.radius * 0.75;
      if (a == parentDisk.angle)
        svg << "<path d = \"M " << p[0] << " " << height - p[1];
      else
        svg << " L " << p[0] << " " << height - p[1];
      found = true;
    }
    if (found)
    {
      svg << "\" fill=\"none\"";
      svg << " stroke-width = \"2\" stroke=\"darkgreen\" />\n";

      // now the arrow head!
      Vector2d cs[] = { Vector2d(0.7, 1.0), Vector2d(0.7, -1.0), Vector2d(1.0, 0.0) };
      for (int i = 0; i < 3; i++)
      {
        double a = parentDisk.angle + (disk.angle - parentDisk.angle)*cs[i][0];
        double r = disk.radius * 0.75 + 4.0*cs[i][1];
        Vector2d p = disk.pos + r * Vector2d(-sin(a), cos(a));
        if (i == 0)
          svg << "<path d = \"M " << p[0] << " " << height - p[1];
        else
          svg << " L " << p[0] << " " << height - p[1];
      }
      svg << "\" fill=\"darkgreen\"";
      svg << " stroke-width = \"2\" stroke=\"darkgreen\" />\n";
    }
    parentDisk = disk;
  }
  svg << "</svg>" << endl;
  svg.close();
}
static double speed = 85.0;

void recurse(vector<Disk> &disks, Disk &disk, double dir)
{
  disks.push_back(disk);
  if (disk.radius < 8.0)
    return;
  Disk child;
  child.radius = disk.radius / 2.0;
  child.angle = disk.angle + speed*dir / child.radius;
  double ang = disk.angle + (child.angle - disk.angle) / 3.0;
  child.pos = disk.pos + (disk.radius + child.radius)*Vector2d(cos(ang), sin(ang));
  recurse(disks, child, -dir);
}
int _tmain(int argc, _TCHAR* argv[])
{
  vector<Disk> disks;
  Disk base;
  base.pos = Vector2d(width / 5.0, height / 2.0);
  base.angle = 0.0;
  base.radius = 200.0;
  recurse(disks, base, 1.0);

  saveSVG("rolling.svg", disks);
}
