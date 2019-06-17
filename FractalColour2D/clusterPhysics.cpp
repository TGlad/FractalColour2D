#include "stdafx.h"
#include "bmp.h"
#include <set>
#include <sstream>
#include <fstream>
static int width = 3200;
static int height = 3200;

void addpixel(vector<BYTE> &out, const Vector2i &pos, double col, int comp)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= height)
    return;
  int ind = 3 * (pos[0] + width * (height - 1 - pos[1]));
  out[ind + comp] = max(0, out[ind + comp] - (int)(255.0*col));
//  out[ind + 1] = max(0, out[ind + 1] - (int)(255.0*col));
//  out[ind + 2] = max(0, out[ind + 2] - (int)(255.0*col));
}

struct Disk
{
  double mass;
  double radius;
  double x, v, time;
};

struct Collision
{
  int disk1;
  int disk1Size, disk2Size;
  double time;
};

void drawHorizLine(vector<BYTE> &out, double x1, double x2, int y, double blend, int comp)
{
  int X1 = (int)x1;
  int X2 = (int)x2;
  for (int x = X1; x <= X2; x++)
  {
    double shade = 1.0;
    if (x == X1)
      shade -= x1 - (double)X1;
    if (x == X2)
      shade -= (double)(X2+1) - x2;
    addpixel(out, Vector2i(x, y), shade*blend, 0);// comp);
    addpixel(out, Vector2i(x, y), shade*blend, 1);
    addpixel(out, Vector2i(x, y), shade*blend, 2);
  }
}

void drawParallelogram(vector<BYTE> &out, double x1, double y1, double x2, double y2, double radius, int comp)
{
  int startY = (int)(y1*(double)height);
  int endY = (int)(y2*(double)height);
  for (int y = startY; y <= endY; y++)
  {
    double blend = 1.0;
    if (y == startY)
      blend -= y1 - (double)(startY);
    if (y == endY)
      blend -= (double)(endY+1)-y2;
    double centre = x1 + (x2 - x1)*max(0.0, min(((((double)y+0.5)/(double)height) - y1)/(y2 - y1), 1.0));
    double xa = ((centre - radius) - 0.18)*4.0;  // 0.63
    double xb = ((centre + radius) - 0.18)*4.0;
    drawHorizLine(out, xa*(double)width, xb*(double)width, y, blend, comp);
  }
}

struct LessThan
{
  bool operator()(const Collision &collision1, const Collision &collision2)
  {
    return collision1.time < collision2.time;
  }
};

void addDisk(vector<Disk> &disks, double x0, double x1, double vel, double radius, int level)
{
  if (level == 0)
    return;
  addDisk(disks, x0, x0 + (x1-x0) / 3.0, -vel, 0, level-1);
  Disk disk;
  disk.x = (x0 + x1) / 2.0;
  disk.time = 0.0;
  disk.radius = (x1 - x0)*0.125;
  disk.mass = disk.radius;
  disk.v = vel;
  disks.push_back(disk);
  addDisk(disks, x1 - (x1 - x0)/3.0, x1, -vel, 0, level-1);

/*  Disk disk;
  disk.radius = radius;
  double mid = (x0 + x1) / 2.0;
  addDisk(disks, x0, mid - disk.radius, -vel, radius/3.0, level - 1);
  disk.x = (x0 + x1) / 2.0;
  disk.time = 0.0;
  disk.mass = disk.radius;
  disk.v = vel;
  disks.push_back(disk);
  addDisk(disks, mid+disk.radius, x1, -vel, radius/3.0, level - 1);*/
}

void applyPhysics(vector<BYTE> &out, int comp, int level)
{
  // 1. construct the pattern
  vector<vector<Disk> > diskTrajs;
#if 0
  double vel = 1.0;
  int numDisks = 10; // even
  for (int i = 0; i < numDisks; i++)
  {
    double x = pow(0.5, numDisks - 1 - i);
    double radius = 0.2*x;

    Disk disk;
    disk.radius = radius;
    disk.mass = disk.radius;
    disk.x = x;
    disk.v = vel;
    disk.time = 0.0;
    vel = -vel;
    vector<Disk> traj;
    traj.push_back(disk);
    diskTrajs.push_back(traj);
  }
#else
  vector<Disk> disks;
  addDisk(disks, 0.0, 1.0, 0.5/3.0, 0.07/*0.125*/, level);
  for (int i = 0; i < (int)disks.size(); i++)
  {
    vector<Disk> traj;
    traj.push_back(disks[i]);
    diskTrajs.push_back(traj);
  }
#endif
  // 2. insert collisions in order based on some comparisons:
  vector<double> collisionTimes(diskTrajs.size() - 1);
  for (int i = 0; i < (int)diskTrajs.size() - 1; i++)
  {
    Disk &disk1 = diskTrajs[i].back();
    Disk &disk2 = diskTrajs[i + 1].back();
    if (disk2.v >= disk1.v)
      collisionTimes[i] = -1.0;
    else
      collisionTimes[i] = ((disk2.x - disk2.radius) - (disk1.x + disk1.radius)) / abs(disk1.v - disk2.v);
  }

  // now go through the collisions in order and extend the trajectories
  for (;;)
  {
    double collisionTime = 1.0;
    int disk1 = 0;
    for (int i = 0; i < (int)collisionTimes.size(); i++)
    {
      if (collisionTimes[i] < collisionTime && collisionTimes[i] >= 0)
      {
        collisionTime = collisionTimes[i];
        disk1 = i;
      }
    }
    collisionTimes[disk1] = -1;
    if (collisionTime >= 1.0)
      break;

    Disk d1, d2;
    d1 = diskTrajs[disk1].back();
    d2 = diskTrajs[disk1+1].back();
    double x = d1.x + d1.radius + (collisionTime - d1.time) * d1.v;
    d1.x = x-d1.radius;
    d2.x = x+d2.radius;
    d1.time = collisionTime;
    d2.time = collisionTime;
    
    // elastic bounce calculation:
#define RELATIVISTIC
#if defined RELATIVISTIC
    double Z = sqrt((1.0-sqr(d1.v))*(1.0-sqr(d2.v)));
    double m2 = sqr(d1.mass) + sqr(d2.mass);
    double m2min = sqr(d1.mass) - sqr(d2.mass);
    double v1n = 2.0*d1.mass*d2.mass*d2.v*Z + 2.0*sqr(d2.mass)*d2.v - m2*d1.v*sqr(d2.v) + m2min*d1.v;
    double v1d = 2.0*d1.mass*d2.mass*Z - 2.0*sqr(d2.mass)*d1.v*d2.v - m2min*sqr(d2.v) + m2;

    double v2n = 2.0*d1.mass*d2.mass*d1.v*Z + 2.0*sqr(d1.mass)*d1.v - m2*d2.v*sqr(d1.v) - m2min*d2.v;
    double v2d = 2.0*d1.mass*d2.mass*Z - 2.0*sqr(d1.mass)*d1.v*d2.v + m2min*sqr(d1.v) + m2;

    double v1 = v1n / v1d;
    double v2 = v2n / v2d;
#else
    double v1 = (d1.mass - d2.mass)*d1.v / (d1.mass + d2.mass) + 2.0*d2.mass*d2.v / (d1.mass + d2.mass);
    double v2 = 2.0*d1.mass*d1.v / (d1.mass + d2.mass) + (d2.mass - d1.mass)*d2.v / (d1.mass + d2.mass);
#endif
    d1.v = v1;
    d2.v = v2;
    diskTrajs[disk1].push_back(d1);
    diskTrajs[disk1 + 1].push_back(d2);

    // now add the collision data to neighbours...
    if (disk1 > 0 && d1.v < diskTrajs[disk1 - 1].back().v)
    {
      Disk &diskA = diskTrajs[disk1 - 1].back();
      Disk &diskB = diskTrajs[disk1].back();
      collisionTimes[disk1 - 1] = ((diskB.x - diskB.radius) - (diskA.x + diskA.radius) + (diskA.v*diskA.time - diskB.v*diskB.time)) / (diskA.v - diskB.v);
    }
    if (disk1+1 < (int)diskTrajs.size()-1 && d2.v > diskTrajs[disk1 + 2].back().v)
    {
      Disk &diskA = diskTrajs[disk1+1].back();
      Disk &diskB = diskTrajs[disk1+2].back();
      collisionTimes[disk1 + 1] = (((diskB.x - diskB.radius) - (diskA.x + diskA.radius)) + (diskA.v*diskA.time - diskB.v*diskB.time)) / (diskA.v - diskB.v);
    }
  }

  // 2. now render the results quite nicely!
  for (auto &traj : diskTrajs)
  {
    for (int i = 1; i < (int)traj.size(); i++)
    {
      double error = abs(traj[i - 1].x + traj[i - 1].v*(traj[i].time - traj[i - 1].time) - traj[i].x);
      if (error > 1e-3)
        cout << "error" << error << endl;
      drawParallelogram(out, traj[i - 1].x, traj[i - 1].time, traj[i].x, traj[i].time, traj[i].radius, comp);
    }
    drawParallelogram(out, traj.back().x, traj.back().time, traj.back().x + traj.back().v*(1.0 - traj.back().time), 1.0, traj.back().radius, comp);
  }
}

int _tmain(int argc, _TCHAR* argv[])
{
  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is white

  for (int comp = 2; comp < 3; comp++)
  {
    applyPhysics(out, comp, 9 + comp);
  }
  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"ssClusterPhysics2Relativistic.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}
