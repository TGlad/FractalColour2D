#include "stdafx.h"
#include "bmp.h"
#include <sstream>
#include <fstream>
#include <iostream>
#define BIG
#if defined BIG
static const int width = 1280;
static const int height = 1280;
#else
static const int width = 160;
static const int height = 160;
#endif

void setPixel(vector<BYTE> &out, int x, int y, const Vector3d &col)
{
  if (x < 0 || x >= width || y < 0 || y >= height)
    return;
  int ind = 3 * (x + width * (height - 1 - y));
  out[ind + 0] = (BYTE)max(0.0, min(col[0]*255.0, 255.0));
  out[ind + 1] = (BYTE)max(0.0, min(col[1] * 255.0, 255.0));
  out[ind + 2] = (BYTE)max(0.0, min(col[2] * 255.0, 255.0));
}

void drawDisk(vector<BYTE> &out, int x, int y, int radius, const Vector3d &col)
{
  for (int i = x - radius; i <= x + radius; i++)
  {
    for (int j = y - radius; j <= y + radius; j++)
    {
      double r2 = sqr(i - x) + sqr(j - y);
      double R2 = sqr(radius);
      if (r2 < R2)
        setPixel(out, i, j, col);
    }
  }
}

void drawCircle(vector<BYTE> &out, int x, int y, int radius, const Vector3d &col)
{
  for (int i = x - radius; i <= x + radius; i++)
  {
    for (int j = y - radius; j <= y + radius; j++)
    {
      double r2 = sqr(i - x) + sqr(j - y);
      double R2 = sqr(radius);
      if (r2 > R2*0.95 && r2 < R2*1.05)
        setPixel(out, i, j, col);
    }
  }
}

static Vector3d bends[3];
static Vector3d flips[3];
static double scale = 2.0;
static Vector3d diag = Vector3d(1, 1, 1) / sqrt(3.0);
static Vector3d xAxis = diag.cross(Vector3d(0, 0, 1)).normalized();
static Vector3d yAxis = diag.cross(xAxis);

Vector3d distort(Vector3d pos, double &radius, Vector3d &curvature, int i)
{
  pos -= bends[i];
  double mult = bends[i].squaredNorm() / pos.squaredNorm();
 // Vector3d curv = pos / bends[i].squaredNorm();
  pos *= mult;
  pos += bends[i];
  radius *= mult*scale;
  pos -= 2.0 * flips[i] * flips[i].dot(pos);
 // curv -= 2.0*flips[i] * flips[i].dot(curv);
  pos *= scale;
 // curv /= scale;
//  curvature += curv;
  return pos;
}

void testThreshold(double offset)
{
  double minAlong = 1e10;
  double maxAlong = 0.0;
  for (int i = 0; i<10000; i++)
  {
    Vector3d pos(random(-0.5, 0.5), random(-0.5, 0.5), random(-0.5, 0.5));
    pos -= diag * pos.dot(diag);
    pos += diag * offset;
    Vector3d c = pos;
    for (int j = 0; j<200; j++)
    {
      int id = rand() % 3;
      double rad = 1.0;
      Vector3d curv(0, 0, 0);
      pos = distort(pos, rad, curv, id) + c;
    }
    double along = pos.dot(diag);
    minAlong = min(minAlong, along);
    maxAlong = max(maxAlong, along);
  }
  cout << "offset: "<< offset << ", min along: " << minAlong << ", maxAlong: " << maxAlong << endl;
}
static Vector3d boundMin, boundMax;

double recurse(Vector3d pos, double radius, Vector3d curvature, int level, const Vector3d &c)
{
 // double xx = c.dot(diag);
  const double thresh = 0.55; // 0.65 for other one // about perfect unless we want to get into paraboloids
  double z = pos.dot(diag);
  if (z - radius > thresh)
  {
    /*     Vector3d pole = pos - curvature / curvature.squaredNorm();
    double rad = (diag.dot(pole) + 0.3)*0.5;
    Vector3d mid = pole - diag * rad;
    double someScale = ..;
    Vector3d inverted = pole + someScale * (pos - pole) / (pos - pole).squaredNorm();
    double distance = inverted.norm() - rad;
    return distance;
    */
    return 100.0;
  }
  if (level == 0)
    return 0;


  double minDist = 100;
  for (int i = 0; i<3; i++)
  {
    double rad = radius;
    Vector3d curv = curvature;
    Vector3d p = distort(pos, rad, curv, i) + c;
    double dist = recurse(p, rad, curv, level - 1, c);
    if (dist == 0.0)
      return 0.0;
 //   minDist = min(minDist, dist);
  }
  return minDist;
}

double getDistanceEstimate(Vector3d pos, double startRad)
{
#if defined BIG
  const int numIterations = 32;// 27;
#else
  const int numIterations = 16;
#endif
  Vector3d curv(0, 0, 0);

  return recurse(pos, /*4.0**/startRad, curv, numIterations, pos);
}

int _tmain(int argc, _TCHAR* argv[])
{
  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey
  vector<double> depths(width*height);
  vector<bool> shadowed(width*height);

  for (int i = 0; i<3; i++)
  {
    bends[i] = Vector3d(0, 0, 0);
    bends[i][i] = 1.0;
    flips[i] = (Vector3d(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0) - bends[i]).normalized();

  //  bends[i] = diag * sqrt(0.5) - flips[i] * sqrt(0.5);
  }
  // testThreshold(-0.2);
  Vector3d cameraPos = Vector3d(1, 0, 0)*2.5 - Vector3d(0, 1, 0)*2.0;
  cameraPos += diag*(1.8 - cameraPos.dot(diag));
  Vector3d cameraTo = 1.6*Vector3d(1, 0, 0);
  cameraTo += diag*(-0.4 - cameraTo.dot(diag));
//#define CLOSEIN
#if defined CLOSEIN
  Vector3d origin = -diag*0.6;
  Vector3d diff = origin - cameraTo;
  cameraTo += diff;
  cameraPos += diff;
  cameraPos = origin + (cameraPos - origin)*0.5;
#endif
  Vector3d cameraDir = (cameraTo-cameraPos).normalized();
  Vector3d cameraUp = diag;
  Vector3d cameraSide = cameraDir.cross(cameraUp).normalized();
  cameraUp = cameraSide.cross(cameraDir);
  Vector3d sunDir = (cameraUp + 0.3*diag).normalized();

  const double fov = 0.4;
  double minDepth = 1e10;
  double maxDepth = -1e10;
  double minOut = 1e10;
  double maxOut = -1e10;
  double avOut = 0.0;
  double avC = 0.0;
  double maxOutX, maxOutY;
#if defined BIG
  ifstream ifile("cached.bin", ios::binary);
  vector<double> ds(160*160);
  ifile.read((char *)&ds[0], sizeof(double)*ds.size()); 
  ifile.close();
#endif
  for (int i = 0; i<width; i++)
  {
    cout << "i: " << i << endl;
    double x = 2.0 * (double)i / (double)(width - 1) - 1.0;
    for (int j = 0; j<height; j++)
    {
      double y = 2.0 * (double)j / (double)(height - 1) - 1.0;
      Vector3d rayDir = (fov*(cameraSide*x + cameraUp*y) + cameraDir).normalized();

#if defined BIG
      int I = i*160/width;
      int J = j*160/width;
      double depth = ds[I * 160 + J];
#else
      double depth = (0.15 + cameraPos.dot(diag)) / -rayDir.dot(diag); // start nearby
#endif
      double distance = 1e10;
      const double minDistance = 0.0;
      double step;
      double oldDepth = depth;
#if defined BIG
#if defined CLOSEIN
      for (step = 0.00125; step >= 0.000125; step /= 2.0)
#else
      for (step = 0.0025; step >= 0.00025; step /= 2.0)
#endif
#else
      for (step = 0.02; step >= 0.001; step /= 2.0)
#endif
      {
        distance = 1e10;
        while (distance > minDistance)
        {
          double pixWidth = depth * 2.0 * fov / (double)(height - 1);
          Vector3d p = cameraPos + depth*rayDir;
          distance = getDistanceEstimate(p, pixWidth);
          oldDepth = depth - step;
          depth += step;
        }
        depth -= 1.75*step;
      }
      Vector3d hitPoint = cameraPos + oldDepth*rayDir;
      minOut = min(minOut, hitPoint.dot(diag));
      if (hitPoint.dot(diag) > maxOut)
      {
        maxOutX = hitPoint.dot(xAxis);
        maxOutY = hitPoint.dot(yAxis);
      }
      maxOut = max(maxOut, hitPoint.dot(diag));
      avOut += hitPoint.dot(diag);
      avC++;
      depths[width*i + j] = oldDepth;
      minDepth = min(minDepth, depth);
      maxDepth = max(maxDepth, depth);
#if 1 // defined BIG
#define SHADOWS
#if defined SHADOWS
      hitPoint -= rayDir*step;
      double pixWidth = oldDepth * 2.0 * fov / (double)(height - 1);
      depth = 0.0;
      shadowed[width*i + j] = false;
      while (depth < 0.2)
      {
        depth += step; 
#if defined CLOSEIN
        if (step < 0.005)
#else
        if (step < 0.01)
#endif
          step *= 1.5;
        distance = getDistanceEstimate(hitPoint + depth*sunDir, pixWidth);
        if (distance == 0.0)
        {
          shadowed[width*i + j] = true;
          break;
        }
      }
#endif
#endif
    }
  }
  cout << "minDepth: " << minDepth << ", maxDepth: " << maxDepth << ", min out: " << minOut << ", av out: " << avOut/avC << ", maxOut: " << maxOut << " x: " << maxOutX << ", y: " << maxOutY << endl;

  // next: proper diffuse light...
  double specular = 0.25;
  for (int i = 1; i<width-1; i++)
  {
    double x = 2.0 * (double)i / (double)(width - 1) - 1.0;
    for (int j = 1; j<height-1; j++)
    {
      double y = 2.0 * (double)j / (double)(height - 1) - 1.0; 
      Vector3d rayDir = (fov*(cameraSide*x + cameraUp*y) + cameraDir).normalized();
      Vector3d upDir = cameraSide.cross(rayDir).normalized();

      double pixWidth = depths[width*i + j] * 2.0 * fov / (double)(height - 1);
      double dif1 = depths[width*i + j] - depths[width*i + j - 1];
      Vector3d normal1 = (-rayDir*pixWidth + upDir*dif1).normalized();
      Vector3d normalS = sunDir - 2.0*(sunDir - normal1*normal1.dot(sunDir));
      double d = max(0.0, -normalS.dot(rayDir));
      d = d*d;
      d = d*d;
      double shade1 = max(0.0, 0.5 + 0.5*normal1.dot(sunDir));
      double spec = d*d*d*d;
      double dif2 = depths[width*i + j + 1] - depths[width*i + j];
      Vector3d normal2 = (-rayDir*pixWidth + upDir*dif2).normalized();
      normalS = sunDir - 2.0*(sunDir - normal2*normal2.dot(sunDir));
      d = max(0.0, -normalS.dot(rayDir));
      d = d*d;
      d = d*d;
      double shade2 = max(0.0, 0.5 + 0.5*normal2.dot(sunDir));
      spec += d*d*d*d;
      spec = specular * spec / 2.0;
      double shade = (shade1 + shade2) / 2.0;

      // next: screen-space ambient-occlusion
      double convex = 0.0;
      double sum = 0.0;
      int ambientRange = width / 80;
      for (int ii = max(0, i - 2 * ambientRange); ii <= min(i + 2 * ambientRange, width - 1); ii += 2)
      {
        for (int jj = max(0, j - 2 * ambientRange); jj <= min(j + 2 * ambientRange, height - 1); jj += 2)
        {
          if (ii == i && jj == j)
            continue;
          double dif = depths[width*ii + jj] - depths[width*i + j];
          double s = 1.0 / (sqr(ii - i) + sqr(jj - j));
          sum += s;
          convex += dif * s;
        }
      }
      double ambientOcclusion = 0.4;
      shade += ambientOcclusion * max(-0.5, min((convex / sum) / ((double)ambientRange * pixWidth), 0.5));
      Vector3d col(shade, shade, shade);

      Vector3d sunColour(1.0, 0.85, 0.3);
      col += spec*sunColour;

#if defined SHADOWS
      if (shadowed[width*i + j])
        col *= 0.6;
#endif

      double blend = min(pow(1.0/1.4, depths[width*i + j] - 2.8), 1.0);
      Vector3d fogColour(0.4, 0.5, 0.7);
      col = col*(blend) + fogColour*(1.0 - blend);
      setPixel(out, i, j, /*0.9**/col);
    }
  }

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  LPCTSTR file = L"mob3d_out.bmp";
  SaveBMP(c, width, height, s2, file);
  delete[] c;
#if !defined BIG
  vector<double> ds = depths;
  for (int i = 0; i < width; i++)
  {
    for (int j = 0; j < height; j++)
    {
      double minDepth = 1e10;
      for (int ii = max(0, i - 1); ii <= min(i + 1, width - 1); ii++)
      {
        for (int jj = max(0, j - 1); jj <= min(j + 1, height - 1); jj++)
        {
          minDepth = min(minDepth, depths[ii*width + jj]);
        }
      }
      ds[i*width + j] = minDepth - 0.025;
    }
  }

  ofstream ofile("cached.bin", ios::binary);
  ofile.write((char *)&ds[0], sizeof(double)*ds.size());
#endif
  return 1;
}
 