#include "stdafx.h"
#include "bmp.h"
#include <set>
#include <sstream>
#include <fstream>
#include "Core\EigenBase.h"
static int width = 1000;
static int height = width; // keep the same

void setPixel(vector<BYTE> &out, const Vector2i &pos, const Vector3d &col)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= height)
    return;
  int ind = 3 * (pos[0] + width * (height - 1 - pos[1]));
  out[ind + 0] = (int)(255.0*min(col[0], 1.0));
  out[ind + 1] = (int)(255.0*min(col[1], 1.0));
  out[ind + 2] = (int)(255.0*min(col[2], 1.0));
}

Vector3d getCol(const Vector3d &pos, const Vector3d &dir, double blend, bool frontSide)
{
  Vector3d sunDir(0, 0, -1);
  Vector3d sunCol(1.0, 1.0, 0.8);
  sunCol *= 0.5;
  Vector3d col(1, 1, 1);
  // a: basic sphere
  Vector3d closest = pos + dir*-pos.dot(dir);
  double perpDist = closest.norm();
  if (perpDist > 1.0)
    return col;
  Vector3d normal = closest;
  if (frontSide)
    normal -= dir*sqrt(1.0 - sqr(perpDist));
  else
    normal += dir*sqrt(1.0 - sqr(perpDist));

  Vector3d sphereColour(1, 0, 0);
  double ambient = frontSide ? 0.1 : 0.5;
  double dotDiffuse = max(0.0, normal.dot(-sunDir));
  Vector3d x = sunDir - normal*normal.dot(sunDir);
  Vector3d y = sunDir - 2.0*x;
  double dot = -y.dot(dir);

  // b: spherical tetrahedron
  Vector3d v0(0.0, 0.0, -1.0);
  Vector3d v1(sqrt(8.0 / 9.0), 0, 1.0 / 3.0);
  Vector3d v2(-sqrt(2.0 / 9.0), sqrt(2.0 / 3.0), 1.0 / 3.0);
  Vector3d v3(-sqrt(2.0 / 9.0), -sqrt(2.0 / 3.0), 1.0 / 3.0);

  Vector3d v12 = (v1 + v2).normalized();
  Vector3d v23 = (v2 + v3).normalized();
  Vector3d v31 = (v3 + v1).normalized();

  Vector3d w1(0.0, 0.0, -1.0);
  Vector3d w2(0.0, 0.0, -1.0);
  Vector3d w3(0.0, 0.0, -1.0);
  Vector3d w12(sqrt(2.0 / 9.0), sqrt(2.0 / 3.0), 1.0 / 3.0);
  Vector3d w23(-sqrt(8.0 / 9.0), 0, 1.0 / 3.0); 
  Vector3d w31(sqrt(2.0 / 9.0), -sqrt(2.0 / 3.0), 1.0 / 3.0);

  Vector3d vs[] = { v0, v1, v2, v3, v12, v23, v31 };
  Vector3d ws[] = { w1, w1, w2, w3, w12, w23, w31 };
  Vector2i edges[] = { Vector2i(1, 4), Vector2i(4, 2), Vector2i(2, 5), Vector2i(5, 3), Vector2i(3, 6), Vector2i(6, 1), Vector2i(4, 5), Vector2i(5, 6), Vector2i(6, 4) };
  Vector2i triEdges[4][3] = { { Vector2i(1, 4), Vector2i(4, 6), Vector2i(6, 1) }, { Vector2i(4, 2), Vector2i(2, 5), Vector2i(5, 4) }, { Vector2i(5, 3), Vector2i(3, 6), Vector2i(6, 5) }, { Vector2i(4, 5), Vector2i(5, 6), Vector2i(6, 4) } };
  int numEdges = 9;
  bool isInside = false;
  int insideI = 0;
  for (int i = 0; i < 4; i++)
  {
    bool inside = true;
    for (int j = 0; j < 3; j++)
    {
      Vector3d va = vs[triEdges[i][j][0]];
      Vector3d vb = vs[triEdges[i][j][1]];
      va += (ws[triEdges[i][j][0]] - va)*blend;
      vb += (ws[triEdges[i][j][1]] - vb)*blend;
      Vector3d sid = va.cross(vb);
      if (normal.dot(sid) < 0.0)
        inside = false;
    }
    if (inside)
      insideI = i+1;
  }
  Vector3d surfaceCols[] = { Vector3d(0.5, 0.5, 0.5), Vector3d(0, 0.8, 0.8), Vector3d(0.1, 0.9, 0.1), Vector3d(0.8, 0.8, 0), Vector3d(1,0,0) };
  col = surfaceCols[insideI] * (dotDiffuse + ambient);
  if (frontSide)
    col += sunCol*dot*dot*dot*dot;

  for (int i = 0; i < numEdges; i++)
  {
    Vector3d va = vs[edges[i][0]];
    Vector3d vb = vs[edges[i][1]];
    va += (ws[edges[i][0]] - va)*blend;
    vb += (ws[edges[i][1]] - vb)*blend;
    va.normalize();
    vb.normalize();
    Vector3d sid = va.cross(vb).normalized();
    double rad = 0.02;
    if (i > 5)
      rad = 0.01;
    if (abs(normal.dot(sid)) > rad)
      continue;
    Vector3d va2 = va.cross(sid);
    Vector3d vb2 = vb.cross(sid);
    if (normal.dot(va2) * vb.dot(va2) < 0.0)
      continue;
    if (normal.dot(vb2) * va.dot(vb2) < 0.0)
      continue;
    col *= 0.25;
  }

  return col;
}

int _tmain(int argc, _TCHAR* argv[])
{
  long s2;
  vector<BYTE> out(width*width * 3); // .bmp pixel buffer

  Vector3d camdir(0, 1, -1.0);
  camdir.normalize();
  Vector3d pos = -camdir * 15.0;
  Vector3d up(0, 0, 1);
  Vector3d side = camdir.cross(up).normalized();
  up = side.cross(camdir);
  double fov = 0.136;

  LPCTSTR files[] = { L"spherical_tetrahedronb1.bmp", L"spherical_tetrahedronb2.bmp", L"spherical_tetrahedronb3.bmp" };
  for (int k = 0; k < 3; k++)
  {
    double blend = (double)k / 2.0;
    memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is white
    for (int i = 0; i < width; i++)
    {
      double x = fov*(double)(i - width / 2) / (double)width;
      for (int j = 0; j < height; j++)
      {
        double y = fov*(double)(j - width / 2) / (double)width;
        Vector3d dir = (side*x + up*y + camdir).normalized();
        Vector3d frontCol = getCol(pos, dir, blend, true);
        Vector3d backCol = getCol(pos, dir, blend, false);
        setPixel(out, Vector2i(i, j), backCol + (frontCol - backCol)*0.8);
      }
    }

    BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    SaveBMP(c, width, height, s2, files[k]);
    delete[] c;
  }
}