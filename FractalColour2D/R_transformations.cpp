#include "stdafx.h"
#include "bmp.h"
#include <fstream>
#include <complex>
typedef complex<double> Complex;

double width = 700;
double height = 700;

void saveSVG(const string &fileName, vector<Vector2d> list)
{
  for (auto &v : list)
    v = width*(Vector2d(3.5, 3.5) + v)/7.0;
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)width << "\" height = \"" << (int)height << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;

  svg << "<path d = \"M " << list[0][0] << " " << height - list[0][1];
  for (int i = 0; i < list.size(); i++)
    svg << " L " << list[i][0] << " " << height - list[i][1];
  svg << "\" fill=\"lightgrey\" stroke-width = \"2\" stroke=\"black\" />\n";

  svg << "<circle cx = \"" << width/2.0 << "\" cy = \"" << height - height/2.0 << "\" r = \"4\" stroke = \"black\" stroke-width = \"1\" fill = \"black\" />" << endl;

  svg << "</svg>" << endl;
  svg.close();
}

vector<Vector2d> makeR()
{
  vector<Vector2d> list;
  vector<Vector2d> corners(6);
  corners[0] = Vector2d(-1.5, -2.5);
  corners[1] = Vector2d(-1.5, 2.5);
  corners[2] = Vector2d(-0.5, 2.5);
  corners[3] = Vector2d(-0.5, -1.5);
  corners[4] = Vector2d(1.5, -1.5);
  corners[5] = Vector2d(1.5, -2.5);
  for (int i = 0; i < corners.size(); i++)
  {
    int j = (i + 1) % 6;
    double dist = (corners[i] - corners[j]).norm();
    for (double x = 0; x < dist; x += 1.0 / 20.0)
    {
      list.push_back(corners[i] + (corners[j] - corners[i])*x / dist);
    }
  }
  return list;
}

int _tmain(int argc, _TCHAR* argv[])
{
  vector<Vector2d> list = makeR();

  saveSVG("R_normal.svg", list);

  vector<Vector2d> reflection = list;
  for (auto &v : reflection)
    v[0] = -v[0];
  saveSVG("R_reflection.svg", reflection);

  vector<Vector2d> translation = list;
  for (auto &v : translation)
    v += Vector2d(1.8, 0.0);
  saveSVG("R_translation.svg", translation);

  vector<Vector2d> rotation = list;
  Matrix2d rot;
  double theta = 1.0;
  rot << cos(theta), sin(theta),
        -sin(theta), cos(theta);
  for (auto &v : rotation)
    v = rot*v;
  saveSVG("R_rotation.svg", rotation);

  vector<Vector2d> scale = list;
  for (auto &v : scale)
    v *= 1.4;
  saveSVG("R_scale.svg", scale);

  vector<Vector2d> euclidean = list;
  for (auto &v : euclidean)
  {
    v = rot*v;
    v *= 1.2;
    v += Vector2d(1.2, 0.0);
  }
  saveSVG("R_Euclidean.svg", euclidean);

  vector<Vector2d> inversion = list;
  for (auto &v : inversion)
  {
    v *= 1.4 / v.squaredNorm();
    v[0] = -v[0];
  }
  saveSVG("R_inversion.svg", inversion);

  vector<Vector2d> mobius = list;
  for (auto &v : mobius)
  {
    v[0] -= 1.0;
    v = rot * v;
    v *= 5.0 / v.squaredNorm();
    v[0] = -v[0];
  }
  saveSVG("R_mobius.svg", mobius);

  vector<Vector2d> compl = list;
  for (auto &v : compl)
  {
    Complex z(v[0], v[1]);
    z = 0.3*pow(z, 2.0) - 0.4*z;
    v[0] = z.real();
    v[1] = z.imag();
  }
  saveSVG("R_complex.svg", compl);

  vector<Vector2d> conformal = list;
  for (auto &v : conformal)
  {
/*    v[0] -= 1.0;
    v = rot * v;
    v *= 5.0 / v.squaredNorm();
    v[0] = -v[0];
    Complex z(v[0], v[1]);
    z += 0.4*pow(z, 2.0);
    z /= 2.0;
    v[0] = z.real();
    v[1] = z.imag();*/
    v[0] -= 2.0;
    v = rot * v;
    v *= 10.0 / v.squaredNorm();
    v[0] = -v[0];
    v[0] -= 2.7;
    v[1] -= 1.0;
  }
  saveSVG("R_conformal.svg", conformal);

  vector<Vector2d> shapepreserving = list;
  for (auto &v : shapepreserving)
  {
    v[0] -= 2.0;
    v = rot * v;
    v *= 10.0 / v.squaredNorm();
    v[0] = -v[0];
    v[0] -= 2.7;
    v[1] -= 1.0;
    if (v[1] < -0.8)
      v[1] = -0.8 + (-0.8 - v[1]);
  }
  saveSVG("R_shapepreserving.svg", shapepreserving);
}
