#include "stdafx.h"
#include "bmp.h"
#include <fstream>
struct Section
{
  Vector2d pos;
  Vector2d up;
  void expand();
  void draw(ofstream &svg, const Vector2d &origin);
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

double scale = 760.0;
Vector2d offset(1.0, 0.5);
static vector<Vector2d> tri;

void saveSVG(const string &fileName, const vector<Section> &list)
{
  static ofstream svg;
  svg.open(fileName.c_str());
  svg << "<svg width = \"" << (int)(2.0*scale) << "\" height = \"" << (int)scale << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;

  for (auto &sec : list)
  {
    Vector2d side(sec.up[1], -sec.up[0]);
    for (int i = 0; i < (int)tri.size(); i++)
    {
      Vector2d end = sec.pos + tri[i][0] * side + tri[i][1] * sec.up;
      svg << "<path d = \"M " << scale*(sec.pos[0] + offset[0]) << " " << scale - scale*(sec.pos[1] + offset[1]);
      svg << " L " << scale*(end[0] + offset[0]) << " " << scale - scale*(end[1] + offset[1]);
      svg << "\" stroke=\"black\" stroke-width=\"1\" />\n"; // stroke=\"black\"
    }
  }

  svg << "</svg>" << endl;
  svg.close();
}

vector<Section> transform(vector<Section> &list, const Vector2d &translation)
{
  vector<Section> newChild = list;
  for (int i = 0; i < (int)newChild.size(); i++)
  {
    newChild[i].pos = -newChild[i].pos + translation;
    newChild[i].up = -newChild[i].up;
  }
  return newChild;
}

int _tmain(int argc, _TCHAR* argv[])
{
  int sides = 5;
  for (int i = 0; i < sides; i++)
  {
    double angle = (double)i * 2.0*pi / (double)sides;
    tri.push_back(2.0 * Vector2d(sin(angle), cos(angle)));
  }
  vector<Section> list(1);
  list[0].pos = Vector2d(0, 0);
  double scl = 0.005;// 1;
  list[0].up = Vector2d(0, scl);
  vector<Vector2d> ends(sides);
  for (int i = 0; i<sides; i++)
    ends[i] = tri[i] * scl;
  string type = "cross"; // "dendrite"; // or "cross"
  double s = type == "dendrite" ? 2.0 : 1.0; 
  for (int j = 0; j < 7; j++)
  {
    vector<Section> newList;
    vector<Vector2d> newEnds(sides);
    for (int i = 0; i < sides; i++)
    {
      vector<Section> newChild = transform(list, s*ends[i]);
      newList.insert(newList.end(), newChild.begin(), newChild.end());
      newEnds[i] = -ends[(i + (sides/2)) % sides] + s*ends[i];
    }
    if (type == "dendrite")
      list.insert(list.end(), newList.begin(), newList.end());
    else
      list = newList;
    for (int i = 0; i < sides; i++)
      ends[i] = newEnds[i];
  }
  saveSVG("pentacross.svg", list);
}
