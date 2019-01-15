#include "stdafx.h"
#include "bmp.h"
#include <fstream>
struct Segment
{
  Vector2d start;
  double angle;
  double length;
  Vector2d end() const { return start + length * Vector2d(sin(angle), cos(angle)); }
};
static double minLength = 0.01;

void recurse(vector<Segment> &segments, const Segment &segment)
{
  segments.push_back(segment);
  if (segment.length < minLength)
    return;
  Segment left;
  double scale = 12.0/5.0;
  left.start = segment.end();
  left.angle = segment.angle - pi / 4.0;
  left.length = segment.length / scale;
  recurse(segments, left);

  Segment right;
  right.start = segment.end();
  right.angle = segment.angle + pi / 4.0;
  right.length = segment.length / scale;
  recurse(segments, right);
}

int _tmain(int argc, _TCHAR* argv[])
{
  vector<Segment> segments;
  Segment first;
  first.start = Vector2d(0,0); 
  first.angle = 0;
  first.length = 1.0;
  recurse(segments, first);

  double minX = 1e10;
  double maxX = -1e10;
  for (auto &segment : segments)
  {
    Vector2d end = segment.end();
    minX = min(minX, end[0]);
    maxX = max(maxX, end[0]);
  }
  cout << "length: " << maxX - minX << endl;
}
