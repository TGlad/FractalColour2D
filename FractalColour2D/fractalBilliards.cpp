// experiments with structural coloration of fractal curves, in particular Koch, Levy, Dragon and random curves.
#include "stdafx.h"
#include "bmp.h"
#include "spectrumToRGB.h"
#include <set>

static int width = 512;
static int height = 512;

void putpixel(vector<BYTE> &out, const Vector2i &pos, int shade)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= width)
    return;
  int ind = 3 * (pos[0] + width*pos[1]);
  out[ind + 0] = out[ind + 1] = out[ind + 2] = shade;
}
int getpixel(vector<BYTE> &out, const Vector2i &pos)
{
  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= width)
    return 0;
  int ind = 3 * (pos[0] + width*pos[1]);
  return out[ind + 0];
}

Vector2i toVector2i(const Vector2d &pos)
{
  Vector2i vec;
  vec[0] = (int)pos[0];
  if ((double)vec[0] > pos[0])
    vec[0]--;
  vec[1] = (int)pos[1];
  if ((double)vec[1] > pos[1])
    vec[1]--;
  return vec;
}

static double ang = atan(4.0);
static double l = sin(ang);
static double w = cos(ang);
static double speed = 100.0;
static double radius = 24.0;
struct Ball
{
  Vector2d pos0;
  Vector2d fwd;
  bool flip;
  double scale;
  double time0;
};
#define SYMMETRIC
#define FLIP
void addBalls(vector<Ball> &balls, Ball &ball)
{
  balls.push_back(ball);
  if (ball.scale < 1.0/8.0)
    return;
  // we need to add new balls that will collide with this in the time period 0 to 1
  Vector2d side(-ball.fwd[1], ball.fwd[0]);
  if (ball.flip)
    side = -side;
  double a = 2.0*(pi / 2.0 - ang);
  Ball ball0;
  ball0.pos0 = ball.pos0 + ball.fwd*speed*ball.scale;
  ball0.scale = ball.scale / 2.0;
#if defined FLIP
  ball0.flip = !ball.flip;
#else
  ball0.flip = ball.flip;
#endif
  ball0.time0 = ball.time0 + 1.0*ball.scale;
  double b = (pi + a) / 2.0;
  ball0.fwd = -ball.fwd*cos(b) - side*sin(b);
 
  Ball ball1;
  ball1.pos0 = ball0.pos0 + 2.0*(speed + radius)*ball.scale*(ball.fwd*cos(a) + side*sin(a));
  ball1.scale = ball.scale / 2.0;
#if defined FLIP
  ball1.flip = ball.flip;
#else
  ball1.flip = !ball.flip;
#endif
  ball1.time0 = ball.time0 + 1.0*ball.scale;
#if defined SYMMETRIC
  double c = b + a;
  ball1.fwd = -ball.fwd*cos(c) - side*sin(c);
#else
  ball1.fwd = ball.fwd*cos(b) + side*sin(b);
#endif

  ball0.pos0 += ball0.fwd * (ball.scale + ball0.scale)*radius;
  ball1.pos0 += ball1.fwd * (ball.scale + ball1.scale)*radius;

  addBalls(balls, ball0);
  addBalls(balls, ball1);
}
void drawDisk(const Vector2d &pos, vector<BYTE> &out, double rad, int shade)
{
  for (int x = (int)(pos[0] - rad); x <= (int)(pos[0] + rad); x++)
    for (int y = (int)(pos[1] - rad); y <= (int)(pos[1] + rad); y++)
      if (sqr(x - pos[0]) + sqr(y - pos[1]) <= sqr(rad))
        putpixel(out, Vector2i(x, y), shade);
}
int _tmain(int argc, _TCHAR* argv[])
{
  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  memset(&out[0], 0, out.size() * sizeof(BYTE)); // background is grey

  vector<Ball> balls;
  Ball ball;
  ball.pos0 = Vector2d(255.0, 25.0);
  ball.fwd = Vector2d(0, 1);
  ball.scale = 1.0;
  ball.time0 = 0;
  ball.flip = false;

  addBalls(balls, ball);

  // now draw the balls
  for (double time = 0; time < 2.5; time += 0.08)
  {
    for (auto &ball : balls)
    {
      double t = (time - ball.time0) / ball.scale;
      Vector2d side(-ball.fwd[1], ball.fwd[0]);
      if (ball.flip)
        side = -side;

      Vector2d p0 = ball.pos0;
      Vector2d p1 = ball.pos0 + ball.fwd * ball.scale * speed;
      double a = 2.0*(pi / 2.0 - ang);
      Vector2d p2 = p1 + (ball.fwd*cos(a) + side*sin(a)) * ball.scale * speed;

      t = fmod(t + 4000.0, 4.0);
      Vector2d p;
      if (t < 1.0)
        p = p0 + (p1 - p0)*t;
      else if (t < 2.0)
        p = p1 + (p2 - p1)*(t - 1.0);
      else if (t < 3.0)
        p = p2 + (p1 - p2)*(t - 2.0);
      else
        p = p1 + (p0 - p1)*(t - 3.0);
      drawDisk(p, out, ball.scale * radius, (int)(time * 63));

      p2 += (ball.fwd*cos(a) + side*sin(a)) * ball.scale * 2.0 * radius;
      p1 = p2 + (ball.fwd*cos(a) + side*sin(a)) * ball.scale * speed;
#if defined SYMMETRIC
      p0 = p1 + (ball.fwd*cos(2.0*a) + side*sin(2.0*a)) * ball.scale * speed;
#else
      p0 = p1 + ball.fwd * ball.scale * speed;
#endif

      if (t < 1.0)
        p = p0 + (p1 - p0)*t;
      else if (t < 2.0)
        p = p1 + (p2 - p1)*(t - 1.0);
      else if (t < 3.0)
        p = p2 + (p1 - p2)*(t - 2.0);
      else
        p = p1 + (p0 - p1)*(t - 3.0);
      drawDisk(p, out, ball.scale * radius, (int)(time * 63));

      p = p0 + (ball.fwd*cos(2.0*a) + side*sin(2.0*a)) * ball.scale * (radius + radius/2.0);
      drawDisk(p, out, ball.scale * radius/2.0, 255);
    }
  }
//  putpixel(out, Vector2i(x, y)*w, start[y][x]);

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
#if defined SYMMETRIC
  LPCTSTR file = L"billiards_symmetric3.bmp";
#else
  LPCTSTR file = L"billiards2.bmp";
#endif
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}

