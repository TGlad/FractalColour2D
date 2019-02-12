#include "stdafx.h"
#include "bmp.h"
#include <sstream>
#include <fstream>
static int width;
static int height;

void setPixel(vector<BYTE> &out, int x, int y, const Vector3i &col)
{
  if (x < 0 || x >= width || y < 0 || y >= height)
    return;
  int ind = 3 * (x + width * (height - 1 - y));
  out[ind + 0] = col[0];
  out[ind + 1] = col[1];
  out[ind + 2] = col[2];
}

void incrementPixel(vector<BYTE> &out, int x, int y, int red, int green, int blue)
{
  if (x < 0 || x >= width || y < 0 || y >= height)
    return;
  int ind = 3 * (x + width * (height - 1 - y));
  out[ind + 0] = min(255, out[ind + 0] + red);
  out[ind + 1] = min(255, out[ind + 1] + green);
  out[ind + 2] = min(255, out[ind + 2] + blue);
}

static const int NX = 28*2;
static const int NY = 7*2;
static const int M = 20;
class BlueGreen
{
public:
  double blues[M][M];
  double greens[M][M];
  BlueGreen operator+(const BlueGreen &b)
  {
    BlueGreen res;
    for (int i = 0; i < M; i++)
      for (int j = 0; j < M; j++)
      {
        res.blues[i][j] = blues[i][j] + b.blues[i][j];
        res.greens[i][j] = greens[i][j] + b.greens[i][j];
      }
    return res;
  }
  BlueGreen operator*(double s)
  {
    BlueGreen res;
    for (int i = 0; i < M; i++)
      for (int j = 0; j < M; j++)
      {
        res.blues[i][j] = blues[i][j] * s;
        res.greens[i][j] = greens[i][j] * s;
      }
    return res;
  }
};


const int specSize = 1000;
static Vector3i spectrum[specSize];
static int choices[16];
static double minX = 0.0;
static double maxX = 0.875;
static double minY = -0.5;
static double maxY = -0.25;
static double threshold = 3.0;
static int K = 0;
static double scale = 2.0;
static Vector2d flips[2];
static Vector2d bends[2];
static Vector2d shifts[2];
static double windowRadius = 1.0;
static BlueGreen bluegreen;
static BlueGreen bluegreens[NX][NY];
static double cx = 0.0;
static double cy = 0.0;
static int strength = 18;

void init()
{
  for (int i = 0; i < specSize; i++)
  {
    double f = (log((double)(i + 8)) - log(8.0)) / log((double)(specSize + 7));
    double f2 = 1.0 - (1.0 - f)*(1.0 - f);

    double red = 0.5*(cos(f * pi) + 1.0);
    double green = 0.375*(1.0f - cos(f * 2.0 * pi));
    double blue = 0.5*(cos((1.0 - f) * pi) + 1.0);

    double scale = 255.0 * f2;
    spectrum[i][0] = (int)(red*scale);
    spectrum[i][1] = (int)(green*scale);
    spectrum[i][2] = (int)(blue*scale);
  }
  flips[0] = Vector2d(1.0, 0);
  flips[1] = Vector2d(1.0, 0.0);
  bends[0] = Vector2d(1, 1) / sqrt(2.0);
  bends[1] = Vector2d(-1, 1) / sqrt(2.0);

  //  shifts[0] = Vector3(0, -0.15f, 0); //  // -0.15 or -0.2, maybe -0.2
  //  shifts[1] = Vector3(0, -0.15f, 0);
  shifts[0] = Vector2d(0, 0); //  // -0.15 or -0.2, maybe -0.2
  shifts[1] = Vector2d(0, 0);
}


Vector2d distort(const Vector2d &point, int i)
{
  Vector2d pos = point - bends[i];
  pos *= bends[i].squaredNorm() / pos.squaredNorm();
  pos += bends[i];
  pos = (shifts[i] + pos - (2.0 * flips[i] * pos.dot(flips[i]))) * scale;
  return pos;
}

void recurseDraw(vector<BYTE> &out, const Vector2d &point, const Vector2d &offset, int depth)
{
  int i = choices[16 - depth];

  int X = round(((double)width - 1.0)*(point[0] + windowRadius) / (2.0*windowRadius));
  int Y = round(((double)height - 1.0)*(point[1] + windowRadius) / (2.0*windowRadius));
  int I = (double)M * 0.5*(1.0 + point[0] / windowRadius);
  int J = (double)M * 0.5*(1.0 + point[1] / windowRadius);
  if (i == 0)
    bluegreen.greens[max(0, min(I, M - 1))][max(0, min(J, M - 1))]++;
  else
    bluegreen.blues[max(0, min(I, M - 1))][max(0, min(J, M - 1))]++;
  if (X >= 0 && X < width && Y >= 0 && Y < height)
  {
    incrementPixel(out, X, Y, 0, i == 0 ? 1 : 0, i == 1 ? 1 : 0);
  }
  Vector2d pos = distort(point, i) + offset;
  if (depth > 5)
    recurseDraw(out, pos, offset, depth - 1);
}

void recurse(vector<BYTE> &out, int &count, const Vector2d &point, const Vector2d &offset, int depth, const Vector2d &point0, bool draw)
{
  for (int i = 0; i<2; i++)
  {
    Vector2d pos = distort(point, i) + offset;
    double dist2 = pos.squaredNorm();
    if (dist2 > threshold)
      continue;
    choices[16 - depth] = i;
    if (depth == 0)
    {
      if (draw)
        recurseDraw(out, point0, offset, 16);
      count++;
    }
    else
      recurse(out, count, pos, offset, depth - 1, point0, draw);
  }
}

void generateFullMap(vector<BYTE> &out)
{
  for (int Y = 0; Y<height; Y++)
  {
    double y = minY + (maxY - minY)*(double)Y / (double)(height - 1);
    for (int X = 0; X<width; X++)
    {
      double x = minX + (maxX - minX)*(double)X / (double)(width - 1);
      int count = 0;
      recurse(out, count, Vector2d(x, y), Vector2d(x, y), 16, Vector2d(x, y), false);
      if (count == 0)
        setPixel(out, X, Y, Vector3i(0,0,0));
      else
        setPixel(out, X, Y, spectrum[min((count * strength) / 10, specSize - 1)]);
    }
  }
}

void recurseFast(vector<BYTE> &out, int &count, const Vector2d &point, const Vector2d &offset, int depth)
{
  double a = (double)(M - 1)*(windowRadius + point[1]) / (2.0*windowRadius);
  double b = (double)M - 1.0 - 1e-10;
  double y = max(0.0, min(a, b));
  int j = (int)y;
  double blendY = y - (double)j;

  a = (double)(M - 1)*(windowRadius + point[0]) / (2.0*windowRadius);
  double x = max(0.0, min(a, b));
  int i = (int)x;
  double blendX = x - (double)i;
  double green = bluegreen.greens[i][j] * (1.0 - blendX)*(1.0 - blendY) + bluegreen.greens[i][j + 1] * (1.0 - blendX)*(blendY)+bluegreen.greens[i + 1][j] * (blendX)*(1.0 - blendY) + bluegreen.greens[i + 1][j + 1] * (blendX)*(blendY);
  double blue = bluegreen.blues[i][j] * (1.0 - blendX)*(1.0 - blendY) + bluegreen.blues[i][j + 1] * (1.0 - blendX)*(blendY)+bluegreen.blues[i + 1][j] * (blendX)*(1.0 - blendY) + bluegreen.blues[i + 1][j + 1] * (blendX)*(blendY);

  i = green > blue ? 0 : 1;

  Vector2d pos = distort(point, i) + offset;
  double dist2 = pos.squaredNorm();
  if (dist2 > threshold)
    return;
  if (depth == 0)
    count++;
  else
    recurseFast(out, count, pos, offset, depth - 1);
}


void generateJuliaMultiset(vector<BYTE> &out)
{
  cout << "generating Julia set. cx: " << cx << ", cy: " << cy << endl;
  for (int Y = 0; Y<height; Y++)
  {
    double y = windowRadius * (-1.0f + 2.0*(double)Y / (double)(height - 1));
    //   double y = minY + (maxY - minY)*(double)Y / (double)(size - 1);
    for (int X = 0; X<width; X++)
    {
      double x = windowRadius * (-1.0f + 2.0f*(double)X / (double)(width - 1));
      //     double x = minX + (maxX - minX)*(double)X / (double)(size - 1);
      int count = 0;
      recurse(out, count, Vector2d(x, y), Vector2d(cx, cy), 16, Vector2d(x, y), false);
      if (count == 0)
        setPixel(out, X, Y, Vector3i(0, 0, 0));
      else
        setPixel(out, X, Y, spectrum[min((count * strength) / 10, specSize - 1)]);
    }
  }
  cout << "finished generating Julia set" << endl;
}

void generateBuddhaSet(vector<BYTE> &out, int type)
{
  cout << "generating Buddha set" << endl;
  if (type == 0 || type == 2 || type == 3)
  {
    for (int i = 0; i < M; i++)
    {
      for (int j = 0; j < M; j++)
      {
        bluegreen.blues[i][j] = 0;
        bluegreen.greens[i][j] = 0;
      }
    }
    for (int Y = 0; Y < height; Y++)
    {
      double y = windowRadius * (-1.0f + 2.0*(double)Y / (double)(height - 1));
      //   double y = minY + (maxY-minY)*(double)Y / (double)(size - 1);
      for (int X = 0; X < width; X++)
      {
        double x = windowRadius * (-1.0f + 2.0f*(double)X / (double)(width - 1));
        //      double x = minX + (maxX - minX)*(double)X / (double)(size - 1);
        int count = 0;
        recurse(out, count, Vector2d(x, y), Vector2d(cx, cy), 16, Vector2d(x, y), true);
      }
    }
    // disperse
    double dispersalRate = 0.1;
    for (int it = 0; it < M; it++)
    {
      BlueGreen temp = bluegreen;
      for (int i = 0; i < M; i++)
      {
        for (int j = 0; j < M; j++)
        {
          double count = 0;
          double averageBlue = 0;
          double averageGreen = 0;
          for (int x = -1; x < 2; x++)
          {
            for (int y = -1; y < 2; y++)
            {
              int X = i + x;
              int Y = j + y;
              if (X >= 0 && X < M && Y >= 0 && Y < M)
              {
                averageBlue += bluegreen.blues[X][Y];
                averageGreen += bluegreen.greens[X][Y];
                count++;
              }
            }
          }
          temp.blues[i][j] += ((averageBlue / count) - temp.blues[i][j])*dispersalRate;
          temp.greens[i][j] += ((averageGreen / count) - temp.greens[i][j])*dispersalRate;
        }
      }
      bluegreen = temp;
    }
  }
  if (type == 0 || type == 3)
    return;
  for (int Y = 0; Y < height; Y++)
  {
    double y = max(0.0, min(-0.5 + ((double)Y * (double)M / (double)height), (double)M - 1.0 - 1e-10));
    int j = (int)y;
    double blendY = y - (double)j;
    for (int X = 0; X < width; X++)
    {
      double x = max(0.0, min(-0.5 + ((double)X * (double)M / (double)width), (double)M - 1.0 - 1e-10));
      int i = (int)x;
      double blendX = x - (double)i;
      double green = bluegreen.greens[i][j] * (1.0 - blendX)*(1.0 - blendY) + bluegreen.greens[i][j + 1] * (1.0 - blendX)*(blendY)+bluegreen.greens[i + 1][j] * (blendX)*(1.0 - blendY) + bluegreen.greens[i + 1][j + 1] * (blendX)*(blendY);
      double blue = bluegreen.blues[i][j] * (1.0 - blendX)*(1.0 - blendY) + bluegreen.blues[i][j + 1] * (1.0 - blendX)*(blendY)+bluegreen.blues[i + 1][j] * (blendX)*(1.0 - blendY) + bluegreen.blues[i + 1][j + 1] * (blendX)*(blendY);

      int inc = 40;
      if (type == 1)
        incrementPixel(out, X, Y, 0, green / 180.0, blue / 180.0);
      else
        incrementPixel(out, X, Y, 0, green > blue ? inc : 0, blue > green ? inc : 0);
    }
  }

  cout << "finished generating buddha set" << endl;
}
void generateJuliaSet(vector<BYTE> &out)
{
  generateBuddhaSet(out, 3);
  for (int Y = 0; Y<height; Y++)
  {
    double y = windowRadius * (-1.0f + 2.0*(double)Y / (double)(height - 1));
    for (int X = 0; X<width; X++)
    {
      double x = windowRadius * (-1.0f + 2.0f*(double)X / (double)(width - 1));
      // so iterate fractal for this position
      int count = 0;
      recurseFast(out, count, Vector2d(x, y), Vector2d(cx, cy), 16);
      if (count > 0)
        setPixel(out, X, Y, Vector3i(0, 0, 0));
    }
  }
}
void generateMetaSet(vector<BYTE> &out)
{
  int oldheight = height;
  int oldwidth = width;
  height = width = 64;
  // precalculate separation data
  for (int I = 0; I < NX; I++)
  {
    cx = minX + (maxX - minX)*(double)I / (double)(NX - 1);
    for (int J = 0; J < NY; J++)
    {
      cy = minY + (maxY - minY)*(double)J / (double)(NY - 1);
      cout << "generating set matrix " << I << ", " << J << ", cx: " << cx << ", cy: " << cy << endl;
      generateBuddhaSet(out, 3);
      bluegreens[I][J] = bluegreen;
    }
  }
  height = oldheight;
  width = oldwidth;

  int numSet = 0;
  for (int Y = 0; Y<height; Y++)
  {
    double y = minY + (maxY - minY)*(double)Y / (double)(height - 1);
    double y2 = (double)(NY - 1)*(double)Y / (double)height;
    int j = y2;
    double blendY = y2 - (double)j;

    for (int X = 0; X<width; X++)
    {
      double x = minX + (maxX - minX)*(double)X / (double)(width - 1);

      // interpolate the bluegreens....
      double x2 = (double)(NX - 1)*(double)X / (double)width;
      int i = x2;
      double blendX = x2 - (double)i;
      bluegreen = bluegreens[i][j] * (1.0 - blendX)*(1.0 - blendY) + bluegreens[i][j + 1] * (1.0 - blendX)*(blendY)+bluegreens[i + 1][j] * (blendX)*(1.0 - blendY) + bluegreens[i + 1][j + 1] * (blendX)*(blendY);

      // so iterate fractal for this position
      int count = 0;
      recurseFast(out, count, Vector2d(x, y), Vector2d(x, y), 16);
      if (count > 0)
      {
        setPixel(out, X, Y, Vector3i(0,0,0));
        numSet++;
      }
    }
  }
  cout << "number of pixels in final output: " << numSet << endl;
}


int _tmain(int argc, _TCHAR* argv[])
{
  init();
  long s2;
  bool on = false;
  if (on)
  {
    width = NX * 64;
    height = NY * 64;
    vector<BYTE> out(width*height * 3); // .bmp pixel buffer
    memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey

    generateFullMap(out);

    BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    LPCTSTR file = L"mobius1.bmp";
    SaveBMP(c, width, height, s2, file);
    delete[] c;
  }
  double cxs[] = { 0.05, 0.325, 0.6 };
  double cys[] = { -0.5, -0.35 };
  double windowSizes[] { 1,0.5, 1.2, 0.7, 1, 1};
  if (on)
  {
    int count = 0;
    for (int i = 0; i < 3; i++)
    {
      cx = cxs[i];
      for (int j = 0; j < 2; j++)
      {
        cy = cys[j];
        width = 2 * 256;
        height = 2 * 256;
        windowRadius = windowSizes[count];
        vector<BYTE> out(width*height * 3); // .bmp pixel buffer
        memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey

        generateJuliaMultiset(out);

        BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
        LPCTSTR file = L"";
        if (count == 0)
          file = L"mobius_0.bmp";
        else if (count == 1)
          file = L"mobius_1.bmp";
        else if (count == 2)
          file = L"mobius_2.bmp";
        else if (count == 3)
          file = L"mobius_3.bmp";
        else if (count == 4)
          file = L"mobius_4.bmp";
        else if (count == 5)
          file = L"mobius_5.bmp";
        count++;
        SaveBMP(c, width, height, s2, file);
        delete[] c;
      }
    }
  }
  windowRadius = windowSizes[1];
  if (on)
  {
    cx = cxs[0];
    cy = cys[1];
    width = 2 * 256;
    height = 2 * 256;
    vector<BYTE> out(width*height * 3); // .bmp pixel buffer
    memset(&out[0], 0, out.size() * sizeof(BYTE)); // background is grey

    generateBuddhaSet(out, 0);

    BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    LPCTSTR file = L"mobius2.bmp";
    SaveBMP(c, width, height, s2, file);
    delete[] c;
  }
  if (on)
  {
    cx = cxs[0];
    cy = cys[1];
    width = 2 * 256;
    height = 2 * 256;
    vector<BYTE> out(width*height * 3); // .bmp pixel buffer
    memset(&out[0], 0, out.size() * sizeof(BYTE)); // background is grey

    generateBuddhaSet(out, 1);

    BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    LPCTSTR file = L"mobius3.bmp";
    SaveBMP(c, width, height, s2, file);
    delete[] c;
  }
  if (on)
  {
    cx = cxs[0];
    cy = cys[1];
    width = 2 * 256;
    height = 2 * 256;
    vector<BYTE> out(width*height * 3); // .bmp pixel buffer
    memset(&out[0], 0, out.size() * sizeof(BYTE)); // background is grey

    generateBuddhaSet(out, 2);

    BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    LPCTSTR file = L"mobius4.bmp";
    SaveBMP(c, width, height, s2, file);
    delete[] c;
  }
  if (on)  // julia set
  {
    cx = cxs[0];
    cy = cys[1];
    width = 2 * 256;
    height = 2 * 256;
    vector<BYTE> out(width*height * 3); // .bmp pixel buffer
    memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey

    generateJuliaSet(out);

    BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    LPCTSTR file = L"mobius5.bmp";
    SaveBMP(c, width, height, s2, file);
    delete[] c;
  }
  if (1)
  {
    width = NX * 32;
    height = NY * 32;
    vector<BYTE> out(width*height * 3); // .bmp pixel buffer
    memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey

    generateMetaSet(out);

    BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    LPCTSTR file = L"mobius6.bmp";
    SaveBMP(c, width, height, s2, file);
    delete[] c;
  }
}
