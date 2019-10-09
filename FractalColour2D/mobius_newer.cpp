#include "stdafx.h"
#include "bmp.h"
#include <sstream>
#include <fstream>
static int width;
static int height;
#include <iostream>

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

static const int NX = 28;// *2;
static const int NY = 7;// *2;
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
static const int maxIterations = 32;
static double startRadius = 1.0 / 1000.0;
static int choices[maxIterations];
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
static int strength = 6;

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


Vector2d distort(const Vector2d &point, int i, double &radius)
{
  Vector2d pos = point - bends[i];
  double mult = bends[i].squaredNorm() / pos.squaredNorm();
  pos *= mult;
  radius *= mult*scale;
  pos += bends[i];
  pos = (shifts[i] + pos - (2.0 * flips[i] * pos.dot(flips[i]))) * scale;
  return pos;
}

void recurseDraw(vector<BYTE> &out, const Vector2d &point, const Vector2d &offset, int depth, double radius)
{
  int i = choices[maxIterations - depth];

//  int X = round((point[0] - minX)* (double)(width - 1) / (maxX - minX));
//  int Y = round((point[1] - minY)* (double)(height - 1) / (maxY - minY));
  int X = round(((double)height - 1.0)*(point[0] + windowRadius) / (2.0*windowRadius));
  int Y = round(((double)height - 1.0)*(point[1] + windowRadius) / (2.0*windowRadius));
  int I = (double)M * 0.5*(1.0 + point[0] / windowRadius);
  int J = (double)M * 0.5*(1.0 + point[1] / windowRadius);
  if (I >= 0 && I < M && J >= 0 && J < M)
  {
    if (i==0)
      bluegreen.greens[max(0, min(I, M - 1))][max(0, min(J, M - 1))]++;
    else
      bluegreen.blues[max(0, min(I, M - 1))][max(0, min(J, M - 1))]++;
  }
  if (X >= 0 && X < width && Y >= 0 && Y < height)
  {
    incrementPixel(out, X, Y, 0, i == 0 ? 1 : 0, i == 1 ? 1 : 0);
  }
  Vector2d pos = distort(point, i, radius) + offset;
  if (depth > 5)
    recurseDraw(out, pos, offset, depth - 1, radius);
}

void recurse(vector<BYTE> &out, int &count, const Vector2d &point, const Vector2d &offset, int depth, const Vector2d &point0, bool draw, double radius)
{
  for (int i = 0; i<2; i++)
  {
    double rad = radius;
    Vector2d pos = distort(point, i, rad) + offset;
    double dist2 = pos.squaredNorm();
//    if (dist2 > threshold)
    if (pos[1]-radius >= 0.9)
      continue;
    choices[maxIterations - depth] = i;
    if (depth == 0)
    {
      if (draw)
        recurseDraw(out, point0, offset, maxIterations, rad);
      count++;
    }
    else
      recurse(out, count, pos, offset, depth - 1, point0, draw, rad);
  }
}

void generateFullMap(vector<BYTE> &out)
{
  for (int Y = 0; Y<height; Y++)
  {
    cout << "y: " << Y << endl;
    double y = minY + (maxY - minY)*(double)Y / (double)(height - 1);
    for (int X = 0; X<width; X++)
    {
      double x = minX + (maxX - minX)*(double)X / (double)(width - 1);
      int count = 0;
      recurse(out, count, Vector2d(x, y), Vector2d(x, y), maxIterations, Vector2d(x, y), false, startRadius);
      if (count == 0)
        setPixel(out, X, Y, Vector3i(0,0,0));
      else
        setPixel(out, X, Y, spectrum[min((count * strength) / 10, specSize - 1)]);
    }
  }
}

void recurseFast(vector<BYTE> &out, int &count, const Vector2d &point, const Vector2d &offset, int depth, double radius)
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

  Vector2d pos = distort(point, i, radius) + offset;
  double dist2 = pos.squaredNorm();
//  if (dist2 > threshold)
  if (pos[1]-radius > 0.9)
    return;
  if (depth == 0)
    count++;
  else
    recurseFast(out, count, pos, offset, depth - 1, radius);
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
      recurse(out, count, Vector2d(x, y), Vector2d(cx, cy), maxIterations, Vector2d(x, y), false, startRadius);
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
        recurse(out, count, Vector2d(x, y), Vector2d(cx, cy), maxIterations, Vector2d(x, y), true, startRadius);
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
      recurseFast(out, count, Vector2d(x, y), Vector2d(cx, cy), maxIterations, startRadius);
      if (count > 0)
        setPixel(out, X, Y, Vector3i(0, 0, 0));
    }
  }
}
void generateMetaSet(vector<BYTE> &out)
{
  double yy = minY + (maxY - minY)*(double)62 / (double)(height - 1);
  double xx = minX + (maxX - minX)*(double)1102 / (double)(width - 1);

  for (int Y = 0; Y < height; Y++)
  {
    double y = windowRadius * (-1.0f + 2.0*(double)Y / (double)(height - 1));
    //   double y = minY + (maxY-minY)*(double)Y / (double)(size - 1);
    for (int X = 0; X < height; X++)
    {
      double x = windowRadius * (-1.0f + 2.0f*(double)X / (double)(height - 1));
      //      double x = minX + (maxX - minX)*(double)X / (double)(size - 1);
      int count = 0;
      recurse(out, count, Vector2d(x, y), Vector2d(xx, yy), maxIterations, Vector2d(x, y), true, startRadius);
    }
  }
  if (0)
  {
    for (int i = 0; i < M; i++)
    {
      for (int j = 0; j < M; j++)
      {
        Vector3i col(0, 0, 255);
        col[1] = min(255, bluegreen.greens[i][j]);
        col[2] = min(255, bluegreen.blues[i][j]);
        setPixel(out, 150 + i, j, col);
      }
    }
  }
  /*  for (int Y = 0; Y < height; Y++)
  {
    double y = minY + (maxY - minY)*(double)Y / (double)(height - 1);
    // for (int X = 0; X<width; X++)
    for (int X = 0; X<width; X++)
    {
      double x = minX + (maxX - minX)*(double)X / (double)(width - 1);
      int count = 0;
      recurse(out, count, Vector2d(x, y), Vector2d(xx, yy), maxIterations, Vector2d(x, y), true);
   //   cout << "count: " << count << endl;
    }
  }*/
  
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

  if (0)
  {
    for (int I = 0; I < NX; I++)
    {
      for (int J = 0; J < NY; J++)
      {
        bluegreen = bluegreens[I][J];
        for (int i = 0; i < M; i++)
        {
          for (int j = 0; j < M; j++)
          {
            Vector3i col(0, 0, 255);
            if (bluegreen.greens[i][j] > bluegreen.blues[i][j])
              col = Vector3i(0, 255, 0);
            col[1] = min(255, 4 * bluegreen.greens[i][j]);
            col[2] = min(255, 4 * bluegreen.blues[i][j]);
            setPixel(out, I*(M + 2) + i, J*(M + 2) + j, col);
          }
        }
      }
    }
  }
  // return;
  int numSet = 0;
  for (int Y = 0; Y<height; Y++)
//  int Y = 62;
  {
    double y = minY + (maxY - minY)*(double)Y / (double)(height - 1);
    double y2 = (double)(NY - 1)*(double)Y / (double)height;
    int j = y2;
    double blendY = y2 - (double)j;

    for (int X = 0; X<width; X++)
 //   int X = 1102;
    {
      double x = minX + (maxX - minX)*(double)X / (double)(width - 1);

      // interpolate the bluegreens....
      double x2 = (double)(NX - 1)*(double)X / (double)width;
      int i = x2;
      double blendX = x2 - (double)i;
      bluegreen = bluegreens[i][j] * (1.0 - blendX)*(1.0 - blendY) + bluegreens[i][j + 1] * (1.0 - blendX)*(blendY)+bluegreens[i + 1][j] * (blendX)*(1.0 - blendY) + bluegreens[i + 1][j + 1] * (blendX)*(blendY);

      // so iterate fractal for this position
      int count = 0;
      recurseFast(out, count, Vector2d(x, y), Vector2d(x, y), maxIterations, startRadius);
      if (count > 0)
      {
        setPixel(out, X, Y, Vector3i(0, 0, 0));
        numSet++;
      }
      if (0)
      {
        //  cout << "count: " << count << endl;
        for (int I = i; I < i + 2; I++)
        {
          for (int J = j; J < j+2; J++)
          {
            for (int ii = 0; ii < M; ii++)
            {
              for (int jj = 0; jj < M; jj++)
              {
                Vector3i col(0, 0, 255);
                if (bluegreens[I][J].greens[ii][jj] > bluegreens[I][J].blues[ii][jj])
                  col = Vector3i(0, 255, 0);
                col[1] = min(255, bluegreens[I][J].greens[ii][jj]);
                col[2] = min(255, bluegreens[I][J].blues[ii][jj]);
                setPixel(out, (M + 2)*I + ii, 10 + (M + 2)*J + jj, col);
              }
            }
          }
        }
        for (int i = 0; i < M; i++)
        {
          for (int j = 0; j < M; j++)
          {
            Vector3i col(0, 0, 255);
            if (bluegreen.greens[i][j] > bluegreen.blues[i][j])
              col = Vector3i(0, 255, 0);
            col[1] = min(255, bluegreen.greens[i][j]);
            col[2] = min(255, bluegreen.blues[i][j]);
            setPixel(out, 200+i, j, col);
          }
        }
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
  if (true)
  {
    width = NX * 128;
    height = NY * 128;
    startRadius = (maxX - minX) / width;
    vector<BYTE> out(width*height * 3); // .bmp pixel buffer
    memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey

    generateFullMap(out);

    BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    LPCTSTR file = L"mobius1e.bmp";
    SaveBMP(c, width, height, s2, file);
    delete[] c;
    return 1;
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
  if (0)
  {
    width = NX * 32;
    height = NY * 32;
    vector<BYTE> out(width*height * 3); // .bmp pixel buffer
    memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey

    generateMetaSet(out);

    BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    LPCTSTR file = L"mobius11.bmp";
    SaveBMP(c, width, height, s2, file);
    delete[] c;
  }
  if (0)
  {
    width = 400;
    height = 400;
    vector<BYTE> out(width*height * 3); // .bmp pixel buffer
    memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey

    double rad = 4.0;
    for (int j = 0; j < 4; j++)
    {
      Vector2d c;
      if (j == 0)
        c = Vector2d(minX, minY);
      else if (j == 1)
        c = Vector2d(minX, maxY);
      else if (j == 2)
        c = Vector2d(maxX, minY);
      else
        c = Vector2d(maxX, maxY);
      for (int i = 0; i < 4096; i++)
      {
        Vector2d point = c;
        for (int it = 0; it < 200; it++)
        {
          int on = (i >> (it % 12)) & 1;
          double radius = 1.0;
          point = distort(point, on, radius) + c;
        }
        Vector2d pos = width * (Vector2d(1, 1) + point / rad) / 2.0;
        setPixel(out, pos[0], pos[1], Vector3i(255.0*(double)j / 2.0, 255.0*(double)(j % 2), 0));
      }
    }
    double h = width*(1.0 + 1.0 / rad) / 2.0;
    for (int x = 0; x < width; x++)
      setPixel(out, x, h, Vector3i(127, 127, 127));

    Vector2d minn = width * (Vector2d(1, 1) + Vector2d(minX, minY) / rad) / 2.0;
    Vector2d maxn = width * (Vector2d(1, 1) + Vector2d(maxX, maxY) / rad) / 2.0;
    for (int x = minn[0]; x <= maxn[0]; x++)
    {
      for (int y = minn[1]; y <= maxn[1]; y++)
      {
        setPixel(out, x, y, Vector3i(127, 127, 127));
      }
    }
    BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    LPCTSTR file = L"endPoints.bmp";
    SaveBMP(c, width, height, s2, file);
    delete[] c;
  }
  if (1)
  {
    width = 400;
    height = 400;
    vector<BYTE> out(width*height * 3); // .bmp pixel buffer
    memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is grey

    double rad = 4.0;
    for (int i = 0; i < 4096*16; i++)
    {
      Vector2d c(random(minX, maxX), random(minY, maxY));
      Vector2d point = c;
      for (int it = 0; it < 200; it++)
      {
        int on = rand()%2;
        double radius = 1.0;
        point = distort(point, on, radius) + c;
      }
      Vector2d pos = width * (Vector2d(1, 1) + point / rad) / 2.0;
      setPixel(out, pos[0], pos[1], Vector3i(0,0,0));
    }
    double h = width*(1.0 + 1.0 / rad) / 2.0;
    for (int x = 0; x < width; x++)
      setPixel(out, x, h, Vector3i(127, 127, 127));

    Vector2d minn = width * (Vector2d(1, 1) + Vector2d(minX, minY) / rad) / 2.0;
    Vector2d maxn = width * (Vector2d(1, 1) + Vector2d(maxX, maxY) / rad) / 2.0;
    for (int x = minn[0]; x <= maxn[0]; x++)
    {
      for (int y = minn[1]; y <= maxn[1]; y++)
      {
        setPixel(out, x, y, Vector3i(127, 127, 127));
      }
    }
    BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    LPCTSTR file = L"endPointsRandom2.bmp";
    SaveBMP(c, width, height, s2, file);
    delete[] c;
  }
}
