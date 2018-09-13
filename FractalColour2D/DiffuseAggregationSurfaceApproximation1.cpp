// Folding.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <vector>
#include "/Code/Eigen/Eigen"
#include "/Code/Eigen/StdVector"
using namespace std;
using namespace Eigen;
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Vector2d);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Vector3d);

/*************************************************************************

file created on:	2002/08/30   19:33
filename: 			Bmp.cpp
author:				Andreas Hartl

visit http://www.runicsoft.com for updates and more information

purpose:	functions to load raw bmp data,
to save raw bmp data,
to convert RGB data to raw bmp data,
to convert raw bmp data to RGB data
and to use the WinAPI to select
a bitmap into a device context

file updated on 2010/09/13
in the 8 years since i first wrote this the windows file functions
have changed their input from char* to LPCTSTR.
Updated this in the code here

**************************************************************************/

#include <windows.h>
#include <stdio.h>       // for memset


/*******************************************************************
BYTE* ConvertRGBToBMPBuffer ( BYTE* Buffer, int width,
int height, long* newsize )


This function takes as input an array of RGB values, it's width
and height.
The buffer gets then transformed to an array that can be used
to write to a windows bitmap file. The size of the array
is returned in newsize, the array itself is the
return value of the function.
Both input and output buffers must be deleted by the
calling function.

The input buffer is expected to consist of width * height
RGB triplets. Thus the total size of the buffer is taken as
width * height * 3.

The function then transforms this buffer so that it can be written
to a windows bitmap file:
First the RGB triplets are converted to BGR.
Then the buffer is swapped around since .bmps store
images uside-down.
Finally the buffer gets DWORD ( 32bit ) aligned,
meaning that each scanline ( 3 * width bytes ) gets
padded with 0x00 bytes up to the next DWORD boundary


*******************************************************************/

BYTE* ConvertRGBToBMPBuffer(BYTE* Buffer, int width, int height, long* newsize)
{

  // first make sure the parameters are valid
  if ((NULL == Buffer) || (width == 0) || (height == 0))
    return NULL;

  // now we have to find with how many bytes
  // we have to pad for the next DWORD boundary	

  int padding = 0;
  int scanlinebytes = width * 3;
  while ((scanlinebytes + padding) % 4 != 0)     // DWORD = 4 bytes
    padding++;
  // get the padded scanline width
  int psw = scanlinebytes + padding;

  // we can already store the size of the new padded buffer
  *newsize = height * psw;

  // and create new buffer
  BYTE* newbuf = new BYTE[*newsize];

  // fill the buffer with zero bytes then we dont have to add
  // extra padding zero bytes later on
  memset(newbuf, 0, *newsize);

  // now we loop trough all bytes of the original buffer, 
  // swap the R and B bytes and the scanlines
  long bufpos = 0;
  long newpos = 0;
  for (int y = 0; y < height; y++)
    for (int x = 0; x < 3 * width; x += 3)
    {
      bufpos = y * 3 * width + x;     // position in original buffer
      newpos = (height - y - 1) * psw + x;           // position in padded buffer

      newbuf[newpos] = Buffer[bufpos + 2];       // swap r and b
      newbuf[newpos + 1] = Buffer[bufpos + 1]; // g stays
      newbuf[newpos + 2] = Buffer[bufpos];     // swap b and r
    }

  return newbuf;
}

/***************************************************************
BYTE* ConvertBMPToRGBBuffer ( BYTE* Buffer,
int width, int height )

This function takes as input the data array
from a bitmap and its width and height.
It then converts the bmp data into an RGB array.
The calling function must delete both the input
and output arrays.
The size of the returned array will be
width * height * 3
On error the returb value is NULL, else the
RGB array.


The Buffer is expected to be the exact data read out
from a .bmp file.
The function will swap it around, since images
are stored upside-down in bmps.
The BGR triplets from the image data will
be converted to RGB.
And finally the function removes padding bytes.
The returned arraay consits then of
width * height RGB triplets.

*****************************************************************/

BYTE* ConvertBMPToRGBBuffer(BYTE* Buffer, int width, int height)
{
  // first make sure the parameters are valid
  if ((NULL == Buffer) || (width == 0) || (height == 0))
    return NULL;

  // find the number of padding bytes

  int padding = 0;
  int scanlinebytes = width * 3;
  while ((scanlinebytes + padding) % 4 != 0)     // DWORD = 4 bytes
    padding++;
  // get the padded scanline width
  int psw = scanlinebytes + padding;

  // create new buffer
  BYTE* newbuf = new BYTE[width*height * 3];

  // now we loop trough all bytes of the original buffer, 
  // swap the R and B bytes and the scanlines
  long bufpos = 0;
  long newpos = 0;
  for (int y = 0; y < height; y++)
    for (int x = 0; x < 3 * width; x += 3)
    {
      newpos = y * 3 * width + x;
      bufpos = (height - y - 1) * psw + x;

      newbuf[newpos] = Buffer[bufpos + 2];
      newbuf[newpos + 1] = Buffer[bufpos + 1];
      newbuf[newpos + 2] = Buffer[bufpos];
    }

  return newbuf;
}


/***********************************************
bool LoadBMPIntoDC ( HDC hDC, LPCTSTR bmpfile )

Takes in a device context and the name of a
bitmap to load. If an error occurs the function
returns false, else the contents of the bmp
are blitted to the HDC

************************************************/

bool LoadBMPIntoDC(HDC hDC, LPCTSTR bmpfile)
{
  // check if params are valid 
  if ((NULL == hDC) || (NULL == bmpfile))
    return false;

  // load bitmap into a bitmap handle
  HANDLE hBmp = LoadImage(NULL, bmpfile, IMAGE_BITMAP, 0, 0,
    LR_LOADFROMFILE);

  if (NULL == hBmp)
    return false;        // failed to load image

  // bitmaps can only be selected into memory dcs:
  HDC dcmem = CreateCompatibleDC(NULL);

  // now select bitmap into the memory dc
  if (NULL == SelectObject(dcmem, hBmp))
  {	// failed to load bitmap into device context
    DeleteDC(dcmem);
    return false;
  }


  // now get the bmp size
  BITMAP bm;
  GetObject(hBmp, sizeof(bm), &bm);
  // and blit it to the visible dc
  if (BitBlt(hDC, 0, 0, bm.bmWidth, bm.bmHeight, dcmem,
    0, 0, SRCCOPY) == 0)
  {	// failed the blit
    DeleteDC(dcmem);
    return false;
  }

  DeleteDC(dcmem);  // clear up the memory dc

  return true;
}

/***************************************************************
bool SaveBMP ( BYTE* Buffer, int width, int height,
long paddedsize, LPCTSTR bmpfile )

Function takes a buffer of size <paddedsize>
and saves it as a <width> * <height> sized bitmap
under the supplied filename.
On error the return value is false.

***************************************************************/

bool SaveBMP(BYTE* Buffer, int width, int height, long paddedsize, LPCTSTR bmpfile)
{
  // declare bmp structures 
  BITMAPFILEHEADER bmfh;
  BITMAPINFOHEADER info;

  // andinitialize them to zero
  memset(&bmfh, 0, sizeof(BITMAPFILEHEADER));
  memset(&info, 0, sizeof(BITMAPINFOHEADER));

  // fill the fileheader with data
  bmfh.bfType = 0x4d42;       // 0x4d42 = 'BM'
  bmfh.bfReserved1 = 0;
  bmfh.bfReserved2 = 0;
  bmfh.bfSize = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER) + paddedsize;
  bmfh.bfOffBits = 0x36;		// number of bytes to start of bitmap bits

  // fill the infoheader

  info.biSize = sizeof(BITMAPINFOHEADER);
  info.biWidth = width;
  info.biHeight = height;
  info.biPlanes = 1;			// we only have one bitplane
  info.biBitCount = 24;		// RGB mode is 24 bits
  info.biCompression = BI_RGB;
  info.biSizeImage = 0;		// can be 0 for 24 bit images
  info.biXPelsPerMeter = 0x0ec4;     // paint and PSP use this values
  info.biYPelsPerMeter = 0x0ec4;
  info.biClrUsed = 0;			// we are in RGB mode and have no palette
  info.biClrImportant = 0;    // all colors are important

  // now we open the file to write to
  HANDLE file = CreateFile(bmpfile, GENERIC_WRITE, FILE_SHARE_READ,
    NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
  if (file == NULL)
  {
    CloseHandle(file);
    return false;
  }

  // write file header
  unsigned long bwritten;
  if (WriteFile(file, &bmfh, sizeof(BITMAPFILEHEADER), &bwritten, NULL) == false)
  {
    CloseHandle(file);
    return false;
  }
  // write infoheader
  if (WriteFile(file, &info, sizeof(BITMAPINFOHEADER), &bwritten, NULL) == false)
  {
    CloseHandle(file);
    return false;
  }
  // write image data
  if (WriteFile(file, Buffer, paddedsize, &bwritten, NULL) == false)
  {
    CloseHandle(file);
    return false;
  }

  // and close file
  CloseHandle(file);

  return true;
}

/*******************************************************************
BYTE* LoadBMP ( int* width, int* height, long* size
LPCTSTR bmpfile )

The function loads a 24 bit bitmap from bmpfile,
stores it's width and height in the supplied variables
and the whole size of the data (padded) in <size>
and returns a buffer of the image data

On error the return value is NULL.

NOTE: make sure you [] delete the returned array at end of
program!!!
*******************************************************************/

BYTE* LoadBMP(int* width, int* height, long* size, LPCTSTR bmpfile)
{
  // declare bitmap structures
  BITMAPFILEHEADER bmpheader;
  BITMAPINFOHEADER bmpinfo;
  // value to be used in ReadFile funcs
  DWORD bytesread;
  // open file to read from
  HANDLE file = CreateFile(bmpfile, GENERIC_READ, FILE_SHARE_READ,
    NULL, OPEN_EXISTING, FILE_FLAG_SEQUENTIAL_SCAN, NULL);
  if (NULL == file)
    return NULL; // coudn't open file


  // read file header
  if (ReadFile(file, &bmpheader, sizeof(BITMAPFILEHEADER), &bytesread, NULL) == false)
  {
    CloseHandle(file);
    return NULL;
  }

  //read bitmap info

  if (ReadFile(file, &bmpinfo, sizeof(BITMAPINFOHEADER), &bytesread, NULL) == false)
  {
    CloseHandle(file);
    return NULL;
  }

  // check if file is actually a bmp
  if (bmpheader.bfType != 'MB')
  {
    CloseHandle(file);
    return NULL;
  }

  // get image measurements
  *width = bmpinfo.biWidth;
  *height = abs(bmpinfo.biHeight);

  // check if bmp is uncompressed
  if (bmpinfo.biCompression != BI_RGB)
  {
    CloseHandle(file);
    return NULL;
  }

  // check if we have 24 bit bmp
  if (bmpinfo.biBitCount != 24)
  {
    CloseHandle(file);
    return NULL;
  }


  // create buffer to hold the data
  *size = bmpheader.bfSize - bmpheader.bfOffBits;
  BYTE* Buffer = new BYTE[*size];
  // move file pointer to start of bitmap data
  SetFilePointer(file, bmpheader.bfOffBits, NULL, FILE_BEGIN);
  // read bmp data
  if (ReadFile(file, Buffer, *size, &bytesread, NULL) == false)
  {
    delete[] Buffer;
    CloseHandle(file);
    return NULL;
  }

  // everything successful here: close file and return buffer

  CloseHandle(file);

  return Buffer;
}

void TestBMPCopy(LPCTSTR input, LPCTSTR output)
{
  int x, y;
  long s;
  BYTE* b = LoadBMP(&x, &y, &s, input);
  SaveBMP(b, x, y, s, output);
  delete[] b;
}
#include <vector>
void drawShape(const vector<Vector2d> &ps, const vector<double> &colRed, const vector<double> &colBlue)
{
  LPCTSTR file = L"c:/Users/Tom/Downloads/levyreflectbFlip.bmp";
  long s2;

  int width = 400;
  int height = 400;
  std::vector<BYTE> out(width*height*3);
  for (int x = 0; x < width; x++)
  {
    for (int y = 0; y < height; y++)
    {
      out[(x + (y*width)) * 3 + 0] = 127;
      out[(x + (y*width)) * 3 + 1] = 220;
      out[(x + (y*width)) * 3 + 2] = 127;
    }
  }
  // now actually draw it
  for (int i = 0; i < (int)ps.size()-1; i++)
  {
    Vector2d pos = (ps[i] + Vector2d(0.5, 0.5)) * (double)width;
    int x = (int)max(0, min((int)pos[0], width - 1));
    int y = (int)max(0, min((int)pos[1], height - 1));
    if (out[(x + (y*width)) * 3 + 1] == 220)
    {
      double blue = (colBlue[2 * i] + colBlue[2 * i + 1]) * 0.5;
      out[(x + (y*width)) * 3 + 0] = (BYTE)(255.0 * colRed[i]);
      out[(x + (y*width)) * 3 + 1] = (BYTE)(255.0 * sqrt(colRed[i] * blue));
      out[(x + (y*width)) * 3 + 2] = (BYTE)(255.0 * blue);
    }
  }

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}
void drawPixelsVsReflectivity(const vector<double> &colRed, const vector<double> &colBlue)
{
  LPCTSTR file = L"c:/Users/Tom/Downloads/pixelvsreflectivity.bmp";
  long s2;

  int width = colRed.size()-1;
  int height = width;
  std::vector<BYTE> out(width*height * 3);

  // now actually draw it
  for (int x = 0; x < width; x++)
  {
    double red = colRed[x];
    double blue = (colBlue[2 * x] + colBlue[2 * x + 1])*0.5;
    double green = sqrt(red*blue);
    for (int y = 0; y < height; y++)
    {
      double shade = (double)y / (double)(height - 1);
      out[(x + (y*width)) * 3 + 0] = (BYTE)(255.0 * (1.0-shade)*(1.0 - (1.0-red)*shade));
      out[(x + (y*width)) * 3 + 1] = (BYTE)(255.0 * (1.0 - shade)*(1.0 - (1.0 - green)*shade));
      out[(x + (y*width)) * 3 + 2] = (BYTE)(255.0 * (1.0 - shade)*(1.0 - (1.0 - blue)*shade));
    }
  }

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}
static double lift = 0.5;
static int type = 0; // 0: Levy, 1: Koch, 2: Dragon, 3: random

void addKochChild(vector<Vector2d> &ps, int order, const Vector2d &p02, const Vector2d &p12, bool flip = false)
{
  const Vector2d p0 = p02;
  const Vector2d p1 = p12;
  Vector2d dir(p1[1] - p0[1], p0[0] - p1[0]);
  if (flip)
    dir = -dir;
  if (type == 3)
    dir *= ((double)(rand() % 1000)) / 500.0 - 1.0;
  Vector2d mid = (p0 + p1)*0.5 + dir * lift;
  bool f1, f2;
  if (type == 0)
    f1 = f2 = flip;
  else if (type == 1)
    f1 = f2 = !flip;
  else if (type == 2)
  {
    f1 = true; f2 = false;
  }
  else
    f1 = f2 = true;
  if (order > 0)
    addKochChild(ps, order - 1, p0, mid, f1);
  ps.push_back(mid);
  if (order > 0)
    addKochChild(ps, order - 1, mid, p1, f2);
}
const double pi = 3.14159265;
double sqr(double x)
{
  return x*x;
}

double angdif(double a1, double a2)
{
  double a = a1 - a2;
  while (a < -pi)
    a += 2.0*pi;
  while (a > pi)
    a -= 2.0*pi;
  return a;
}

void colourReflected(const vector<Vector2d> &ps, int i, vector<double> &cols)
{
  bool useReflectionModel = true;
  // 1. go from i backwards and add on each colour, if it can see it...
  Vector2d r1 = ps[i] - ps[i - 1];
  r1.normalize();
  Vector2d r2 = ps[i+1] - ps[i];
  r2.normalize();
  Vector2d dir1(r1[1], -r1[0]);
  Vector2d dir2(r2[1], -r2[0]);
  Vector2d normal = dir1 + dir2;
  normal.normalize();
  Vector2d dir(r1 + r2);
  dir.normalize();
  double largestAngle;
  double smallestAngle;
  if (useReflectionModel)
  {
    double a1 = atan2((ps[i - 1] - ps[i]).dot(normal), (ps[i] - ps[i - 1]).dot(dir)) * 2.0;
    double a2 = atan2((ps[i] - ps[i + 1]).dot(normal), (ps[i + 1] - ps[i]).dot(dir)) * 2.0; // reflect
    largestAngle = max(a1, a2);
    smallestAngle = min(a1, a2);
    if (largestAngle == smallestAngle)
    {
      largestAngle += 0.05;
      smallestAngle -= 0.05;
    }
  }
  else
  {
    smallestAngle = atan2((ps[i - 1] - ps[i]).dot(dir), (ps[i - 1] - ps[i]).dot(normal));
    largestAngle = atan2((ps[i + 1] - ps[i]).dot(dir), (ps[i + 1] - ps[i]).dot(normal));
  }
  double l1 = largestAngle;
  double s1 = smallestAngle;

  double angle = smallestAngle;
  bool onSmallSide = true;
  for (int j = i - 2; j >= 0; j--)
  {
    double X = (ps[j] - ps[i]).dot(dir);
    double Y = (ps[j] - ps[i]).dot(normal);
    double newAng = atan2(X, Y);
    if (newAng > 1.57 && angle < -1.57)
      onSmallSide = false;
    if (newAng < -1.57 && angle > 1.57)
      onSmallSide = true;
    angle = newAng;
    if (onSmallSide)
      smallestAngle = max(smallestAngle, angle);
    else
      largestAngle = min(largestAngle, angle);
  }

  // 2. go from i forwards and find largest angle
  angle = largestAngle;
  onSmallSide = false;
  for (int j = i + 2; j < (int)ps.size(); j++)
  {
    double X = (ps[j] - ps[i]).dot(dir);
    double Y = (ps[j] - ps[i]).dot(normal);
    double newAng = atan2(X, Y);
    if (newAng > 1.57 && angle < -1.57)
      onSmallSide = false;
    if (newAng < -1.57 && angle > 1.57)
      onSmallSide = true;
    angle = newAng;
    if (onSmallSide)
      smallestAngle = max(smallestAngle, angle);
    else
      largestAngle = min(largestAngle, angle);
  }
  if (useReflectionModel)
    cols[i] = max(0.0, largestAngle - smallestAngle) / (l1 - s1);
  else
    cols[i] = max(0.0, largestAngle - smallestAngle) / pi;
}
const bool FLIP = false; // false for Flip (confusingly)
void drawPixelsVsDimension()
{
  LPCTSTR file = L"c:/Users/Tom/Downloads/levyDimensionReflectbFlip.bmp";
  long s2;
  int depth = 9; // 10 for Bl
  int width = 2 << depth;
  int height = 400;
  std::vector<BYTE> out(width*height * 3);
  //  for (dimension = 0.0; dimension < 1.9; dimension += 0.05)
  cout << "extra red %: " << endl;
  for (int d = 1; d < height / 10; d++)
  {
    double dimension = 1.0 + (double)d / 80.0;
    lift = sqrt(sqr(pow(2.0, -1.0 / dimension)) - sqr(0.5));
    double avg[2];
    vector<double> seenCols[2];
    vector<Vector2d> pss[2];
    int order = 1 + depth;
    pss[1].push_back(Vector2d(-0.5, 0));
    Vector2d p2 = Vector2d(0.5, 0);
    srand(8);
    addKochChild(pss[1], order, pss[1][0], p2, FLIP);
    pss[1].push_back(p2);
    for (int j = 0; j < pss[1].size(); j += 2)
      pss[0].push_back(pss[1][j]);
    for (int k = 0; k < 2; k++)
    {
      int order = k + depth;
      vector<Vector2d> ps = pss[k];

      // 2. colour curve based on proportion of sky visible and reflectance (radiosity) off rest of curve
      vector<double> cols(ps.size());
      for (int i = 0; i < (int)cols.size(); i++)
        cols[i] = 0;
      for (int i = 1; i < (int)cols.size() - 1; i++)
        colourReflected(ps, i, cols);
      cols[0] = cols[1];
      cols.back() = cols[cols.size() - 2];


      vector<double> vals(ps.size() + 1);
      vector<double> zs(ps.size() + 1);
      for (int i = 0; i < (int)vals.size(); i++)
        zs[i] = 1e10;
      for (int i = 1; i < (int)ps.size() - 1; i++)
      {
        int index1 = (int)(ps.size() * (ps[i][0] - ps[0][0]) / (ps.back()[0] - ps[0][0]));
        int index2 = (int)(ps.size() * (ps[i + 1][0] - ps[0][0]) / (ps.back()[0] - ps[0][0]));
        index1 = max(0, min(index1, (int)ps.size()));
        index2 = max(0, min(index2, (int)ps.size()));
        for (int j = index1; j < index2; j++)
        {
          double dist = (ps[i][1] + ps[i + 1][1])*0.5;
          if (dist < zs[j])
          {
            zs[j] = dist;
            vals[j] = (cols[i] + cols[i + 1])*0.5;
          }
        }
      }
      seenCols[k] = vals;

      avg[k] = 0.0;
      double num = 0.0;
      for (int i = 1; i < (int)vals.size(); i++)
        avg[k] += vals[i];
      avg[k] /= (double)(vals.size() - 1);
    }
    cout << 100.0 * (-1.0 + avg[0] / avg[1]) << ", ";
    for (int x = 0; x < width; x++)
    {
      double red = seenCols[0][x];
      double blue = seenCols[1][max(0, 2*x - 1)]*0.25 + seenCols[1][2 * x]*0.5 + seenCols[1][2 * x + 1]*0.25;
      double green = sqrt(red*blue);
      for (int y = d * 10; y < (d + 1) * 10; y++)
      {
        out[(x + y*width) * 3 + 0] = (BYTE)(255.0*red);
        out[(x + y*width) * 3 + 1] = (BYTE)(255.0*green);
        out[(x + y*width) * 3 + 2] = (BYTE)(255.0*blue);
      }
    }
  }

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}

int _tmain(int argc, _TCHAR* argv[])
{
  drawPixelsVsDimension();
  double dimension = 1.4;
  lift = sqrt(sqr(pow(2.0, -1.0 / dimension)) - sqr(0.5));
  double lastav = 1.0;
  double lastAvg = 1.0;
  double avg[6];
  double av[6];
  vector<Vector2d> p1;
  vector<double> colours[2];
  vector<double> seenCols[2];
  for (int k = 0; k < 6; k++)
  {
    int order = k + 5;
    // 1. generate Koch curve into array:
    vector<Vector2d> ps;
    ps.push_back(Vector2d(-0.5, 0));
    Vector2d p2 = Vector2d(0.5, 0);
    srand(10);
    addKochChild(ps, order, ps[0], p2, FLIP);
    ps.push_back(p2);
    if (k == 4)
      p1 = ps;

    // 2. colour curve based on proportion of sky visible and reflectance (radiosity) off rest of curve
    vector<double> cols(ps.size());
    for (int i = 0; i < (int)cols.size(); i++)
      cols[i] = 0;
    for (int i = 1; i < (int)cols.size() - 1; i++)
      colourReflected(ps, i, cols);
    cols[0] = cols[1];
    cols.back() = cols[cols.size() - 2];
    if (k >= 4)
      colours[k-4] = cols;

    avg[k] = 0.0;
    for (int i = 1; i < (int)cols.size() - 1; i++)
      avg[k] += cols[i];
    avg[k] /= (double)(cols.size() - 2);

    vector<double> vals(ps.size() + 1);
    vector<double> zs(ps.size() + 1);
    for (int i = 0; i < (int)vals.size(); i++)
      zs[i] = 1e10;
    for (int i = 1; i < (int)ps.size() - 1; i++)
    {
      int index1 = (int)(ps.size() * (ps[i][0] - ps[0][0]) / (ps.back()[0] - ps[0][0]));
      int index2 = (int)(ps.size() * (ps[i + 1][0] - ps[0][0]) / (ps.back()[0] - ps[0][0]));
      index1 = max(0, min(index1, (int)ps.size()));
      index2 = max(0, min(index2, (int)ps.size()));
      for (int j = index1; j < index2; j++)
      {
        double dist = (ps[i][1] + ps[i + 1][1])*0.5;
        if (dist < zs[j])
        {
          zs[j] = dist;
          vals[j] = (cols[i] + cols[i+1])*0.5;
        }
      }
    }
    if (k >= 4)
      seenCols[k - 4] = vals;
    double num = 0.0;
    av[k] = 0.0;
    for (int i = 0; i < (int)vals.size(); i++)
    {
      if (zs[i] < 1e10)
      {
        av[k] += vals[i];
        num++;
      }
    }
    av[k] /= num;
  }
  for (int k = 1; k < 6; k++)
  {
    cout << "gain: " << avg[k - 1] / avg[k] << ", " << av[k - 1] / av[k] << ", vals:  " << avg[k] << ", " << av[k] << endl;
  }
  cout << "p1 size: " << p1.size() << endl;
//  drawShape(p1, colours[0], colours[1]);
//  drawPixelsVsReflectivity(seenCols[0], seenCols[1]);


  // experiments:
  // 1. diffuse sky, exact Koch, reflectance 0.5
  //   b. overall colour (adjacent orders e.g. 10, 11) vs dimension  
  //   e. overall colour vs randomness for 1 dimension
  // 2. specular sun, exact Koch:
  //   a. overall colour vs sun angle for 2 dimensions
  //   b. colour range ordered for 1 sun angle, 2 dimensions
  //   c. 5th,95th percentile colour vs dimension for 1 sun angle
  //   d. overall colour vs dimension for 2 sun angles
  //   e. colour vs x position for 2 sun angles and 2 dimensions and 2 randomnesses
  //   f. 5th,95th percentile colour vs dimension for +ve Levy, -ve Levy and Dragon versions for 2 sun angles
  //   g. 5th,95th percentile colour vs dimension for 2 randomnesses for 2 sun angles
  //   h. use histogram of colours to approximate specularity vs dimension 
}
