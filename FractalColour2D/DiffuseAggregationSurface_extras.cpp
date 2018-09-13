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

static double lift = 0.5;
static int type = 1; // 0: Levy, 1: Koch, 2: Dragon, 3: random
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

void drawCurve(const vector<Vector2d> &points, const vector<int> &hits)
{
  int width = 400;
  int height = width;
  long s2;
  LPCTSTR file = L"c:/Users/Tom/Downloads/curve.bmp";
  vector<BYTE> out(width*height * 3);
  for (int i = 0; i < (int)out.size(); i++)
  {
    out[i] = 127;
  }
  for (int i = 0; i < (int)points.size(); i++)
  {
    double x = (points[i][0] + 0.5) * (double)width;
    double y = (points[i][1] + 0.5) * (double)height;
    int X = max(0, min((int)x, width - 1));
    int Y = max(0, min((int)y, height - 1));
    BYTE col = hits[i]==1 ? 255 : 0;
    out[3 * (X + Y*width) + 0] = col;
    out[3 * (X + Y*width) + 1] = col;
    out[3 * (X + Y*width) + 2] = col;
  }
  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}

int _tmain_basic(int argc, _TCHAR* argv[])
{
  // procedure, axes: scale, dimension, incident angle, view angle
  // 100x100 incident vs view, 4x4 scale vs dimension
  int order = 11;
  bool flip = true; // true = concave
  int gap = 4;
  int maxDimensions = 8;
  int maxScales = 8;
  int maxAngles = 2*maxScales * 2;
  int width = (maxAngles + gap) * maxScales + 2;
  int height = (maxAngles+gap) * maxDimensions - gap;
  long s2;
  vector<BYTE> out(width*height * 3);
  double rv = 0.0;
  double bv = 0.0;
  vector<double> rr(maxScales);
  vector<double> bb(maxScales);
// #define MATRIX_TABLE
#if defined(MATRIX_TABLE)
  LPCTSTR file = L"c:/Users/Tom/Downloads/randomDiscrete.bmp";
  for (int d = 0; d < maxDimensions; d++)
  {
    double dimension = 1.0 + (0.2*(double)d) / ((double)maxDimensions);
    lift = sqrt(sqr(pow(2.0, -1.0 / dimension)) - sqr(0.5));
    for (int s = 0; s < maxScales; s++)
    {
      double rval = 0.0;
      double gval = 0.0;
      double bval = 0.0;
      int scale = (1 + s) * (maxAngles / (2 * maxScales));
      int K = width*(maxAngles + gap)*d + (maxAngles + gap)*s;
      for (int anglev = -scale; anglev < scale; anglev++) // view angles
      {
        double viewAngle = asin(abs((double)anglev) / (double)scale);
        if (anglev < 0)
          viewAngle = -viewAngle;
#else
  LPCTSTR file = L"c:/Users/Tom/Downloads/kochDimPower2flip.bmp";
  for (int d = 0; d < 200; d++)
  {
    double dimension = 1.0 + 0.9*(double)d / 200.0;
    lift = sqrt(sqr(pow(2.0, -1.0 / dimension)) - sqr(0.5));
    for (int s = 0; s<maxScales; s++)
    {
      rr[s] = 0;
      bb[s] = 0;
      double rval = 0.0;
      double gval = 0.0;
      double bval = 0.0;
      int angles = (1 + s) * (maxAngles / (2 * maxScales));
      if (angles % 1)
        cout << "woops" << endl;
      int scale = angles / 2;
      int K = width*d + (maxAngles + gap)*s;
      for (int anglev = -angles; anglev <= angles; anglev++) // view angles
      {
        int anglei = anglev;
        double viewAngle = asin(abs((double)anglev) / (double)angles);
#endif
        // generate curve
        vector<Vector2d> ps;
        ps.push_back(0.5*Vector2d(-cos(viewAngle), sin(viewAngle)));
        Vector2d p2 = -ps.back();
        srand(20);
        addKochChild(ps, order, ps[0], p2, flip);
        ps.push_back(p2);

        // generate depth map of curve
        int maxSteps = 400; // integrate over how many
        vector<double> vLengths(maxSteps);
        double maxLength = 100.0;
        for (int i = 0; i < (int)vLengths.size(); i++)
          vLengths[i] = maxLength;
        for (int i = 0; i < (int)ps.size()-1; i++)
        {
          double x1 = (double)(maxSteps-1) * (ps[i][0]   + 0.5);
          double x2 = (double)(maxSteps-1) * (ps[i+1][0] + 0.5);
          for (int X = max(0, (int)ceil(x1)); X < min(maxSteps-1, (int)ceil(x2)); X++)
          {
            double blend = ((double)X - x1) / (x2 - x1);
            if (blend > 1.0 || blend < 0.0)
              cout << "unusual blend " << blend << " at " << i << ", " << s << ", " << d << endl;
            double dist = ps[i][1] * (1.0 - blend) + ps[i + 1][1] * blend;
            if (dist < vLengths[X])
              vLengths[X] = dist;
          }
        }
#if !defined(MATRIX_TABLE)
        {
          double iAngle = viewAngle;
#else
        for (int anglei = -scale; anglei < scale; anglei++) // incidence angles
        {
          double iAngle = asin(abs((double)anglei) / (double)scale);
          if (anglei < 0)
            iAngle = -iAngle;
#endif
          double incidenceAngle = iAngle - viewAngle;
          Vector2d angleDir(-sin(incidenceAngle), cos(incidenceAngle));
          Vector2d angleNorm(cos(incidenceAngle), sin(incidenceAngle));
          // now work out total lengths for each view index
          Vector2d red(0, 0);
          Vector2d green(0, 0);
          Vector2d blue(0, 0);
          vector<double> totalLengths(vLengths.size());
          vector<Vector2d> points;
          vector<int> hits;
          for (int i = 0; i < (int)vLengths.size(); i++)
          {
            if (vLengths[i] >= maxLength)
              continue;
            Vector2d p1(-0.5 + (double)i / ((double)maxSteps - 1.0), vLengths[i]);
            // now check for occlusion
            bool blocked = false;
            if (incidenceAngle > 0.0)
            {
              for (int j = i + 1; j < (int)vLengths.size(); j++)
              {
                Vector2d p2(-0.5 + (double)j / ((double)maxSteps - 1.0), vLengths[j]);
                if ((p2 - p1).dot(angleNorm) < 0.0) // has been blocked
                {
                  if (vLengths[j] == 100.0)
                    cout << "blocked by background at: " << d << ", " << s << ", " << anglev << ", " << anglei << ", something is wrong" << endl;
                  blocked = true;
                  break;
                }
              }
            }
            else
            {
              for (int j = i - 1; j >= 0; j--)
              {
                Vector2d p2(-0.5 + (double)j / ((double)maxSteps - 1.0), vLengths[j]);
                if ((p2 - p1).dot(angleNorm) > 0.0) // has been blocked
                {
                  if (vLengths[j] == 100.0)
                    cout << "blocked by background at: " << d << ", " << s << ", " << anglev << ", " << anglei << ", something is wrong" << endl;
                  blocked = true;
                  break;
                }
              }
            }
            points.push_back(p1);
            hits.push_back(blocked ? 0 : 1);
            if (!blocked)
            {
              double incidenceLength = p1.dot(angleDir);
              double totalLength = incidenceLength + vLengths[i];

              Vector3d frequencies = (double)scale * 2.0*pi * Vector3d(1, 1.5, 2);// 3, 4, 5);

              red += Vector2d(sin(totalLength*frequencies[0]), cos(totalLength*frequencies[0]));
              green += Vector2d(sin(totalLength*frequencies[1]), cos(totalLength*frequencies[1]));
              blue += Vector2d(sin(totalLength*frequencies[2]), cos(totalLength*frequencies[2]));
            }
          }
//          if (s == maxScales-1 && d == maxDimensions-1 && anglev == -scale/2 && anglei==scale/2)//-scale / 4 && anglei == scale / 2)
//            drawCurve(points, hits);
          double r = red.norm()   / (double)vLengths.size();
          double g = green.norm() / (double)vLengths.size();
          double b = blue.norm()  / (double)vLengths.size();
          g = sqrt(r*b);
          rval += r;
          gval += g;
          bval += b;
#if defined(MATRIX_TABLE) 
          out[(K + anglei + maxAngles / 2 + (anglev + maxAngles / 2)*width) * 3 + 0] = (BYTE)(255.0*min(1.0, r));
          out[(K + anglei + maxAngles / 2 + (anglev + maxAngles / 2)*width) * 3 + 1] = (BYTE)(255.0*min(1.0, g));
          out[(K + anglei + maxAngles / 2 + (anglev + maxAngles / 2)*width) * 3 + 2] = (BYTE)(255.0*min(1.0, b));
#else
          out[(K + anglei + maxAngles / 2) * 3 + 0] = (BYTE)(255.0*min(1.0, r));
          out[(K + anglei + maxAngles / 2) * 3 + 1] = (BYTE)(255.0*min(1.0, g));
          out[(K + anglei + maxAngles / 2) * 3 + 2] = (BYTE)(255.0*min(1.0, b));
#endif
        }
        ps.clear();
      }
#if defined(MATRIX_TABLE) 
      cout << " r/b: " << rval / bval << endl;
#else 
      rv += rval;
      bv += bval;
      rr[s] += rval;
      bb[s] += bval;
      rval /= (double)angles;
      gval /= (double)angles;
      bval /= (double)angles;
      out[(K + angles + 4 + maxAngles / 2) * 3 + 0] = (BYTE)(255.0*min(1.0, rval));
      out[(K + angles + 4 + maxAngles / 2) * 3 + 1] = (BYTE)(255.0*min(1.0, gval));
      out[(K + angles + 4 + maxAngles / 2) * 3 + 2] = (BYTE)(255.0*min(1.0, bval));
      out[(K + angles + 3 + maxAngles / 2) * 3 + 0] = (BYTE)(255.0*min(1.0, rval));
      out[(K + angles + 3 + maxAngles / 2) * 3 + 1] = (BYTE)(255.0*min(1.0, gval));
      out[(K + angles + 3 + maxAngles / 2) * 3 + 2] = (BYTE)(255.0*min(1.0, bval));
#endif
    }
  }
  for (int s = 0; s < maxScales; s++)
    cout << "scale " << s << " red shift: " << rr[s] / bb[s] << endl;

  cout << "r/b " << rv / bv << endl;
  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  SaveBMP(c, width, height, s2, file);
  delete[] c;
  return 1;
}

int _tmain_interferencewdimension(int argc, _TCHAR* argv[])
{
  int order = 11;
  bool flip = false; // true = concave
  int width = 100;
  int height = 100;
  long s2;
  vector<BYTE> out(width*height * 3);
  LPCTSTR file = L"c:/Users/Tom/Downloads/colourInterferenceN2.bmp";
  for (int d = 0; d < height; d++)
  {
    double dimension = 1.0 + (double)d / (double)height;
    lift = sqrt(sqr(pow(2.0, -1.0 / dimension)) - sqr(0.5));
    // generate curve
    double viewAngle = 0.0;
    vector<Vector2d> ps;
    ps.push_back(0.5*Vector2d(-cos(viewAngle), sin(viewAngle)));
    Vector2d p2 = -ps.back();
    srand(20);
    addKochChild(ps, order, ps[0], p2, flip);
    ps.push_back(p2);
    for (int s = 0; s < width; s++)
    {
      // generate depth map of curve
      int maxSteps = 400; // integrate over how many
      vector<double> vLengths(maxSteps);
      double maxLength = 100.0;
      for (int i = 0; i < (int)vLengths.size(); i++)
        vLengths[i] = maxLength;
      for (int i = 0; i < (int)ps.size() - 1; i++)
      {
        double x1 = (double)(maxSteps - 1) * (ps[i][0] + 0.5);
        double x2 = (double)(maxSteps - 1) * (ps[i + 1][0] + 0.5);
        for (int X = max(0, (int)ceil(x1)); X < min(maxSteps - 1, (int)ceil(x2)); X++)
        {
          double blend = ((double)X - x1) / (x2 - x1);
          double dist = ps[i][1] * (1.0 - blend) + ps[i + 1][1] * blend;
          if (dist < vLengths[X])
            vLengths[X] = dist;
        }
      }
      double frequency = 2.0*pi * (double)s / ((double)width * lift);
      Vector2d intensity(0, 0);
      Vector2d green(0, 0);
      Vector2d blue(0, 0);
      for (int i = 0; i < (int)vLengths.size(); i++)
      {
        double totalLength = vLengths[i] * 2.0;
        double frequencyg = frequency * 1.25;
        double frequencyb = frequency * 1.5;
        intensity += Vector2d(sin(totalLength*frequency), cos(totalLength*frequency));
        green += Vector2d(sin(totalLength*frequencyg), cos(totalLength*frequencyg));
        blue += Vector2d(sin(totalLength*frequencyb), cos(totalLength*frequencyb));
      }
      double intr = intensity.norm() / (double)vLengths.size();
      double intg = green.norm() / (double)vLengths.size();
      double intb = blue.norm() / (double)vLengths.size();
      out[(s + d*width) * 3 + 0] = (BYTE)(255.0*min(1.0, intr));
      out[(s + d*width) * 3 + 1] = (BYTE)(255.0*min(1.0, intg));
      out[(s + d*width) * 3 + 2] = (BYTE)(255.0*min(1.0, intb));
    }
  }
  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  SaveBMP(c, width, height, s2, file);
  delete[] c;
  return 1;
}

int _tmain(int argc, _TCHAR* argv[])
{
  // this one plots reflected intensity against angle (x) for each frequency (y) for several dimensions (coarse y)
  // note angle is normalised to twice bend angle for dimension, and frequency is normalised to lift height
  int order = 11;
  bool flip = false; // true = concave
  int numAngles = 200;
  int width = numAngles;
  int numDimensions = 5;
  int numFrequencies = 60;
  int gap = 3;
  int height = numDimensions * (numFrequencies + gap);
  long s2;
  vector<BYTE> out(width*height * 3);
//#define UNNORMAL
#if defined(UNNORMAL)
  LPCTSTR file = L"c:/Users/Tom/Downloads/unnormalInterferenceMatrix4.bmp";
#else
  LPCTSTR file = L"c:/Users/Tom/Downloads/interferenceMatrix4.bmp";
#endif
  for (int d = 0; d < numDimensions; d++)
  {
    double dimension = 1.0 + sqr((double)(d+1) / (double)(numDimensions + 1));
    lift = sqrt(sqr(pow(2.0, -1.0 / dimension)) - sqr(0.5));
    double bendAngle = atan2(lift, 0.5);
    for (int ang = 0; ang < numAngles; ang++)
    {
#if defined(UNNORMAL)
      double viewAngle = 1.57 * (-1.0 + 2.0*(double)ang / (double)(numAngles - 1)); // -2 to 2 bend angles
#else
      double viewAngle = 2.0 * bendAngle * (-1.0 + 2.0*(double)ang / (double)(numAngles - 1)); // -2 to 2 bend angles
#endif
      // generate curve
      vector<Vector2d> ps;
      ps.push_back(0.5*Vector2d(-cos(viewAngle), sin(viewAngle)));
      Vector2d p2 = -ps.back();
      srand(20);
      addKochChild(ps, order, ps[0], p2, flip);
      ps.push_back(p2);
      bool firstS = false;
      for (int s = 0; s < numFrequencies; s++)
      {
        // generate depth map of curve
        int maxSteps = 400; // integrate over how many
        vector<double> vLengths(maxSteps);
        double maxLength = 100.0;
        for (int i = 0; i < (int)vLengths.size(); i++)
          vLengths[i] = maxLength;
        for (int i = 0; i < (int)ps.size() - 1; i++)
        {
          double x1 = (double)(maxSteps - 1) * (ps[i][0] + 0.5);
          double x2 = (double)(maxSteps - 1) * (ps[i + 1][0] + 0.5);
          for (int X = max(0, (int)ceil(x1)); X < min(maxSteps - 1, (int)ceil(x2)); X++)
          {
            double blend = ((double)X - x1) / (x2 - x1);
            double dist = ps[i][1] * (1.0 - blend) + ps[i + 1][1] * blend;
            if (dist < vLengths[X])
              vLengths[X] = dist;
          }
        }
        double freqScale = 20; // 5.0;
#if defined(UNNORMAL)
        double frequency = 2.0*pi * freqScale * (double)s / (double)(numFrequencies - 1);
#else
        double frequency = 2.0*pi * freqScale * (double)s / ((double)(numFrequencies - 1) * lift);
#endif
        double frequencyg = frequency * 1.25;
        double frequencyb = frequency* 1.5; // 2.0 for basic
        Vector2d intensity(0, 0);
        Vector2d green(0, 0);
        Vector2d blue(0, 0);
        for (int i = 0; i < (int)vLengths.size(); i++)
        {
          if (vLengths[i] == 100.0)
            continue;
          double totalLength = vLengths[i] * 2.0;
          intensity += Vector2d(sin(totalLength*frequency), cos(totalLength*frequency));
          green += Vector2d(sin(totalLength*frequencyg), cos(totalLength*frequencyg));
          blue += Vector2d(sin(totalLength*frequencyb), cos(totalLength*frequencyb));
        }
        double intr = intensity.norm() / (double)vLengths.size();
        double intg = green.norm() / (double)vLengths.size();
        double intb = blue.norm() / (double)vLengths.size();
    //    intg = sqrt(intr*intb);
        int index = ((d*(numFrequencies + gap) + s)*width + ang) * 3;
        out[index + 0] = (BYTE)(255.0*min(1.0, intr));
        out[index + 1] = (BYTE)(255.0*min(1.0, intr));
        out[index + 2] = (BYTE)(255.0*min(1.0, intr));
      }
    }
  }
  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}
