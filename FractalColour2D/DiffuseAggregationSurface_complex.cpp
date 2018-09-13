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

double angdif(double a1, double a2)
{
  double a = a1 - a2;
  while (a < -pi)
    a += 2.0*pi;
  while (a > pi)
    a -= 2.0*pi;
  return a;
}

void colourReflected(const vector<Vector2d> &ps, int i, Vector2d &angleRange, double &viewAngle)
{
  // find and return incident angle range and view angle:

  Vector2d r1 = ps[i] - ps[i - 1];
  r1.normalize();
  Vector2d r2 = ps[i + 1] - ps[i];
  r2.normalize();
  Vector2d dir1(r1[1], -r1[0]);
  Vector2d dir2(r2[1], -r2[0]);
  Vector2d normal = dir1 + dir2;
  normal.normalize();
  Vector2d dir(r1 + r2);
  dir.normalize();
  double largestAngle;
  double smallestAngle;
  viewAngle = atan2(normal[0], normal[1]);

  double a1 = atan2(-dir[0], -dir[1]);
  double a2 = atan2( dir[0],  dir[1]);
  largestAngle = max(a1, a2);
  smallestAngle = min(a1, a2);
  if (largestAngle == smallestAngle)
  {
    largestAngle += 0.05;
    smallestAngle -= 0.05;
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
  angleRange = Vector2d(smallestAngle, largestAngle);
}


int _tmain(int argc, _TCHAR* argv[])
{
  // TODO: 
  // do once for each dimension:
  // calculate curve, for each row, calculate colour matrix for incidence vs reflection angle, quadratic weight function
  // then for big curve, for each vertex, calculate incidence angle range and reflection angle, integrate above matrix to produce colour.
  int order = 11;
  bool flip = false; // true = concave
  const int numAngles = 200;
  int width = numAngles;
  int numDimensions = 5;
  int numFrequencies = 60;
  int gap = 3;
  int height = numDimensions * (numFrequencies + gap);
  long s2;
  vector<BYTE> out(width*height * 3);
  LPCTSTR file = L"c:/Users/Tom/Downloads/realColours.bmp";
  for (int d = 0; d < numDimensions; d++)
  {
    double dimension = 1.0 + sqr((double)(d) / (double)(numDimensions + 1));
    lift = sqrt(sqr(pow(2.0, -1.0 / dimension)) - sqr(0.5));

    // generate big curve
    double bigCurveAngle = 0.0;
    vector<Vector2d> pbs;
    pbs.push_back(0.5*Vector2d(-cos(bigCurveAngle), sin(bigCurveAngle)));
    Vector2d p2 = -pbs.back();
    srand(20);
    addKochChild(pbs, order, pbs[0], p2, flip);
    pbs.push_back(p2);

    // now we want to go through each vertex in turn, and calculate both the view angle and the incidence angle range
    vector<Vector2d> angleRanges(pbs.size());
    vector<double> viewAngles(pbs.size());
    for (int i = 0; i < pbs.size(); i++)
    {
      colourReflected(pbs, i, angleRanges[i], viewAngles[i]);
    }
    // now calculate just the seen values:
    vector<double> indices(pbs.size() + 1);
    vector<double> zs(pbs.size() + 1);
    for (int i = 0; i < (int)indices.size(); i++)
      zs[i] = 1e10;
    for (int i = 1; i < (int)pbs.size() - 1; i++)
    {
      int index1 = (int)(pbs.size() * (pbs[i][0] - pbs[0][0]) / (pbs.back()[0] - pbs[0][0]));
      int index2 = (int)(pbs.size() * (pbs[i + 1][0] - pbs[0][0]) / (pbs.back()[0] - pbs[0][0]));
      index1 = max(0, min(index1, (int)pbs.size()));
      index2 = max(0, min(index2, (int)pbs.size()));
      for (int j = index1; j < index2; j++)
      {
        double dist = (pbs[i][1] + pbs[i + 1][1])*0.5;
        if (dist < zs[j])
        {
          zs[j] = dist;
          indices[j] = i;
        }
      }
    }

    for (int s = 0; s < numFrequencies; s++)
    {
      Vector3d matrix[numAngles][numAngles];
      for (int anglev = 0; anglev < numAngles; anglev++)
      {
        double viewAngle = 1.57 * (-1.0 + 2.0*(double)anglev / (double)(numAngles - 1));
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

        for (int anglei = 0; anglei < numAngles; anglei++)
        {
          double iAngle = 1.57 * (-1.0 + 2.0*(double)anglei / (double)(numAngles - 1));
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
              double weight = (double)i / (double)vLengths.size();
              weight *= 4.0*(1.0 - weight);

              Vector3d frequencies = (double)scale * 2.0*pi * Vector3d(1, 1.25, 1.5);

              red += weight * Vector2d(sin(totalLength*frequencies[0]), cos(totalLength*frequencies[0]));
              green += weight * Vector2d(sin(totalLength*frequencies[1]), cos(totalLength*frequencies[1]));
              blue += weight * Vector2d(sin(totalLength*frequencies[2]), cos(totalLength*frequencies[2]));
            }
          }
          double r = red.norm() / (double)vLengths.size();
          double g = green.norm() / (double)vLengths.size();
          double b = blue.norm() / (double)vLengths.size();
          matrix[anglev][anglei] = Vector3d(r, g, b);
        }
        ps.clear();
      }
      // now go through all the columns and draw the big curve using the known angle ranges... 
      for (int i = 0; i < pbs.size(); i++)
      {
        int index = indices[i];
        int startIndex = (int)((double)numAngles*0.5*(1.0 + angleRanges[index][0] / 1.57));
        int endIndex = (int)((double)numAngles*0.5*(1.0 + angleRanges[index][1] / 1.57));
        int viewAng = (int)((double)numAngles*0.5*(1.0 + viewAngles[index] / 1.57)); // TODO: is this right?
        Vector3d colour(0, 0, 0);
        for (int col = startIndex; col < endIndex; col++)
        {
          colour += matrix[viewAng][col];
        }
        int ind = ((d*(numFrequencies + gap) + s)*width + i) * 3;
        out[ind + 0] = (BYTE)(255.0*min(1.0, colour[0]));
        out[ind + 1] = (BYTE)(255.0*min(1.0, colour[1]));
        out[ind + 2] = (BYTE)(255.0*min(1.0, colour[2]));
      }
    }
  }
  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
  SaveBMP(c, width, height, s2, file);
  delete[] c;
}
