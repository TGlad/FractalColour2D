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
void TestBMPCopy2()
{
  LPCTSTR files[3] = { L"c:/Users/Tom/Downloads/red1cm8x8.bmp", L"c:/Users/Tom/Downloads/green1cm8x8-2.bmp", L"c:/Users/Tom/Downloads/blue1cm8x8.bmp" };
  LPCTSTR output = L"c:/Users/Tom/Downloads/col1cmdirectq.bmp";
  int width, height;
  long s, s2;

  int window = 1;
  std::vector<BYTE> out;
  for (int color = 0; color < 3; color++)
  {
//    window = 3 + color;
    BYTE* pix1 = LoadBMP(&width, &height, &s, files[color]);
    BYTE* pix = ConvertBMPToRGBBuffer(pix1, width, height);
    if (out.size() == 0)
      out.resize((width/window)*(height/window) * 3);
    for (int x = 0; x < width / window; x++)
    {
      for (int y = 0; y < height / window; y++)
      {
        double shade = 0;
        for (int i = 0; i < window; i++)
        {
          for (int j = 0; j < window; j++)
          {
            BYTE *col = pix + 3 * (window*x + i + width*(window*y + j));
//            double r = col[0];
//            double g = (double)(col[1] - 127) / 127.0;
//            double b = (double)(col[2] - 127) / 127.0;
//            shade += r*sqrt(b*b + g*g);
            double r = col[0];
            double g = (double)col[1];// / 255.0;
            shade += g;
          }
        }
        out[3 * (x + (width / window)*y) + color] = (BYTE)min(shade / double(window*window), 255.0);
      }
    }
    delete[] pix1;
    delete[] pix;
  }

  BYTE* c = ConvertRGBToBMPBuffer(&out[0], width/window, height/window, &s2);
  SaveBMP(c, width/window, height/window, s2, output);
  delete[] c;
}

static double lift = 0.5;

void addKochChild(vector<Vector2d> &ps, int order, const Vector2d &p0, const Vector2d &p1, bool flip = false)
{
  Vector2d dir(p1[1] - p0[1], p1[0] - p0[0]);
  if (flip)
    dir = -dir;
  Vector2d mid = (p0 + p1)*0.5 + dir * lift;
  if (order > 0)
    addKochChild(ps, order - 1, p0, mid, flip);
  ps.push_back(mid);
  if (order > 0)
    addKochChild(ps, order - 1, mid, p1, flip);
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

void colourReflectedOld(const vector<Vector2d> &ps, int i, vector<double> &cols)
{
  vector<Vector3d> angles(ps.size());
  // 1. go from i backwards and add on each colour, if it can see it...
  double smallestAngle = atan2(ps[i-1][0] - ps[i][0], ps[i-1][1] - ps[i][1]);
  cols[i] = 0.0;
  double reflectance = 1.0;
  int index = i - 1;
  angles[index][0] = smallestAngle;
  angles[index][1] = 0.0;
  int minIndex = index;
  double sa = smallestAngle;
  for (int j = i - 2; j >= 0; j--)
  {
    double angle = atan2(ps[j][0] - ps[i][0], ps[j][1] - ps[i][1]);
    sa = min(sa, angle);
    double distance = sqr(ps[i][0] - ps[j][0]) + sqr(ps[i][1] - ps[j][1]);
    if (angdif(angle, angles[index][0])<0.0)
    {
      while (index >= minIndex && angdif(angle, angles[index][0])<0.0)
        index--;
    }
    else
    {
      while (index < i - 1 && angdif(angle, angles[index][0])>0.0)
        index++;
    }
    if (index < minIndex || distance < angles[index][1])
    {
      angles[index][0] = angle;
      angles[index][1] = distance;
      angles[index][2] = cols[j];
    }
    minIndex = min(minIndex, index);
  }
  if (sa != angles[minIndex][0])
    int h = 0;
  for (int j = i - 2; j >= minIndex; j--)
  {
    double addition = angdif(angles[j + 1][0], angles[j][0]) / (2.0*pi);
    cols[i] += angles[j][2] * addition * reflectance;
  }

  // 2. go from i forwards and find largest angle
  double largestAngle = atan2(ps[i+1][0] - ps[i][0], ps[i+1][1] - ps[i][1]);
  index = i + 1;
  angles[index][0] = largestAngle;
  angles[index][1] = 0.0;
  int maxIndex = index;
  double la = largestAngle;
  for (int j = i + 2; j < (int)ps.size(); j++)
  {
    double angle = atan2(ps[j][0] - ps[i][0], ps[j][1] - ps[i][1]);
    la = max(la, angle);
    double distance = sqr(ps[i][0] - ps[j][0]) + sqr(ps[i][1] - ps[j][1]);
    
    if (angdif(angle, angles[index][0])>0.0)
    {
      while (index <= maxIndex && angdif(angle, angles[index][0])>0.0)
        index++;
    }
    else
    {
      while (index > i + 1 && angdif(angle, angles[index][0])<0.0)
        index--;
    }
    if (index > maxIndex || distance < angles[index][1])
    {
      angles[index][0] = angle;
      angles[index][1] = distance;
      angles[index][2] = cols[j];
    }
    maxIndex = max(maxIndex, index);
  }
  if (la != angles[maxIndex][0])
    int h = 3;
  for (int j = i + 2; j < maxIndex; j++)
  {
    double addition = angdif(angles[j][0], angles[j-1][0]) / (2.0*pi);
    cols[i] += angles[j][2] * addition * reflectance;
  }

  cols[i] += reflectance * angdif(angles[minIndex][0], angles[maxIndex][0]) / (2.0*pi);
}

void colourReflected(const vector<Vector2d> &ps, int i, vector<double> &cols)
{
  vector<Vector3d> angles(ps.size());
  // 1. go from i backwards and add on each colour, if it can see it...
  double smallestAngle = atan2(ps[i - 1][0] - ps[i][0], ps[i - 1][1] - ps[i][1]);
  double largestAngle = atan2(ps[i + 1][0] - ps[i][0], ps[i + 1][1] - ps[i][1]);
  double angleDif = largestAngle - smallestAngle;
  while (angleDif < 0.0)
    angleDif += 2.0*pi;
  while (angleDif > 2.0*pi)
    angleDif -= 2.0*pi;
  angleDif = 2.0*pi - angleDif;
  for (int j = i - 2; j >= 0; j--)
  {
    double angle = atan2(ps[j][0] - ps[i][0], ps[j][1] - ps[i][1]);
    double dif = angdif(angle, smallestAngle);
    if (dif < 0.0)
    {
      smallestAngle += dif;
      angleDif += dif;
    }
  }

  // 2. go from i forwards and find largest angle
  for (int j = i + 2; j < (int)ps.size(); j++)
  {
    double angle = atan2(ps[j][0] - ps[i][0], ps[j][1] - ps[i][1]);
    double dif = angdif(angle, largestAngle);
    if (dif > 0.0)
    {
      largestAngle += dif;
      angleDif -= dif;
    }
  }

  cols[i] = max(0.0, angleDif) / (2.0*pi);
}
int _tmain(int argc, _TCHAR* argv[])
{
//  TestBMPCopy2(); // combine images together

  double dimension = 1.2;
  lift = sqrt(sqr(pow(2.0, -1.0/dimension)) - sqr(0.5));
  // 1. generate Koch curve into array:
  vector<Vector2d> ps;
  ps.push_back(Vector2d(-0.5, 0));
  Vector2d p2 = Vector2d(0.5, 0);
  int order = 9;
  addKochChild(ps, order, ps[0], p2, true);
  ps.push_back(p2);

  // 2. colour curve based on proportion of sky visible and reflectance (radiosity) off rest of curve
  vector<double> cols(ps.size());
  for (int i = 0; i < (int)cols.size(); i++)
    cols[i] = 0;
  int numberOfReflections = 1;
  for (int r = 0; r < numberOfReflections; r++)
  {
    for (int i = 1; i < (int)cols.size() - 1; i++)
      colourReflected(ps, i, cols);
    cols[0] = cols[1];
    cols.back() = cols[cols.size() - 2];

    double avg = 0.0;
    for (int i = 1; i < cols.size()-1; i++)
      avg += cols[i];
    avg /= (double)(cols.size()-2);
    cout << "colour: " << avg << endl;
  }
  
  double col = 0.0;
  double num = 0.0;
  for (int i = 1; i < ps.size()-1; i++)
  {
    bool obscured = false;
    for (int j = 1; j < ps.size() - 2; j++)
    {
      if (ps[j][1] < ps[i][1])
      {
        if ((ps[j][0] < ps[i][0] && ps[j + 1][0] > ps[i][0]) || (ps[j+1][0] < ps[i][0] && ps[j][0] > ps[i][0]))
        {
          obscured = true;
          break;
        }
      }
    }
    if (!obscured)
    {
      col += cols[i];
      num++;
    }
  }
  double res = col / num;
  cout << "res: " << res;

  // 3. extract viewed colour per x position
  vector<Vector3d> data(ps.size());
  int index = 0;
  data[index][0] = ps[0][0];
  data[index][1] = ps[0][1];
  data[index][2] = cols[0];
  for (int i = 1; i < (int)ps.size(); i++)
  {
    double xPos = ps[i][0];
    double distance = ps[i][1];
    while (index < i && xPos > data[index][0])
      index++;
    if (index < i)
      while (xPos < data[index][0] && index > 0)
        index--;
    if (index == i || distance < data[index][1])
    {
      data[index][0] = xPos;
      data[index][1] = distance;
      data[index][2] = cols[i];
    }
  }
  vector<double> xCol(ps.size());
  int j = 0;
  double avgColour = 0.0;
  for (int i = 0; i < (int)ps.size(); i++)
  {
    avgColour += data[j][2];
    double xPos = ps[0][0] + (ps.back()[0] - ps[0][0])*(double)i / (double)(xCol.size() - 1);
    while (data[j][0] < xPos && j < (int)ps.size()) // average
      j++;
    double blend = (data[j][0] - xPos) / (data[j][0] - data[j - 1][0]);
    xCol[i] = cols[j - 1] * blend + cols[j] * (1.0 - blend); // linear interpolation
  }
  avgColour /= (double)ps.size();

  // experiments:
  // 1. diffuse sky, exact Koch, reflectance 0.5
  //   a. overall shade vs order
  //   b. overall colour (adjacent orders e.g. 10, 11) vs dimension  
  //   c. colour vs x position for 2 dimensions
  //   d. overall colour vs dimension for +ve Levy, -ve Levy and Dragon versions
  //   e. overall colour vs randomness for 1 dimension
  //   f. overall colour vs dimension for 2 randomnesses
  //   g. overall colour vs reflectance for 2 dimensions
  //   h. overall colour vs dimension for 2 reflectances
  //   i. 5th,95th percentile colour vs dimension for 2 reflectances
  //   j. find strongest 90th percentile colour (max - min) across all: dimension, reflectance, type, randomness
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
