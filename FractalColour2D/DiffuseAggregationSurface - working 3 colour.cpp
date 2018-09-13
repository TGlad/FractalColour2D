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
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Vector4d);

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
#include <strstream>

/**
* Convert a wavelength in the visible light spectrum to a RGB color value that is suitable to be displayed on a
* monitor
*
* @param wavelength wavelength in nm
* @return RGB color encoded in int. each color is represented with 8 bits and has a layout of
* 00000000RRRRRRRRGGGGGGGGBBBBBBBB where MSB is at the leftmost
*/
Vector3d wavelengthToRGB(double wavelength)
{
  return srgbXYZ2RGB(cie1931WavelengthToXYZFit(wavelength));
}

/**
* Convert XYZ to RGB in the sRGB color space
* <p>
* The conversion matrix and color component transfer function is taken from http://www.color.org/srgb.pdf, which
* follows the International Electrotechnical Commission standard IEC 61966-2-1 "Multimedia systems and equipment -
* Colour measurement and management - Part 2-1: Colour management - Default RGB colour space - sRGB"
*
* @param xyz XYZ values in a double array in the order of X, Y, Z. each value in the range of [0.0, 1.0]
* @return RGB values in a double array, in the order of R, G, B. each value in the range of [0.0, 1.0]
*/
Vector3d srgbXYZ2RGB(const Vector3d &xyz) 
{
  double x = xyz[0];
  double y = xyz[1];
  double z = xyz[2];
  double rl = 3.2406255 * x + -1.537208  * y + -0.4986286 * z;
  double gl = -0.9689307 * x + 1.8757561 * y + 0.0415175 * z;
  double bl = 0.0557101 * x + -0.2040211 * y + 1.0569959 * z;

  return Vector3d(
    srgbXYZ2RGBPostprocess(rl),
      srgbXYZ2RGBPostprocess(gl),
      srgbXYZ2RGBPostprocess(bl)
  );
}

/**
* helper function for {@link #srgbXYZ2RGB(double[])}
*/
double srgbXYZ2RGBPostprocess(double c) 
{
  // clip if c is out of range
  c = c > 1 ? 1 : (c < 0 ? 0 : c);

  // apply the color component transfer function
  c = c <= 0.0031308 ? c * 12.92 : 1.055 * pow(c, 1. / 2.4) - 0.055;

  return c;
}

/**
* A multi-lobe, piecewise Gaussian fit of CIE 1931 XYZ Color Matching Functions by Wyman el al. from Nvidia. The
* code here is adopted from the Listing 1 of the paper authored by Wyman et al.
* <p>
* Reference: Chris Wyman, Peter-Pike Sloan, and Peter Shirley, Simple Analytic Approximations to the CIE XYZ Color
* Matching Functions, Journal of Computer Graphics Techniques (JCGT), vol. 2, no. 2, 1-11, 2013.
*
* @param wavelength wavelength in nm
* @return XYZ in a double array in the order of X, Y, Z. each value in the range of [0.0, 1.0]
*/
Vector3d cie1931WavelengthToXYZFit(double wavelength) 
{
  double wave = wavelength;

  double x;
  {
    double t1 = (wave - 442.0) * ((wave < 442.0) ? 0.0624 : 0.0374);
    double t2 = (wave - 599.8) * ((wave < 599.8) ? 0.0264 : 0.0323);
    double t3 = (wave - 501.1) * ((wave < 501.1) ? 0.0490 : 0.0382);

    x = 0.362 * exp(-0.5 * t1 * t1)
      + 1.056 * exp(-0.5 * t2 * t2)
      - 0.065 * exp(-0.5 * t3 * t3);
  }

  double y;
  {
    double t1 = (wave - 568.8) * ((wave < 568.8) ? 0.0213 : 0.0247);
    double t2 = (wave - 530.9) * ((wave < 530.9) ? 0.0613 : 0.0322);

    y = 0.821 * exp(-0.5 * t1 * t1)
      + 0.286 * exp(-0.5 * t2 * t2);
  }

  double z;
  {
    double t1 = (wave - 437.0) * ((wave < 437.0) ? 0.0845 : 0.0278);
    double t2 = (wave - 459.0) * ((wave < 459.0) ? 0.0385 : 0.0725);

    z = 1.217 * exp(-0.5 * t1 * t1)
      + 0.681 * exp(-0.5 * t2 * t2);
  }

  return Vector3d( x, y, z );
}
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

struct Bit
{
  Vector2d p, p2;
  int index;
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};
bool lessThanBit(const Bit &a, const Bit &b)
{
  return a.p[0] < b.p[0];
}
int _tmain(int argc, _TCHAR* argv[])
{
  // TODO: 
  // do once for each dimension:
  // calculate curve, for each row, calculate colour matrix for incidence vs reflection angle, quadratic weight function
  // then for big curve, for each vertex, calculate incidence angle range and reflection angle, integrate above matrix to produce colour.
  int order = 12;
  bool flip = false; // true = concave

#define DOKOCH
#if defined(DOKOCH)
  const int width = 330;
  int totalWidth = width;
  const int numFrequencies = 110;
  int numDimensions = 1;
#else
  const int width = 400;
  int totalWidth = width + 4;
  const int numFrequencies = 80;
  int numDimensions = 5;
#endif
  int gap = 3;
  long s2;
#if defined(DOKOCH)
  int height = numFrequencies;
  int maxRots = 60;
  for (int rot = 0; rot < 1/*maxRots*/; rot++)
  {
    double bigCurveAngle = 0.03 + 60.0 * (pi / 180.0) * (double)rot / (double)maxRots;
    vector<BYTE> out(totalWidth*height * 3);
    memset(&out[0], 255, out.size() * sizeof(BYTE));
    int d = 0;
    {
    lift = 0.5*tan(30.0 * pi / 180.0);
    Vector2d a(-0.5, -lift);
    Vector2d b(0.5, -lift);
    Vector2d c(0, 0.5 / sin(60.0 * pi / 180.0));
    double cosa = cos(bigCurveAngle);
    double sina = -sin(bigCurveAngle);
    vector<Vector2d> ps;
    ps.push_back(Vector2d(a[0] * cosa - a[1] * sina, a[0] * sina + a[1] * cosa));
    Vector2d p2(b[0] * cosa - b[1] * sina, b[0] * sina + b[1] * cosa);
    addKochChild(ps, order + 2, ps[0], p2, flip);
    ps.push_back(p2);
    Vector2d p3(c[0] * cosa - c[1] * sina, c[0] * sina + c[1] * cosa);
    addKochChild(ps, order + 2, ps.back(), p3, flip);
    ps.push_back(p3);
#else
  int height = numDimensions * (numFrequencies + gap);
  double bigCurveAngle = 0.0;
  vector<BYTE> out(totalWidth*height * 3);
  memset(&out[0], 0, out.size() * sizeof(BYTE));
  for (int d = 0; d < numDimensions; d++)
  {
    vector<Vector2d> ps;
    double dimension = 1.0 + sqr((double)(d) / (double)(numDimensions + 1));
    lift = sqrt(sqr(pow(2.0, -1.0 / dimension)) - sqr(0.5));
    ps.push_back(0.5*Vector2d(-cos(bigCurveAngle), sin(bigCurveAngle)));
    Vector2d p2 = -ps.back();
    srand(20);
    int add = max(0, 4 + d - numDimensions);
    addKochChild(ps, order + add, ps[0], p2, flip);
    ps.push_back(p2);
#endif
    // now find the visible indices
    vector<Bit, Eigen::aligned_allocator<Bit> > fwd(ps.size()); // stays in order
    vector<bool> visible(ps.size());
    for (int i = 0; i < ps.size(); i++)
    {
      int next = min(i + 1, (int)ps.size() - 1);
      Bit b;
      b.p = ps[i];
      b.p2 = ps[next];
      b.index = i;
      fwd[i] = b;
      visible[i] = true;
    }
    sort(fwd.begin(), fwd.end(), lessThanBit);
    int i = 0;
    for (int i = 0; i < fwd.size(); i++)
    {
      int j = i + 1;
      while (j < fwd.size() && fwd[j].p[0] < fwd[i].p2[0])
      {
        double blend = (fwd[j].p[0] - fwd[i].p[0]) / (fwd[i].p2[0] - fwd[i].p[0]);
        double dist = fwd[i].p2[1] * blend + fwd[i].p[1] * (1.0 - blend);
        if (fwd[j].p[1] > dist) // i to i+1 line is in front of this vertex
          visible[fwd[j].index] = false;
        j++;
      }
      j = i - 1;
      while (j >= 0 && fwd[j].p[0] > fwd[i].p2[0])
      {
        double blend = (fwd[j].p[0] - fwd[i].p[0]) / (fwd[i].p2[0] - fwd[i].p[0]);
        double dist = fwd[i].p2[1] * blend + fwd[i].p[1] * (1.0 - blend);
        if (fwd[j].p[1] > dist) // i to i+1 line is in front of this vertex
          visible[fwd[j].index] = false;
        j--;
      }
    }
    vector<int> pstops3(ps.size());
    vector<Vector2d> ps3;
    vector<double> origis;
    for (int i = 0; i < visible.size(); i++)
    {
      if (visible[i])
      {
        pstops3[i] = ps3.size();
        ps3.push_back(ps[i]);
        origis.push_back(i);
      }
    }
    int numAngles = 40;
    Vector3d pixels[numFrequencies][width];
    memset(pixels, 0, sizeof(Vector3d)*numFrequencies*width);
    for (int iAngle = 0; iAngle < numAngles; iAngle++)
    {
      double incidentAngle = 1.5 * (-1.0 + 2.0*(double)iAngle / (double)(numAngles - 1));
      cout << "angle: " << iAngle << endl;

      // next, filter to just the indices visible from light angle
      //    double incidentAngle = 0.0;
      Vector2d normal(sin(incidentAngle), cos(incidentAngle));
      Vector2d sideways(normal[1], -normal[0]);
      vector<bool> shadowed(ps3.size());
      for (int i = 0; i < shadowed.size(); i++)
        shadowed[i] = false;

      // now find the visible indices
      vector<Bit, Eigen::aligned_allocator<Bit> > fwd; // stays in order
      int shadowedBackFace = 0;
      for (int i = 0; i < (int)ps3.size(); i++)
      {
        int prev = max(0, i - 1);
        int next = min(i + 1, (int)ps3.size() - 1);
        double side0 = ps3[prev].dot(sideways);
        double side1 = ps3[i].dot(sideways);
        double side2 = ps3[next].dot(sideways);
        Bit nb;
        nb.p = Vector2d(side1, ps3[i].dot(normal));
        nb.p2 = Vector2d(side2, ps3[next].dot(normal));
        nb.index = i; 
        fwd.push_back(nb);
        if (side1 < side0 && side2 < side1)
        {
          shadowed[i] = true; // shadowed
          shadowedBackFace++;
        }
      }
      sort(fwd.begin(), fwd.end(), lessThanBit);
      int shadowedOccluded = 0;
      int shadowedOccludedBackFace = 0;
      for (int i = 0; i < fwd.size(); i++)
      {
        int j = i + 1;
        while (j < fwd.size() && fwd[j].p[0] < fwd[i].p2[0])
        {
          double blend = (fwd[j].p[0] - fwd[i].p[0]) / (fwd[i].p2[0] - fwd[i].p[0]);
          double dist = fwd[i].p2[1] * blend + fwd[i].p[1] * (1.0 - blend);
          if (fwd[j].p[1] > dist)
          {
            shadowed[fwd[j].index] = true;
            shadowedOccluded++;
          }
          j++;
        }
        j = i - 1; // for case of back face going backwards, we still want to occlude stuff behind it (for Levy curves)
        while (j >= 0 && fwd[j].p[0] > fwd[i].p2[0])
        {
          double blend = (fwd[j].p[0] - fwd[i].p[0]) / (fwd[i].p2[0] - fwd[i].p[0]);
          double dist = fwd[i].p2[1] * blend + fwd[i].p[1] * (1.0 - blend);
          if (fwd[j].p[1] > dist)
          {
            shadowed[fwd[j].index] = true;
            shadowedOccludedBackFace++;
          }
          j--;
        }
      }
      //    cout << "shadowed back face: " << shadowedBackFace << ", occluded: " << shadowedOccluded << ", occluded by back face: " << shadowedOccludedBackFace << endl;

      vector<Vector2d> vs2 = ps3;

      // now, for each index that is closest to a pixel centre, find the pixel radius around it, to get an average length.
      for (int s = 0; s < numFrequencies; s++)
      {
        // now, for each visible index, find ray lengths in +- one pixel radius of this...
        double avgcount = 0;
        int ind = 0;
        double rr = 0;
        double gg = 0;
        double bb = 0;
        // next I want curve size to grow but nothing else, 
        // this means coherence width reduces with the wavelength...
        // or curve actually grows and coherence width and wavelength stays the same
        int indexStart = (numFrequencies - 1) - s;
        int indexEnd = width - indexStart;
        int diff = indexEnd - indexStart;
        for (int i = 0; i < (int)ps3.size(); i++)
          vs2[i] = ps3[i] * (double)diff / (double)width;
        for (int index = indexStart; index < indexEnd; index++)
          //      for (int index = 0; index < width; index++)
        {
          double x = -0.5 + (double)index / (double)(width - 1);
          while (vs2[ind][0] < x && ind < (int)vs2.size() - 1)
            ind++;
          if (ind > 0 && abs(vs2[ind - 1][0] - x) < abs(vs2[ind][0] - x)) // pick closest
            ind--;
          Vector2d red(0, 0);
          Vector2d green(0, 0);
          Vector2d blue(0, 0);
          int i;
          Vector3d freqs(1, 1.25, 1.5);
          double offset = (freqs[2] / (double)width);
          // instead we want offset to grow with s, starting at half its normal value...
          //        double offset = (1.0 + 7.0*((double)s / (double)numFrequencies)) * 0.5*(freqs[2] / (double)width);
          int count = 0;
          for (i = ind; vs2[i][0] > x - offset && i > 0; i--);
          i++;
          Vector3d totalWeight(1e-10, 1e-10, 1e-10);

          //        double scale = 0.5*pi * (double)width * (1.0 + 10.0*(double)s / (double)numFrequencies);
          double scale = 0.5*pi * (double)width * 5.0;// 4.0;
          ////        double scale = 0.5*pi * (double)width * 2.0;
          for (; vs2[i][0] < x + offset && i < (int)vs2.size() - 1; i++)
          {
            double depth = vs2[i][1] + vs2[i].dot(normal);
            
            Vector3d freq = scale * freqs;
            Vector3d weight;
            weight[0] = max(0.0, 1.0 - sqr((vs2[i][0] - x) * freqs[0] / offset));
            weight[1] = max(0.0, 1.0 - sqr((vs2[i][0] - x) * freqs[1] / offset));
            weight[2] = max(0.0, 1.0 - sqr((vs2[i][0] - x) * freqs[2] / offset));
            count++;
            totalWeight += weight;
            if (!shadowed[i])
            {
              red += weight[0] * Vector2d(sin(depth*freq[0]), cos(depth*freq[0]));
              green += weight[1] * Vector2d(sin(depth*freq[1]), cos(depth*freq[1]));
              blue += weight[2] * Vector2d(sin(depth*freq[2]), cos(depth*freq[2]));
            }
          }
//          cout << "total weight: " << totalWeight[0];
          avgcount += count;
          double r = red.norm() / totalWeight[0];
          double g = green.norm() / totalWeight[1];
          double b = blue.norm() / totalWeight[2];
          if (numAngles >= 10)
          {
            r *= 5.0;
            g *= 5.0;
            b *= 5.0;
      /*      if (iAngle == numAngles / 6)
            {
              double sunStrength = 8.0;
              r *= sunStrength;
              g *= sunStrength;
              b *= sunStrength;
            }*/
          }
          rr += r;
          gg += g;
          bb += b;
          pixels[s][index] += Vector3d(r, g, b);
        }
#if !defined(DOKOCH)
        int inde = ((d*(numFrequencies + gap) + s)*totalWidth + width + 1) * 3;
        out[inde + 0] = out[inde + 3] += (BYTE)(255.0*min(rr / (double)(width * numAngles), 1.0));
        out[inde + 1] = out[inde + 4] += (BYTE)(255.0*min(gg / (double)(width * numAngles), 1.0));
        out[inde + 2] = out[inde + 5] += (BYTE)(255.0*min(bb / (double)(width * numAngles), 1.0));
#endif
//        cout << "average count: " << (double)avgcount / (double)width << endl;
      }
    }
    for (int s = 0; s < numFrequencies; s++)
    {
      for (int i = 0; i < width; i++)
      {
        pixels[s][i] /= (double)numAngles;
        double mult = 1.0;
        double m = max(pixels[s][i][0], max(pixels[s][i][1], pixels[s][i][2]));
        if (m > 1.0) // should never happen
          mult = 1.0 / m;

        int inde = ((d*(numFrequencies + gap) + s)*totalWidth + i) * 3;
        if (pixels[s][i] == Vector3d(0, 0, 0))
        {
          out[inde + 0] = (BYTE)255;
          out[inde + 1] = (BYTE)255;
          out[inde + 2] = (BYTE)255;
        }
        else
        {
          out[inde + 0] = (BYTE)(255.0*mult*pixels[s][i][0]);
          out[inde + 1] = (BYTE)(255.0*mult*pixels[s][i][1]);
          out[inde + 2] = (BYTE)(255.0*mult*pixels[s][i][2]);
        }
      }
    }
  }
#if defined(DOKOCH)
  BYTE* c = ConvertRGBToBMPBuffer(&out[0], totalWidth, height, &s2);
  wstringstream str;
  str << L"c:/Users/Tom/Downloads/kochs2/koch" << rot << L".bmp";
  wstring strng = str.str();
  const TCHAR * bll = strng.c_str();
  LPCTSTR file = bll;// str.str().c_str();
  SaveBMP(c, totalWidth, height, s2, file);
  delete[] c;
  }
#else
  BYTE* c = ConvertRGBToBMPBuffer(&out[0], totalWidth, height, &s2);
  LPCTSTR file = L"c:/Users/Tom/Downloads/dragonAmbient3.bmp";
  SaveBMP(c, totalWidth, height, s2, file);
  delete[] c;
#endif
}
