// experiments with structural coloration of fractal curves, in particular Koch, Levy, Dragon and random curves.
#include "stdafx.h"
#include "bmp.h"
#include "spectrumToRGB.h"

static double lift;
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
  int order = 12; // iteration depth of fractal curve
  bool flip = false; // false: curve bends towards camera, true: concave, curve bends away

#define DOKOCH // animation of spinning Koch snowflake
#if defined(DOKOCH)
  const int width = 330; // pixels along botto medge
  int totalWidth = width;
  const int numRows = 110; // this is the height  
  int numDimensions = 1;
#else
  const int width = 400;
  int totalWidth = width + 4;
  const int numRows = 80; // this is height per dimension 
  int numDimensions = 5; // number of stacked fractals of increasing dimension
#endif
  int gap = 3;
  long s2;
#if defined(DOKOCH)
  int height = numRows;
  int maxRots = 60;
  for (int rot = 0; rot < maxRots; rot++) // one per view angle in the gif
  {
    double bigCurveAngle = 0.03 + 60.0 * (pi / 180.0) * (double)rot / (double)maxRots; // 0 to 60 degrees is sufficient as the snowflake has hexagonal symmetry
    vector<BYTE> out(totalWidth*height * 3); // .bmp pixel buffer
    memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is white
    int d = 0;
    {
    lift = 0.5*tan(30.0 * pi / 180.0); // lift is the height of the mid point when the width is unit length. This is the correct lift for the Koch snowflake
    Vector2d a(-0.5, -lift);
    Vector2d b(0.5, -lift);
    Vector2d c(0, 0.5 / sin(60.0 * pi / 180.0));
    double cosa = cos(bigCurveAngle);
    double sina = -sin(bigCurveAngle);
    vector<Vector2d> ps;
    ps.push_back(Vector2d(a[0] * cosa - a[1] * sina, a[0] * sina + a[1] * cosa));
    Vector2d p2(b[0] * cosa - b[1] * sina, b[0] * sina + b[1] * cosa);
    addKochChild(ps, order + 2, ps[0], p2, flip); // recursive function creates the first edge of the snowflake
    ps.push_back(p2);
    Vector2d p3(c[0] * cosa - c[1] * sina, c[0] * sina + c[1] * cosa);
    addKochChild(ps, order + 2, ps.back(), p3, flip); // creates the second edge of the snowflake (we only need two for this animation)
    ps.push_back(p3);
#else
  int height = numDimensions * (numRows + gap);
  double bigCurveAngle = 0.0; // look at these curves head-on
  vector<BYTE> out(totalWidth*height * 3);
  memset(&out[0], 0, out.size() * sizeof(BYTE));
  for (int d = 0; d < numDimensions; d++) 
  {
    vector<Vector2d> ps;
    double dimension = 1.0 + sqr((double)(d) / (double)(numDimensions + 1)); // dimension above 1 increases in a squared rather than linear way, just more useful
    lift = sqrt(sqr(pow(2.0, -1.0 / dimension)) - sqr(0.5));  // lift is the height of the mid point when the width is unit length.
    ps.push_back(0.5*Vector2d(-cos(bigCurveAngle), sin(bigCurveAngle)));
    Vector2d p2 = -ps.back();
    srand(20); 
    int add = max(0, 4 + d - numDimensions); // this adds some more iterations to the higher dimension fractals as otherwise they would be too sparse
    addKochChild(ps, order + add, ps[0], p2, flip); // recursively generate curve
    ps.push_back(p2);
#endif
    // now find the visible indices of the curve
    vector<Bit, Eigen::aligned_allocator<Bit> > fwd(ps.size()); // stays in order
    vector<bool> visible(ps.size());
    for (int i = 0; i < (int)ps.size(); i++)
    {
      int next = min(i + 1, (int)ps.size() - 1);
      Bit b;
      b.p = ps[i];
      b.p2 = ps[next];
      b.index = i;
      fwd[i] = b;
      visible[i] = true;
    }
    sort(fwd.begin(), fwd.end(), lessThanBit); // so all points go 'forward' from left to right 
    int i = 0;
    for (int i = 0; i < (int)fwd.size(); i++) // basically, for each edge you cull out the vertices behind the edge
    {
      int j = i + 1;
      while (j < (int)fwd.size() && fwd[j].p[0] < fwd[i].p2[0]) // check forwards
      {
        double blend = (fwd[j].p[0] - fwd[i].p[0]) / (fwd[i].p2[0] - fwd[i].p[0]);
        double dist = fwd[i].p2[1] * blend + fwd[i].p[1] * (1.0 - blend);
        if (fwd[j].p[1] > dist) // i to i+1 line is in front of this vertex
          visible[fwd[j].index] = false;
        j++;
      }
      j = i - 1;
      while (j >= 0 && fwd[j].p[0] > fwd[i].p2[0]) // the edge could also go backwards
      {
        double blend = (fwd[j].p[0] - fwd[i].p[0]) / (fwd[i].p2[0] - fwd[i].p[0]);
        double dist = fwd[i].p2[1] * blend + fwd[i].p[1] * (1.0 - blend);
        if (fwd[j].p[1] > dist) // i to i+1 line is in front of this vertex
          visible[fwd[j].index] = false;
        j--;
      }
    }
    vector<Vector2d> ps3; // the set of visible points on the curve (not in left to right order)
    for (int i = 0; i < (int)visible.size(); i++)
      if (visible[i])
        ps3.push_back(ps[i]);
 
    int numAngles = 20; // integration over all incident light angles in 180 degree sky
    Vector3d pixels[numRows][width];
    memset(pixels, 0, sizeof(Vector3d)*numRows*width);
    for (int iAngle = 0; iAngle < numAngles; iAngle++)
    {
      double incidentAngle = 1.57 * (-1.0 + 2.0*(double)iAngle / (double)(numAngles - 1)); // -90 to 90 degrees
      cout << "angle: " << iAngle << endl;

      // next, filter to just the indices visible from light angle
      Vector2d normal(sin(incidentAngle), cos(incidentAngle));
      Vector2d sideways(normal[1], -normal[0]);
      vector<bool> shadowed(ps3.size());
      for (int i = 0; i < (int)shadowed.size(); i++)
        shadowed[i] = false;

      // now find the visible indices
      vector<Bit, Eigen::aligned_allocator<Bit> > fwd; // stays in order
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
        if (side1 < side0 && side2 < side1) // back faces..
          shadowed[i] = true; // ..are shadowed
      }
      sort(fwd.begin(), fwd.end(), lessThanBit);
      for (int i = 0; i < (int)fwd.size(); i++) // now iterate over each edge and cull any vertices behind it (just like above)
      {
        int j = i + 1;
        while (j < (int)fwd.size() && fwd[j].p[0] < fwd[i].p2[0]) // edge going left to right
        {
          double blend = (fwd[j].p[0] - fwd[i].p[0]) / (fwd[i].p2[0] - fwd[i].p[0]);
          double dist = fwd[i].p2[1] * blend + fwd[i].p[1] * (1.0 - blend);
          if (fwd[j].p[1] > dist)
            shadowed[fwd[j].index] = true;
          j++;
        }
        j = i - 1; // for case of back face going backwards, we still want to occlude stuff behind it (for Levy curves)
        while (j >= 0 && fwd[j].p[0] > fwd[i].p2[0]) // edge going right to left
        {
          double blend = (fwd[j].p[0] - fwd[i].p[0]) / (fwd[i].p2[0] - fwd[i].p[0]);
          double dist = fwd[i].p2[1] * blend + fwd[i].p[1] * (1.0 - blend);
          if (fwd[j].p[1] > dist)
            shadowed[fwd[j].index] = true;
          j--;
        }
      }
      vector<Vector2d> vs2 = ps3; // vs2 are the visible vertices

      for (int s = 0; s < numRows; s++) // vertical (each curve of increasing scale)
      {
        // now, for each visible index, find ray lengths in +- one pixel radius of this...
        double avgcount = 0;
        int ind = 0;
        double rr = 0;
        double gg = 0;
        double bb = 0;
        int indexStart = (numRows - 1) - s; // trapezoid shape
        int indexEnd = width - indexStart;
        int diff = indexEnd - indexStart;
        for (int i = 0; i < (int)ps3.size(); i++)
          vs2[i] = ps3[i] * (double)diff / (double)width; // scale the curve
        int numColours = 10; // colour spectrum resolution (number of frequency samples to use)
        vector<double> freqs(numColours); 
        for (int i = 0; i < numColours; i++)
          freqs[i] = 1.0 + 0.7*(double)i / (double)(numColours - 1); // linear ramp in frequency from 1 to 1.7 (times frequency of red light)
        for (int index = indexStart; index < indexEnd; index++) // for each pixel along the row
        {
          double x = -0.5 + (double)index / (double)(width - 1);
          while (vs2[ind][0] < x && ind < (int)vs2.size() - 1) // find vertex closest to pixel centre
            ind++;
          if (ind > 0 && abs(vs2[ind - 1][0] - x) < abs(vs2[ind][0] - x)) // pick closest
            ind--;
          vector<Vector2d> vals(numRows); // complex amplitude accumulators (I used 2d vector rather than actual complex number)
          for (int i = 0; i < numColours; i++)
            vals[i] = Vector2d(0, 0); 
          int i;
          double offset = (freqs.back() / (double)width); // max spot size radius
          // if we want to grow spot size with s then...
          //        double offset = (1.0 + 7.0*((double)s / (double)numRows)) * 0.5*(freqs[2] / (double)width);
          int count = 0;
          for (i = ind; vs2[i][0] > x - offset && i > 0; i--); // go to left-most vertex on spot
          i++;
          vector<double> totalWeight(numColours);
          for (int j = 0; j < numColours; j++)
            totalWeight[j] = 1e-10; // nonzero to avoid a division by 0 later

          //        double scale = 0.5*pi * (double)width * (1.0 + 10.0*(double)s / (double)numRows);
          double scale = 0.5*pi * (double)width * 5.0; // controls how many wavelengths you get in the spot size (this doesn't affect the colour much, unless too low).
          ////        double scale = 0.5*pi * (double)width * 2.0;
          for (; vs2[i][0] < x + offset && i < (int)vs2.size() - 1; i++) // for all vertices within the spot size
          {
            double depth = vs2[i][1] + vs2[i].dot(normal); // total path length of the ray
            
            for (int j = 0; j < numColours; j++)
            {
              double freq = scale * freqs[j];
              // just a parabolic weighting to diminish intensity on ends of spot. Note that this makes the spot width grow with light wavelength (though results are similar without this)
              double weight = max(0.0, 1.0 - sqr((vs2[i][0] - x) * freqs[j] / offset)); 
              totalWeight[j] += weight;
              if (!shadowed[i])
                vals[j] += weight * Vector2d(sin(depth*freq), cos(depth*freq)); // here's the complex number part
            }
            count++;
          }
          avgcount += count;
          Vector3d rgb(0, 0, 0);
          for (int j = 0; j < numColours; j++)
            rgb += wavelengthToRGB(690.0 / freqs[j]) *vals[j].norm() / totalWeight[j];
          rgb *= 3.0/(double)numColours; // make it brighter
          if (numAngles >= 10)
          {
            rgb *= 6.0; // and again
      /*      if (iAngle == numAngles / 6) // if we want an additional point light source in the sky
            {
              double sunStrength = 8.0;
              rgb *= sunStrength; 
            }*/
          }
          pixels[s][index] += rgb;
        }
#if !defined(DOKOCH) // this just shows average colour down the right hand side
        int inde = ((d*(numRows + gap) + s)*totalWidth + width + 1) * 3;
        out[inde + 0] = out[inde + 3] += (BYTE)(255.0*min(rr / (double)(width * numAngles), 1.0));
        out[inde + 1] = out[inde + 4] += (BYTE)(255.0*min(gg / (double)(width * numAngles), 1.0));
        out[inde + 2] = out[inde + 5] += (BYTE)(255.0*min(bb / (double)(width * numAngles), 1.0));
#endif
//        cout << "average count: " << (double)avgcount / (double)width << endl; // enable to see how many vertices are used in the interference for each pixel
      }
    }
    for (int s = 0; s < numRows; s++)
    {
      for (int i = 0; i < width; i++)
      {
        pixels[s][i] /= (double)numAngles;
        double mult = 1.0;
        double m = max(pixels[s][i][0], max(pixels[s][i][1], pixels[s][i][2]));
        if (m > 1.0) 
          mult = 1.0 / m; // saturate too bright reflections without losing the colour

        int inde = ((d*(numRows + gap) + s)*totalWidth + i) * 3;
        if (pixels[s][i] == Vector3d(0, 0, 0)) // background
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
  str << L"kochs/koch" << rot << L".bmp";
  wstring strng = str.str();
  const TCHAR * bll = strng.c_str();
  LPCTSTR file = bll;
  SaveBMP(c, totalWidth, height, s2, file);
  delete[] c;
  }
#else
  BYTE* c = ConvertRGBToBMPBuffer(&out[0], totalWidth, height, &s2);
  LPCTSTR file = L"koch.bmp";
  SaveBMP(c, totalWidth, height, s2, file);
  delete[] c;
#endif
}
