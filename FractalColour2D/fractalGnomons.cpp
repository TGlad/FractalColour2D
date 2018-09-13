// experiments with structural coloration of fractal curves, in particular Koch, Levy, Dragon and random curves.
#include "stdafx.h"
#include "bmp.h"
#include "spectrumToRGB.h"
#include <set>
#include<conio.h>

static int width = 1600;// 800;
static int height = 1200;// 600;

void putpixel(vector<BYTE> &out, const Vector2d &point, double green, double blue)
{
  Vector2i pos;
  pos[0] = int((point[0] + 1.0) * (double)height*0.8);
  pos[1] = int((point[1] + 0.7) * (double)height*0.8);

  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= height)
    return;
  int ind = 3 * (pos[0] + width*pos[1]);
  out[ind + 0] = 255;
  out[ind + 1] = (int)(green * 255.0);
  out[ind + 2] = (int)(blue * 255.0);
}
Vector3i getPixel(vector<BYTE> &out, const Vector2d &point)
{
  Vector2i pos;
  pos[0] = int((point[0] + 1.0) * (double)height*0.8);
  pos[1] = int((point[1] + 0.7) * (double)height*0.8);

  if (pos[0] < 0 || pos[0] >= width || pos[1] < 0 || pos[1] >= height)
    return Vector3i(0,0,0);
  int ind = 3 * (pos[0] + width*pos[1]);
  return Vector3i(out[ind + 0], out[ind + 1], out[ind + 2]);
}
double random(double start, double end)
{
  return start + (end - start)*(double)(rand() % 1000) / 1000.0;
}

static int fatness = 8;
static int fatwidth = width / fatness;
static int fatheight = height / fatness;
vector<BYTE> fat[2];
void putfatpixel(const Vector2d &point, double green, double blue)
{
  for (int i = 0; i < 2; i++)
  {
    double d = i == 0 ? -0.25 : 0.25;
    Vector2i pos;
    pos[0] = int(d + (point[0] + 1.0) * (double)fatheight / 2.0);
    pos[1] = int(d + (point[1] + 1.0) * (double)fatheight / 2.0);

    if (pos[0] < 0 || pos[0] >= fatwidth || pos[1] < 0 || pos[1] >= fatheight)
      return;

    int ind = 3 * (pos[0] + fatwidth*pos[1]);

    fat[i][ind + 0] = 255;
    fat[i][ind + 1] = (int)(green * 255.0);
    fat[i][ind + 2] = (int)(blue * 255.0);
  }
}
Vector3i getFatPixel(int i, const Vector2d &point)
{
  double d = i == 0 ? -0.25 : 0.25;
  Vector2i pos;
  pos[0] = int(d + (point[0] + 1.0) * (double)fatheight / 2.0);
  pos[1] = int(d + (point[1] + 1.0) * (double)fatheight / 2.0);

  if (pos[0] < 0 || pos[0] >= fatwidth || pos[1] < 0 || pos[1] >= fatheight)
    return Vector3i(0, 0, 0);
  int ind = 3 * (pos[0] + fatwidth*pos[1]);
  return Vector3i(fat[i][ind + 0], fat[i][ind + 1], fat[i][ind + 2]);
}

int _tmain(int argc, _TCHAR* argv[])
{
  long s2;
  vector<BYTE> out(width*height * 3); // .bmp pixel buffer
  for (int i = 0; i < 2; i++)
    fat[i].resize(fatwidth*fatheight * 3);
  int seed = 5;
  srand(seed); 
  // maxang1.5
  // srand2, ii=111 is good one!
  // srand5, ii=196 filled out
  // srand10 ii=47, maxang3.14
  // srand12 ii=199, "  closest to connected
  // srand13 ii=29 - beautiful wreathes
  // srand15 ii=11 - very nice leaves
  int stopii = 1521728;
  int iis[] = { 160, 369, 659, 875, 1170, 1257, 2919, 2936, 4322, 5690, 8466, 9717, 9720, 10112, 11144, 11691, 12378, 14081, 15198, 15699, 15723, 15974, 16295, 17962, 19673, 20346, 21195, 21304, 21849, 22437, 22991, 25349, 26529, 28668, 29510, 29632, 31153, 31865, 32262, 32421, 38214, 38389, 38553, 38583 };
  int iis2[] = { 1148124, 1521728, 1930700, 2129150, 2536776, 2930738, 3372686, 3472683 };
  int iic = 0;
/*  for (int ii = 0; ii < stopii; ii++)
    for (int i = 0; i < 10; i++)
      random(1, 2);*/
  double maxang = 3.14;
  Vector2d pos_[2];
  Vector2d posi_[2];
/*  double angle2_ = random(-maxang, maxang);
  pos_[1] = Vector2d(random(-1.0, 1.0), random(-1.0, 1.0));
  double anglei1_ = random(-maxang, maxang);
  double scalei_ = random(0.6, 0.8);
  posi_[0] = Vector2d(random(-0.5, 0.5), random(-0.5, 0.5));
  double anglei2_ = random(-maxang, maxang);
  posi_[1] = Vector2d(random(-0.5, 0.5), random(-0.5, 0.5));
  */
/*  double angle1_ = 0.05625;
  double angle2_ = -1.754;
  pos_[0] = Vector2d(0.5, 0.0);
  pos_[1] = Vector2d(-0.617, 0.5914);
  double anglei1_ = 1.9755;// random(-maxang, maxang);
  double scalei_ = sqrt(0.5);
  posi_[0] = Vector2d(-0.2467, -0.36195);// random(-0.5, 0.5), random(-0.5, 0.5));
  double anglei2_ = 0.63635;// random(-maxang, maxang);
  posi_[1] = Vector2d(-0.1895, -0.0965);// random(-0.5, 0.5), random(-0.5, 0.5));
  */
  double angle1_ = 0.0825;
  double angle2_ = -1.832775;
  pos_[0] = Vector2d(0.4493, 0.0);
  pos_[1] = Vector2d(-0.5925, 0.53515);
  double anglei1_ = 2.00338;// random(-maxang, maxang);
  double scalei_ = 0.6933;// sqrt(0.5);
  posi_[0] = Vector2d(-0.2417, -0.3413);// random(-0.5, 0.5), random(-0.5, 0.5));
  double anglei2_ = 0.6476;// random(-maxang, maxang);
  posi_[1] = Vector2d(-0.1895, -0.0965);// random(-0.5, 0.5), random(-0.5, 0.5));

  /*
  double angle2_ = 0;// random(-maxang, maxang);
  pos_[1] = Vector2d(0,0);// random(-1.0, 1.0), random(-1.0, 1.0));
  double anglei1_ = 0;// random(-maxang, maxang);
  double scalei_ = 0;// random(0.6, 0.8);
  posi_[0] = Vector2d(0,0);// random(-0.5, 0.5), random(-0.5, 0.5));
  double anglei2_ = 0;// random(-maxang, maxang);
  posi_[1] = Vector2d(0,0);// random(-0.5, 0.5), random(-0.5, 0.5)); */

  double oldGap = 1e10;
  double speed = 0.005;
  for (int ii = 0; ii < 4000000; ii++) // 10000000
  {
    {
      char c = _getch();
      if (c == 'l')
        pos_[1][0] += speed;
      if (c == 'j')
        pos_[1][0] -= speed;
      if (c == 'i')
        pos_[1][1] += speed;
      if (c == 'k')
        pos_[1][1] -= speed;

      if (c == 'h')
        pos_[0][0] += speed;
      if (c == 'f')
        pos_[0][0] -= speed;
      if (c == 't')
        pos_[0][1] += speed;
      if (c == 'g')
        pos_[0][1] -= speed;
      if (c == 'd')
        posi_[0][0] += speed;
      if (c == 'a')
        posi_[0][0] -= speed;
      if (c == 'w')
        posi_[0][1] += speed;
      if (c == 's')
        posi_[0][1] -= speed;
      if (c == '=')
        scalei_ += speed;
      if (c == '-')
        scalei_ -= speed;
      if (c == '0')
        anglei2_ += speed;
      if (c == '9')
        anglei2_ -= speed;
      if (c == '8')
        anglei1_ += speed;
      if (c == '7')
        anglei1_ -= speed;
      if (c == '6')
        angle2_ += speed;
      if (c == '5')
        angle2_ -= speed;
      if (c == '4')
        angle1_ += speed;
      if (c == '3')
        angle1_ -= speed;
      if (c == 'x')
        speed *= 1.5;
      if (c == 'z')
        speed /= 1.5;
//      cout << "new input" << posi_[0].transpose() << posi_[1].transpose() << endl;
    }

    // now how do I search this 14D space?
    // I either need realtime or output loads of pictures, 
    // or optimise for some criteria about proximity, e.g maximum distance as a percent of overall area
    Matrix2d mat[2];  // mat + pos i for shape 1 equals shape i in orientation defined by rmat + rpos i (which can have any scale) 
    Vector2d pos[2];
    Matrix2d mati[2];
    Vector2d posi[2];
    double angle1 = angle1_;
    pos[0] = pos_[0];
    double len = 0.0;// ii == 0 ? 0 : 0.02;
    double angle2 = angle2_ + random(-maxang*len, maxang*len);
    pos[1] = pos_[1] + Vector2d(random(-len, len), random(-len, len));

    double anglei1 = anglei1_ + random(-maxang*len, maxang*len);
    double scalei = scalei_ + random(-0.5*len, 0.5*len);
    posi[0] = posi_[0] + Vector2d(random(-0.5*len, 0.5*len), random(-0.5*len, 0.5*len));

    double anglei2 = anglei2_ + random(-maxang*len, maxang*len);
    posi[1] = posi_[1] + Vector2d(random(-0.5*len, 0.5*len), random(-0.5*len, 0.5*len));

    cout << "index: " << ii << ", angle 1: " << angle1 << ", pos[0]: " << pos[0].transpose() << ", angle2: " << angle2 << ", pos[1]: " << pos[1].transpose() << ", anglei1: " << anglei1 << ", scalei: " << scalei << ", posi[0]: " << posi[0].transpose() << ", anglei2: " << anglei2 << ", posi[1]: " << posi[1].transpose() << endl;

    /*
    if (ii == iis2[iic])
    {
      iic++;
      if (iic == 8)
        break;
    }
      //  else
    continue;*/
 /*   if (ii == iis2[iic] && iic < 44)
    {
      iic++;
    }
//    else
      continue;
      */
    Vector2d v1(cos(angle1), sin(angle1));
    Vector2d v2(cos(angle2), sin(angle2));
    mat[0] << v1[0], v1[1],
      -v1[1], v1[0];
    mat[1] << v2[0], v2[1],
      -v2[1], v2[0];

    mati[0] << cos(anglei1), sin(anglei1),
      -sin(anglei1), cos(anglei1);
    mati[0] *= scalei;

    mati[1] << cos(anglei2), sin(anglei2),
      -sin(anglei2), cos(anglei2);
    mati[1] *= scalei;

    int collided = 0;
    int wideHit[2] = { 0,0 };
    memset(&out[0], 0, out.size() * sizeof(BYTE));
    for (int i = 0; i < 2; i++)
      memset(&fat[i][0], 0, fat[i].size() * sizeof(BYTE));
    vector<Vector2d> points[2], newPoints[2];
    newPoints[1].push_back(Vector2d(0, 0));
    int top = 28;
    for (int i = 0; i < top; i++) // 22
    {
      int p = i % 2;
      // place latest points in 0, 1 onto layout p to update shape p
      vector<Vector2d> newps;
      for (auto &point : newPoints[0])
        newps.push_back(mati[p] * point + posi[p]);
      for (auto &point : newPoints[1])
        newps.push_back(mati[p] * (mat[p] * point + pos[p]) + posi[p]);
      
      points[p].insert(points[p].end(), newPoints[p].begin(), newPoints[p].end());
      newPoints[p] = newps;
      if (p == 0)
      {
        for (auto &point : newPoints[p])
        {
          Vector3i col = getPixel(out, point);
          if (col == Vector3i(255,0,255) || col == Vector3i(255,255,0))
          {
            collided ++;
          }
          Vector3i fcol1 = getFatPixel(0, point);
          Vector3i fcol2 = getFatPixel(1, point);
          if (fcol1 == Vector3i(255, 255, 0) || fcol2 == Vector3i(255, 255, 0))
            wideHit[0] ++;
          if (fcol1 == Vector3i(255, 0, 255) || fcol2 == Vector3i(255, 0, 255))
            wideHit[1] ++;
          if (i >= top - 2)
            putpixel(out, point, 1.0, 1.0);
          putfatpixel(point, 1.0, 1.0);
        }
      }
      else
      {
        for (auto &point : newPoints[p])
        {
          Vector2d p[2] = { mat[0] * point + pos[0], mat[1] * point + pos[1] };
          if (getPixel(out, p[0]) == Vector3i(255, 255, 255) || getPixel(out, p[1]) == Vector3i(255, 255, 255))
          {
            collided ++;
          }
          if (i >= top - 1)
          {
            putpixel(out, p[0], 1, 0);
            putpixel(out, p[1], 0, 1);
          }
          for (int j = 0; j < 2; j++)
          {
            Vector3i fcol1 = getFatPixel(0, p[j]);
            Vector3i fcol2 = getFatPixel(1, p[j]);
            if (fcol1 == Vector3i(255, 255, 255) || fcol2 == Vector3i(255, 255, 255))
              wideHit[j] ++;
          }
          putfatpixel(p[0], 1, 0);
          putfatpixel(p[1], 0, 1);
        }
      }
 //     if (collided>3 && ii != 0)
 //       break;
    }
 //   if (ii != 0 && collided>3 || wideHit[0]<=10 || wideHit[1]<=10)
 //     continue;

/*    Vector2d means[2];
    double difs[2];
    for (int x = 0; x < 2; x++)
    {
      means[x] = Vector2d(0, 0);
      for (auto &p : points[x])
        means[x] += p;
      means[x] /= (double)points[x].size();
      difs[x] = 0;
      for (auto &p : points[x])
        difs[x] += (p-means[x]).squaredNorm();
      difs[x] /= (double)points[x].size();
      difs[x] = sqrt(difs[x]);
    }
    Vector2d p[2] = { mat[0] * means[1] + pos[0], mat[1] * means[1] + pos[1] };
    double gap = ((p[0] - means[0]).norm() + (p[1] - means[0]).norm()) / (difs[0] + difs[1]);
    cout << "num: " << ((p[0] - means[0]).norm() + (p[1] - means[0]).norm()) << ", den: " << difs[0] << ", " << difs[1] << endl;
    if (gap >= oldGap)
      continue;
    oldGap = gap;
    angle2_ = angle2;
    pos_[1] = pos[1];
    anglei1_ = anglei1;
    scalei_ = scalei;
    posi_[0] = posi[0];
    anglei2_ = anglei2;
    posi_[1] = posi[1]; */

    BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    wstringstream str;
    str << L"gnomon" << seed << "/shape" << stopii << L".bmp";
    wstring strng = str.str();
    const TCHAR * bll = strng.c_str();
    LPCTSTR file = bll;
    SaveBMP(c, width, height, s2, file);
    delete[] c;
    /*
    for (int i = 0; i < 2; i++)
    {
      BYTE* c = ConvertRGBToBMPBuffer(&fat[i][0], fatwidth, fatheight, &s2);
      wstringstream str;
      str << L"gnomon/fattest" << i << "_" << ii << L".bmp";
      wstring strng = str.str();
      const TCHAR * bll = strng.c_str();
      LPCTSTR file = bll;
      SaveBMP(c, fatwidth, fatheight, s2, file);
      delete[] c;
    }*/
  }
}

