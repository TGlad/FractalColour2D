#include "stdafx.h"
#include "bmp.h"
#include <sstream>
#include <fstream>
#include <complex>
typedef complex<double> Complex;

int _tmain(int argc, _TCHAR* argv[])
{
  const int res = 10;
  double radii[res];
  memset(radii, 0, res*sizeof(double));
  Vector2d a(0, 1);
  Vector2d b(1, 1);
  vector<Vector2d> farey;
  farey.push_back(a);
  farey.push_back(b);
  radii[1] += 1.0/2.0;
  for (int i = 0; i < 10; i++)
  {
    for (int j = farey.size() - 2; j >= 0; j--)
    {
      Vector2d x = farey[j] + farey[j + 1];
      if (x[1] < res)
        radii[x[1]] += 1.0 / (2.0*x[1] * x[1]);
      farey.insert(farey.begin() + j+1, x);
    }
    for (auto x : farey)
      cout << x[0] << "/" << x[1] << " ";
    cout << endl;
  }
  double total = 0.0;
  for (int i = 0; i < res; i++)
  {

  }
}
