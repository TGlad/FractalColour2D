#include "stdafx.h"
#include "bmp.h"
#include "spectrumToRGB.h"
// spacetime pixel
static const int width = 256;
static const int height = width;

double func(double x)
{
  return abs(x);
}
/*
static ofstream svg;
double scale = 900.0;
void PolyMesh::openSVG(const string &fileName, int number)
{
svg.open(fileName.c_str());

svg << "<svg width = \"" << (int)(scale * (double)number*1.05) << "\" height = \"" << (int)scale << "\" xmlns = \"http://www.w3.org/2000/svg\">" << endl;
}

void PolyMesh::saveSVG(const Vector3d &offset, double shade)
{
for (int i = 0; i < (int)faces.size(); i++)
{
//    cout << "i:" << i << ", head node: " << faces[i].head->nodeID << ", uv0: " << nodes[faces[i].head->nodeID].uv[0] << ", uv1: " << nodes[faces[i].head->nodeID].uv[1]  << endl;
svg << "<path d = \"M " << scale*(nodes[faces[i].head->nodeID].uv[0] + offset[0]) << " " << scale*(nodes[faces[i].head->nodeID].uv[1] + offset[1]);
Face::FaceNode *fn = faces[i].head->next;
do
{
svg << " L " << scale*(nodes[fn->nodeID].uv[0] + offset[0]) << " " << scale*(nodes[fn->nodeID].uv[1] + offset[1]);
fn = fn->next;
} while (fn != faces[i].head->next);
svg << "\" fill=\"transparent\" stroke=\"black\" style=\"stroke-opacity:" << shade << "\" />\n";
}
}

void PolyMesh::closeSVG()
{
svg << "</svg>" << endl;
svg.close();
}*/
int _tmain(int argc, _TCHAR* argv[])
{
  Vector2d NE(sqrt(0.5), sqrt(0.5));
  Vector2d SE(sqrt(0.5), -sqrt(0.5));
  long s2;
  int maxFrames = 40;
  for (int i = 0; i < maxFrames; i++)
  {
    double scale = pow(2.0 + sqrt(3.0), (double)i / (double)(maxFrames - 1));
    cout << "scale: " << scale << endl;
    vector<BYTE> out(width*height * 3); // .bmp pixel buffer
    memset(&out[0], 255, out.size() * sizeof(BYTE)); // background is white

    for (int x = 0; x < width; x++)
    {
      for (int y = 0; y < height; y++)
      {
        Vector2d e = 3.0 * Vector2d(0.5 + (double)(x-width/2), 0.5 + (double)(y-width/2)) / (double)width;
        double depth = func(sqr(e[0]) - sqr(e[1]));
        int maxX = 40;
        int maxY = 40;
        bool hit = true;
        for (int X = -maxX; X <= maxX; X++)
        {
          for (int Y = -maxY; Y <= maxY; Y++)
          {
            if (X == 0 && Y == 0)
              continue;
            Vector2d e2 = Vector2d(X, sqrt(3)*Y);
            Vector2d e3 = NE*(NE.dot(e2) * scale) + SE*(SE.dot(e2) / scale);
            e3 -= e;
            double depth2 = func(sqr(e3[0]) - sqr(e3[1]));
            if (depth2 < depth)
              hit = false;
          }
        }

        int shade = hit ? 255 : 0;
        int ind = 3 * (x + y*width);
        out[ind + 0] = shade;
        out[ind + 1] = shade;
        out[ind + 2] = shade;
      }
    }

    BYTE* c = ConvertRGBToBMPBuffer(&out[0], width, height, &s2);
    wstringstream str;
    str << L"pixel/test" << i << L".bmp";
    wstring strng = str.str();
    const TCHAR * bll = strng.c_str();
    LPCTSTR file = bll;
    SaveBMP(c, width, height, s2, file);
    delete[] c;
  }
}
