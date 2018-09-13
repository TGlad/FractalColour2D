
main()
{
  double angle[2] = { 0.0825, -1.832775 };
  Vector2d pos[0] = Vector2d(0.4493, 0.0);
  Vector2d pos[1] = Vector2d(-0.5925, 0.53515);
  double anglei[2] = { 2.00338, 0.6476 };
  Vector2d posi[0] = Vector2d(-0.2417, -0.3413);
  Vector2d posi[1] = Vector2d(-0.1895, -0.0965);

  double scalei = 0.6933;

  Matrix2d mat[2];  // mat + pos i for shape 1 equals shape i in orientation defined by rmat + rpos i (which can have any scale) 
  Matrix2d mati[2];
  Vector2d v1(cos(angle[0]), sin(angle[0]));
  Vector2d v2(cos(angle[1]), sin(angle[1]));
  mat[0] << v1[0], v1[1],
           -v1[1], v1[0];
  mat[1] << v2[0], v2[1],
           -v2[1], v2[0];

  mati[0] << cos(anglei[0]), sin(anglei[0]),
            -sin(anglei[0]), cos(anglei[0]);
  mati[0] *= scalei;

  mati[1] << cos(anglei[1]), sin(anglei[1]),
            -sin(anglei[1]), cos(anglei[1]);
  mati[1] *= scalei;

  vector<Vector2d> points[2], newPoints[2];
  newPoints[1].push_back(Vector2d(0, 0));
  int iterations = 28; // must be even 
  for (int i = 0; i < top; i++)
  {
    int p = i % 2;
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
        if (i >= iterations - 2) // earlier iterations haven't converged yet
          putpixel(point, "white");
    }
    else
    {
      for (auto &point : newPoints[p])
      {
        if (i >= iterations - 1) // earlier iterations haven't converged yet
        {
          putpixel(mat[0] * point + pos[0], "yellow");
          putpixel(mat[1] * point + pos[1], "magenta");
        }
      }
    }
  }
}

