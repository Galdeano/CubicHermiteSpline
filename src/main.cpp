
#include <cstdio>  /* printf, fgets */
#include <cstdlib> /* atoi */
#include <iostream>
#include <vector>
#include "CHSpline.h"

int main(int argc, char* argv[])
{

  if ((argc != 1) && (argc != 8))
  {
    printf("Usage: %s t0 t1 p0 p1 v0 v1 Np \n", argv[0]);
    printf("As :  %s 0.0 5.0 -1.0 5.0 -1.0 0.0 100 \n", argv[0]);
    exit(1);
  }

  if (argc == 1)
  {
    printf("Demo mode, real usage: %s t0 t1 p0 p1 v0 v1 Np \n", argv[0]);
    printf("As :  %s 0.0 5.0 -1.0 5.0 -1.0 0.0 100 (as here)\n", argv[0]);
  }

  Spline Splines;
  double t0 = 0.0, t1 = 5.0, p0 = -1.0, p1 = 5.0, v0 = -1.0, v1 = 0.0;
  int Np = 100;

  if (argc == 8)
  {
    t0 = atof(argv[1]);
    t1 = atof(argv[2]);
    p0 = atof(argv[3]);
    p1 = atof(argv[4]);
    v0 = atof(argv[5]);
    v1 = atof(argv[6]);
    Np = atoi(argv[7]);
  }

  std::vector<double> ti, pi, vi;
  ti.push_back(t0);
  ti.push_back(t1);
  pi.push_back(p0);
  pi.push_back(p1);
  vi.push_back(v0);
  vi.push_back(v1);
  Splines.initSpline(ti, pi, vi);

  Splines.initSpline(ti, pi, vi);
  std::cout << "Spline=[";
  for (int i = 0; i < Np; i++)
  {
    if (i == 0)
    {
      std::cout << Splines.evalSpline(t0 + ((double)i) * (t1 - t0) / Np);
    }
    else
    {
      std::cout << "," << Splines.evalSpline(t0 + ((double)i) * (t1 - t0) / Np);
    }
  }
  std::cout << "]" << std::endl;

  return 0;
}
