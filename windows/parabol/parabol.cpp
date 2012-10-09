// parabol.cpp: главный файл проекта.

#include "stdafx.h"

#include "stdio.h"
#include "omp.h"
#define _USE_MATH_DEFINES
#include "math.h"
#include "time.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>

using namespace System;
using namespace std;

double L, N, h, delta, T, K, koeff_teploprov, lambda, c, ro, gam;
int i,j;

inline double right_side(double x, double t)
{
	double y = sin(x*t);
	return y;
}

inline double initial_conditions(double x)
{
	double y = exp(0.15 * x);
	return y;
}

inline double left_boundary(double t)
{
	double y = 1;
	return y;
}

inline double right_boundary(double t)
{
	double y = 2.117;
    return y;
}

double parabol()
{
    double  **v = new  double*[N+1];
      for (i = 0; i<=N; i++)
            v[i]=new  double[K];

    double *x = new double[N+1];
    double *t = new double[K+1];

    h=L/N;
    delta=T/K;
 
    gam=koeff_teploprov * koeff_teploprov * delta/h*h;

    double start = omp_get_wtime ();
   

      for(i = 0; i <= N-1; i++)
        {
          x[i]=(i-1) * h;
          for (j = 0;  j <= K-1; j++)
          {
                    t[j]=(j-1) * delta;
            if(j == 0)
            {
                 v[i][j] = initial_conditions(x[i]);
            }
            else if (i == 0)
            {
                v[i][j] = left_boundary(t[j]);
                
            }
            else if (i == N-1)
            {
               v[i][j] = right_boundary(t[j]);
            }
            else
            {
                v[i][j+1] = gam * v[i-1][j] + (1-2*gam) * v[i][j] + gam * v[i+1][j] + delta * right_side(x[i],t[j]);
            }
            
            
        }
    }
    double finish = omp_get_wtime ();
    double duration = (finish - start);
	
	for ( i = 0; i <= N; i++)
        delete [] v[i];
	delete [] v;
    
    delete [] x;
    delete [] t;
	return duration;

}
int main(array<System::String ^> ^args)
{
    delta = 0.001;
    L = 1;
    T = 100;
    N = 100;
    K = 100;
    lambda = 46,5;
    c = 460;
    ro = 7,7*10^3;

    koeff_teploprov = lambda/c*ro;

    double duration = parabol();
}
