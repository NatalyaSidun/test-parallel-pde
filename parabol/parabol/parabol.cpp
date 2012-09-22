// parabol.cpp: главный файл проекта.
#include "stdafx.h"
#include "stdio.h"
#include "omp.h"
#define _USE_MATH_DEFINES
#include "math.h"
#include <time.h>

#include "stdlib.h"

#include <iostream>
#include <fstream>

using namespace System;
using namespace std;

double duration, duration_parallel, acceleration;
int N, T, i, j, num_threads, Tstart, Tfinish, Nstart, Nfinish, Th, Nh;

inline double parabol_linear(int N, int T)
{
	double  **v = new  double*[N+1];
      for (i = 0; i<=N; i++)
            v[i]=new  double[T];

	double *f =new  double [T];



	double *r =new  double [N];

}
int main(array<System::String ^> ^args)
{

}
