// wave.cpp: главный файл проекта.
#include "stdafx.h"
#include "stdio.h"
#include "omp.h"
#define _USE_MATH_DEFINES
#include "math.h"
#include "time.h"
#include "stdlib.h"

using namespace System;

int main(int argc, char **argv)
{
	double start, finish;
    double duration;

    int N = 10; //количество разбиений по координате
	int T = 10; //количество разбиений по времени

	double t0 = 0;
	double tk = 1;
	
	double x0 = 0;
	double xk = 1;

	double tau = (tk-t0)/T;
	double h   = (xk-x0)/N;

	double v[10][10];

	double vj0 = 0;
	double vjN = 0;

	double x, t;
	int i, j;
	double alpha = 0.001;
	
	int rank, size;
    int nr;

	printf("tau = %4.10f \n", tau);
	printf("h = %4.10f \n", h);

	
	getchar();

	if (alpha*tau >= h)
	{
		printf("—хема неустойчива!");
	}
	start = time(0);

	for (j = 0; j <= T; j++)
	{
		v[0][j]   = 0;
		v[N][j]   = 0;
	}
	for (i = 1; i < N; i++)
	{
		x = x0 + (h*i);

		v[i][0] = 0.1 * sin(M_PI*x); 
		
			
	}
	 omp_set_num_threads(1);
#pragma omp parallel private (rank, size)
	{
	 size = omp_get_num_threads();
     rank = omp_get_thread_num();
	  nr = N / size;
#pragma omp for schedule (static, nr) private(i,j)
	   
	for( i = 1; i <= N-1; i++ )
	{
        for(j = 1; j < T; j++ ) 
		{
			x = x0 + (h*i);
			t  = t0 + (tau*j);
			v[i][j+1] = 2*v[i][j] - v[i][j-1] + (alpha*tau/h)*(v[i+1][j] - 2*v[i][j] + v[i-1][j]);
        }
	}

	  }
	  for( i = 0; i < N; i++ )
	{
        for(j = 0; j < T; j++ ) 
		{
			printf("%f   ", v[i][j]);
        }
	}
	
	
	/*for (i = 1; i < N; i++)
	{
		
		for (j = 0; j < T-1; j++)
		{
			x = x0 + (h*i);
			t  = t0 + (tau*j);
			v[i][j+1] = 2*v[i][j] - v[i][j-1] + (alpha*tau/h)*(v[i+1][j] - 2*v[i][j] + v[i-1][j]);
		}
	}*/
	finish = time(0);
	/*for( i = 0; i < N; i++ )
	{
        for(j = 0; j < T; j++ ) 
		{
			printf("%f   ", v[i][j]);
        }
		getchar();
	}*/

	duration = (finish - start);		
	printf(" Time = %f \n   ", duration);
	//getchar();
	return 0;
}
