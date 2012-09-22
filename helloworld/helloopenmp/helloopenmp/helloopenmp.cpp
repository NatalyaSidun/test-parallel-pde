// helloopenmp.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "stdio.h"
#include "omp.h"
#define _USE_MATH_DEFINES
#include "math.h"
#include "time.h"
#include "stdlib.h"

#include "cylinder.cpp"

#define M 100 
int i, j, k, n, m, t;
double R, p, h; 

double u[M][M], a[M], c[M], conc[k][n];

int _tmain(int argc, _TCHAR* argv[])
{   
	double start_time;
	double finish_time;

	srand ( time(0) );
	int i,j;
    int rank, size;
    int nr;
	int numth;
	numth = omp_get_num_procs();
	printf("%d procs \n",numth);
	for( i = 0; i < M; i++ ) {
        for( j = 0; j < M; j++ ) {
            u[i][j] = rand();
			printf("%f   ", u[i][j]);
        }
		printf("\n");
        a[i] = rand();
        c[i] = 0.0;
    }
	 omp_set_num_threads(4);
	 	start_time = time(0);
#pragma omp parallel private (rank, size)
	{
	 size = omp_get_num_threads();
     rank = omp_get_thread_num();
     nr = M / size;
	 #pragma omp for schedule (static, nr) private(j)

        for( i = 0; i < M; i++ ) {
            for( j = 0; j < M; j++ ) {
                c[i] += (u[i][j] * a[j]);
            }
        }
	
    }
	finish_time = time(0);
	for( i = 0; i < M; i++ ) {
        printf("C[%d] = %f\n", i, c[i] );
	}
	printf("%f   ", finish_time-start_time);
	getchar();
    

	printf("Enter start concentration");
	p = getchar();

	printf("Enter number H");
	double H = getchar();

	printf("Enter number Q1");
	double Q1 = getchar();

	printf("Enter number V");
	double V = getchar();

	h = 0.5;
	R = 500;
	m = ceil(R/h);
	double Kn=0;

	double Qst = 1.2;
	double pi  = double(M_PI);

	double Um = GetRegress(V, H);
	double d   = 0.032+21.8m;
	k = 0;
	double beta1=abs(d-1.2)/2*pi*H;

	for (i <= m; i++;) 
	{
		conc[i][0] = p;
	}
	i = 0;
	for (i <= m; i++;)
	{
		 j = 0;
		 for(j <= t; j++;)
		 {
			 if (i == 0 )
			 {
				 conc[i][j] = ((2*Qst/2*pi*H)*Q1+abs(1-((Qst/2*pi*H)/h^2-(2*d/h^2)+Kn*1))*conc(i,j)+abs((d/h^2)-(Qst/2*pi*H/(h^2)))*conc(i+1, j));
			 }
			 if (i == m)
			 {
				 conc[i][j] = ((d+beta1)/((2*j-1)*h^2)*p+abs(1-2*d*(1/h^2)+Kn*1)*conc(j,i)+abs(d-beta1)/((2*j-1)*h^2)*conc(j-1,i));
			 }
		 }
	}
}
double GetRegress(double V, double H)
{
	/*double a0 = 0.003613; 
    double a1 = -0.0002751;
    double a2 = 0.001180;
    double a3 = 0.0001461;
    double a4 = 0.000009729;
    double a5 = -0.0007189;
    double a6 = 0.00009925;
    double a7 = 0.0000003875;

    double Um = a0 + a1 * V + a2 * H + a3 * V^2 + a4 * H^2 + a5 * V * H + a6 * V^2 * H + a7 * V * H^2;
	return Um;*/
}
