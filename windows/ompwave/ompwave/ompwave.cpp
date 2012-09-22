// ompwave.cpp: главный файл проекта.

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

double t0 = 0;//
double tk = 2;

int i, j, N, Nstart, Nfinish, Tstart, Tfinish, T, Th, Nh, num_threads;
	
double x0 = 0;
double xk = 1;
double tau, h, y, x, alpha, e;

std::fstream resFile;



inline double linear_wave(double x0, double h, int N, int T, double alpha, double tau)
{

		double  **f = new  double*[N+1];
      for (int i = 0; i<=N; i++)
            f[i]=new double[T];


double start = omp_get_wtime ();
for( int i = 0; i <= N-1; i++ )
	{
		double y = x0 + (h*i);
		for(int j = 0; j <= T-1; j++ ) 
			{
				
				if (i == 0 || i == N-1)
					{
						f[i][j] = 0;
					}
				else if ( j == 0 || j == T-1)
					{
						f[i][j] = 0.1 * sin(M_PI*y); 
					}
				else
					{
						f[i][j+1] = 2*f[i][j] - f[i][j-1] + (alpha*tau/h)*(f[i+1][j] - 2*f[i][j] + f[i-1][j]);
					}
			}
	}
	double finish = omp_get_wtime ();
	double wtick = omp_get_wtick();

	double duration = (finish - start);
	
	for ( i = 0; i <= N; i++)
        delete [] f[i];
	delete [] f;

	return duration;
}

inline double parallel_wave(int num_threads, double x0, double h, int N, int T, double alpha, double tau)
{
	double  **v = new  double*[N+1];
      for (int i = 0; i<=N; i++)
            v[i]=new  double[T];

omp_set_num_threads(num_threads);

double start = omp_get_wtime ();
#pragma omp parallel  shared (tau,h, alpha, x)
{
	
	#pragma omp for private(i,j)
	   
		for(int i = 0; i <= N-1; i++ )
			{

//			#pragma omp critical
//				{
					double x = x0 + (h*i);
	//			}
				

				for( int j = 0; j <= T-1; j++ ) 
				{
					if (i == 0 || i == N-1)
					{
						v[i][j] = 0;
					}
					else if ( j == 0 || j == T-1)
					{
						v[i][j] = 0.1 * sin(M_PI*x); 
					}
					else
					{
						v[i][j+1] = 2*v[i][j] - v[i][j-1] + (alpha*tau/h)*(v[i+1][j] - 2*v[i][j] + v[i-1][j]);
					}
        }
	}
 
}

	double finish = omp_get_wtime ();
	double wtick = omp_get_wtick( );
	double duration = (finish - start);
	 
	for ( i = 0; i <= N; i++)
        delete [] v[i];
	delete [] v;

return duration;
}

inline double linear_siedel(double x0, double xk, double h, int N, int T)
{
	double dmax, temp, dm, eps;
	 srand ( time(NULL) );


	double  **u = new  double*[N];
      for (int i = 0; i<=N; i++)
            u[i]=new double[T];
	

	  getchar();
	
	double start = omp_get_wtime ();
	for ( int i = 1; i < N-1; i++ ) 
			{
					for (int j = 1; j < T-1; j++ ) 
					{ 
						u[i][j] = rand() % 100 + 100;
					}
			}
	do { 
			dmax = 0;
		
			u[0][0] = x0;
			u[N-1][T-1] = x0;
			u[0][T-1] = xk;
			u[N-1][0] = xk;
			for ( int i = 1; i < N-1; i++ ) 
			{
					for (int j = 1; j < T-1; j++ ) 
					{ 
							temp = u[i][j]; 

								u[i][j] = 0.25 * (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]); 
							
							
							dm = fabs(temp-u[i][j]); 
							if ( dmax < dm ) dmax = dm; 
					} 
				
					
			}
	}
	while (dmax>eps);

	double finish = omp_get_wtime ();
	double wtick = omp_get_wtick( );
	double duration = (finish - start);
	return duration;
}

int main(array<System::String ^> ^args)
{
	double duration, parallel_duration, linear_duration, linear_siedel_duration, parallel_seidel_duration;

	Console::WriteLine(L"¬ведите минимальный размер сетки по координате: ");
	scanf ("%d",&Nstart);

	Console::WriteLine(L"¬ведите максимальный размер сетки по координате: ");
	scanf ("%d",&Nfinish);


	Console::WriteLine(L"¬ведите минимальный размер сетки по времени: ");
	scanf ("%d",&Tstart);

	Console::WriteLine(L"¬ведите максимальный размер сетки по времени: ");
	scanf ("%d",&Tfinish);

	Console::WriteLine(L"¬ведите шаг изменений по координате: ");
	scanf ("%d",&Nh);

	Console::WriteLine(L"¬ведите шаг изменений по времени: ");
	scanf ("%d",&Th);
	
	resFile.open("result.txt", std::ios_base::out);


	for (N = Nstart; N < Nfinish; N+=Nh)
	{
		for (T = Tstart; T <= Tfinish; T+=Th)
		{
			 tau = (tk-t0)/T;
			 h   = (xk-x0)/N;
			 e = 0.1;
			 alpha = 1;

			 printf(" T = %d \n   ", T);

			  resFile << "T =  ", 
			  resFile << T;
		      resFile << "\n";

			 printf(" N = %d \n   ", N);

			 
			  resFile << "N =  ", 
			  resFile << N;
		      resFile << "\n";


			 if (alpha*tau >= h)
				{
					Console::WriteLine(L"—хема неустойчива!");
					continue;
				}

			 for ( num_threads = 1; num_threads <= omp_get_num_procs(); num_threads *= 2 )
			 {
					printf(" Num threads = %d \n   ", num_threads);
					
					 resFile << "Num threads =  ", 
					 resFile << num_threads;
					 resFile << "\n";

					linear_duration = linear_wave(x0, h, N, T, alpha, tau);

					parallel_duration = parallel_wave(num_threads, x0, h, N, T, alpha, tau);
	        
					 double acceleration = linear_duration - parallel_duration;
					 printf(" Acceleration = %10.10f \n   ", acceleration);

					 resFile << "Acceleration =   ", 
					 resFile << acceleration;
					 resFile << "\n";

					 Console::WriteLine(L"................................... \n");
					 resFile << "...................... \n";

			 }
			 
			 
			 Console::WriteLine(L"------//------------- \n");
			  resFile << "----------//------------- \n";
			
		}
	}
	
	resFile.close();
	//linear_siedel_duration = linear_siedel(x0, xk, h, N, T);	
	//printf(" Time linear_seidel = %f \n   ", linear_siedel_duration);
	
	getchar();
    return 0;
}

