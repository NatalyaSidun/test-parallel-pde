// matrices.cpp: главный файл проекта.

#include "stdafx.h"
#include "stdio.h"
#include "omp.h"
#define _USE_MATH_DEFINES
#include "math.h"
#include <time.h>

#include "stdlib.h"

#include <iostream>
#include <fstream>

using namespace std;

double duration, duration_parallel, acceleration;
int N, T, i, j, num_threads, Tstart, Tfinish, Nstart, Nfinish, Th, Nh;


double s;

std::fstream resFile;

inline double matrix_vector(int N, int T)
{
	double  **v = new  double*[N+1];
      for (i = 0; i<=N; i++)
            v[i]=new  double[T];

	double *f =new  double [T];

	double *r =new  double [N];


	   for (j = 0; j <= T-1; j++ ) 
		{ 
			f[j] = rand();
		}

	  srand ( time(0) );

	  for (i = 0; i <= N-1; i++ ) 
			{
					for (int j = 0; j <= T-1; j++ ) 
					{ 
						v[i][j] = rand();
					}
			}
	
	/*
		Вывод массивов на экран
	*/
	/*printf("Вектор: ");
	getchar();
	for (int j = 0; j <= T-1; j++ ) 
	{
		{
			printf("%4.10f \n", f[j]);
		}
	}
	getchar();
	printf("Матрица: ");
	getchar();
	for ( int i = 0; i <= N-1; i++ ) 
			{
					for (int j = 0; j <= T-1; j++ ) 
					{ 
						printf("%4.10f   ", v[i][j]);
					}
					
				printf("\n");
			}
	getchar();*/
double start = omp_get_wtime ();


	for ( int i = 0; i <= N-1; i++ ) 
			{
				
				r[i] = 0;

					for ( j = 0; j <= T-1; j++ ) 
					{ 
						r[i]  += v[i][j] * f[j];
					
					}
					
				
			}




	double finish = omp_get_wtime ();
	double wtick = omp_get_wtick( );

	for ( i = 0; i <= N; i++)
        delete [] v[i];
	delete [] v;
	delete [] f;
	delete [] r;

	double duration = (finish - start);
	return duration;
}

inline double matrix_vector_parallel(int num_threads, int N, int T)
{
	double  **v = new  double*[N+1];
      for ( i = 0; i<=N; i++)
            v[i]=new  double[T];

	double *f =new  double [T];

	double *r =new  double [N];


	   for (j = 0; j <= T-1; j++ ) 
		{ 
			f[j] =  rand();
		}

	  srand ( time(0) );

	  for (i = 0; i <= N-1; i++ ) 
			{
					for (j = 0; j <= T-1; j++ ) 
					{ 
						v[i][j] = rand()%100;
					}
			}
	
	/*
		Вывод массивов на экран
	*/
	/*printf("Вектор: ");
	getchar();
	for (int j = 0; j <= T-1; j++ ) 
	{
		{
			printf("%4.10f \n", f[j]);
		}
	}
	getchar();
	printf("Матрица: ");
	getchar();
	for ( int i = 0; i <= N-1; i++ ) 
			{
					for (int j = 0; j <= T-1; j++ ) 
					{ 
						printf("%4.10f   ", v[i][j]);
					}
					
				printf("\n");
			}
	getchar();*/
	//double start = clock();
	  double start=omp_get_wtime ();
	  omp_set_num_threads(num_threads);

#pragma omp parallel  shared(r, v)
	{
	#pragma omp for private(i, j)  reduction(+:s)

		for (i = 0; i <= N-1; i++ ) 
			{
			 s = 0;
				r[i] = 0;
					for (j = 0; j <= T-1; j++ ) 
					{ 
						s += v[i][j] * f[j];
					
					}
					r[i] = s;
			}


	}
	double finish = omp_get_wtime();
	double wtick = omp_get_wtick();

	for ( i = 0; i <= N; i++)
        delete [] v[i];
	delete [] v;
	delete [] f;
	delete [] r;

	double duration = (finish - start);
	return duration;

}

int main()
{
   // printf("Здравствуй, мир!");

	double duration, parallel_duration, linear_duration, linear_siedel_duration, parallel_seidel_duration;

	printf("Type min matrix size: ");
	scanf ("%d",&Nstart);

	printf("Type max matrix size: ");
	scanf ("%d",&Nfinish);


	printf("Type min vector size: ");
	scanf ("%d",&Tstart);

	printf("Type max vector size: ");
	scanf ("%d",&Tfinish);

	printf("Type matrix change step: ");
	scanf ("%d",&Nh);

	printf("Type matrix change step:  ");
	scanf ("%d",&Th);

	resFile.open("result.txt", std::ios_base::out);


	for (N = Nstart; N < Nfinish; N+=Nh)
	{
		for (T = Tstart; T <= Tfinish; T+=Th)
		{
		      if (T == N)
			{
			duration = matrix_vector(N, T);
			//printf("Linear Time = %4.10f \n", duration);

			 printf(" T = %d \n   ", T);

			  resFile << "T =  ", 
			  resFile << T;
		      resFile << "\n";

			 printf(" N = %d \n   ", N);

			 
			  resFile << "N =  ", 
			  resFile << N;
		      resFile << "\n";

			for ( num_threads = 1; num_threads <= omp_get_num_procs(); num_threads *= 2 )
			{
				printf(" Num threads = %d \n   ", num_threads);
					
				resFile << "Num threads =  ", 
				resFile << num_threads;
				resFile << "\n";

				printf("Linear Time = %4.10f \n", duration);

				resFile << "Linear Time =  ", 
				resFile << duration;
				resFile << "\n";

				duration_parallel = matrix_vector_parallel(num_threads ,N, T);
				printf("Parallel Time = %4.10f \n", duration_parallel);


				resFile << "Parallel Time =  ", 
				resFile << duration_parallel;
				resFile << "\n";

				acceleration = duration - duration_parallel;

				printf(" Acceleration = %10.10f \n   ", acceleration);

				resFile << "Acceleration =   ", 
			    resFile << acceleration;
			    resFile << "\n";

				printf("................................... \n");
			    resFile << "...................... \n";
			}
		      }
			
		}
	}

	getchar();
    return 0;
}

