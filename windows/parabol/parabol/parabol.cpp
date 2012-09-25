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
double x0, xk, t0, tk, gam, a, x, t, Th, Nh;
int N, T, i, j, num_threads, Tstart, Tfinish, Nstart, Nfinish;

std::fstream resFile;
clock_t startlinear, finishlinear;

inline double parabol_linear(int N, int T)
{
	double  **v = new  double*[N+1];
      for ( i = 0; i<=N; i++)
            v[i]=new  double[T];

	//начальные услови€  
	for (i = 1; i < N; i++)
	{
		x = x0 + (Nh * i);

		v[i][0] = exp(0.15 * x); 
			
	}
	//граничные услови€
	for (j = 0; j <= T-1; j++)
	{
		v[0][j]   = 1;
		v[N][j]   = 2.117;
	}

	gam=a*a*Th/Nh*Nh;

	startlinear = clock();

	for (i = 1; i < N-1; i++)
	{
		for(j = 1; j < T-1; j++)
		{
			x = x0 + (Nh * i);
			t = t0 + (Th * j); 
			v[i][j+1] = gam*v[i-1][j] + (1 - 2*gam)*v[i][j] + gam*v[i+1][j] - Th*sin(x*t); 
		}
	}
 finishlinear = clock();

 	for ( i = 0; i <= N; i++)
        delete [] v[i];
	delete [] v;

 double duration = ( finishlinear - startlinear ) / (double) CLOCKS_PER_SEC;
  
 return duration;

}

inline double parabol_parallel(int num_threads, int N, int T)
{
	double  **v = new  double*[N+1];
      for (i = 0; i <= N; i++)
            v[i] = new  double[T];

	//начальные услови€  
	for (i = 1; i < N; i++)
	{
		x = x0 + (Nh * i);

		v[i][0] = exp(0.15 * x); 
			
	}
	//граничные услови€
	for (j = 0; j <= T; j++)
	{
		v[0][j]   = 1;
		v[N][j]   = 2.117;
	}

	gam=a*a*Th/Nh*Nh;
	omp_set_num_threads(num_threads);
	double start = omp_get_wtime ();
#pragma omp parallel 
	{	
	#pragma omp for
		for (i = 1; i <= N-1; i++)
			{
			for(j = 1; j <= T-1; j++)
				{
					x = x0 + (Nh * i);
					t = t0 + (Th * j); 
					v[i][j+1] = gam*v[i-1][j] + (1 - 2*gam)*v[i][j] + gam*v[i+1][j] - Th*sin(x*t); 
				}		
			}
	}
	double finish = omp_get_wtime ();

	for ( i = 0; i <= N; i++)
        delete [] v[i];
	delete [] v;

	double duration = (finish - start);
	return duration;
}
int main(array<System::String ^> ^args)
{
	double duration, parallel_duration, linear_duration;

	Console::WriteLine(L"¬ведите начальное значение по координате \n");
	scanf ("%d",&x0);

	Console::WriteLine(L"¬ведите конечное значение по координате \n");
	scanf ("%d",&xk);

	Console::WriteLine(L"¬ведите начальное значение по времени \n");
	scanf ("%d",&t0);

	Console::WriteLine(L"¬ведите конечное значение по времени \n");
	scanf ("%d",&tk);
	/*x0 = 0;
	xk = 1;
	t0 = 0;
	tk = 10; */

	a = 0.4;
	
	Console::WriteLine(L"¬ведите минимальное число разбиений по координате \n");
	scanf ("%d",&Nstart);

	Console::WriteLine(L"¬ведите максимальное число разбиений по координате: \n");
	scanf ("%d",&Nfinish);

	Console::WriteLine(L"¬ведите минимальное число разбиений по времени: \n");
	scanf ("%d",&Tstart);

	Console::WriteLine(L"¬ведите максимальное число разбиений по времени: \n");
	scanf ("%d",&Tfinish);


	Console::WriteLine(L"¬ведите шаг изменений разбиений по координате: \n");
	scanf ("%d",&Nh);

	Console::WriteLine(L"¬ведите шаг изменений разбиений по времени: \n");
	scanf ("%d",&Th);

	resFile.open("result.txt", std::ios_base::out);

	getchar();

	printf(" Tstart = %d \n   ", Tstart);
	resFile << "Tstart =  ", 
	resFile << Tstart;
	resFile << "\n";

	printf(" Tfinish = %d \n   ", Tfinish);
	resFile << "Tfinish =  ", 
	resFile << Tfinish;
	resFile << "\n";

	printf(" Nstart = %d \n   ", Nstart);
	resFile << "Nstart =  ", 
	resFile << Nstart;
	resFile << "\n";

	printf(" Nfinish = %d \n   ", Nfinish);
	resFile << "Nfinish =  ", 
	resFile << Nfinish;
	resFile << "\n";

	printf(" Nh = %d \n   ", Nh);
	resFile << "Nh =  ", 
	resFile << Nstart;
	resFile << "\n";

	printf(" Th = %d \n   ", Th);
	resFile << "Th =  ", 
	resFile << Nstart;
	resFile << "\n";

	printf("................................... \n");
	resFile << "...................... \n";

	printf(" x0 = %d \n   ", x0);
	resFile << "x0 =  ", 
	resFile << x0;
	resFile << "\n";

	printf(" xk = %d \n   ", xk);
	resFile << "x0 =  ", 
	resFile << xk;
	resFile << "\n";

	printf(" t0 = %d \n   ", t0);
	resFile << "t0 =  ", 
	resFile << t0;
	resFile << "\n";

	printf(" tk = %d \n   ", tk);
	resFile << "tk =  ", 
	resFile << tk;
	resFile << "\n";

	getchar();

	for (N = Nstart; N <= Nfinish; N += Nh)
		{
		for (T = Tstart; T <= Tfinish; T += Th)
		{
			Nh = (xk - x0) / N;
			Th = (tk - t0) / T;

			if(Th <= (Nh*Nh)/2*(a*a) )
			{
				
				printf(" T = %d \n   ", T);

				  resFile << "T =  ", 
				  resFile << T;
			      resFile << "\n";

				 printf(" N = %d \n   ", N);

			 
				  resFile << "N =  ", 
			      resFile << N;
				  resFile << "\n";


				duration = parabol_linear(N, T);
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

					duration_parallel = parabol_parallel(num_threads ,N, T);
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
					getchar();
				}
			}
			else
			{

				Console::WriteLine(L"—хема неустойчива!");
				printf(" N = %d \n   ", N);
				resFile << "N =  ", 
				resFile << N;
			    resFile << "\n";


				printf(" Nh = %10.4f \n   ", Nh);
				resFile << "Nh =  ", 
				resFile << Nh;
			    resFile << "\n";

				printf(" T = %d \n   ", T);
				resFile << "T =  ", 
				resFile << T;
			    resFile << "\n";

				printf(" Th = %10.4f \n   ", Th);
				resFile << "Th =  ", 
				resFile << Th;
			    resFile << "\n";


				printf(" (Nh*Nh)/2*(a*a) = %10.4f \n", (Nh*Nh)/2*(a*a));
				resFile << "(Nh*Nh)/2*(a*a) =  ", 
				resFile << (Nh*Nh)/2*(a*a);
			    resFile << "\n";

				printf("................................... \n");
				resFile << "...................... \n";
				//getchar();
			}
		}
	}
}
