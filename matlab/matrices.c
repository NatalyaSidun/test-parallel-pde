#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <omp.h>

double duration, duration_parallel, acceleration;
int N, T, i, j, num_threads, Tstart, Tfinish, Nstart, Nfinish, Th, Nh;

//clock_t startlinear, finishlinear;
double s;

double matrix_vector(int N, int T)
{
	double  **v = new  double*[N+1];
      for (i = 0; i<=N; i++)
            v[i]=new  double[T];

	double *f = new  double [T];
	double *r = new  double [N];


	   for (j = 0; j <= T-1; j++ ) 
		{ 
			f[j] = rand()%100;
		}

	  srand ( time(0) );

	  for (i = 0; i <= N-1; i++ ) 
			{
					for (int j = 0; j <= T-1; j++ ) 
					{ 
						v[i][j] = rand()%100;
					}
			}
      double startlinear = omp_get_wtime();
	for ( int i = 0; i <= N-1; i++ ) 
			{
				
				r[i] = 0;

					for ( j = 0; j <= T-1; j++ ) 
					{ 
						r[i] += sin(exp(v[i][j]) * cos(f[j]));
					
					}
			}

 double finishlinear = omp_get_wtime();

 double duration = ( finishlinear - startlinear ) / (double) CLOCKS_PER_SEC;

	for ( i = 0; i <= N; i++)
        delete [] v[i];
	delete [] v;
	delete [] f;
	delete [] r;
  return duration;
}
double matrix_vector_parallel(int num_threads, int N, int T)
{
	double  **v = new  double*[N+1];
      for ( i = 0; i<=N; i++)
            v[i]=new  double[T];

	double *f =new  double [T];
	double *r =new  double [N];

	 srand ( time(0) );

	   for (j = 0; j <= T-1; j++ ) 
		{ 
			f[j] =  rand()%100;
		}
	 

	  for (i = 0; i <= N-1; i++ ) 
			{
					for (j = 0; j <= T-1; j++ ) 
					{ 
						v[i][j] = rand()%100;
					}
			}
	
	  omp_set_num_threads(num_threads);
	double start=omp_get_wtime ();
double k;
#pragma omp parallel  
	{

	#pragma omp for firstprivate(j) lastprivate(i) reduction(+: k)

		for (i = 0; i <= N-1; i++ ) 
			{
				r[i] = 0;
				 k = 0;
					for (j = 0; j <= T-1; j++ ) 
					{ 
						k += sin(exp(v[i][j]) * cos(f[j]));
					
					}
				r[i] = k;
			}


	}
	double finish = omp_get_wtime();

	for ( i = 0; i <= N; i++)
        delete [] v[i];
	delete [] v;
	delete [] f;
	delete [] r;

	double duration = (finish - start);
	return duration;
}
void mexFunction()
{
    double duration, parallel_duration, linear_duration;
    Nstart = 100;
    Nfinish = 1000;
    Nh = 100;
    Tstart = 100;
    Tfinish = 1000;
    Th = 100;
    for (N = Nstart; N < Nfinish; N+=Nh)
	{
		for (T = Tstart; T <= Tfinish; T+=Th)
		{
			if(N == T) 
			{
                mexPrintf(" T = %d \n   ", T);
                
                duration = matrix_vector(N, T);
                mexPrintf(" N = %d \n   ", N);
                
                for ( num_threads = 1; num_threads <= omp_get_num_procs(); num_threads *= 2 )
                {
                    mexPrintf(" Num threads = %d \n   ", num_threads);
                    mexPrintf("Linear Time = %4.10f \n", duration);
                    duration_parallel = matrix_vector_parallel(num_threads ,N, T);
                    
                    mexPrintf("Parallel Time = %4.10f \n", duration_parallel);
                    
                    acceleration = duration - duration_parallel;

                    mexPrintf(" Acceleration = %10.10f \n   ", acceleration);
                    mexPrintf(" ................................... \n");
                    
                }
            }
}
