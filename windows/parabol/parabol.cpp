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

double L, h, delta, T, koeff_teploprov, lambda, c, ro, gam;
int i, j, N, K, Nh, Kh, Nstart, Nfinish, Kstart, Kfinish, num_threads;

std::fstream resFile;

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
    double  **v = new  double *[N+1];
      for (i = 0; i<=N; i++)
            v[i]=new  double[K+1];

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

double parallel_parabol(int num_threads)
{
    double  **v = new  double *[N+1];
      for (i = 0; i <= N; i++)
            v[i] = new  double[K+1];

    double *x = new double[N+1];
    double *t = new double[K+1];

    h=L/N;
    delta=T/K;
 
    gam=koeff_teploprov * koeff_teploprov * delta/h*h;
omp_set_num_threads(num_threads);
    double start = omp_get_wtime ();
   
#pragma omp parallel shared (h, delta, x, t, K, N, v, gam) private (i,j)
    {
   #pragma omp for 
      for(i= 0; i <= N-1; i++)
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
    L = 1;
    T = 100;
    N = 100;
    K = 100;
    lambda = 46,5;
    c = 460;
    ro = 7,7*10^3;

    printf("Type length, м: \n");
	scanf_s ("%d",&L);

    printf("Type time, с: \n");
	scanf_s ("%d",&T);



    printf("Type min coordinates mesh size: \n");
	scanf_s ("%d",&Nstart);

	printf("Type max coordinates mesh size: \n");
	scanf_s ("%d",&Nfinish);


	printf("Type min time mesh size: \n");
	scanf_s ("%d",&Kstart);

	printf("Type max time mesh size: ");
	scanf_s ("%d",&Kfinish);

	printf("Type coordinates increment step: \n");
	scanf_s ("%d",&Nh);

	printf("Type time increment step: \n");
	scanf_s ("%d",&Kh);

    resFile.open("result.txt", std::ios_base::out);

    printf("Type lambda: \n");
    printf("lambda = %5.5", lambda);
   // scanf_s ("%5.5",&lambda);
    resFile << "lambda =  ", 
    resFile << lambda;
    resFile << "\n";

    printf("Type density: \n");
    printf("density = %5.5", ro);
    //scanf_s ("%5.5",&ro);
    resFile << "density =  ", 
    resFile << ro;
    resFile << "\n";


    printf("Type specific heat: \n");
    printf("specific heat = %5.5", c);
  //  scanf_s ("%10.10",&c);

    resFile << "specific heat =  ", 
    resFile << c;
    resFile << "\n";


    koeff_teploprov = lambda/c*ro;
    printf("Type thermal conductivity coefficient= %5.5", koeff_teploprov);

    resFile << "thermal conductivity coefficient=  ", 
    resFile << koeff_teploprov;
    resFile << "\n";

    
    printf(" Time T = %d \n   ", T);

    resFile << "Time T =  ", 
	resFile << T;
    resFile << "\n";

    printf("Length L = %d \n   ", T);

    resFile << "Length L =  ", 
	resFile << L;
    resFile << "\n";


     printf("Min coordinates mesh size Nstart = %d \n   ", N);

    resFile << "Min coordinates mesh size Nstart =  ", 
	resFile << Nstart;
    resFile << "\n";


    printf("min time mesh size Kstart = %d \n   ", K);

    resFile << "min time mesh size Kstart =  ", 
	resFile << Kstart;
    resFile << "\n";


    printf("max coordinates mesh size Nstart = %d \n   ", N);

    resFile << "max coordinates mesh size Nfinish =  ", 
	resFile << Nfinish;
    resFile << "\n";

    printf("coordinates increment step Nh = %d \n   ", Nh);

    resFile << "coordinates increment step Nh =  ", 
	resFile << Nh;
    resFile << "\n";


    printf("time increment step Kh = %d \n   ", Kh);

    resFile << "time increment step Kh =  ", 
	resFile << Kh;
    resFile << "\n";

    printf("------Start Computing!------------- \n");
    resFile << "----------Start Computing!------------- \n";
    for (N = Nstart; N < Nfinish; N+=Nh)
	{
		for (K = Kstart; K <= Kfinish; K+=Kh)
        {
             printf(" K = %d \n   ", K);

			  resFile << "K =  ", 
			  resFile << K;
		      resFile << "\n";

              printf(" N = %d \n   ", N);

			  resFile << "N =  ", 
			  resFile << N;
		      resFile << "\n";

              double duration = parabol();

               printf(" Linear duration = %10.10f \n   ", duration);

			   resFile << "Linear Duration =   ", 
			   resFile << duration;
			   resFile << "\n";

              for ( num_threads = 1; num_threads <= omp_get_num_procs(); num_threads++ )
              {
                  

                 printf("Num threads = %d \n", num_threads);

                 resFile << "Num threads =   ", 
				 resFile << num_threads;
				 resFile << "\n";

                 double parallel_duration = parallel_parabol(num_threads);
                 
                 printf("Parallel duration = %10.10f \n", parallel_duration);

                  resFile << "Parallel duration =   ", 
				  resFile << parallel_duration;
				  resFile << "\n";

                  double acceleration = duration - parallel_duration;

					 printf(" Acceleration = %10.10f \n   ", acceleration);

					 resFile << "Acceleration =   ", 
					 resFile << acceleration;
					 resFile << "\n";


                     printf("................................... \n");
					 resFile << "...................... \n";
              }

               printf(".............../Finish iteration!/.................... \n");
			   resFile << "........../Finish iteration!/............ \n";

        }
    }
 resFile.close();
}
