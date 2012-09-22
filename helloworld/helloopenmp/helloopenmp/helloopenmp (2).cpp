// helloopenmp.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "stdio.h"
#include "omp.h"


int _tmain(int argc, _TCHAR* argv[])
{
	int numth;
	numth = omp_get_num_procs();
	int a; 
	int b;
	printf("%d procs \n",numth);
#pragma omp parallel num_threads(numth)
	b = omp_get_thread_num();
	printf ("current thread %d", b);
	a=omp_get_num_threads();
	printf("%d threads \n", a );
	getchar;
	return 0;
}

