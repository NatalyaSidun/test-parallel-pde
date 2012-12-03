// matches.cpp: главный файл проекта.

#include "stdafx.h"
#include "math.h"
#include <conio.h>
#include <iostream>
#include <iomanip>
using namespace System;
using namespace std;
int main()
{
   int N, K;
   cout << "Enter cubes' quantity: \n";
   cin>>N;
   K = (N-1)*8+12;
   cout << "You need " << setw(3) << K <<"  matches"<< endl; 
   getch();

}
