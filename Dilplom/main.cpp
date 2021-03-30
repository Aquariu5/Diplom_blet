#include <iostream>
#include <cmath>
#include <chrono>

#include "temptest.h"
#include "temptest2d.h"
using namespace std;

int main()
{
	srand(time(NULL));
	adept::Stack stack;
	/*temptest temp(501, 50000, 1e-2); // кол-во узлов, временных шагов, тау
	/*temp.analytical(); 
	//temp.numerical_classic_with_k();
	temp.numerical_classic_with_sigma();
	//temp.one_equation();

	//-ku
	/*temptest temp(201, 20000, 1e-4); // кол-во узлов, временных шагов, тау
	temp.analytical(); 
	temp.numerical_classic_with_k();
	cout << temp.max << endl;
	*/

	auto begin = chrono::system_clock::now();
	temptest2d temp(501, 50000, 1e-4);
	//temp.numerical_classic_with_sigma();
	temp.pyroman();
	//temp.norm_with_sigma();
	//temp.init_fractures();
	auto end = chrono::system_clock::now();
	std::cout << "time : " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;

	system("pause");
	return 0;
}