#include "temptest.h"
#include "mymath.h"
#include <iostream>
temptest::temptest()
{
}

temptest::temptest(int n, int T, double tau)
{

	a = 0;
	l = 100.;
	d = (n - 1) * 0.2; // узлы для граничек
	this->n = n;
	this->T = T;
	h = (l - a) / (n - 1);
	this->tau = tau;

	uold.resize(n);
	unew.resize(n);
	pold.resize(n);
	pnew.resize(n);

	//левая полуволна
	uold_l.resize(n);
	unew_l.resize(n);
	pold_l.resize(n);
	pnew_l.resize(n);

	//правая полуволна
	uold_r.resize(n);
	unew_r.resize(n);
	pold_r.resize(n);
	pnew_r.resize(n);

	uanal.resize(n);
	panal.resize(n);

	r1.resize(n);
	r2.resize(n);
	C1.resize(n);
	C2.resize(n);

	//вектора для метода ньютона
	uold_ad.resize(n);
	unew_ad.resize(n);
	pold_ad.resize(n);
	pnew_ad.resize(n);

	xnew_ad.resize(2 * n);
	xold_ad.resize(2 * n);
	xsol_ad.resize(2 * n);

	f1_ad.resize(n);
	f2_ad.resize(n);
	f_ad.resize(2 * n);
	
	jacobian.resize(2 * n);
	for (int i = 0; i < 2 * n; ++i)
	{
		jacobian[i].resize(2 * n);
	}
	//параметры среды
	c = 15;
	rho = 10;

	//параметры затухания
	m = 1.5;

	W = 15;

	//обычное затухание
	k1 = 0;
	k2 = 0;

	//параметры начальной волны
	a1 = 1e-2; // пологость
	a2 = 1e-3; // пологость
	alpha1 = 0.1; // амплитуда
	alpha2 = 0.1; // амплитуда

	left.open("leftders.txt");
	right.open("rightders.txt");

	left << "numleft" << "\t" << "manleft" << std::endl;
	right << "numright" << "\t" << "manright" << std::endl;
}
double temptest::u0(double x)
{

	return alpha1 * exp(-a1 * (x - l / 2) * (x - l / 2));
}

double temptest::p0(double x)
{
	return alpha2 * exp(-a2 * (x - l / 2) * (x - l / 2));
}



void temptest::analytical()
{
	int k = 0;
	//h = l / (n - 1);
	rho = 583;
	c = 0.1456;
	a1 = 1e-1; // пологость
	a2 = 1e-1; // пологость
	alpha1 = 1; // амплитуда
	alpha2 = 1; // амплитуда

	double P0divU0 = 4.28 * 1e4;
	double U0 = 3.5 * 1e3;
	double P0 = 1.5 * 1e8;
	double prc = 1.7 * 1e7;

	for (int it = 0; it < T; ++it)
	{
#pragma omp parallel for
		for (int i = 0; i < n; ++i)
		{
			double u_minus = u0(i * h - c * it * tau);
			double u_plus = u0(i * h + c * it * tau);
			double p_plus = P0divU0 * p0(i * h + c * it * tau);
			double p_minus = P0divU0 * p0(i * h - c * it * tau);
			uanal[i] = 0.5 * (u_minus + u_plus + (p_minus \
				- p_plus) / (rho * c)) \
			  * exp(-k1 * it * tau);

			panal[i] = 0.5 * ((rho * c) * (u_minus - u_plus) + p_minus \
				+ p_plus) \
			  * exp(-k2 * it * tau);

		}

		if (it % 100 == 0)
		{
			std::ofstream foutU("anU" + std::to_string(k) + ".txt");
			std::ofstream foutP("anP" + std::to_string(k++) + ".txt");
			//for (int i = 0; i < n; ++i)
			for (int i = d; i < n - d; ++i)
			{

				foutU << i * h << '\t' << uanal[i] << std::endl;
				foutP << i * h << '\t' << panal[i] << std::endl;
			}
		}
	}
}

void temptest::numerical_tvd()
{
	//начальное условие
	for (int i = 0; i < n; ++i)
	{
		uold[i] = u0(i * h);
		pold[i] = p0(i * h);
	}
	std::ofstream foutU0("ch_tvdU0.txt");
	std::ofstream foutP0("ch_tvdP0.txt");

	// rj
	for (int i = 1; i < n - 1; ++i)
	{
		r1[i] = rj(i, pold);
		r2[i] = rj(i, uold);
	}

	//Cj на основе rj
	for (int i = 0; i < n; ++i)
	{
		C1[i] = Cj(i, r1);
		C2[i] = Cj(i, r2);

		foutU0 << i * h << '\t' << uold[i] << std::endl;
		foutP0 << i * h << '\t' << pold[i] << std::endl;

	}

	mu1 = tau / (rho * h);
	mu2 = rho * c * c * tau / h;

	foutU0.close();
	foutP0.close();
	for (int it = 1; it < T; ++it)
	{
		std::ofstream foutU("ch_tvdU" + std::to_string(it) + ".txt");
		std::ofstream foutP("ch_tvdP" + std::to_string(it) + ".txt");

		for (int i = 1; i < n - 1; ++i)
		{
			r1[i] = rj(i, pold);
			r2[i] = rj(i, uold);
		}
		for (int i = 0; i < n; ++i)
		{
			C1[i] = Cj(i, r1);
			C2[i] = Cj(i, r2);
		}

		for (int i = 1; i < n - 1; ++i) // края не трогаются
		{
			double mult1 = (\
				mu1 + mu1 * (1 - mu1) / 2 * \
				(C1[i] * (pold[i + 1] - pold[i]) / (pold[i] - pold[i - 1]) - C1[i - 1])\
				);
			if (!(mult1 >= 0 && mult1 <= 1))
				printf("Mult1 is bad, it's %f\n",mult1);

			check_tvd(i, pold, C1); // проверка условия 2

			unew[i] = uold[i] - (pold[i] - pold[i - 1]) * mult1;// -sigma_left(i * h) / rho - sigma_right(i * h) / rho;
				
			double mult2 = (\
				mu2 + mu2 * (1 - mu2) / 2 * \
				(C2[i] * (uold[i + 1] - uold[i]) / (uold[i] - uold[i - 1]) - C2[i - 1])\
				);

			if (!(mult2 >= 0 && mult2 <= 1))
				printf("Mult2 is bad, it's %f\n", mult2);

			check_tvd(i, uold, C2); // проверка условия 2

			pnew[i] = pold[i] - (uold[i] - uold[i - 1]) * mult2; // -sigma_left(i * h) / rho - sigma_right(i * h) / rho;
				

		} // dx cycle

		//unew[0] = unew[1];
		//unew[n - 1] = unew[n - 2];


		//for (int i = 0; i < d; ++i)
		//{
		//	unew[i] -= sigma_left(i * h) * uold[i] / rho;
		//}
		//
		//for (int i = n - d; i < n; ++i)
		//{
		//	unew[i] -= sigma_right(i * h) * uold[i] / rho;
		//}

		//отрисовка всего перед сменой слоя
		for (int i = 0; i < n; ++i)
		{
			foutU << i * h << '\t' << unew[i] << std::endl;
			foutP << i * h << '\t' << pnew[i] << std::endl;
		}

		pold = pnew;
		uold = unew;
	} // time cycle
}

void temptest::numerical_classic()
{
	//начальное условие
	for (int i = 0; i < n; ++i)
	{
		uold[i] = u0(i * h);
		pold[i] = p0(i * h);
	}
	std::ofstream foutU0("ch_clU0.txt");
	std::ofstream foutP0("ch_clP0.txt");


	for (int i = 0; i < n; ++i)
	{
		foutU0 << i * h << '\t' << uold[i] << std::endl;
		foutP0 << i * h << '\t' << pold[i] << std::endl;
	}

	for (int it = 1; it < T; ++it)
	{
		std::ofstream foutU("ch_clU" + std::to_string(it) + ".txt");
		std::ofstream foutP("ch_clP" + std::to_string(it) + ".txt");
		
		for (int i = 1; i < n - 1; ++i) // края не трогаются
		{
			mu1 = tau / (rho * h);
			mu2 = rho * c * c * tau / h;

			//unew[i] = -mu1 * (pold[i + 1] - pold[i]) + uold[i]  + sigma_left(i * h) / rho * uold[i] - sigma_right(i * h) / rho * uold[i]; // bad
			unew[i] = mu1 * mu2 * 0.5 * (uold[i + 1] - 2 * uold[i] + uold[i - 1]) - mu1 * 0.5 * (pold[i + 1] - pold[i - 1]) + uold[i] + sigma_left(i * h) / rho * uold[i] - sigma_right(i * h) / rho * uold[i];// -sigma_left(i * h) / rho * uold[i] - sigma_right(i * h) / rho; // best
			//pnew[i] = -mu2 * (uold[i + 1] - uold[i]) + pold[i] -sigma_left(i * h) / rho * pold[i] - sigma_right(i * h) / rho * pold[i];
			pnew[i] = mu1 * mu2 * 0.5 * (pold[i + 1] - 2 * pold[i] + pold[i - 1]) - mu2 * 0.5 * (uold[i + 1] - uold[i - 1]) + pold[i] - sigma_left(i * h) * rho * c * c * pold[i] - sigma_right(i * h) * rho * c * c * pold[i]; // best

		}
		//unew[0] -= sigma_left(0) / rho * uold[0];
		//unew[n-1] -= sigma_left((n-1)*h) / rho * uold[n-1];

		//unew[0] = unew[1];
		//pnew[0] = pnew[1];
		//unew[n - 1] = unew[n - 2];
		//pnew[n - 1] = pnew[n - 2];

		//for (int i = 0; i < n; ++i) // с ГУ
		for (int i = d; i < n - d; ++i)
		{
			foutU << i * h << '\t' << unew[i] << std::endl;
			foutP << i * h << '\t' << pnew[i] << std::endl;
		}

		pold = pnew;
		uold = unew;

	}
}

void temptest::numerical_classic_with_k()
{
	//начальное условие
	for (int i = 0; i < n; ++i)
	{
		uold[i] = u0(i * h);
		pold[i] = p0(i * h);

		// для полуволн также
		uold_l[i] = u0(i * h);
		pold_l[i] = p0(i * h);

		uold_r[i] = u0(i * h);
		pold_r[i] = p0(i * h);
	}
	std::ofstream foutU0("ch_cl_kU0.txt");
	std::ofstream foutP0("ch_cl_kP0.txt");


	for (int i = 0; i < n; ++i)
	{
		foutU0 << i * h << '\t' << uold[i] << std::endl;
		foutP0 << i * h << '\t' << pold[i] << std::endl;
	}

	int k = 1;

	//h = l / (n - 1);
	rho = 583;
	c = 0.1456;
	a1 = 1e2; // пологость
	a2 = 1e2; // пологость
	alpha1 = 1; // амплитуда
	alpha2 = 1; // амплитуда

	//k1 = 0;
	//k2 = 0;
	for (int it = 1; it < T; ++it)
	{


		// цикл с одной большой волной в одну сторону 

#define SINGLE
#ifdef SINGLE
		//printf("Single\n");
#pragma omp parallel for
		for (int i = 1; i < n - 1; ++i) // края не трогаются
		{
			mu1 = tau / (rho * h);
			mu2 = rho * c * c * tau / h;

			unew[i] = mu1 * mu2 * 0.5 * (uold[i + 1] - 2 * uold[i] + uold[i - 1]) \
				- mu1 * 0.5 * (pold[i + 1] - pold[i - 1]) + uold[i] \
/*-k1u*/		- tau * k1 * uold[i] + mu1 * tau * 0.25 * (k1 + k2) * (pold[i + 1] - pold[i - 1]) + tau * tau * k1 * k1 * 0.5 * uold[i];																																																   //pnew[i] = -mu2 * (uold[i + 1] - uold[i]) + pold[i] -sigma_left(i * h) / rho * pold[i] - sigma_right(i * h) / rho * pold[i];
			pnew[i] = mu1 * mu2 * 0.5 * (pold[i + 1] - 2 * pold[i] + pold[i - 1]) \
				- mu2 * 0.5 * (uold[i + 1] - uold[i - 1]) + pold[i] \
/*-k2u*/        - tau * k2 * pold[i] + mu2 * tau * 0.25 * (k1 + k2) * (uold[i + 1] - uold[i - 1]) + tau * tau * k2 * k2 * 0.5 * pold[i];

		} // end for cycle

		//for (int i = 0; i < n; ++i) // с ГУ
		if (it % 100 == 0)
		{
			std::ofstream foutU("ch_cl_kU" + std::to_string(k) + ".txt");
			std::ofstream foutP("ch_cl_kP" + std::to_string(k++) + ".txt");
			for (int i = 0; i < n; ++i)
			{
				foutU << i * h << '\t' << unew[i] << std::endl;
				foutP << i * h << '\t' << pnew[i] << std::endl;
			}
		}
		pold = pnew;
		uold = unew;

#else
		// здесь две волны в разные стороны
		for (int i = 1; i < n - 1; ++i) // края не трогаются
		{
			// + для правой, - для левая
			mu1 = tau / (rho * h);
			mu2 = rho * c * c * tau / h;

			//инит для правой волны
			unew_r[i] = mu1 * mu2 * 0.5 * (uold_r[i + 1] - 2 * uold_r[i] + uold_r[i - 1]) \
				- mu1 * 0.5 * (pold_r[i + 1] - pold_r[i - 1]) + uold_r[i] \
				/*- tau * k1 * uold_r[i] + mu1 * tau * 0.25 * (k1 + k2) * (pold_r[i + 1] - pold_r[i - 1]) + tau * tau * k1 * k1 * 0.5 * uold_r[i]*/;

			pnew_r[i] = mu1 * mu2 * 0.5 * (pold_r[i + 1] - 2 * pold_r[i] + pold_r[i - 1]) \
				- mu2 * 0.5 * (uold_r[i + 1] - uold_r[i - 1]) + pold_r[i] \
				/*- tau * k2 * pold_r[i] + mu2 * tau * 0.25 * (k1 + k2) * (uold_r[i + 1] - uold_r[i - 1]) + tau * tau * k2 * k2 * 0.5 * pold_r[i]*/;
			
			//инит для левой волны
			mu1 = -mu1;
			mu2 = -mu2;

			unew_l[i] = mu1 * mu2 * 0.5 * (uold_l[i + 1] - 2 * uold_l[i] + uold_l[i - 1]) \
				- mu1 * 0.5 * (pold_l[i + 1] - pold_l[i - 1]) + uold_l[i] \
				/*- tau * k1 * uold_l[i] + mu1 * tau * 0.25 * (k1 + k2) * (pold_l[i + 1] - pold_l[i - 1]) + tau * tau * k1 * k1 * 0.5 * uold_l[i]*/;

			pnew_l[i] = mu1 * mu2 * 0.5 * (pold_l[i + 1] - 2 * pold_l[i] + pold_l[i - 1]) \
				- mu2 * 0.5 * (uold_l[i + 1] - uold_l[i - 1]) + pold_l[i] \
				/*- tau * k2 * pold_l[i] + mu2 * tau * 0.25 * (k1 + k2) * (uold_l[i + 1] - uold_l[i - 1]) + tau * tau * k2 * k2 * 0.5 * pold_l[i]*/;

			unew[i] = unew_r[i];//(unew_l[i] + unew_r[i]) / 2; 
			pnew[i] = pnew_r[i];//fabs(pnew_l[i] + pnew_r[i]) / 2;

		} // end for cycle


		for (int i = 0; i < n; ++i)
		{
			foutU << i * h << '\t' << unew[i] << std::endl;
			foutP << i * h << '\t' << pnew[i] << std::endl;
		}

		//для скорости левой и правой
		uold_l = unew_l;
		uold_r = unew_r;

		//для давления левого и правого
		pold_l = pnew_l;
		pold_r = pnew_r;
#endif
	} // end while time cycle

	max_wrong();
}

void temptest::numerical_classic_with_sigma()
{
	rho = 583;
	c = 0.1456;
	a1 = 1e-1; // пологость
	a2 = 1e-1; // пологость
	alpha1 = 1; // амплитуда
	alpha2 = 1; // амплитуда

	double P0divU0 = 4.28 * 1e4;
	double U0 = 3.5 * 1e3;
	double P0 = 1.5 * 1e8;
	double prc = 1.7 * 1e7;

	//начальное условие
	for (int i = 0; i < n; ++i)
	{
		uold[i] =  u0(i * h);
		pold[i] = P0divU0 * p0(i * h);

		// для полуволн также
		uold_l[i] = u0(i * h);
		pold_l[i] = p0(i * h);

		uold_r[i] = u0(i * h);
		pold_r[i] = p0(i * h);
	}
	std::ofstream foutU0("ch_cl_sigmaU0.txt");
	std::ofstream foutP0("ch_cl_sigmaP0.txt");


	for (int i = d; i < n - d; ++i)
	//for (int i = 0; i < n; ++i) // с ГУ
	{
		foutU0 << i * h << '\t' << uold[i] << std::endl;
		foutP0 << i * h << '\t' << pold[i] << std::endl;
	}

	int k = 1;
	std::ofstream up("up.txt");


	for (int it = 1; it <= T; ++it)
	{


		// цикл с одной большой волной в одну сторону 
		double c1 = 1 / rho;
		double c2 = rho * c * c;
#define SINGLE
#ifdef SINGLE
		//printf("Single\n");
//#pragma omp parallel for
		for (int i = 2; i < n - 2; ++i) // края не трогаются
		{
			mu1 = tau / (rho * h);
			mu2 = rho * c * c * tau / h;
			k1_r = sigma_right(i * h); // обновление коэффа
			k1_l = sigma_left(i * h);

			k2_r = sigma_right(i * h); // обновление коэффа
			k2_l = sigma_left(i * h);
			//k2, k2_l, k2_r = 0;

			double f_minus_u = std::min(1., ksi_j_minus(i, uold));
			f_minus_u = std::max(0., f_minus_u);

			double f_plus_u = std::min(1., ksi_j_plus(i, uold));
			f_plus_u = std::max(0., f_plus_u);

			double f_minus_p = std::min(1., ksi_j_minus(i, pold));
			f_minus_p = std::max(0., f_minus_p);

			double f_plus_p = std::min(1., ksi_j_plus(i, pold));
			f_plus_p = std::max(0., f_plus_p);

			/* Слаемое с производной сi * (sigma'2L + sigma'2R) * pold  для unew */
			/*sigma'2L, sigma'2R -- затухания для давления*/
//			unew[i] = mu1 * mu2 * 0.5 * h * (u_or_p_j_plus(i, uold) - u_or_p_j_minus(i, uold)) \
				- mu1 * 0.5 * h * (u_or_p_j_plus(i, pold) + u_or_p_j_minus(i, pold)) + uold[i] \
/*-k1u*/		- tau * uold[i] * (k1_l + k1_r) + mu1 * tau * 0.25 * (k1_l + k1_r + k2_l + k2_r) * (pold[i + 1] - pold[i - 1]) \
				+ tau * tau * 0.5 * (k1_l + k1_r) * (k1_l + k1_r) * uold[i] \
/*if k2!=0*/	- tau * tau * 0.5 * c1 * (sigma_right_dx(i * h) + sigma_left_dx(i * h))  * pold[i];																																																   //pnew[i] = -mu2 * (uold[i + 1] - uold[i]) + pold[i] -sigma_left(i * h) / rho * pold[i] - sigma_right(i * h) / rho * pold[i];
//			pnew[i] = mu1 * mu2 * 0.5 * h * (u_or_p_j_plus(i, pold) - u_or_p_j_minus(i, pold)) \
				- mu2 * 0.5 * h * (u_or_p_j_plus(i, uold) + u_or_p_j_minus(i, uold)) + pold[i] \
				- tau * pold[i] * (k2_l + k2_r) + mu2 * tau * 0.25 * (k1_l + k1_r + k2_l + k2_r) * (uold[i + 1] - uold[i - 1]) \
				+ tau * tau * 0.5 * (k2_l + k2_r) * (k2_l + k2_r) * pold[i] \
/*if k2!=0*/	- tau * tau * 0.5 * c2 * (sigma_right_dx(i * h) + sigma_left_dx(i * h))  * uold[i];

			double u1 = 1 - tau * (k1_l + k1_r) + tau * tau * 0.5 * (k1_l + k1_r) * (k1_l + k1_r); //un
			double u2 = mu1 * mu2 * 0.5; // uxxn
			double u3 = tau * 0.25 * mu1 * (k1_l + k1_r + k2_l + k2_r) - (mu1 * 0.5); // pi+1 - pi-1
			double u4 = tau * tau * 0.5 * c1 * (sigma_right_dx(i * h) + sigma_left_dx(i * h)); // pn
			//u3 = 0;
			//u1 - если убрать - вообще нет отрисовки (нельзя убирать)
			//u2 - оказывает хорошее влияние - наоборот сглаживает решение
			//u3 - если убрать - картинка стоит ()
			//u4 - если убрать - изменений нет от слова вообще

			unew[i] = u1 * uold[i] + u2 * (uold[i + 1] - 2 * uold[i] + uold[i - 1]) \
				+ u3 * (pold[i + 1] - pold[i - 1]) + u4 * pold[i];
			//	+ u3 * h * (u_or_p_j_plus(i, pold) * ksi_j_plus(i, pold) + u_or_p_j_minus(i, pold) * ksi_j_minus(i, pold) ) + u4 * pold[i]; // регуряризация
			//	+ u3 * h * (u_or_p_j_plus(i, pold) + u_or_p_j_minus(i, pold)) + u4 * pold[i];
			double p1 = 1 - tau * (k2_l + k2_r) + tau * tau * 0.5 * (k2_l + k2_r) * (k2_l + k2_r); //pn
			double p2 = mu1 * mu2 * 0.5; // pxxn
			double p3 = tau * 0.25 * mu2 * (k1_l + k1_r + k2_l + k2_r) - (mu2 * 0.5); // ui+1 - ui-1
			double p4 = tau * tau * 0.5 * c2 * (sigma_right_dx(i * h) + sigma_left_dx(i * h)); // un
			pnew[i] = p1 * pold[i] + p2 * (pold[i + 1] - 2 * pold[i] + pold[i - 1]) \
				+ p3 * (uold[i + 1] - uold[i - 1]) + p4 * uold[i];
			//	+ p3 * h * (u_or_p_j_plus(i, uold) * ksi_j_plus(i, uold) + u_or_p_j_minus(i, uold) * ksi_j_minus(i, uold)) + p4 * uold[i]; // регуряризация
			//	+ p3 * h * (u_or_p_j_plus(i, uold) + u_or_p_j_minus(i, uold)) + p4 * uold[i];
		
		} // end for cycle

		  //for (int i = 0; i < n; ++i) // с ГУ

		if (it % 100 == 0)
		{
			std::ofstream foutU("ch_cl_sigmaU" + std::to_string(k) + ".txt");
			std::ofstream foutP("ch_cl_sigmaP" + std::to_string(k++) + ".txt");

			for (int i = d; i < n - d; ++i) // без поглощающих слоев
			//for (int i = 0; i < n; ++i) // без ГУ
			{
				foutU << i * h << '\t' << unew[i] << std::endl;
				foutP << i * h << '\t' << pnew[i] << std::endl;
			}
		}
		pold = pnew;
		uold = unew;
#else
		// здесь две волны в разные стороны
		for (int i = 1; i < n - 1; ++i) // края не трогаются
		{
			// + для правой, - для левая
			mu1 = tau / (rho * h);
			mu2 = rho * c * c * tau / h;

			//инит для правой волны
			unew_r[i] = mu1 * mu2 * 0.5 * (uold_r[i + 1] - 2 * uold_r[i] + uold_r[i - 1]) \
				- mu1 * 0.5 * (pold_r[i + 1] - pold_r[i - 1]) + uold_r[i] \
				/*- tau * k1 * uold_r[i] + mu1 * tau * 0.25 * (k1 + k2) * (pold_r[i + 1] - pold_r[i - 1]) + tau * tau * k1 * k1 * 0.5 * uold_r[i]*/;

			pnew_r[i] = mu1 * mu2 * 0.5 * (pold_r[i + 1] - 2 * pold_r[i] + pold_r[i - 1]) \
				- mu2 * 0.5 * (uold_r[i + 1] - uold_r[i - 1]) + pold_r[i] \
				/*- tau * k2 * pold_r[i] + mu2 * tau * 0.25 * (k1 + k2) * (uold_r[i + 1] - uold_r[i - 1]) + tau * tau * k2 * k2 * 0.5 * pold_r[i]*/;

			//инит для левой волны
			mu1 = -mu1;
			mu2 = -mu2;

			unew_l[i] = mu1 * mu2 * 0.5 * (uold_l[i + 1] - 2 * uold_l[i] + uold_l[i - 1]) \
				- mu1 * 0.5 * (pold_l[i + 1] - pold_l[i - 1]) + uold_l[i] \
				/*- tau * k1 * uold_l[i] + mu1 * tau * 0.25 * (k1 + k2) * (pold_l[i + 1] - pold_l[i - 1]) + tau * tau * k1 * k1 * 0.5 * uold_l[i]*/;

			pnew_l[i] = mu1 * mu2 * 0.5 * (pold_l[i + 1] - 2 * pold_l[i] + pold_l[i - 1]) \
				- mu2 * 0.5 * (uold_l[i + 1] - uold_l[i - 1]) + pold_l[i] \
				/*- tau * k2 * pold_l[i] + mu2 * tau * 0.25 * (k1 + k2) * (uold_l[i + 1] - uold_l[i - 1]) + tau * tau * k2 * k2 * 0.5 * pold_l[i]*/;

			unew[i] = unew_r[i];//(unew_l[i] + unew_r[i]) / 2; 
			pnew[i] = pnew_r[i];//fabs(pnew_l[i] + pnew_r[i]) / 2;

		} // end for cycle


		for (int i = 0; i < n; ++i)
		{
			foutU << i * h << '\t' << unew[i] << std::endl;
			foutP << i * h << '\t' << pnew[i] << std::endl;
		}

		//для скорости левой и правой
		uold_l = unew_l;
		uold_r = unew_r;

		//для давления левого и правого
		pold_l = pnew_l;
		pold_r = pnew_r;
#endif
	} // end while time cycle
}

void temptest::one_equation()
{
	std::ofstream fout("1eq0.txt");
	std::ofstream fouta("1eqan0.txt");
	for (int i = 1; i < n - 1; ++i)
	{
		uold[i] = u0(i * h);

	}

	for (int i = 0; i < n; ++i)
	{
		fout << i * h << '\t' << uold[i] << std::endl;
		fouta << i * h << '\t' << uold[i] << std::endl;
	}

	mu1 = tau / (rho * h);
	mu2 = rho * c * c * tau / h;
	double a = 1 / rho;
	fout.close();
	fouta.close();
	int k = 1;
	for (int it = 1; it < T; ++it)
	{
		for (int i = 2; i < n - 2; ++i)
		{
			double f_minus = std::min(1., ksi_j_minus(i, uold));
			f_minus = std::max(0., f_minus);

			double f_plus = std::min(1., ksi_j_plus(i, uold));
			f_plus = std::max(0., f_plus);

			//unew[i] = -mu1 * 0.5 * (uold[i + 1] - uold[i - 1]) + (uold[i + 1] + uold[i - 1]) / 2; // гасится
			//unew[i] = mu1 * mu1 * 0.5 * (uold[i + 1] - 2 * uold[i] + uold[i - 1]) - mu1 * 0.5 * (uold[i + 1] - uold[i - 1]) + uold[i]; // Лакс-Вердорф
			unew[i] = mu1 * mu1 * 0.5 * (uold[i + 1] - 2 * uold[i] + uold[i - 1]) - mu1 * 0.5 * h * ( u_or_p_j_plus(i, uold) * f_plus + u_or_p_j_minus(i, uold) * f_minus) + uold[i]; // Лакс-Вердорф с кси

			//unew[i] = mu1 * mu1 * 0.5 * (uold[i] - 2 * uold[i - 1] + uold[i - 2]) - mu1 * 0.5 * (3 * uold[i] - 4 * uold[i - 1] + uold[i - 2]) + uold[i]; // tvd



			//unew[i] = uold[i] - a * tau * u_or_p_j_minus(i, uold) \
			//- a * tau * 0.5 * (1 - a * tau / h) * (f_plus * u_or_p_j_plus(i, uold) - f_minus * u_or_p_j_minus(i, uold));
		}


		if (it % 100 == 0)
		{
			std::ofstream foutU("1eq" + std::to_string(k) + ".txt");
			std::ofstream foutUa("1eqan" + std::to_string(k++) + ".txt");

			for (int i = 0; i < n; ++i) // без поглощающих слоев
											//for (int i = 1; i < n - 1; ++i) // без ГУ
			{
				foutU << i * h << '\t' << unew[i] << std::endl;
				foutUa << i * h << '\t' << u0(i * h - it * tau / rho) << std::endl;
			}
		}
		uold = unew;

	}
}

void temptest::newton_gu1(adept::Stack &stack)
{
	// н.у.
	for (int i = 0; i < n; ++i)
	{
		pold_ad[i] = p0(i * h);
		uold_ad[i] = u0(i * h);
	}
	std::ofstream foutU0("ch_newtonU0.txt");
	std::ofstream foutP0("ch_newtonP0.txt");

	for (int i = 0; i < n; ++i)
	{
		foutU0 << i * h << '\t' << uold_ad[i] << std::endl;
		foutP0 << i * h << '\t' << pold_ad[i] << std::endl;
	}
	// г.у
	uold_ad[0] = 0; // u0 = 0
	pold_ad[0] = 0; // p0 = 0

	uold_ad[n - 1] = 0; // un = 0
	pold_ad[n - 1] = 0; // pn = 0



	for (int i = 0; i < n; ++i)
	{
		pnew_ad[i] = p0(i * h) + 1e-5 * (i + 1); // начальное приближение корня
		unew_ad[i] = u0(i * h) + 1e-5 * (i + 1); // начальное приближение корня
	}


	for (int i = 1; i < T; ++i)
	{
		int iter = 1;
		std::ofstream foutU("ch_newtonU" + std::to_string(i) + ".txt");
		std::ofstream foutP("ch_newtonP" + std::to_string(i) + ".txt");
		encode_from_shorts_to_long(xnew_ad, unew_ad, pnew_ad); // unew_ad, pnew_ad в xnew_ad ПЕРЕД рекордингом
		while (true)
		{
			stack.new_recording();
			decode_from_long_to_shorts(xnew_ad, unew_ad, pnew_ad);
			// невязки 
			f1_ad[0] = unew_ad[0];
			f1_ad[n - 1] = unew_ad[n - 1];
			f2_ad[0] = pnew_ad[0];
			f2_ad[n - 1] = pnew_ad[n - 1];

			residual1(unew_ad, uold_ad, pnew_ad, pold_ad); // инит невязок (МОГУТ БЫТЬ ГУ ПОГЛОЩЕНИЯ, ПРОВЕРЬ!)
			residual2(unew_ad, uold_ad, pnew_ad, pold_ad); // инит невязок

			encode_from_shorts_to_long(f_ad, f1_ad, f2_ad); // f1_ad, f2_ad в f

			getJacobian(stack, jacobian, 2 * n); // cоставление якоби из f и xnew_ad
			findRoot(2 * n); // заполнить xsol
			decode_from_long_to_shorts(xsol_ad, unew_ad, pnew_ad); // заполнить unew_ad, pnew_ad из xsol_ad
			if (checkresidual1() && checkresidual2())
			{
				printf("Solution with %d iter\n", iter);

				for (int i = 0; i < n; ++i) // c ГУ
				//for (int i = 1; i < n - 1; ++i) // без ГУ
				{
					foutU << i * h << '\t' << unew_ad[i] << std::endl;
					foutP << i * h << '\t' << pnew_ad[i] << std::endl;
				}

				//обновить решение
				uold_ad = unew_ad;
				pold_ad = pnew_ad;

				break;
			}
			else
			{
				++iter;
				printf("iter: %d\n", iter);
				xnew_ad = xsol_ad; // сместить решение
			}
		}
	}
}

void temptest::newton_gu3(adept::Stack & stack)
{
	// н.у.
	for (int i = 0; i < n; ++i)
	{
		pold_ad[i] = p0(i * h);
		uold_ad[i] = u0(i * h);
	}
	std::ofstream foutU0("ch_newton_gu3U0.txt");
	std::ofstream foutP0("ch_newton_gu3P0.txt");

	for (int i = 0; i < n; ++i)
	{
		foutU0 << i * h << '\t' << uold_ad[i] << std::endl;
		foutP0 << i * h << '\t' << pold_ad[i] << std::endl;
	}
	// г.у
	//uold_ad[0] = 0; // u0 = 0

	pold_ad[0] = 1; // p0 = 0

	//uold_ad[n - 1] = 0; // un = 0

	pold_ad[n - 1] = 1; // pn = 0


	for (int i = 0; i < n; ++i)
	{
		pnew_ad[i] = p0(i * h) + 1e-5 * (i + 1); // начальное приближение корня
		unew_ad[i] = u0(i * h) + 1e-5 * (i + 1); // начальное приближение корня
	}


	for (int i = 1; i < T; ++i)
	{
		int iter = 1;
		std::ofstream foutU("ch_newton_gu3U" + std::to_string(i) + ".txt");
		std::ofstream foutP("ch_newton_gu3P" + std::to_string(i) + ".txt");
		encode_from_shorts_to_long(xnew_ad, unew_ad, pnew_ad); // unew_ad, pnew_ad в xnew_ad ПЕРЕД рекордингом
		while (true)
		{
			stack.new_recording();
			decode_from_long_to_shorts(xnew_ad, unew_ad, pnew_ad);
			// невязки 

			//f1_ad[0] = unew_ad[0];
			//f1_ad[n - 1] = unew_ad[n - 1];
			residual1_left(unew_ad, uold_ad, pnew_ad, pold_ad); // затухание слева
			residual1_right(unew_ad, uold_ad, pnew_ad, pold_ad); // затухание справа
			residual1_gu3(unew_ad, uold_ad, pnew_ad, pold_ad); // здесь затухания отсутствуют
			// ГУ для давления
			f2_ad[0] = pnew_ad[0] - 1;
			f2_ad[n - 1] = pnew_ad[n - 1] - 1;
			residual2_gu3(unew_ad, uold_ad, pnew_ad, pold_ad); // здесь все кроме граничек для давления


			encode_from_shorts_to_long(f_ad, f1_ad, f2_ad); // f1_ad, f2_ad в f

			getJacobian(stack, jacobian, 2 * n); // cоставление якоби из f и xnew_ad
			findRoot(2 * n); // заполнить xsol
			decode_from_long_to_shorts(xsol_ad, unew_ad, pnew_ad); // заполнить unew_ad, pnew_ad из xsol_ad
			if (checkresidual1_gu3() && checkresidual2_gu3())
			{
				printf("Solution with %d iter\n", iter);

				for (int i = 0; i < n; ++i)
				{
					foutU << i * h << '\t' << unew_ad[i] << std::endl;
					foutP << i * h << '\t' << pnew_ad[i] << std::endl;
				}

				//обновить решение
				uold_ad = unew_ad;
				pold_ad = pnew_ad;

				break;
			}
			else
			{
				++iter;
				//decode_from_long_to_shorts();// из xsol_ad в unew_ad, pnew_ad
				printf("iter: %d\n", iter);
				xnew_ad = xsol_ad; // сместить решение
			}
		}
	}
}

void temptest::residual1(const std::vector<adouble> &unew_ad, const std::vector<adouble> &uold_ad, const std::vector<adouble> &pnew_ad, const std::vector<adouble> &pold_ad)
{
	for (int i = 1; i < n - 1; ++i)
	{
		f1_ad[i] = (unew_ad[i] - uold_ad[i]) / tau + (pnew_ad[i + 1] - pnew[i]) / (rho * h); \
		//	+ sigma_left(i * h) / rho * unew_ad[i] - sigma_right(i * h) / rho * unew_ad[i];
		//f1_ad[i] = unew_ad[i] - i;
	}
}

void temptest::residual2(const std::vector<adouble> &unew_ad, const std::vector<adouble> &uold_ad, const std::vector<adouble> &pnew_ad, const std::vector<adouble> &pold_ad)
{
	for (int i = 1; i < n - 1; ++i)
	{
		f2_ad[i] = (pnew_ad[i] - pold_ad[i]) / tau + (unew_ad[i + 1] - unew[i]) * (rho * c * c) / h; \
		//	+ sigma_left(i * h) * (rho * c * c) * pnew_ad[i] - sigma_right(i * h) * (rho * c * c) * pnew_ad[i];
		//f2_ad[i] = pnew_ad[i] - (i + pnew_ad.size());
	}
}

void temptest::residual1_left(const std::vector<adouble>& unew_ad, const std::vector<adouble>& uold_ad, const std::vector<adouble>& pnew_ad, const std::vector<adouble>& pold_ad)
{
	for (int i = 0; i < d; ++i)
	{
		f1_ad[i] = (unew_ad[i] - uold_ad[i]) / tau + (pnew_ad[i + 1] - pnew[i]) / (rho * h) + sigma_left(i * h) * unew_ad[i] / rho;
	}
}

// правое затухание
void temptest::residual1_right(const std::vector<adouble>& unew_ad, const std::vector<adouble>& uold_ad, const std::vector<adouble>& pnew_ad, const std::vector<adouble>& pold_ad)
{
	f1_ad[n - 1] = unew_ad[n - 1] - unew_ad[n - 2];
	for (int i = n - d; i < n - 1; ++i)
	{
		f1_ad[i] = (unew_ad[i] - uold_ad[i]) / tau + (pnew_ad[i + 1] - pnew[i]) / (rho * h) + sigma_right(i * h) * unew_ad[i] / rho;
	}

}

// слева и справа отсекаем для поглощения
void temptest::residual1_gu3(const std::vector<adouble>& unew_ad, const std::vector<adouble>& uold_ad, const std::vector<adouble>& pnew_ad, const std::vector<adouble>& pold_ad)
{
	for (int i = d; i < n - d; ++i)
	{
		f1_ad[i] = (unew_ad[i] - uold_ad[i]) / tau + (pnew_ad[i + 1] - pnew[i]) / (rho * h);
	}
}

//берем все, кроме граничек
void temptest::residual2_gu3(const std::vector<adouble>& unew_ad, const std::vector<adouble>& uold_ad, const std::vector<adouble>& pnew_ad, const std::vector<adouble>& pold_ad)
{
	for (int i = 1; i < n - 1; ++i)
	{
		f2_ad[i] = (pnew_ad[i] - pold_ad[i]) / tau + (unew_ad[i + 1] - unew[i]) * (rho * c * c) / h;
	}
}

bool temptest::checkresidual1_gu3()
{
	residual1_gu3(unew_ad, uold_ad, pnew_ad, pold_ad);
	for (int i = d; i < n - d; ++i)
	{
		if (fabs(f1_ad[i]) > 1e-8)
			return false;
	}
	return true;
}

bool temptest::checkresidual2_gu3()
{
	residual2_gu3(unew_ad, uold_ad, pnew_ad, pold_ad);
	for (int i = 1; i < n - 1; ++i)
	{
		if (fabs(f2_ad[i]) > 1e-8)
			return false;
	}
	return true;
}

bool temptest::checkresidual1()
{
	residual1(unew_ad, uold_ad, pnew_ad, pold_ad);
	for (int i = 1; i < n - 1; ++i)
	{
		if (fabs(f1_ad[i]) > 1e-8)
			return false;
	}
	return true;
}

bool temptest::checkresidual2()
{
	residual2(unew_ad, uold_ad, pnew_ad, pold_ad);
	for (int i = 1; i < n - 1; ++i)
	{
		if (fabs(f2_ad[i]) > 1e-8)
			return false;
	}
	return true;
}


//из длинного в короткие
void temptest::decode_from_long_to_shorts(const std::vector<adouble> &longvec, std::vector<adouble> &short1, std::vector<adouble> &short2)
{
	// есть xsol длинный вектор
	// есть unew_ad, pnew_ad - короткие вектора
	for (int i = 0; i < short1.size(); ++i)
	{
		short1[i] = longvec[i]; // unew_ad заполнение
	}
	int k = 0;
	for (int i = short1.size(); i < longvec.size(); ++i)
	{
		short2[k++] = longvec[i]; // pnew_ad заполнение
	}
}

// из коротких сделать длинный
void temptest::encode_from_shorts_to_long(std::vector<adouble> &longvec, const std::vector<adouble> &short1,const std::vector<adouble> &short2)
{
	//первая половинка
	for (int i = 0; i < short1.size(); ++i)
	{
		longvec[i] = short1[i];
	}

	//вторая половинка
	int k = 0;
	for (int i = short1.size(); i < longvec.size(); ++i)
	{
		longvec[i] = short2[k++];
	}
}
// указывем какой вектор участвует для вычисления (j не могут быть крайними)
double temptest::rj(int j, const std::vector<double>& f)
{
	return (f[j] - f[j - 1]) / (f[j + 1] - f[j]);
}

// указывем какой вектор участвует для вычисления
double temptest::Cj(int j, std::vector<double>& r)
{
	std::vector<double> mins = { 2 * r[j], 0.5 * (1 + r[j]), 2 };
	double tmin = *std::min_element(mins.begin(), mins.end());
	return std::max(0., tmin);
}

double temptest::sigma_right(double x)
{

	if (x >= (n - d - 1) * h && x <= l)
	{
		double res = (m + 1) * W * log(10) * pow((x - ((n - d - 1) * h)), m) / pow(d * h, m + 1);
		return res;
	}
	else
		return 0;
}

double temptest::sigma_left(double x)
{

	if (x >= 0 && x <= d * h)
		return (m + 1) * W * log(10) * pow(d * h - x, m) / pow(d * h, m + 1);
	else
		return 0;
}

double temptest::sigma_right_dx(double x)
{
	double eps = 1e-6;
	if (x >= (n - d - 1) * h && x <= l)
	{
		double num = (sigma_right(x + eps) - sigma_right(x)) / eps;
		double man = (m + 1) * m * W * log(10) * pow((x - ((n - d - 1) * h)), m - 1) / pow(d * h, m + 1);
		//right << num << '\t' << man << std::endl;
		return man;
	}
	return 0;
}

double temptest::sigma_left_dx(double x)
{
	double eps = 1e-6;
	if (x >= 0 && x <= d * h)
	{
		double num = (sigma_left(x + eps) - sigma_left(x)) / eps;
		double man = (m + 1) * m * W * log(10) * pow(d * h - x, m - 1) / pow(d * h, m + 1);
		//left << num << '\t' << man << std::endl;
		return man;
	}
	return 0;
}

double temptest::ksi_j_plus(int j, std::vector<double> &vec)
{
	double chisl = vec[j] - vec[j - 1];
	double znam = vec[j + 1] - vec[j];
	if (znam != 0)
	{
		return chisl / znam;
	}

	return 0;
}

double temptest::ksi_j_minus(int j, std::vector<double> &vec)
{
	double chisl = vec[j - 1] - vec[j - 2];
	double znam = vec[j] - vec[j - 1];
	if (znam != 0)
	{
		return chisl / znam;
	}

	return 0;
}

double temptest::u_or_p_j_plus(int j, std::vector<double>& vec)
{
	return (vec[j + 1] - vec[j]) / h;
}
double temptest::u_or_p_j_minus(int j, std::vector<double>& vec)
{
	return (vec[j] - vec[j - 1]) / h;
}

void temptest::plot_sigma()
{
	//double m_i = 0.1;
	double W_i = 11.2;
	for (int i = 0; i <= d; ++i)
	{
		std::ofstream fout("left_sigma" + std::to_string(i) + ".txt");
		//while (m_i <= 10)
		while (W_i <= 15)
		{
			//double func = (m_i + 1) * W * log(10) * pow(i * h, m_i) / pow(d * h, m_i + 1);
			double func = (m + 1) * W_i * log(10) * pow(i * h, m) / pow(d * h, m + 1);
			//fout << m_i << '\t' << func << std::endl;
			fout << W_i << '\t' << func << std::endl;
			//m_i += 0.1;
			W_i += 0.1;
		}
		//m_i = 0.1;
		W_i = 0.1;
	}
}

void temptest::max_wrong()
{
	for (int i = 6; i < n - 6; ++i)
	{
		//if (unew[i] < 0)
		//	continue;
		if (fabs(unew[i] - uanal[i]) > max)
		{
			max = fabs(uanal[i] - unew[i]);
			//std::cout << unew[i] << '\t' << uanal[i] << std::endl;
			std::cout << "abc" << '\t' << max << std::endl;
			std::cout << "rel" << '\t' << max / uanal[i] << std::endl;
		}
	}

}

void temptest::check_tvd(int j, const std::vector<double> &f, const std::vector<double> &C)
{
	double frac = (f[j] - f[j - 1]) / (f[j + 1] - f[j]);

	if (frac >= 0)
	{
		if (!(C[j] >= 0 && C[j] <= 2))
			printf("Bad Cj. It's %f instead of [0;2]\n", C[j]);
	}
	else
	{
		if (!(C[j] >= 0 && C[j] <= 0.01))
			printf("Bad Cj. It's %f instead of [0;0.01]\n", C[j]);
	}
}

temptest::~temptest()
{
}

void temptest::getJacobian(adept::Stack &stack, std::vector<std::vector<double>> &Jacobian, int n)
{

	std::vector <double> jac(n*n);   // вектор для получения матрицы Якоби
	for (int i = 0; i < n; i++)
		f_ad[i].set_gradient(1.0); // определение вектора в качестве целевой функции

	stack.compute_adjoint();  // запуск записи в стэк
	stack.independent(&xnew_ad[0], n); // определение вектора независимых переменных
	stack.dependent(&f_ad[0], n);  // определение вектора невязок
	stack.jacobian(&jac[0]); // вычисление матрицы Якоби
	int k = 0;

	// запись матрицы Якоби в двумерный вектор из одномермного 
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; ++j) {
			Jacobian[j][i] = jac[k++];
		}
	}
}

void temptest::findRoot(int n) 
{
	//printf("Det: %f\n", det(jacobian, n));
	if (abs(det(jacobian, n)) < 1e-5 || isnan(det(jacobian, n))) {
		printf("Singular matrix\n");
		Print2D(jacobian, n);
		return;
	}

	std::vector<std::vector<double>> invjac(n);
	for (int i = 0; i < n; i++) {
		invjac[i].resize(n);
	}

	Inversion(jacobian, invjac, n);
	//Print2D(invjac, n);
	for (int i = 0; i < n; i++) {
		adouble sum = 0;
		for (int j = 0; j < n; j++) {
			sum += invjac[i][j] * f_ad[j];
		}
		xsol_ad[i] = xnew_ad[i] - sum;
	}

}