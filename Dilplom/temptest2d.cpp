#include "temptest2d.h"
#include <chrono>
#include <omp.h>

#include <iomanip>
#include <iostream>
temptest2d::temptest2d()
{
}

temptest2d::~temptest2d()
{

}

temptest2d::temptest2d(int n, int T, double tau)
{
	// размеры области 100 см х 400 см
	// длина трещины 3 см
	// 15% - концентрация трещин
	// lw 25 см (в 8 раз больше длины трещины
	// d1 25 см (= lw)
	// d2 50 см (= 2 *lw)


	K = 2.25;
	rho = 1;
	d = (n - 1) * 0.1;// 50; // узлы для граничек
	
	//отжиг
	//fp
	lw = d / 2.5; // длина для высчитывания fp
	Nw = (n - 2 * d) / lw; // кол-во блоков в одном направлении 

	//fd
	d1 = 2 * lw; // fd
	d2 = 4 * lw; //fd
	T_pyro = 1e-3;
	//D2 = -5;
	//
	a = 0;
	l = 1;
	this->n = n;
	this->T = T;
	h = (l - a) / (n - 1);
	this->tau = tau;

	uold.resize(n);
	unew.resize(n);
	pold.resize(n);
	pnew.resize(n);
	vold.resize(n);
	vnew.resize(n);

	fractures.resize(n);
	fractures1.resize(n);
	fractures2.resize(n);
	fracture_centers.resize(n);
	fracture_centers1.resize(n);
	fracture_centers2.resize(n);
	for (int i = 0; i < n; ++i)
	{
		uold[i].resize(n);
		unew[i].resize(n);
		pold[i].resize(n);
		pnew[i].resize(n);
		vold[i].resize(n);
		vnew[i].resize(n);
		fractures[i].resize(n);
		fractures1[i].resize(n);
		fractures2[i].resize(n);
		fracture_centers[i].resize(n);
		fracture_centers1[i].resize(n);
		fracture_centers2[i].resize(n);
	}

	//параметры начальной волны
	a1_p = 0.5 * 1e-2; // пологость
	a2_p = 0.5 * 1e-2; // пологость

	a1_u = 0.5 * 1e-1; // пологость
	a2_u = 0.5 * 1e-1; // пологость

	a1_v = 0.5 * 1e-1; // пологость
	a2_v = 0.5 * 1e-1; // пологость

	alpha_u = 1e-4; // амплитуда
	alpha_v = 1e-4; // амплитуда
	alpha_p = 1e-2; // амплитуда

	// параметры затухания
	m = 1.5;
	W = 14;

	this->lod10 = (m + 1) *  W * log(10);

}

void temptest2d::numerical_classic()
{
	//начальное условие
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			uold[i][j] = u0(i * h, j * h);
			pold[i][j] = p0(i * h, j * h);
			vold[i][j] = v0(i * h, j * h);
		}

	}
	std::ofstream foutU0("ch_2d_U0.txt");
	std::ofstream foutP0("ch_2d_P0.txt");
	std::ofstream foutV0("ch_2d_V0.txt");

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			foutU0 << uold[i][j] << '\t';
			foutP0 << pold[i][j] << '\t';
			foutV0 << vold[i][j] << '\t';
		}
		foutU0 << std::endl;
		foutP0 << std::endl;
		foutV0 << std::endl;
	}

	int k = 1;
	std::ofstream up("up.txt");
	for (int it = 1; it <= T; ++it)
	{


		// цикл с одной большой волной в одну сторону 
		double p1 = tau * K / (2 * h);
		double uv1 = tau / (2 * h * rho);
		double uvp2 = (tau * tau * K) / (2 * rho * h * h);
		//printf("Single\n");

		//for (int i = 0; i < n; ++i) // с ГУ
		for (int i = 1; i < n - 1; ++i) // края не трогаются
		{
			for (int j = 1; j < n - 1; ++j)
			{
				pnew[i][j] = pold[i][j] - p1 * (uold[i + 1][j] - uold[i - 1][j] + vold[i][j + 1] - vold[i][j - 1]) \
					+ uvp2 * (pold[i + 1][j] - 4 * pold[i][j] + pold[i - 1][j] + pold[i][j + 1] + pold[i][j - 1]);

				unew[i][j] = uold[i][j] - uv1 * (pold[i + 1][j] - pold[i - 1][j])\
					+ uvp2 * (uold[i + 1][j] - 2 * uold[i][j] + uold[i - 1][j]);

				vnew[i][j] = vold[i][j] - uv1 * (pold[i][j + 1] - pold[i][j - 1])\
					+ uvp2 * (vold[i][j + 1] - 2 * vold[i][j] + uold[i][j - 1]);
			}
		} // end for cycle

		if (it % 100 == 0)
		{
			std::ofstream foutU("ch_2d_U" + std::to_string(k) + ".txt");
			std::ofstream foutP("ch_2d_P" + std::to_string(k) + ".txt");
			std::ofstream foutV("ch_2d_V" + std::to_string(k++) + ".txt");

			//for (int i = d; i < n - d; ++i) // без поглощающих слоев
			for (int i = 0; i < n; ++i) // без ГУ
			{
				for (int j = 0; j < n; ++j)
				{
					foutU << unew[i][j] << '\t';
					foutP << pnew[i][j] << '\t';
					foutV << vnew[i][j] << '\t';
				}
				foutU << std::endl;
				foutP << std::endl;
				foutV << std::endl;
			}
		}
		pold = pnew;
		uold = unew;
		vold = vnew;
	} // end while time cycle
}

void temptest2d::numerical_classic_with_sigma() 
{
	//начальное условие

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{ 
			uold[i][j] = u0(i * h, j * h);
			pold[i][j] = p0(i * h, j * h);
			vold[i][j] = v0(i * h, j * h);
			/*if (i >= 3 * d && i <= n - 3 * d && j >= 3 * d && j <= n - 3 * d)
			{
				pold[i][j] += alpha_p;
			}*/
			/*else
				pold[i][j] += alpha_p / 2.;*/
		}
		
	}

	std::ofstream foutU0("ch_2d_sigma_U0.txt");
	std::ofstream foutP0("ch_2d_sigma_P0.txt");
	std::ofstream foutV0("ch_2d_sigma_V0.txt");

	for (int i = d; i < n - d; ++i)
	{
		for (int j = d; j < n - d; ++j)
		{
			foutU0 << uold[i][j] << '\t';
			foutP0 << pold[i][j] << '\t';
			foutV0 << vold[i][j] << '\t';
		}
		foutU0 << std::endl;
		foutP0 << std::endl;
		foutV0 << std::endl;
	}

	int k = 1;
	std::ofstream up("up.txt");

	auto begin = std::chrono::system_clock::now();

	std::vector<std::vector<double>> sigmaSum(n, std::vector<double>(n));
	//std::vector<double> sum2jArray(n);
	std::vector<double> sum2ijArray(n);
#pragma omp parallel for
	for (int i = 1; i < n - 1; ++i) // края не трогаются
	{
		for (int j = 1; j < n - 1; ++j)
		{
			sigmaSum[i][j] = sigma_summ(i * h, j * h);
		}

		sum2ijArray[i] = sigma_left_down_der(i * h) + sigma_right_up_der(i * h);
		//sum2jArray[i] = sigma_left_down_der(i * h) + sigma_right_up_der(i * h);
	}


	for (int it = 1; it <= T; ++it)
	{


		// цикл с одной большой волной в одну сторону 
		double p1 = tau * K / (2 * h);
		double uv1 = tau / (2 * h * rho);
		double uvp2 = (tau * tau * K) / (2 * rho * h * h);
		//printf("Single\n");

		//for (int i = 0; i < n; ++i) // с ГУ

		
		

#pragma omp parallel for
		for (int i = 1; i < n - 1; ++i) // края не трогаются
		{
			for (int j = 1; j < n - 1; ++j)
			{
#if 0
				auto sum2j = sigma_left_down_der(j * h) + sigma_right_up_der(j * h);
				auto sum2i = sigma_left_down_der(i * h) + sigma_right_up_der(i * h);
#else
				auto sum2j = sum2ijArray[j]; // сумма производных сигм
				auto sum2i = sum2ijArray[i]; // сумма производных сигм
#endif

#if 0
				auto sigSum = sigma_summ(i * h, j * h);
#else
				auto sigSum = sigmaSum[i][j]; // sigma_summ(x,y)
#endif
				auto poldIJ = pold[i][j];
				auto uoldIJ = uold[i][j];
				auto voidIJ = vold[i][j];
				auto sig2Sum = sigSum * sigSum; // (sigma_summ(x,y))^2

				auto kTau = 0.5 * tau * tau;

				pnew[i][j] = poldIJ - tau * (K * (vec_x(i, j, uold) + vec_y(i, j, vold)) + sigSum * poldIJ) \
					+ kTau *
					(
						K / rho * (vec_xx(i, j, pold) + vec_yy(i, j, pold)) + K * vec_x(i, j, uold) * 2 * sigSum \
						+ K * vec_y(i, j, vold) * 2 * sigSum \
						+ K * uoldIJ * (sum2i) \
						+ K * voidIJ * (sum2j) \
						+ sig2Sum * poldIJ
					);

				unew[i][j] = uoldIJ - tau * (vec_x(i, j, pold) / rho + sigSum * uoldIJ) \
					+ kTau *
					(
						K / rho * vec_xx(i, j, uold) + vec_x(i, j, pold) / rho * 2 * sigSum \
						+ (sum2i) * poldIJ / rho \
						+ sig2Sum * uoldIJ
					);

				vnew[i][j] = voidIJ - tau * (vec_y(i, j, pold) / rho + sigSum * voidIJ) \
					+ kTau *
					(
						K / rho * vec_yy(i, j, vold) + vec_y(i, j, pold) / rho * 2 * sigSum \
						+ (sum2j) * poldIJ / rho \
						+ sig2Sum * voidIJ
					);
			}
		} // end for cycle
		
#if 1
		if (it % 100 == 0)
		{
			std::ofstream foutU("ch_2d_sigma_U" + std::to_string(k) + ".txt");
			std::ofstream foutP("ch_2d_sigma_P" + std::to_string(k) + ".txt");
			std::ofstream foutV("ch_2d_sigma_V" + std::to_string(k++) + ".txt");

			for (int i = d; i < int(n - d); ++i) // без поглощающих слоев
			//for (int i = 0; i < n; ++i) // с ГУ
			{
				for (int j = d; j < n - d; ++j) // без поглощающих слоев
				//for (int j = 0; j < n; ++j) // с ГУ
				{
					foutU << unew[i][j] << '\t';
					foutP << pnew[i][j] << '\t';
					foutV << vnew[i][j] << '\t';
				}
				foutU << std::endl;
				foutP << std::endl;
				foutV << std::endl;
			}
		}
#endif
		pold = pnew;
		uold = unew;
		vold = vnew;
	} // end while time cycle
	auto end = std::chrono::system_clock::now();
	std::cout << "time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;
}

void temptest2d::norm_with_sigma()
{
	//начальное условие

	//параметры начальной волны
	a1_p = 0.5 * 1e2; // пологость
	a2_p = 0.5 * 1e2; // пологость

	a1_u = 0.5 * 1e2; // пологость
	a2_u = 0.5 * 1e2; // пологость

	a1_v = 0.5 * 1e2; // пологость
	a2_v = 0.5 * 1e2; // пологость

	alpha_u = 1; // амплитуда
	alpha_v = 1; // амплитуда
	alpha_p = 1; // амплитуда

	double P0divU0 = 4.28 * 1e4;
	double U0 = 3.5 * 1e3;
	double P0 = 1.5 * 1e8;
	double prc = 1.7 * 1e7;
	double eps = 1e-12;
	//init_fractures();
	//perkolation(2, (n - 1) / 10);
	
	std::ifstream fin("fracture.txt");

	for (int i = d; i < n - d; ++i)
	{
		for (int j = d; j < n - d; ++j)
		{
			fin >> fractures[i][j];
		}
	}
//#pragma omp parallel for
	enum MyEnum
	{
		No,
		Up,
		Left,
		Bottom,
		Right,
		Center
	};
	int string_for_bottom = -1;
	for (int i = 1; i < n - 1; ++i)
	{
		for (int j = 1; j < n - 1; ++j)
		{

			if ( (i != 0 && j != 0 && i != n - 1 && j != n - 1) && (fractures[i][j] == No && (fractures[i + 1][j] == Up || fractures[i - 1][j] == Up || fractures[i][j + 1] == Up || fractures[i][j - 1]) ) )
			{
				continue;
			}
			if (fractures[i][j] == No) // нет трещины обычное решение
			{
				// сверху
				/*if (i < n - 2 && fractures[i + 1][j] == No && fractures[i + 2][j] == Up)
				{
					uold[i + 1][j] = uold[i][j];
					pold[i + 1][j] = pold[i][j];
					vold[i + 1][j] = vold[i][j];
				}
				//слева
				else if (j < n - 2 && fractures[i][j + 1] == No && fractures[i][j + 2] == Up)
				{
					uold[i][j + 1] = uold[i][j];
					pold[i][j + 1] = pold[i][j];
					vold[i][j + 1] = vold[i][j];
				}
				else*/
				{
					uold[i][j] = 0;// u0(i* h, j* h);
					pold[i][j] = P0divU0 * p0(i * h, j * h)/* + prc / P0*/;
					vold[i][j] = 0;// v0(i* h, j* h);
				}
				/*if (fractures[i][j - 1] == Right) // предыдущий столбец трещина
				{
					uold[i][j - 1] = uold[i][j];
					pold[i][j - 1] = pold[i][j];
					vold[i][j - 1] = vold[i][j];
				}
				else if (fractures[i - 1][j] == Bottom) // предыдущая строка с трещинами
				{
					uold[i - 1][j] = uold[i][j];
					pold[i - 1][j] = pold[i][j];
					vold[i - 1][j] = vold[i][j];
				}*/
				
			}
			/*else if (fractures[i][j] == Up)
			{
				uold[i][j] = uold[i - 1][j];
				pold[i][j] = pold[i - 1][j];
				vold[i][j] = vold[i - 1][j];
			}
			else if (fractures[i][j] == Left)
			{
				uold[i][j] = uold[i][j - 1];
				pold[i][j] = pold[i][j - 1];
				vold[i][j] = vold[i][j - 1];
			}
			else if (fractures[i][j] == Center)
			{
				uold[i][j] = std::nan("1");
				pold[i][j] = std::nan("1");
				vold[i][j] = std::nan("1");
			}*/
			else if (fractures[i][j] == Up)
			{
				//cверху
				/*if (i > 2 && fractures[i - 1][j] == No && fractures[i - 2][j] == No)
				{
					uold[i - 1][j] = uold[i - 2][j];
					pold[i - 1][j] = pold[i - 2][j];
					vold[i - 1][j] = vold[i - 2][j];
				}
				//слева
				if (j > 2 && fractures[i][j - 1] == No && fractures[i][j - 2] == No)
				{
					uold[i][j - 1] = uold[i][j - 2];
					pold[i][j - 1] = pold[i][j - 2];
					vold[i][j - 1] = vold[i][j - 2];
				}*/
				/*uold[i][j] = std::nan("1");
				pold[i][j] = std::nan("1");
				vold[i][j] = std::nan("1");*/
				uold[i][j] = eps;
				pold[i][j] = eps;
				vold[i][j] = eps;
			}

			if (j > 2 && fractures[i][j - 2] == Up && fractures[i][j] == No && fractures[i][j - 1] == No)
			{
				//справа
				uold[i][j - 1] = uold[i][j];
				pold[i][j - 1] = pold[i][j];
				vold[i][j - 1] = vold[i][j];
			}
			if (i > 2 && fractures[i - 2][j] == Up && fractures[i][j] == No && fractures[i - 1][j] == No)
			{
				//снизу
				uold[i - 1][j] = uold[i][j];
				pold[i - 1][j] = pold[i][j];
				vold[i - 1][j] = vold[i][j];
			}

			//сверху
			if (i < n - 2 && fractures[i + 1][j] == No && fractures[i + 2][j] == Up)
			{
				uold[i + 1][j] = uold[i][j];
				pold[i + 1][j] = pold[i][j];
				vold[i + 1][j] = vold[i][j];
			}
			//слева
			if (j < n - 2 && fractures[i][j + 1] == No && fractures[i][j + 2] == Up)
			{
				uold[i][j + 1] = uold[i][j];
				pold[i][j + 1] = pold[i][j];
				vold[i][j + 1] = vold[i][j];
			}
		}

	}

	std::ofstream foutU0("ch_norm_2d_sigma_U0.txt");
	std::ofstream foutP0("ch_norm_2d_sigma_P0.txt");
	std::ofstream foutV0("ch_norm_2d_sigma_V0.txt");


	//for (int i = 2 * d; i < n - int(2 * d); ++i)
	for (int i = d; i < n - d; ++i)
	{
		//for (int j = 2 * d; j < n - int(2 * d); ++j)
		for (int j = d; j < n - d; ++j)
		{
			/*foutU0 << std::setw(20) << uold[i][j] << std::setw(20);
			foutP0 << std::setw(20) << pold[i][j] << std::setw(20);
			foutV0 << std::setw(20) << vold[i][j] << std::setw(20);*/

			if (fabs(uold[i][j] - eps) < 1e-8)
				foutU0 << std::setw(20) << std::nan("1") << std::setw(20);
			else
				foutU0 << std::setw(20) << uold[i][j] << std::setw(20);

			if (fabs(pold[i][j] - eps) < 1e-8)
				foutP0 << std::setw(20) << std::nan("1") << std::setw(20);
			else
				foutP0 << std::setw(20) << pold[i][j] << std::setw(20);

			if (fabs(vold[i][j] - eps) < 1e-8)
				foutV0 << std::setw(20) << std::nan("1") << std::setw(20);
			else
				foutV0 << std::setw(20) << vold[i][j] << std::setw(20);
		}
		foutU0 << std::endl;
		foutP0 << std::endl;
		foutV0 << std::endl;
	}

	int k = 1;
	std::ofstream up("up.txt");

	auto begin = std::chrono::system_clock::now();

	std::vector<std::vector<double>> sigmaSum(n, std::vector<double>(n));
	std::vector<double> sum2ijArray(n);
	//std::vector<double> sum2iArray(n);
#pragma omp parallel for
	for (int i = 1; i < n - 1; ++i) // края не трогаются
	{
		for (int j = 1; j < n - 1; ++j)
		{
			sigmaSum[i][j] = sigma_summ(i * h, j * h);
		}

		sum2ijArray[i] = sigma_left_down_der(i * h) + sigma_right_up_der(i * h);
		//sum2jArray[i] = sigma_left_down_der(i * h) + sigma_right_up_der(i * h);
	}


	for (int it = 1; it <= T; ++it)
	{


		// цикл с одной большой волной в одну сторону 
		//double p1 = tau * K / (2 * h);
		//double uv1 = tau / (2 * h * rho);
		//double uvp2 = (tau * tau * K) / (2 * rho * h * h);
		//printf("Single\n");

		//for (int i = 0; i < n; ++i) // с ГУ




//#pragma omp parallel for
		for (int i = 1; i < n - 1; ++i) // края не трогаются
		{

			for (int j = 1; j < n - 1; ++j)
			{
#if 0
				auto sum2j = sigma_left_down_der(j * h) + sigma_right_up_der(j * h);
				auto sum2i = sigma_left_down_der(i * h) + sigma_right_up_der(i * h);
#else
				auto sum2j = sum2ijArray[j];
				auto sum2i = sum2ijArray[i];
#endif

#if 0
				auto sigSum = sigma_summ(i * h, j * h);
#else
				auto sigSum = sigmaSum[i][j];
#endif
				auto poldIJ = pold[i][j];
				auto uoldIJ = uold[i][j];
				auto voidIJ = vold[i][j];
				auto sig2Sum = sigSum * sigSum;

				auto kTau = 0.5 * tau * tau;

				K = 12.36;
				rho = 583.3;

				if ((i != 0 && j != 0 && i != n - 1 && j != n - 1) && (fractures[i][j] == No && (fractures[i + 1][j] == Up || fractures[i - 1][j] == Up || fractures[i][j + 1] == Up || fractures[i][j - 1]) ) )
				{
					continue;
				}
				if (fractures[i][j] == No)
				{
				
					/*if (j > 2 && fractures[i][j - 2] == Up && fractures[i][j - 1] == No)
					{
						//справа
						unew[i][j - 1] = unew[i][j];
						pnew[i][j - 1] = pnew[i][j];
						vnew[i][j - 1] = vnew[i][j];
					}
					else if (i > 2 && fractures[i - 2][j] == Up && fractures[i - 1][j] == No)
					{
						//снизу
						unew[i - 1][j] = unew[i][j];
						pnew[i - 1][j] = pnew[i][j];
						vnew[i - 1][j] = vnew[i][j];
					}

					else if (i < n - 2 && fractures[i + 1][j] == No && fractures[i + 2][j] == Up)
					{ // сверху
						unew[i + 1][j] = unew[i][j];
						pnew[i + 1][j] = pnew[i][j];
						vnew[i + 1][j] = vnew[i][j];
					}
					//слева
					else if (j < n - 2 && fractures[i][j + 1] == No && fractures[i][j + 2] == Up)
					{
						unew[i][j + 1] = unew[i][j];
						pnew[i][j + 1] = pnew[i][j];
						vnew[i][j + 1] = vnew[i][j];
					}
					else*/
					{
						pnew[i][j] = poldIJ - tau * (K * (vec_x(i, j, uold) + vec_y(i, j, vold)) + sigSum * poldIJ) \
							+ kTau *
							(
								K / rho * (vec_xx(i, j, pold) + vec_yy(i, j, pold)) + K * vec_x(i, j, uold) * 2 * sigSum \
								+ K * vec_y(i, j, vold) * 2 * sigSum \
								+ K * uoldIJ * (sum2i) \
								+ K * voidIJ * (sum2j) \
								+ sig2Sum * poldIJ
								);

						unew[i][j] = uoldIJ - tau * (vec_x(i, j, pold) / rho + sigSum * uoldIJ) \
							+ kTau *
							(
								K / rho * vec_xx(i, j, uold) + vec_x(i, j, pold) / rho * 2 * sigSum \
								+ (sum2i)*poldIJ / rho \
								+ sig2Sum * uoldIJ
								);

						vnew[i][j] = voidIJ - tau * (vec_y(i, j, pold) / rho + sigSum * voidIJ) \
							+ kTau *
							(
								K / rho * vec_yy(i, j, vold) + vec_y(i, j, pold) / rho * 2 * sigSum \
								+ (sum2j)*poldIJ / rho \
								+ sig2Sum * voidIJ
								);
					}
					/*if (fractures[i][j - 1] == Right) // предыдущий столбец трещина
					{
						unew[i][j - 1] = unew[i][j];
						pnew[i][j - 1] = pnew[i][j];
						vnew[i][j - 1] = vnew[i][j];
					}
					else if (fractures[i - 1][j] == Bottom) // предыдущая строка с трещинами
					{
						unew[i - 1][j] = unew[i][j];
						pnew[i - 1][j] = pnew[i][j];
						vnew[i - 1][j] = vnew[i][j];
					}*/
				}
				/*else if (fractures[i][j] == Up)
				{
					unew[i][j] = unew[i - 1][j];
					pnew[i][j] = pnew[i - 1][j];
					vnew[i][j] = vnew[i - 1][j];
				}
				else if (fractures[i][j] == Left)
				{
					unew[i][j] = unew[i][j - 1];
					pnew[i][j] = pnew[i][j - 1];
					vnew[i][j] = vnew[i][j - 1];
				}
				else if (fractures[i][j] == Center)
				{
					unew[i][j] = std::nan("1");
					pnew[i][j] = std::nan("1");
					vnew[i][j] = std::nan("1");
				}*/
				else if (fractures[i][j] == Up)
				{
					//cверху
					/*if (i > 2 && fractures[i - 1][j] == No && fractures[i - 2][j] == No)
					{
						unew[i - 1][j] = unew[i - 2][j];
						pnew[i - 1][j] = pnew[i - 2][j];
						vnew[i - 1][j] = vnew[i - 2][j];
					}
					//слева
					if (j > 2 && fractures[i][j - 1] == No && fractures[i][j - 2] == No)
					{
						unew[i][j - 1] = unew[i][j - 2];
						pnew[i][j - 1] = pnew[i][j - 2];
						vnew[i][j - 1] = vnew[i][j - 2];
					}*/

					/*unew[i][j] = std::nan("1");
					pnew[i][j] = std::nan("1");
					vnew[i][j] = std::nan("1");*/

					unew[i][j] = eps;
					pnew[i][j] = eps;
					vnew[i][j] = eps;
				}

				if (j > 2 && fractures[i][j - 2] == Up && fractures[i][j] == No && fractures[i][j - 1] == No)
				{
					//справа
					unew[i][j - 1] = unew[i][j];
					pnew[i][j - 1] = pnew[i][j];
					vnew[i][j - 1] = vnew[i][j];
				}
				if (i > 2 && fractures[i - 2][j] == Up && fractures[i][j] == No && fractures[i - 1][j] == No)
				{
					//снизу
					unew[i - 1][j] = unew[i][j];
					pnew[i - 1][j] = pnew[i][j];
					vnew[i - 1][j] = vnew[i][j];
				}

				if (i < n - 2 && fractures[i + 1][j] == No && fractures[i + 2][j] == Up)
				{ // сверху
					unew[i + 1][j] = unew[i][j];
					pnew[i + 1][j] = pnew[i][j];
					vnew[i + 1][j] = vnew[i][j];
				}
				//слева
				if (j < n - 2 && fractures[i][j + 1] == No && fractures[i][j + 2] == Up)
				{
					unew[i][j + 1] = unew[i][j];
					pnew[i][j + 1] = pnew[i][j];
					vnew[i][j + 1] = vnew[i][j];
				}

			} // j cycle
		} // i cycle
	 

		/*for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{

				if (i >= n / 2 - 1 && i <= n / 2 + 10)
				{
					if (j >= d * 2 + 5 && j <= d * 2 + d * 3)
					{
						pnew[i][j] = 0;
						unew[i][j] = 0;
						vnew[i][j] = 0;
					}
				}

			}
		}*/
#if 1
		if (it % 100 == 0)
		{
			std::ofstream foutU("ch_norm_2d_sigma_U" + std::to_string(k) + ".txt");
			std::ofstream foutP("ch_norm_2d_sigma_P" + std::to_string(k) + ".txt");
			std::ofstream foutV("ch_norm_2d_sigma_V" + std::to_string(k++) + ".txt");

			//for (int i = 2 * d; i < n - (2 * d); ++i) // без поглощающих слоев
			for (int i = d; i < n - d; ++i)
												 //for (int i = 0; i < n; ++i) // с ГУ
			{
				//for (int j = 2 * d; j < n - (2 * d); ++j) // без поглощающих слоев
				for (int j = d; j < n - d; ++j)
												//for (int j = 0; j < n; ++j) // с ГУ
				{
					if (fabs(unew[i][j] - eps) < 1e-8)
						foutU << std::setw(20) << std::nan("1") << std::setw(20);
					else
						foutU << std::setw(20) << unew[i][j] << std::setw(20);

					if (fabs(pnew[i][j] - eps) < 1e-8)
						foutP << std::setw(20) << std::nan("1") << std::setw(20);
					else
						foutP << std::setw(20) << pnew[i][j] << std::setw(20);

					if (fabs(vnew[i][j] - eps) < 1e-8)
						foutV << std::setw(20) << std::nan("1") << std::setw(20);
					else
						foutV << std::setw(20) << vnew[i][j] << std::setw(20);

					/*foutU << std::setw(20) << unew[i][j] << std::setw(20);
					foutP << std::setw(20) << pnew[i][j] << std::setw(20);
					foutV << std::setw(20) << vnew[i][j] << std::setw(20);*/
					
				}
				foutU << std::endl;
				foutP << std::endl;
				foutV << std::endl;
			}
		}
#endif
		pold = pnew;
		uold = unew;
		vold = vnew;
	} // end while time cycle
	auto end = std::chrono::system_clock::now();
	std::cout << "time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;
}



double temptest2d::sigma_x_right_or_y_up(const double& xy)
{
	auto leftBound = (n - d - 1) * h;

	if (xy >= leftBound && xy <= l)
	{
		double res = this->lod10 * pow(xy - leftBound, m) / pow(d * h, m + 1);
		return res;
	}
	else
		return 0;
}

double temptest2d::sigma_x_left_or_y_down(const double& xy)
{
	auto rightBound = d * h;
	if (xy >= 0 && xy <= rightBound)
	{
		double res = this->lod10 * pow(rightBound - xy, m) / pow(rightBound, m + 1);
		return  res;
	}
	else
		return 0;
}

// рандомный разброс трещин
void temptest2d::perkolation(int min_len, int max_len, std::vector<std::vector<int>> &fractures, std::vector<std::vector<int>>& fracture_centers, double proc)
{
	/*std::vector<std::vector<double>> numArrays;
	numArrays.resize(n);
	for (int i = 0; i < n; ++i)
	{
		numArrays[i].resize(n);
	}

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			numArrays[i][j] = rand();
		}
	}*/
	//while (true)
	int all = (n  - 1 - 2 * d) * proc;
	try
	{
		for (int k = 0; k < all; ++k)
		{
			/*if (all * 0.05 <= singles)
			{
				printf("enough\n");
				std::ofstream fout("fracture.txt");

				for (int i = d; i < n - d; ++i)
				{
					for (int j = d; j < n - d; ++j)
					{
						fout << fractures[i][j] << '\t';
					}
					fout << std::endl;
				}
				break;
			}*/
			int x = d + rand() % (n - 2 * d);
			int y = d + rand() % (n - 2 * d);
			int len = 0;
			try {
				if (min_len == max_len)
				{
					len = min_len;
				}
				else
				{
					len = min_len + rand() % (max_len - min_len);
				}
			}
			catch (std::logic_error &e)
			{
				std::cout << e.what() << std::endl;
			}
			fractures[x][y] = 1;
			fracture_centers[x][y] = 1;

			if (rand() % 2 == 0) // vert
			{
				for (int i = 1; i <= len; ++i)
				{
					
					//fractures[x][y - 1] = 1;
					//fractures[x][y + 1] = 1;
					if (x + i > n - d || x - i < d)
						break;

					//fractures[x - i][y - 1] = 1;
					fractures[x - i][y] = 1;
					//fractures[x - i][y + 1] = 1;

					//fractures[x + i][y - 1] = 1;
					fractures[x + i][y] = 1;
					//fractures[x + i][y + 1] = 1;
				}

			}
			else
			{
				for (int i = 1; i <= len; ++i)
				{
					//fractures[x - 1][y] = 1;
					//fractures[x + 1][y] = 1;
					if (y + i > n - d || y - i < d)
						break;

					//fractures[x - 1][y + i] = 1;
					fractures[x][y + i] = 1;
					//fractures[x + 1][y + i] = 1;

					//fractures[x - 1][y - i] = 1;
					fractures[x][y - i] = 1;
					//fractures[x + 1][y - i] = 1;
				}
			}
		}

		std::ofstream fout("fracture.txt");
		for (int i = d; i < n - d; ++i)
		{
			for (int j = d; j < n - d; ++j)
			{
				fout << fractures[i][j] << '\t';
			}
			fout << std::endl;
		}
	}
	catch (std::exception& e)
	{
		std::cout << e.what() << std::endl;
	}
}

void temptest2d::pyroman()
{

	//инициализация стартовая трещинами
	int lmin = (n - 1) / 100;
	int lmax = (n - 1) / 100;
	try
	{
		perkolation(lmin, lmax, fractures1, fracture_centers1, 0.5); // от 5 до 10 длины (при проценте от 0.25 до 0.5 опускается до 0.12 - 0.13)

		/*fracture_centers1[n / 2][n / 2] = 1;
		for (int i = d; i < n - d; ++i)
		{
			fractures1[n / 2][i] = 1;
		}
		write_fractures(fractures1);*/
	}
	catch (std::exception& e)
	{
		std::cout << e.what() << std::endl;
	}
	double f1 = 5;
	double f2 = 5;
	int k = 0;

	fractures2 = fractures1; // инит fractures который будет немножко изменен (нужно чтобы трещины не повернулись на 90 градусов а остались теми же кроме тех, у кого центры поменялись)
	while (true)
	{
		//change_some_centers(fracture_centers1, fracture_centers2, 0.01); // изменить процент центров трещин и поместить в (fracture_centers2) (old ver)
		std::vector<std::pair<int, int>> old_coords;
		std::vector<std::pair<int, int>> new_coords;
		change_some_centers_with_pairs(fracture_centers1, old_coords, new_coords, 0.01); // заполнение old_coords, new_coords на основе fr_cents1

		std::cout << "-----after changing centers-----" << std::endl;
		std::cout << "frcent1count " << fracture_centers_count(fracture_centers1) << std::endl;
		std::cout << "frt1count " << fracture_centers_count(fractures1) << std::endl;
		std::cout << "frcent2count " << fracture_centers_count(fracture_centers2) << std::endl;
		std::cout << "frt2count " << fracture_centers_count(fractures2) << std::endl;
		std::cout << "----------------------" << std::endl;
		// fractures1 хранит сами трещины
		// 

		//fill_lengths_having_centers(fracture_centers2, fractures2, lmin, lmax); // заполнение (fractures2) (n - 1) / 100 = 5, (n - 1) / 50 = 10 (old ver)

		fracture_centers2 = fracture_centers1; // копируем, затем на основе вектора пар заменим без всяких циклов
		fractures2 = fractures1; // видоизменим немного на основе fractures1
		fill_lengths_having_centers_with_pairs(fracture_centers2, fractures2, old_coords, new_coords, lmin, lmax); // меняем fractures2 и fr_centers2

		/*std::cout << "---------------------" << std::endl;
		std::cout << "after filling lengths" << std::endl;
		std::cout << "frcent1count " << fracture_centers_count(fracture_centers1) << std::endl;
		std::cout << "frt1count " << fracture_centers_count(fractures1) << std::endl;
		std::cout << "frcent2count " << fracture_centers_count(fracture_centers2) << std::endl;
		std::cout << "frt2count " << fracture_centers_count(fractures2) << std::endl;
		
		write_fractures(fractures1, -1);
		write_fractures(fractures2, -2);

		write_fractures(fracture_centers1, -10);
		write_fractures(fracture_centers2, -20);*/

		double f1 = f_computing(fractures1, fracture_centers1);
		
		std::cout << "f1 " << f1 << std::endl;
		
		if (f1 < 0.02) // идеально расположены трещины
		{
			fractures = fractures1;
			fracture_centers = fracture_centers1;
			std::cout << "Fractured ready" << std::endl;
			break;
		}
		double f2 = f_computing(fractures2, fracture_centers2);
		std::cout << "f2 " << f2 << std::endl;
		
		std::cout << "delta" << f2 - f1 << std::endl;
		if (f2 - f1 < 0) // новая модель
		{
			fractures1 = fractures2; // теперь fractures1 обновился, через него опять изменится fractures2, т.е. новый f1 будет совпадать со старым f2
			fracture_centers1 = fracture_centers2;

			/*std::cout << "---------------------" << std::endl;
			std::cout << "after delta < 0" << std::endl;
			std::cout << "frcent1count " << fracture_centers_count(fracture_centers1) << std::endl;
			std::cout << "frt1count " << fracture_centers_count(fractures1) << std::endl;
			std::cout << "frcent2count " << fracture_centers_count(fracture_centers2) << std::endl;
			std::cout << "frt2count " << fracture_centers_count(fractures2) << std::endl;*/

			write_fractures(fractures1, k++);
		}
		else
		{
			double pacc = Pacc(f2 - f1) * 100;
			double rnd = rand() % 101; // от 0 до 100 число
			std::cout << "rnd " << rnd << ", pacc " << pacc << std::endl;
			if (rnd <= pacc)
			{
				std::cout << "accepted perkolation" << std::endl;
				fractures1 = fractures2; // теперь fractures1 обновился, через него опять изменится fractures2, т.е. новый f1 будет совпадать со старым f2
				fracture_centers1 = fracture_centers2;

				/*std::cout << "---------------------" << std::endl;
				std::cout << "after accepting perkolation" << std::endl;
				std::cout << "frcent1count " << fracture_centers_count(fracture_centers1) << std::endl;
				std::cout << "frt1count " << fracture_centers_count(fractures1) << std::endl;
				std::cout << "frcent2count " << fracture_centers_count(fracture_centers2) << std::endl;
				std::cout << "frt2count " << fracture_centers_count(fractures2) << std::endl;*/
			}
			else
			{
				// иначе fracture1 остается без изменений и по новому изменится 1% всех трещин
				// на данный момент fractures1 - ниче такой и fractures2 - ужасный
				std::cout << "declined perkolation" << std::endl;
			}
		}
	}
}

double temptest2d::Pacc(double delta_f)
{
	return exp(- delta_f / T_pyro);
}

void temptest2d::change_some_centers(const std::vector<std::vector<int>>& fracture_centers1, std::vector<std::vector<int>>& fracture_centers2, double perc)
{
	int count1 = fracture_centers_count(fracture_centers1);
	if (count1 == 0)
	{
		std::cout << "Null vector of fractures centers1 change some centers" << std::endl;
		return;
	}
	int count2 = count1 * perc;

	fracture_centers2 = fracture_centers1;
	if (count2 == 0)
	{
		std::cout << "count2 is zero!!!!" << std::endl;
		return;
	}
	/*for (int i = d; i < n - d; ++i)
	{
		for (int j = d; j < n - d; ++j)
		{
			i
		}
	}*/
	while (count2 > 0)
	{
		int x = d + rand() % (n - 2 * d);
		int y = d + rand() % (n - 2 * d);

		//std::cout << "x: " << x << std::endl;
		//std::cout << "y: " << y << std::endl;
		//for (int i = d; i < n - d; ++i)
		//{
			int i = d + rand() % (n - 2 * d);
			int chance = rand() % 2;
			if (chance == 0)
			{
				if (fracture_centers2[i][y] == 1) // если нашлась строка с фикс столбцом (i,y)
				{
					int newcol = d + rand() % (n - 2 * d); // уносим на совсем др столбец
					if (fracture_centers2[x][newcol] == 0)
					{
						fracture_centers2[i][y] = 0; // зануляем старое значение
						fracture_centers2[x][newcol] = 1; // переносим на новое место единицу
						--count2;
					}
				}
			}
			else
			{
				if (fracture_centers2[x][i] == 1) // (x,i)
				{
					//x = d + rand() % (n - 2 * d); // уносим совсем на др строку
					int newrow = d + rand() % (n - 2 * d); // уносим на совсем др строку
					if (fracture_centers2[newrow][y] == 0)
					{
						fracture_centers2[x][i] = 0; // зануляем старое значение
						fracture_centers2[newrow][y] = 1; // переносим на новое место единиуц
						--count2;
					}
				}
			}
			
		//} // for

	} // while
}

void temptest2d::change_some_centers_with_pairs(const std::vector<std::vector<int>>& fracture_centers1, std::vector<std::pair<int, int>>& old_coords, std::vector<std::pair<int, int>>& new_coords, double perc)
{
	int count1 = fracture_centers_count(fracture_centers1);
	if (count1 == 0)
	{
		std::cout << "Null vector of fractures centers1 change some centers" << std::endl;
		return;
	}
	int count2 = count1 * perc;

	//fracture_centers2 = fracture_centers1;
	if (count2 == 0)
	{
		std::cout << "count2 is zero!!!!" << std::endl;
		return;
	}
	/*for (int i = d; i < n - d; ++i)
	{
		for (int j = d; j < n - d; ++j)
		{
			i
		}
	}*/
	while (count2 > 0)
	{
		int x = d + rand() % (n - 2 * d);
		int y = d + rand() % (n - 2 * d);

		//std::cout << "x: " << x << std::endl;
		//std::cout << "y: " << y << std::endl;
		//for (int i = d; i < n - d; ++i)
		//{
		int i = d + rand() % (n - 2 * d);
		int chance = rand() % 2;
		if (chance == 0) // проход по строке
		{
			if (fracture_centers1[i][y] == 1) // если нашлась строка с фикс столбцом (i,y)
			{
				int newcol = d + rand() % (n - 2 * d); // уносим на совсем др столбец
				if (fracture_centers1[x][newcol] == 0)
				{
					//fracture_centers1[i][y] = 0; // зануляем старое значение
					old_coords.push_back({ i, y });
					//fracture_centers1[x][newcol] = 1; // переносим на новое место единицу
					new_coords.push_back({ x, newcol });
					--count2;
				}
			}
		}
		else // по столбцу
		{
			if (fracture_centers1[x][i] == 1) // (x,i)
			{
				//x = d + rand() % (n - 2 * d); // уносим совсем на др строку
				int newrow = d + rand() % (n - 2 * d); // уносим на совсем др строку
				if (fracture_centers1[newrow][y] == 0)
				{
					//fracture_centers1[x][i] = 0; // зануляем старое значение
					old_coords.push_back({ x, i });
					//fracture_centers1[newrow][y] = 1; // переносим на новое место единиуц
					new_coords.push_back({ newrow, y });
					--count2;
				}
			}
		}

		//} // for

	} // while
}

void temptest2d::fill_lengths_having_centers(const std::vector<std::vector<int>>& fractures_centers, std::vector<std::vector<int>>& fractures, int min_len, int max_len)
{
	// 
	int count = fracture_centers_count(fractures_centers); // центры

	int len = 0;
	if (min_len == max_len)
	{
		len = min_len;
	}
	else
	{ 
		len = min_len + rand() % (max_len - min_len);
	}
		
	null_vector(fractures); // обнуление вектора, т.к полагаемся на центры
	for (int x = d; x < n - d; ++x)
	{
		for (int y = d; y < n - d; ++y)
		{
			if (fractures_centers[x][y] == 1) // если центр, то от него ответвляем длины 
			{
				fractures[x][y] = 1;
				if (rand() % 2 == 0) // vert
				{
					//std::cout << "vert" << std::endl;
					for (int i = 1; i <= len; ++i)
					{
						//fractures[x][y - 1] = 1;
						//fractures[x][y + 1] = 1;
						if (x + i > n - d || x - i < d)
							break;

						//fractures[x - i][y - 1] = 1;
						fractures[x - i][y] = 1;
						//fractures[x - i][y + 1] = 1;

						//fractures[x + i][y - 1] = 1;
						fractures[x + i][y] = 1;
						//fractures[x + i][y + 1] = 1;
					}

				}
				else
				{
					//std::cout << "horiz" << std::endl;
					for (int i = 1; i <= len; ++i)
					{
						//fractures[x - 1][y] = 1;
						//fractures[x + 1][y] = 1;
						if (y + i > n - d || y - i < d)
							break;

						//fractures[x - 1][y + i] = 1;
						fractures[x][y + i] = 1;
						//fractures[x + 1][y + i] = 1;

						//fractures[x - 1][y - i] = 1;
						fractures[x][y - i] = 1;
						//fractures[x + 1][y - i] = 1;
					}
				}
			}
		}
	}
}

// horiz == true, vert == false
void temptest2d::update_fracture(std::vector<std::vector<int>>& fractures, const std::pair<int, int> &old_point, const std::pair<int, int> &new_point, int len)
{
	int vert_count = 0;
	int horiz_count = 0;
	//int value = old_point.first == std::pair<int, int>().first && old_point.second == std::pair<int, int>().second ? 1 : 0; // если old point пустой, значит new_point задан и заполняется единицами (значение которым будет заполнен вектор)

	//std::pair<int, int> point = value == 1 ? new_point : old_point;

	// случай, если задан old_point, т.е. хотим опеределить тип трещины, чтобы занулить
	for (int i = -len; i <= len; ++i) 
	{
		if (fractures[old_point.first + i][old_point.second] == 1) { // вертикалка
			++vert_count;
		}

		if (fractures[old_point.first][old_point.second + i] == 1) { // горизонталка
			++horiz_count;
		}
	}

	//процесс удоления старой трещИны
	if (horiz_count > vert_count) { // горизонт стопудов

		for (int i = -len; i <= len; ++i)
		{
			if (fractures[old_point.first][old_point.second + i] == 1) { // горизонталка (столб меняется)

				if (fractures[old_point.first + 1][old_point.second + i] == 0 || fractures[old_point.first - 1][old_point.second + i] == 0) // если трещина не пересекает вертикальную
					fractures[old_point.first][old_point.second + i] = 0; // зануляем старую трещИну увы
			}

			fractures[new_point.first][new_point.second + i] = 1; // добавляем новую трещину
		}

	}
	else 
	{
		for (int i = -len; i <= len; ++i)
		{
			if (fractures[old_point.first + i][old_point.second] == 1) { // вертикалка ( строка меняется)

				if (fractures[old_point.first + i][old_point.second + 1] == 0 || fractures[old_point.first + i][old_point.second - 1] == 0) // если трещина не пересекает горизонтальную
					fractures[old_point.first + i][old_point.second] = 0;
			}

			fractures[new_point.first + i][new_point.second] = 1; // добавляем новую трещину
		}
	}
}

void temptest2d::fill_lengths_having_centers_with_pairs(std::vector<std::vector<int>>& fractures_centers, std::vector<std::vector<int>>& fractures, const std::vector<std::pair<int, int>>& old_coords, const std::vector<std::pair<int, int>>& new_coords, int min_len, int max_len)
{
	// fracture_centers заполнен без правок

	int len = 0;
	if (min_len == max_len)
	{
		len = min_len;
	}
	else
	{
		len = min_len + rand() % (max_len - min_len);
	}
	
	//обновим fracture_centers
	//теперь поправим fracture (затронем только лишь трешины центры которых есть в парах
	//занулим

	if (old_coords.size() != new_coords.size()) {
		std::cout << "Different dimentions for pairs!" << std::endl;
		return;
	}

	for (int i = 0; i < old_coords.size(); ++i) // проход по обоим парам
	{
		fractures_centers[old_coords[i].first][old_coords[i].second] = 0;
		fractures_centers[new_coords[i].first][new_coords[i].second] = 1;

		update_fracture(fractures, old_coords[i], new_coords[i], max_len); // обновить трещИну
		
	}

}
int temptest2d::fracture_centers_count(const std::vector<std::vector<int>>& fracture_centers)
{
	int Nfrac = 0;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			if (fracture_centers[i][j] == 1)
				++Nfrac;
		}
	}

	return Nfrac;
}

void temptest2d::null_vector(std::vector<std::vector<int>>& fracture_centers)
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			fracture_centers[i][j] = 0;
		}
	}
}

double temptest2d::f_p_computing(const std::vector<std::vector<int>>& fractures)
{
	double p_x_total = 0;
	double p_y_total = 0;
	int frac_len = (n - 1) / 100;
	for (int i = 0; i < Nw; ++i)
	{
		double px = (1 - p_i_x(i, fractures, frac_len = 0));
		double wi = w_i(i);
		p_x_total += wi * px * px / Nw;
	}

	p_x_total = sqrt(p_x_total);

	for (int i = 0; i < Nw; ++i)
	{
		double py = (1 - p_i_y(i, fractures, frac_len = 0));
		double wi = w_i(i);
		p_y_total += wi * py * py / Nw;
	}

	p_y_total = sqrt(p_y_total);
	double res = (p_x_total + p_y_total) / 2.;
	return res;
}

double temptest2d::f_d_computing(const std::vector<std::vector<int>>& fracture_centers)
{
	double d2 = D2_computing(fracture_centers);
	return fabs(0.5 * (d2 - 2));
}

double temptest2d::f_computing(const std::vector<std::vector<int>>& fractures, const std::vector<std::vector<int>>& fracture_centers)
{
	double f_p = f_p_computing(fractures); // согласно расположению трещин учитывая длины
	//std::cout << "f_p: " <<  f_p << std::endl;
	double f_d = f_d_computing(fracture_centers); // связность согласно центрам
	//std::cout << "f_d: " << f_d << std::endl;
	return 0.8 * f_p + 0.2 * f_d;
}

double temptest2d::w_i(int i)
{
	//double res = 2 * double(Nw - i) / (Nw * Nw + Nw);
	double res = double(rand() % RAND_MAX) / RAND_MAX / 100.;
	return res;
}



double temptest2d::p_i_x(int i, const std::vector<std::vector<int>>& fractures, int frac_len)
{
	int width_len = lw;
	if (frac_len != 0)
	{
		//std::cout << "pix with rectangle windows" << std::endl;
		width_len = int((1 + i) / 2. * frac_len);
	}
		int Nperc = 0;
		for (int j = d; j < n - d; j += width_len) // блоки (номер столбца в блоке i) |_|_|_|  - насечки (ширина блока)
		{
			for (int iy = j; iy < j + lw; ++iy) // внутри блока |||| - один блок между двумя насечками
			{
				if (iy > n - d)
					break;
				int len = 0;												   //___
				for (int ix = d + i * lw; ix < (i + 1) * lw + d; ++ix) // высота блока   ___
				{															   //___
					if (ix > n - d)
						break;
					if (fractures[ix][iy] == 1)
						++len;
				}
				if (len == lw)
				{
					++Nperc;
					break; // из iy, т.к. в этом блоке уже есть вертикальная трещина, значит ячейка в расчет попадает
				}
			}
		}
		double res = double(Nperc) / Nw;
		//std::cout << "p_i_x: Nperc " << Nperc << ", Nw " << Nw << std::endl;
		return res;
}

double temptest2d::p_i_y(int i, const std::vector<std::vector<int>>& fractures, int frac_len)
{
	int width_len = lw;
	if (frac_len != 0)
	{
		//std::cout << "pix with rectangle windows" << std::endl;
		width_len = int((1 + i) / 2. * frac_len);
	}
	int Nperc = 0;
	for (int j = d; j < n - d; j += width_len) // блоки (номер столбца в блоке i) |_|_|_|  - насечки
	{
		for (int ix = j; ix < j + lw; ++ix) // внутри блока |||| - один блок между двумя насечками
		{
			if (ix > n - d)
				break;
			int len = 0;												   //___
			for (int iy = d + i * lw; iy < (i + 1) * lw + d; ++iy)  // высота блока  ___
			{															   //___
				if (iy > n - d)
					break;
				if (fractures[ix][iy] == 1)
					++len;
			}
			if (len == lw)
			{
				++Nperc;
				break; // из ix, т.к. в этом блоке уже есть вертикальная трещина, значит ячейка в расчет попадает
			}
		}
	}
	//std::cout << "p_i_y: Nperc " << Nperc << ", Nw " << Nw << std::endl;
	double res = double(Nperc) / Nw;
	return res;
}


double temptest2d::D2_computing(const std::vector<std::vector<int>>& fracture_centers)
{
	double p1_2_sum = 0;
	p1_2_sum = Pj2sum(d1, fracture_centers);
	p1_2_sum = log(p1_2_sum);
	double p2_2_sum = 0;
	p2_2_sum = Pj2sum(d2, fracture_centers);
	p2_2_sum = log(p2_2_sum);

	
	double module = fabs(p2_2_sum - p1_2_sum);
	double D2 = module / (log(d2) - log(d1));
	return D2;
}

double temptest2d::Pjomega(int idxj, int d1, const std::vector<std::vector<int>>& fracture_centers)
{
	int Nperc = 0;
	int n_b = (n - 2 * d) / d1; // кол-во блоков в одном направлении
	
	/*for (int j = d; j < n - d; j += d1) // блоки (номер столбца в блоке i) |_|_|_|  - насечки
	{
		for (int ix = j; ix < j + 1; ++ix) // внутри блока |||| - один блок между двумя насечками
		{
			int len = 0;												   //___
			for (int iy = idxj * lw; iy < (idxj + 1) * lw; ++iy)  // высота блока  ___
			{															   //___
				if (fractures[ix][iy] == 1)
					++len;
			}
			if (len == lw)
			{
				++Nperc;
				break; // из ix, т.к. в этом блоке уже есть вертикальная трещина, значит ячейка в расчет попадает
			}
		}
	}*/

	/*for (int i = d; i < n - d; ++i)
	{
		for (int j = d; j < n - d; ++j)
		{
			if 
		}
	}*/
	int Nfrac = fracture_centers_count(fracture_centers);

	int k = -1;
	for (int i = d, px = 0; i < n - d; i += d1, ++px)
	{
		for (int j = d, py = 0; j < n - d; j += d1, ++py)
		{
			int Njomega = 0;
			for (int ix = i; ix < i + d1; ++ix)
			{
				if (ix >= (n - 1) - d)
					break;
				for (int iy = j; iy < j + d1; ++iy)
				{
					//cout << "(" << ix << "," <<  iy << ")" << endl;
					if (iy >= (n - 1) - d)
						break;
					if (fracture_centers[ix][iy] == 1)
						++Njomega;
				}
			}
			//cout << "blocky " << j << endl;
			k = px * n_b + py;
			if (k == idxj)
			{
				double res = double(Njomega) / Nfrac;
				//if (Njomega != 0)
				//	std::cout << "i,j = " << i << "," << j << " , px,nb,py: " << px << ", " << n_b << ", " << py << ", njomega: "  << Njomega << std::endl;
				return res;
			}
			else
				Njomega = 0;
		}

		//cout << "blockx " << i << endl;

	}

	return -5;
}

double temptest2d::Pj2sum(int d1, const std::vector<std::vector<int>>& fracture_centers)
{
	int n_b = (n - 2 * d) / d1; // кол-во блоков в одном направлении
	int all = n_b * n_b;

	double sum = 0;
	for (int i = 0; i < all; ++i)
	{
		double pj = Pjomega(i, d1, fracture_centers);
		if (pj < 0)
			std::cout << "Pj is negative" << std::endl;
		sum += (pj * pj);
	}

	return sum;
}

void temptest2d::perkolation_old(int min_len, int max_len, std::vector<std::vector<int>>& fractures, std::vector<std::vector<int>>& fracture_centers)
{
	try
	{
		for (int k = 0; k < 15; ++k)
		{


			int all = (n - 2 * d) * (n - 2 * d);
			int singles = 0;

			for (int i = d; i < n - d; ++i)
			{
				for (int j = d; j < n - d; ++j)
				{
					if (fractures[i][j] == 1)
						singles += 1;
				}
			}

			/*if (all * 0.05 <= singles)
			{
				printf("enough\n");
				std::ofstream fout("fracture.txt");

				for (int i = d; i < n - d; ++i)
				{
					for (int j = d; j < n - d; ++j)
					{
						fout << fractures[i][j] << '\t';
					}
					fout << std::endl;
				}
				break;
			}*/
			int x = d + rand() % (n - 2 * d);
			int y = d + rand() % (n - 2 * d);
			int len = 0;
			try {
				len = min_len + rand() % (max_len - min_len);
			}
			catch (std::logic_error& e)
			{
				std::cout << e.what() << std::endl;
			}
			fractures[x][y] = 1;

			fractures[x][y] = 1;

			fracture_centers[x][y] = 1;
			fracture_centers[x][y] = 1;

			if (rand() % 2 == 0) // vert
			{
				for (int i = 1; i <= len; ++i)
				{
					
					fractures[x][y - 1] = 1;
					fractures[x][y + 1] = 1;
					if (x + i > n - d || x - i < d)
						break;

					//fractures[x - i][y - 1] = 1;
					fractures[x - i][y] = 1;
					//fractures[x - i][y + 1] = 1;

					//fractures[x + i][y - 1] = 1;
					fractures[x + i][y] = 1;
					//fractures[x + i][y + 1] = 1;
				}

			}
			else
			{
				for (int i = 1; i <= len; ++i)
				{
					//fractures[x - 1][y] = 1;
					//fractures[x + 1][y] = 1;
					if (y + i > n - d || y - i < d)
						break;

					//fractures[x - 1][y + i] = 1;
					fractures[x][y + i] = 1;
					//fractures[x + 1][y + i] = 1;

					//fractures[x - 1][y - i] = 1;
					fractures[x][y - i] = 1;
					//fractures[x + 1][y - i] = 1;
					
				}
			}
		}

		std::ofstream fout("fracture.txt");
		for (int i = d; i < n - d; ++i)
		{
			for (int j = d; j < n - d; ++j)
			{
				fout << fractures[i][j] << '\t';
			}
			fout << std::endl;
		}
	}
	catch (std::exception& e)
	{
		std::cout << e.what() << std::endl;
	}
}

double temptest2d::u0(double x, double y)  // задание u0
{
	return alpha_u * exp(- ( a1_u * (x - l / 2)*(x - l / 2) + a2_u * (y - l / 2 ) * (y - l / 2)) );
	//return 1.0;
}

double temptest2d::p0(double x, double y)  // задание u0
{
	return alpha_p * exp(- (a1_p * (x - l / 2)*(x - l / 2) + a2_p * (y - l / 2) * (y - l / 2)) );
}

double temptest2d::v0(double x, double y)  // задание u0
{
	return alpha_v * exp(- (a1_v * (x - l / 2)*(x - l / 2) + a2_v * (y - l / 2) * (y - l / 2)) );
	//return 1.5;
}

double temptest2d::sigma_right_up_der(double xy)
{
	double eps = 1e-6;
	if (xy >= (n - d - 1) * h && xy <= l)
	{
		double man = (m + 1) * m * W * log(10) * pow((xy - ((n - d - 1) * h)), m - 1) / pow(d * h, m + 1);
		//right << num << '\t' << man << std::endl;
		return man;
	}
	return 0;
}

double temptest2d::sigma_left_down_der(double xy)
{
	double eps = 1e-6;
	if (xy >= 0 && xy <= d * h)
	{
		double man = (m + 1) * m * W * log(10) * pow(d * h - xy, m - 1) / pow(d * h, m + 1);
		//left << num << '\t' << man << std::endl;
		return man;
	}
	return 0;
}

double temptest2d::sigma_summ(double x, double y)
{
	return sigma_x_right_or_y_up(x) + sigma_x_right_or_y_up(y) + sigma_x_left_or_y_down(x) + sigma_x_left_or_y_down(y);
}

double temptest2d::vec_xx(int i, int j, const std::vector<std::vector<double>>& vec)
{
	return (vec[i + 1][j] - 2 * vec[i][j] + vec[i - 1][j]) / (h * h);
}

double temptest2d::vec_yy(int i, int j, const std::vector<std::vector<double>>& vec)
{
	return (vec[i][j + 1] - 2 * vec[i][j] + vec[i][j - 1]) / (h * h);
}

double temptest2d::vec_x(int i, int j, const std::vector<std::vector<double>>& vec)
{
	return (vec[i + 1][j] - vec[i - 1][j]) / (2 * h);
}

double temptest2d::vec_y(int i, int j, const std::vector<std::vector<double>>& vec)
{
	return (vec[i][j + 1] - vec[i][j - 1]) / (2 * h);
}

void temptest2d::write_fractures(const std::vector< std::vector<int>>& fractures, int k)
{
	std::ofstream fout("fracture" + std::to_string(k) + ".txt");
	for (int i = d; i < n - d; ++i)
	{
		for (int j = d; j < n - d; ++j)
		{
			fout << fractures[i][j] << '\t';
		}
		fout << std::endl;
	}
}

void temptest2d::init_fractures()
{
	

	/*for (int i = 3 * d ; i < n - 4 * d; ++i)
	{
		fractures[3.5 * d][i] = 1;
		fractures[3.5 * d + 1][i] = 1;
	}

	for (int i = 3.5 * d; i < 5 * d; ++i)
	{
		fractures[i][3 * d] = 1;
		fractures[i][3 * d + 1] = 1;
	}

	for (int i = 3 * d; i < n - 4 * d; ++i)
	{
		fractures[5 * d][i] = 1;
		fractures[5 * d + 1][i] = 1;
	}

	for (int i = 5 * d; i >= 2 * d; --i)
	{
		fractures[i][n - 4 * d] = 1;
		fractures[i][n - 4 * d + 1] = 1;
	}*/

	/*for (int i = 3 * d; i < n - 4 * d; ++i)
	{
		fractures[3.5 * d][i] = 1;
		fractures[3.5 * d + 1][i] = 1;
	}

	for (int i = 3 * d; i < n - 4 * d; ++i)
	{
		fractures[3.7 * d][i] = 1;
		fractures[3.7 * d + 1][i] = 1;
	}*/

	/*for (int i = 3 * d; i < n - 4 * d; ++i)
	{
		fractures[4.5 * d][i] = 1;
		fractures[4.5 * d + 1][i] = 1;
	}*/

	/*for (int i = 5 * d; i >= 2 * d; --i)
	{
		fractures[i][7 * d] = 1;
		fractures[i][7 * d + 1] = 1;
	}

	for (int i = 5 * d; i >= 2 * d; --i)
	{
		fractures[i][7 * d + 6] = 1;
		fractures[i][7 * d + 7] = 1;
	}*/

	//верхний горизонт
	/*for (int i = 3.5 * d; i < n - 3.5 * d; ++i)
	{
		fractures[i][4 * d] = 1;
		fractures[i][4 * d + 1] = 1;
		fractures[i][4 * d + 2] = 1;
	}

	//вверх от левого горизонта
	for (int i = 4 * d; i >= 2.5 * d; --i)
	{
		fractures[3.5 * d][i] = 1;
		fractures[3.5 * d + 1][i] = 1;
		fractures[3.5 * d + 2][i] = 1;
	}

	//влево от верхнего цикла
	for (int i = 3.5 * d; i >= 2.5 * d; --i)
	{
		fractures[i][2.5 * d] = 1;
		fractures[i][2.5 * d + 1] = 1;
		fractures[i][2.5 * d + 2] = 1;
	}

	for (int i = 7 * d; i >= 3.5 * d; --i)
	{
		fractures[n - 4.5 * d][i] = 1;
		fractures[n - 4.5 * d + 1][i] = 1;
		fractures[n - 4.5 * d + 2][i] = 1;
	}

	for (int i = n - 4.5 * d; i < n - 2.5 * d; ++i)
	{
		fractures[i][6 * d] = 1;
		fractures[i][6 * d + 1] = 1;
		fractures[i][6 * d + 2] = 1;
	}*/

	for (int i = d; i < n - d; ++i)
	{
		fractures[i][6 * d] = 1;
		fractures[i][6 * d] = 1;
		fractures[i][6 * d] = 1;
	}
	std::ofstream fout("fracture.txt");

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			fout << fractures[i][j] << '\t';
		}
		fout << std::endl;
	}
}
