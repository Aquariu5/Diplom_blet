#pragma once
#include "temptest.h"
class temptest2d: public temptest
{
public:
	temptest2d();
	~temptest2d();
	temptest2d(int n, int T, double tau);

	void numerical_classic();
	void numerical_classic_with_sigma();
	void norm_with_sigma();
	double sigma_x_left_or_y_down(const double& xy);
	double sigma_x_right_or_y_up(const double& xy);

	
	//отжиг
	double f_p_computing(const std::vector<std::vector<int>>& fractures);
	double f_d_computing(const std::vector<std::vector<int>>& fracture_centers);
	double f_computing(const std::vector<std::vector<int>>& fractures, const std::vector<std::vector<int>>& fracture_centers);
	double w_i(int i);

	double p_i_x(int i, const std::vector<std::vector<int>>& fractures, int frac_len = 0); // для fp с квадратными окнами
	double p_i_y(int i, const std::vector<std::vector<int>>& fractures, int frac_len = 0); // для fp с квадратными окнами

	double D2_computing(const std::vector<std::vector<int>>& fracture_centers);
	double Pjomega(int j, int d1, const std::vector<std::vector<int>> &fracture_centers);
	double Pj2sum(int d1, const std::vector<std::vector<int>>& fracture_centers);
	void perkolation_old(int min_len, int max_len, std::vector<std::vector<int>>& fractures, std::vector<std::vector<int>>& fracture_centers);
	void perkolation(int min_len, int max_len, std::vector<std::vector<int>>& fractures, std::vector<std::vector<int>>& fracture_centers, double proc);
	void pyroman();
	double Pacc(double delta_f);
	void change_some_centers(const std::vector<std::vector<int>>& fracture_centers1, std::vector<std::vector<int>>& fracture_centers2, double perc);
	void fill_lengths_having_centers(const std::vector<std::vector<int>>& fractures_centers, std::vector<std::vector<int>>& fractures, int min_len, int max_len);
	int fracture_centers_count(const std::vector<std::vector<int>>& fracture_centers);
	void null_vector(std::vector<std::vector<int>>& fracture_centers);
	///
	virtual double u0(double x, double y); // задание u0
	virtual double p0(double x, double y); // задание p0
	virtual double v0(double x, double y); // задание v0

	virtual double sigma_right_up_der(double xy);
	virtual double sigma_left_down_der(double xy);

	double sigma_summ(double x, double y);
	//операторы
	double vec_xx(int i, int j, const std::vector < std::vector<double> > &vec);
	double vec_yy(int i, int j, const std::vector < std::vector<double> > &vec);

	double vec_x(int i, int j, const std::vector < std::vector<double> > &vec);
	double vec_y(int i, int j, const std::vector < std::vector<double> > &vec);

	void init_fractures();
	void write_fractures(const std::vector < std::vector<int> >& vec, int k = 0);
private:
	std::vector < std::vector<double> > uold;
	std::vector < std::vector<double> > unew;

	std::vector < std::vector<double> > pold;
	std::vector < std::vector<double> > pnew;

	std::vector < std::vector<double> > vold;
	std::vector < std::vector<double> > vnew;

	std::vector < std::vector<int> > fractures;
	std::vector < std::vector<int>> fracture_centers;

	std::vector < std::vector<int> > fractures1;
	std::vector < std::vector<int>> fracture_centers1;

	std::vector < std::vector<int> > fractures2;
	std::vector < std::vector<int>> fracture_centers2;

	//отжиг
	int Nw;
	int lw;
	//double D2;
	double d1; // для Pjomega - размер квадрата
	double d2; // для Pjomega - размер квадрата
	double T_pyro;
	//альфы  для амплитуд
	double alpha_u;
	double alpha_v;
	double alpha_p;

	//ашки  для пологости
	double a1_u;
	double a2_u;

	double a1_p;
	double a2_p;

	double a1_v;
	double a2_v;

	//упругость среды
	double K;


	double lod10;

	

};
