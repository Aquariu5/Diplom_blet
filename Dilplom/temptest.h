#pragma once
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include "adept.h"

using namespace adept;
class temptest
{
public:
	temptest();
	temptest(int n, int T, double tau);

	virtual double u0(double x); // ������� u0
	virtual double p0(double x); // ������� p0

	void analytical(); // ������������� �������


	void numerical_tvd();
	void numerical_classic(); // ������������ �������
	void numerical_classic_with_k(); // ������� � ���������� ��������� -k1(2)u
	virtual void numerical_classic_with_sigma(); // ������� � ����������� ��������� -sigma1(2)u
	void one_equation();

	//������ ��� ������� ��1
	void newton_gu1(adept::Stack & stack); // ������ � �� 1 ����
	void residual1(const std::vector<adouble> &unew_ad, const std::vector<adouble> &uold_ad, const std::vector<adouble> &pnew_ad, const std::vector<adouble> &pold_ad);
	void residual2(const std::vector<adouble> &unew_ad, const std::vector<adouble> &uold_ad, const std::vector<adouble> &pnew_ad, const std::vector<adouble> &pold_ad);
	bool checkresidual1();
	bool checkresidual2();

	//������ ��� ������� ��3
	void newton_gu3(adept::Stack & stack); // ������ � �� 3 ����

	//��������� �����
	void residual1_left(const std::vector<adouble> &unew_ad, const std::vector<adouble> &uold_ad, const std::vector<adouble> &pnew_ad, const std::vector<adouble> &pold_ad);

	//��������� ������
	void residual1_right(const std::vector<adouble> &unew_ad, const std::vector<adouble> &uold_ad, const std::vector<adouble> &pnew_ad, const std::vector<adouble> &pold_ad);

	//����� ���� ��� 1 ��-�
	void residual1_gu3(const std::vector<adouble> &unew_ad, const std::vector<adouble> &uold_ad, const std::vector<adouble> &pnew_ad, const std::vector<adouble> &pold_ad);

	// ����� ���� ��� 2 ��-�
	void residual2_gu3(const std::vector<adouble> &unew_ad, const std::vector<adouble> &uold_ad, const std::vector<adouble> &pnew_ad, const std::vector<adouble> &pold_ad);

	bool checkresidual1_gu3();
	bool checkresidual2_gu3();

	void decode_from_long_to_shorts(const std::vector<adouble> &longvec, std::vector<adouble> &short1, std::vector<adouble> &short2); // �� �������� ������� � ��������
	void encode_from_shorts_to_long(std::vector<adouble> &longvec, const std::vector<adouble> &short1, const std::vector<adouble> &short2); // �� �������� � ������� ������
	double rj(int j, const std::vector<double> & f);
	double Cj(int j,  std::vector<double> & r);

	//������� ���������
	double sigma_right(double x);
	double sigma_left(double x);

	// ����������� �� ������� ���������
	virtual double sigma_right_dx(double x);
	virtual double sigma_left_dx(double x);
	
	double ksi_j_plus(int j, std::vector<double> &vec);
	double ksi_j_minus(int j, std::vector<double> &vec);

	double u_or_p_j_plus(int j, std::vector<double> &vec);
	double u_or_p_j_minus(int j, std::vector<double> &vec);
	void plot_sigma();
	
	void max_wrong();
	double max = 0;
	void check_tvd(int j, const std::vector<double>& f, const std::vector<double>& C);
	void check(double &mu);
	~temptest();

	void getJacobian(adept::Stack & stack, std::vector< std::vector<double> > &jacobian, int n);

	void findRoot(int n); // ����������� �������

private:
	// ��� ��������� �����
	std::vector<double> uold;
	std::vector<double> unew;
	std::vector<double> pold;
	std::vector<double> pnew;

	// ��� ����� �����
	std::vector<double> uold_l;
	std::vector<double> unew_l;
	std::vector<double> pold_l;
	std::vector<double> pnew_l;

	//��� ������ �����
	std::vector<double> uold_r;
	std::vector<double> unew_r;
	std::vector<double> pold_r;
	std::vector<double> pnew_r;

	
	std::vector<double> uanal;
	std::vector<double> panal;

	//������� ��� TVD �����
	std::vector<double> r1;
	std::vector<double> r2;
	std::vector<double> C1;
	std::vector<double> C2;

	// ��� ������ �������
	std::vector<adouble> uold_ad; // ������ n
	std::vector<adouble> unew_ad; // ������ n
	std::vector<adouble> pold_ad; // ������ n
	std::vector<adouble> pnew_ad; // ������ n

	std::vector<adouble> xnew_ad; // ������ 2n
	std::vector<adouble> xold_ad; // ������ 2n
	std::vector<adouble> xsol_ad; // ������ 2n

	std::vector<adouble> f1_ad; // ������ n
	std::vector<adouble> f2_ad; // ������ n
	std::vector<adouble> f_ad; // ������ 2n
	std::vector< std::vector<double> > jacobian; // (2n)x(2n)

	std::ofstream left;
	std::ofstream right;
	double mu1, mu2; // ������������ ������ c1(2) * tau / h
	//��������� �����, ����� ��� ��

	//������� �������� ���������
	double k1, k2;

	//����� � ������ ���������� ��� ��������
	double k1_l;
	double k1_r;

	//����� � ������ ���������� ��� ��������
	double k2_l;
	double k2_r;


	//��������� ���������� �������
	double a1;
	double a2;
	double alpha1;
	double alpha2;

protected:
	int n;
	int T;
	double h;
	double tau;
	double l;
	double a;
	// ��������� �����
	double c;
	double rho;

	//��� ������� ��������� �����
	double m;
	int d;
	double W;

	
};
