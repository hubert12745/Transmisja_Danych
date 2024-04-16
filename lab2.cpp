#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <complex>
#include <iterator>
#include <chrono>
#include "matplotlibcpp.h"
using namespace std::chrono;

namespace plt = matplotlibcpp;
#define M_PI 3.14159265358979323846
//const std::vector<int> x = { 2, 0 ,-1, 3, 10, 2, 3, 1, 0 , 50, 0 ,-1, 3, 10, 2, 3, 1, 0 , 50 };
//const int N = x.size();
double t = 0;//s
double fs = 10000;//Hz
double Tc = 1;//s

int n = (int)floor(Tc * fs);
int lg2n = (int)log2(n);
const int N = pow(2, lg2n);
const int N_freq = (int)N / 2;
double T = 1 / fs;//
double f = 1 / T;
std::vector<double> t_arr;
std::vector<double> x_arr;
std::vector<double> y_arr;
std::vector<double> z_arr;
std::vector<double> v_arr;
std::vector<double> u_arr;
std::complex<double>* a = new std::complex<double>[N];
std::complex<double>* b = new std::complex<double>[N];

void adaptToComplex(std::vector<double> x, std::complex<double>* a) {
	for (int i = 0; i < N; i++) {
		a[i] = std::complex<double>(x[i], 0);
	}
}
unsigned int bitReverse(unsigned int x, int log2n) {
	int n = 0;
	int mask = 0x1;
	for (int i = 0; i < log2n; i++) {
		n <<= 1;
		n |= (x & 1);
		x >>= 1;
	}
	return n;
}

void fft(std::complex<double>* a, std::complex<double>* b, int log2n) {
	const std::complex<double> J(0, 1);
	int n = 1 << log2n;
	for (unsigned int i = 0; i < n; ++i) {
		b[bitReverse(i, log2n)] = a[i];
	}
	for (int s = 1; s <= log2n; ++s) {
		int m = 1 << s;
		int m2 = m >> 1;
		std::complex<double> w(1, 0);
		std::complex<double> wm = exp(-J * (M_PI / m2));
		for (int j = 0; j < m2; ++j) {
			for (int k = j; k < n; k += m) {
				std::complex<double> t = w * b[k + m2];
				std::complex<double> u = b[k];
				b[k] = u + t;
				b[k + m2] = u - t;
			}
			w *= wm;
		}
	}
}
//template<class Iter_T>
//void fft(Iter_T a, Iter_T b, int log2n)
//{
//	typedef typename iterator_traits<Iter_T>::value_type complex;
//	const complex J(0, 1);
//	int n = 1 << log2n;
//	for (unsigned int i = 0; i < n; ++i) {
//		b[bitReverse(i, log2n)] = a[i];
//	}
//	for (int s = 1; s <= log2n; ++s) {
//		int m = 1 << s;
//		int m2 = m >> 1;
//		complex w(1, 0);
//		complex wm = exp(-J * (PI / m2));
//		for (int j = 0; j < m2; ++j) {
//			for (int k = j; k < n; k += m) {
//				complex t = w * b[k + m2];
//				complex u = b[k];
//				b[k] = u + t;
//				b[k + m2] = u - t;
//			}
//			w *= wm;
//		}
//	}
//}

void func_x()
{
	for (int n = 0; n < N; n++) {
		t = n / fs;
		t_arr.push_back(t);
		double x = cos(2 * M_PI * f * t) * cos(2.5 * pow(t, 0.2) * M_PI);
		x_arr.push_back(x);
	}
}

void func_y(std::vector<double> x_arr, std::vector<double> t_arr)
{
	for (int n = 0; n < N; n++) {
		double y = (x_arr[n] * t_arr[n]) / (3 + cos(20 * M_PI * t_arr[n]));
		y_arr.push_back(y);
	}
}
void func_z(std::vector<double> x_arr, std::vector<double> t_arr, std::vector<double> y_arr) {
	for (int n = 0; n < N; n++) {
		double z = (pow(t_arr[n], 2) * abs(x_arr[n] * y_arr[n] - (2 / (10 + y_arr[n]))));
		z_arr.push_back(z);
	}
}

void func_v(std::vector<double> x_arr, std::vector<double> t_arr, std::vector<double> y_arr, std::vector<double> z_arr) {
	for (int n = 0; n < N; n++) {
		double v = pow(z_arr[n], 3) + 3 * sin(z_arr[n] * y_arr[n]) * abs(y_arr[n] - x_arr[n]);
		v_arr.push_back(v);
	}
}
void func_u()
{
	for (int n = 0; n < N; n++)
	{
		t = n / fs;
		t_arr.push_back(t);
		double u;
		if (t < 0.1)
		{
			u = (sin(6 * M_PI * t) * cos(5 * M_PI * t));
		}
		if (t < 0.4 && t >= 0.1)
		{
			u = -1.1 * t * cos(41 * M_PI * pow(t, 2));
		}
		if (t < 0.72 && t >= 0.4)
		{
			u = t * sin(20 * pow(t, 4));
		}
		if (t < 1 && t >= 0.72)
		{
			u = 3.3 * (t - 0.72) * cos(27 * t + 1.3);
		}
		u_arr.push_back(u);
	}
}

int H_arr[] = { 5, 20, 50 };
std::vector<std::vector<double>> b_arr;
void func()
{
	for (int n = 0; n < N; n++)
	{
		t = n / fs;
		t_arr.push_back(t);
	}
	for (int H = 0; H < 3; H++)
	{
		std::vector<double> bn_arr;
		for (int n = 0; n < N; n++)
		{
			t = t_arr[n];
			double sum = 0;
			for (int h = 1; h <= H_arr[H]; h++)
			{
				sum += (pow(-1, h) / h) * sin(h * M_PI * t);
			}
			double b = (2 / M_PI) * sum;
			bn_arr.push_back(b);
		}
		b_arr.push_back(bn_arr);
	}
}
std::vector<double> fk(N_freq);
std::vector<std::complex<double>> X(N);

void DFT(std::vector<double> x ,std::vector<std::complex<double>> cplx) {
	//for (int k = 0; k < N; k++) {
	//	X[k] = 0;
	//	for (int n = 0; n < N; n++) {
	//		double real = x[n] * cos((2 * M_PI * k * n) / N);
	//		double imag = x[n] * -sin((2 * M_PI * k * n) / N);
	//		X[k] += std::complex<double>(real, imag);
	//	}
	//}
	for (int k = 0; k < N; k++) {
		X[k] = 0;
		for (int n = 0; n < N; n++) {
			double angle = 2.0 * M_PI * k * n / N;
			std::complex<double> exp_term(cos(angle), -sin(angle)); // Complex exponential
			X[k] += x[n] * exp_term;
		}
	}
}

std::vector<double> spectrum(std::vector<std::complex<double>> cplx) {
	std::vector<double> M(N_freq);
	for (int k = 0; k < N_freq; k++) {
		double real = (double)X[k].real();
		double imag = (double)X[k].imag();
		double m = sqrt(pow(real, 2) + pow(imag, 2));
		M[k] = m;
	}
	return M;

}
std::vector<double> amplitude(std::vector<double> spectrum) {
	std::vector<double> M_amp(N_freq);

	for (int k = 0; k < N_freq; k++) {
		double amp = 10 * log10(spectrum[k]);
		M_amp[k] = amp;
		fk[k] = k * fs / N;
	}

	return M_amp;
}

int main() {
	//x
	int log2n = (int)log2(N);
	
	func_x();
	std::vector<double> temp = x_arr;
	DFT(x_arr,X);
	std::cout << "DFT" << std::endl;
	for (int i = 0; i < 10; i++) {
		std::cout << X[i] << std::endl;
	}
	adaptToComplex(x_arr, a);
	fft(a, b, log2n);
	std::cout << "FFT" << std::endl;
	for (int i = 0; i < 10; i++) {
		std::cout << b[i] << std::endl;
	}
	for (int i = 0; i < x_arr.size(); i++) {
		if (abs(x_arr[i] - temp[i]) > 0.0001) {
			std::cout << "Error in index "<< i << std::endl;
			break;
		}
	}
	//std::vector<double> M = spectrum(X);

	//std::vector<double> M_amp = amplitude(M);

	//plt::plot(fk, M_amp);
	//plt::grid();
	//plt::savefig("x.png");
	////plt::show();
	//plt::clf();
	////y
	//func_y(x_arr, t_arr);
	//DFT(y_arr,X);

	//M = spectrum(X);

	//M_amp = amplitude(M);

	//plt::plot(fk, M_amp);
	//plt::grid();
	//plt::savefig("y.png");
	////plt::show();	
	// plt::clf();
	////z
	//func_z(x_arr, t_arr, y_arr);
	//DFT(z_arr,X);

	//M = spectrum(X);

	//M_amp = amplitude(M);

	//plt::plot(fk, M_amp);
	//plt::grid();
	//plt::savefig("z.png");
	////plt::show();
	//plt::clf();
	////v
	//func_v(x_arr, t_arr, y_arr, z_arr);
	//DFT(v_arr,X);

	//M = spectrum(X);

	//M_amp = amplitude(M);

	//plt::plot(fk, M_amp);
	//plt::grid();
	//plt::savefig("v.png");
	////plt::show();	
	//plt::clf();
	//
	////u
	//func_u();
	//DFT(u_arr,X);

	//M = spectrum(X);

	// M_amp = amplitude(M);

	//plt::plot(fk, M_amp);
	//plt::grid();
	//plt::savefig("u.png");
	////plt::show();	
	//plt::clf();

	////b1
	//fs = 22050;
	//func();
	//DFT(b_arr[0], X);

	//M = spectrum(X);

	//M_amp = amplitude(M);

	//plt::plot(fk, M_amp);
	//plt::grid();
	//plt::savefig("b1.png");
	////plt::show();	
	//plt::clf();
	////b2
	//DFT(b_arr[1], X);

	//M = spectrum(X);
	//M_amp = amplitude(M);

	//plt::plot(fk, M_amp);
	//plt::grid();
	//plt::savefig("b2.png");
	////plt::show();	
	//plt::clf();
	////b3
	//DFT(b_arr[2], X);

	//M = spectrum(X);

	//M_amp = amplitude(M);

	//plt::plot(fk, M_amp);
	//plt::grid();
	//plt::savefig("b3.png");
	////plt::show();
	//plt::clf();

}