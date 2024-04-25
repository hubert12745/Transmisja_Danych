#include <vector>
#include <cmath>
#include <iostream>
#include <complex>
#include <algorithm>
#include <iterator>
#include "matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;
#define M_PI 3.14159265358979323846


int fs = 8000;
int Tc = 1;
//int n = (int)floor(Tc * fs);
//int lg2n = (int)log2(n);
//const int N = pow(2, lg2n);
//const int N_freq = (int)N / 2;
const int N = Tc * fs;
int lg2n = (int)log2(N);
const int N_fft = pow(2, lg2n);
const int N_freq = (int)N_fft / 2;
int w = 1;
int W = 2;
vector<int> B = { 1,0,1,1,0,1,1,0, 1,1 };
int M = w * B.size();
double Tb = (double)Tc / M;
double Tbp = (double)Tb * fs;
double fn = W * 1 / Tb;
double A1 = 1.0 / 2.0;
double A2 = 2.0;

complex<double>* a = new complex<double>[N];
complex<double>* b = new complex<double>[N];

//implementacja fft https://www.oreilly.com/library/view/c-cookbook/0596007612/ch11s18.html
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

const double PI = 3.1415926536;

template<class Iter_T>
void fft(Iter_T a, Iter_T b, int log2n)
{
	typedef typename iterator_traits<Iter_T>::value_type complex;
	const complex J(0, 1);
	int n = 1 << log2n;
	for (unsigned int i = 0; i < n; ++i) {
		b[bitReverse(i, log2n)] = a[i];
	}
	for (int s = 1; s <= log2n; ++s) {
		int m = 1 << s;
		int m2 = m >> 1;
		complex w(1, 0);
		complex wm = exp(-J * (PI / m2));
		for (int j = 0; j < m2; ++j) {
			for (int k = j; k < n; k += m) {
				complex t = w * b[k + m2];
				complex u = b[k];
				b[k] = u + t;
				b[k + m2] = u - t;
			}
			w *= wm;
		}
	}
}
//koniec fft

vector<double> translateToDecibels(vector<double> x) {
	vector<double> result;
	for (int i = 0; i < x.size(); i++) {
		result.push_back(10 * log10(x[i]));
	}
	return result;
}
void adaptToComplex(vector<double> x, complex<double>* a) {
	for (int i = 0; i < N; i++) {
		a[i] = complex<double>(x[i], 0);
	}
}

void drawSpectrum(vector<double> x, string title, int mode = 0, int val = 0) {
	adaptToComplex(x, a);
	fft(a, b, lg2n);
	vector<double> freq;
	vector<double> magnitude;
	for (int i = 0; i < N_freq; i++) {
		freq.push_back(i * fs / N);
		magnitude.push_back(abs(b[i]) / N);
	}
	if (mode) magnitude = translateToDecibels(magnitude);
	plt::plot(freq, magnitude);
	plt::grid();
	if (mode) {
		plt::ylim(-val, val);
		plt::show();
	}

	if (!mode)plt::savefig(title);
	plt::clf();
}
vector<double> time() {
	vector<double> t;

	for (int i = 0; i < N; i++) {
		double time = static_cast<double>(i) / fs;
		t.push_back(time);
	}
	return t;
}

vector<double> ASK(vector<double> time) {

	vector<double> s(N);
	int z = (int)Tbp;

	for (int i = 0; i < B.size(); i++) {
		int start = i * z;
		int end = start + z;
		for (int j = start; j < end; j++) {
			if (B[i] == 0) {

				double a1 = A1 * sin(2 * M_PI * fn * time[j]);
				s[j] = a1;
			}
			else {
				double a2 = A2 * sin(2 * M_PI * fn * time[j]);
				s[j] = a2;
			}
		}
	}
	return s;
}vector<double> PSK(vector<double> time) {

	vector<double> s(N);
	int z = (int)Tbp;

	for (int i = 0; i < B.size(); i++) {
		int start = i * z;
		int end = start + z;
		for (int j = start; j < end; j++) {
			if (B[i] == 0) {

				double a1 = sin(2 * M_PI * fn * time[j]);
				s[j] = a1;
			}
			else {
				double a2 = sin(2 * M_PI * fn * time[j]+M_PI);
				s[j] = a2;
			}
		}
	}
	return s;
}vector<double> FSK(vector<double> time) {
	double fn1 = (double)(W+1)/Tb;
	double fn2 = (double)(W+2)/Tb;
	vector<double> s(N);
	int z = (int)Tbp;

	for (int i = 0; i < B.size(); i++) {
		int start = i * z;
		int end = start + z;
		for (int j = start; j < end; j++) {
			if (B[i] == 0) {

				double a1 = sin(2 * M_PI * fn1 * time[j]);
				s[j] = a1;
			}
			else {
				double a2 = sin(2 * M_PI * fn2 * time[j]);
				s[j] = a2;
			}
		}
	}
	return s;
}


int main() {
	vector<double> t = time();
	vector<double> ask = ASK(t);
	cout << "_____________________" << endl;
	cout << ask.size() << endl;
	plt::plot(t, ask);
	plt::grid(true);
	plt::savefig("za.png");
	plt::clf();
	vector<double> psk = PSK(t);
	plt::plot(t, psk);
	plt::grid(true);
	plt::savefig("zp.png");
	plt::clf();
	vector<double> fsk = FSK(t);
	plt::plot(t, fsk);
	plt::grid(true);
	plt::savefig("zf.png");
	plt::clf();

	return 0;
}

