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
int fn = 5000;
int fm = 250;
float t = 0;//s
float fs = 22000;//Hz
float Tc = 1;//s
int n = (int)floor(Tc * fs);
int lg2n = (int)log2(n);
const int N = pow(2, lg2n);
const int N_freq = (int)N / 2;
float T = 1 / fs;//
float f = 1 / T;
string extension = ".png";
string za = "za_";
string zf = "zp_";
string zp = "zf_";
string wid = "_widmo";

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
int checkFirst(vector<double> x, int val) {
	int i = 0, first = 0;
	while (x[i] < val && i < x.size()) {
		first = i;
		i++;
	}
	return first;
}
int checkLast(vector<double> x, int val) {
	int last = 0;
	for (size_t i = 0; i < x.size(); i++)
	{
		if(x[i]>val) last = i;
	}
	return last;
}
void checkRange(vector<double> x, int val) {
	int range = checkLast(x, val)- checkFirst(x, val);
	printf("Range: %d\n", range);
}
vector<double> moveScale(vector<double> x) {
	vector<double> result;
	double max = *max_element(x.begin(), x.end());
	for (int i = 0; i < x.size(); i++) {
		result.push_back(x[i] - max);
	}
	return result;
}
vector<double> translateToDecibels(vector<double> x, int val) {
	vector<double> result;
	for (int i = 0; i < x.size(); i++) {
		result.push_back(10 * log10(x[i]));
	}
	result = moveScale(result);
	checkRange(result,-val);
	return result;
}
void adaptToComplex(vector<double> x, complex<double>* a) {
	for (int i = 0; i < N; i++) {
		a[i] = complex<double>(x[i], 0);
	}
}

vector<double> timeVector() {
	vector<double> t;
	for (int n = 0; n < N; n++) {
		t.push_back(n / fs);
	}
	return t;
}

double m(double t) {
	return sin(2 * M_PI * fm * t);
}
vector<vector<double>> modulateAM(double k) {
	vector<double> t = timeVector();
	vector<double> s;
	//vector<double> m;
	vector<vector<double>> result;
	for (int i = 0; i < N; i++) {
		s.push_back(((k * m(t[i])) + 1) * cos(2 * M_PI * fn * t[i]));
	}
	result.push_back(t);
	result.push_back(s);
	return result;
}
vector<vector<double>> modulatePM(double k) {
	vector<double> t = timeVector();
	vector<double> s;
	//vector<double> m;
	vector<vector<double>> result;
	for (int i = 0; i < N; i++) {
		s.push_back(cos((2 * M_PI * fn * t[i]) + k * m(t[i])));
	}
	result.push_back(t);
	result.push_back(s);
	return result;
}
vector<vector<double>> modulateFM(double k) {
	vector<double> t = timeVector();
	vector<double> s;
	//vector<double> m;
	vector<vector<double>> result;
	for (int i = 0; i < N; i++) {
		s.push_back(cos((2 * M_PI * fn * t[i]) + ((k / fm) * m(t[i]))));
	}
	result.push_back(t);
	result.push_back(s);
	return result;
}

vector<vector<double>> modulate(double k, int mode)
{
	switch (mode) {
	case 1:
		return modulateAM(k);
		break;
	case 2:
		return modulatePM(k);
		break;
	case 3:
		return modulateFM(k);
		break;

	}
}

void drawModulationPlots(double* kA, double* kP) {
	string stringset = "abc";
	for (int i = 0; i < 3; i++) {
		//AM
		string title = za + stringset[i] + extension;
		vector<vector<double>> s = modulate(kA[i], 1);
		plt::plot(s[0], s[1]);
		plt::grid();
		//plt::ylim(-2, 2);
		plt::savefig(title);
		plt::clf();

		//PM
		title = zp + stringset[i] + extension;
		vector<vector<double>> s1 = modulate(kP[i], 2);
		plt::plot(s1[0], s1[1]);
		plt::grid();
		//plt::ylim(-2, 2);
		plt::savefig(title);
		plt::clf();


		//FM
		title = zf + stringset[i] + extension;
		vector<vector<double>> s2 = modulate(kP[i], 3);
		plt::plot(s2[0], s2[1]);
		plt::grid();
		//plt::ylim(-2, 2);
		plt::savefig(title);
		plt::clf();

	}
}

void drawSpectrum(vector<double> x, string title, int mode = 0,int val = 0) {
	adaptToComplex(x, a);
	fft(a, b, lg2n);
	vector<double> freq;
	vector<double> magnitude;
	for (int i = 0; i < N_freq; i++) {
		freq.push_back(i * fs / N);
		magnitude.push_back(abs(b[i]) / N);
	}
	if (mode) magnitude = translateToDecibels(magnitude,val);
	plt::plot(freq, magnitude);
	plt::grid();
	if (mode) {
		plt::ylim(-val, val);
		plt::show();
	}

	if(!mode)plt::savefig(title);
	plt::clf();
}


int main() {
	int b = 12;
	double kA[] = { 0.5, 8, 50 };
	double kP[] = { -5, 3, 10 };
	drawModulationPlots(kA, kP);
	string stringset = "abc";
	for (int i = 0; i < 3; i++) {
		vector<vector<double>> s = modulate(kA[i], 1);
		drawSpectrum(s[1], za + stringset[i] + wid + extension);
		vector<vector<double>> s1 = modulate(kP[i], 2);
		drawSpectrum(s1[1], zp + stringset[i] + wid + extension);
		vector<vector<double>> s2 = modulate(kP[i], 3);
		drawSpectrum(s2[1], zf + stringset[i] + wid + extension);
	}

	for (int i = 0; i < 3; i++) {
		vector<vector<double>> s = modulate(kA[i], 1);
		cout << za + stringset[i] << " b"<<b<< endl;
		drawSpectrum(s[1], za + stringset[i] + wid + extension,1, b);
		vector<vector<double>> s1 = modulate(kP[i], 2);
		cout << zp + stringset[i] << " b" << b << endl;
		drawSpectrum(s1[1], zp + stringset[i] + wid + extension,1,b);
		vector<vector<double>> s2 = modulate(kP[i], 3);
		cout << zf + stringset[i] << " b" << b << endl;
		drawSpectrum(s2[1], zf + stringset[i] + wid + extension,1,b);
	}


	return 0;
}