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


int fs = 1000;
int Tc = 1;
//int n = (int)floor(Tc * fs);
//int lg2n = (int)log2(n);
//const int N = pow(2, lg2n);
//const int N_freq = (int)N / 2;
const int N = (int)floor(Tc * fs);
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
		if (x[i] > val) last = i;
	}
	return last;
}
void checkRange(vector<double> x, int val) {
	int range = checkLast(x, val) - checkFirst(x, val);
	cout << range << endl;
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
	
	if (val != 0) {
		result = moveScale(result);
		checkRange(result, -val);
	}
	return result;
}
void adaptToComplex(vector<double> x, complex<double>* a) {
	for (int i = 0; i < N; i++) {
		a[i] = complex<double>(x[i], 0);
	}
}

void drawSpectrum(vector<double> x, string title, int value, bool save = 0) {

	adaptToComplex(x, a);
	fft(a, b, lg2n);
	vector<double> freq;
	vector<double> magnitude;
	for (int i = 0; i < N_freq; i++) {
		freq.push_back(i * fs / N);
		magnitude.push_back(abs(b[i]) / N);
	}
	magnitude = translateToDecibels(magnitude,value);
	plt::plot(freq, magnitude);
	plt::grid();
	//if (value != 0) {
	//	plt::ylim(-value, value);
	//	plt::show();
	//}
	save ? plt::show() : plt::savefig(title);
	//if (save)plt::savefig(title);
	plt::clf();
}
vector<double> timeVector() {
	vector<double> t;

	for (int i = 0; i < N; i++) {
		double time = static_cast<double>(i) / fs;
		t.push_back(time);
	}
	return t;
}

vector<vector<double>> ASK() {
	vector<vector<double>> result;
	vector<double> time = timeVector();
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
	result.push_back(time);
	result.push_back(s);
	return result;
}
	vector<vector<double>> PSK() {
		vector<vector<double>> result;
		vector<double> time = timeVector();
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
					double a2 = sin(2 * M_PI * fn * time[j] + M_PI);
					s[j] = a2;
				}
			}
		}
		result.push_back(time);
		result.push_back(s);
		return result;
	}
	vector<vector<double>> FSK() {
		vector<vector<double>> result;
		vector<double> time = timeVector();
		double fn1 = (double)(W + 1) / Tb;
		double fn2 = (double)(W + 2) / Tb;
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
		result.push_back(time);
		result.push_back(s);
		return result;
	}

	vector<double> multiplyBySine(vector<double> signal) {
		vector<double> time = timeVector();

		int z = (int)Tbp;
		int size = signal.size();
		vector<double> s(size);
		for (int i = 0; i < size; i++) {
			double a = sin(2*M_PI * fn * time[i]);
			s[i] = a;
		}
		for (int i = 0; i < size; i++) {
			double a = signal[i] * s[i];
			s[i] = a;
		}
		return s;
	}
	void multiplyPlot(vector<vector<double>> (*func)()) {
		vector<vector<double>> result = func();
		result[1] = multiplyBySine(result[1]);
		plt::plot(result[0], result[1]);
		plt::show();
	}

	vector<vector<double>> calculateSum(vector<vector<double>>(*func)()) {
		vector<vector<double>> data = func();
		vector<double> s = multiplyBySine(data[1]);
		vector<double> sumSignal(s.size());
		double sum = 0;
		int z = (int)Tbp;
		for (int i = 0; i < B.size(); i++) {
			int start = i * z;
			int end = start + z;
			sum = 0;
			//cout << "\t bit:"<< i <<endl;
			for (int j = start; j < end; j++) {
				double a = s[j];
				sum += a;
				sumSignal[j] = sum;
				//cout << sum << endl;
			}
		}
		data[1] = sumSignal;
		return data;

	}
	void drawPlots(bool show = false) {
		vector<vector<double>> ask = ASK();
		plt::plot(ask[0], ask[1]);
		plt::grid(true);
		if(show) plt::show();
		plt::savefig("za.png");
		plt::clf();
		vector<vector<double>> psk = PSK();
		plt::plot(psk[0], psk[1]);
		plt::grid(true);
		if (show) plt::show();
		plt::savefig("zp.png");
		plt::clf();
		vector<vector<double>> fsk = FSK();
		plt::plot(psk[0], fsk[1]);
		plt::grid(true);
		if (show) plt::show();
		plt::savefig("zf.png");
		plt::clf();
	}
	void drawSpectrums(int val, bool save = 0) {
		int oldFn = fn;
		fn = 250;
		vector<vector<double>> ask = ASK();
		if (val != 0) cout << "ASK: ";
		drawSpectrum(ask[1], "za_widmo.png",val, save);
		vector<vector<double>> psk = PSK();
		if (val != 0) cout << "PSK: ";
		drawSpectrum(psk[1], "zp_widmo.png",val, save);
		vector<vector<double>> fsk = FSK();
		if (val != 0) cout << "FSK: ";
		drawSpectrum(fsk[1], "zf_widmo.png",val,save);
		fn = oldFn;
	}
	

	int main() {
		vector<int> bd = { 3, 6, 12 };
		//drawPlots();
		//fn = 250;
		vector<vector<double>> ask = calculateSum(ASK);
		
		plt::plot(ask[0], ask[1]);
		plt::show();
		plt::clf();
		//multiplyPlot(ASK);
		//drawSpectrums(0,true);

		//for (int i : bd) {
		//	cout << "B" << i << "D: "<<endl;
		//	drawSpectrums(i);
		//	cout << endl;
		//}
		return 0;
	}

