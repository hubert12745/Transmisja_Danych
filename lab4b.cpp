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
int fs = 8000; //czêstotliwoœæ próbkowania
int Tc = 1;//czas trwania sygna³u
const int N = Tc * fs;//iloœæ próbek
int lg2n = (int)log2(N);//logarytm dwójkowy z iloœci próbek
const int N_fft = pow(2, lg2n);//iloœæ próbek po fft
const int N_freq = (int)N_fft / 2;//dziedzina czêstotliwoœci dla widma
int w = 1;
int W = 2;
vector<int> B = { 1,0,1,1,0,1,1,0, 1,1 };
int M = w * B.size();
double Tb = (double)Tc / M;
double Tbp = (double)Tb * fs;
double fn = W * 1 / Tb;
double A1 = 1.0 / 2.0;
double A2 = 2.0;


 class global {
public:
	static vector<double> time() {
		vector<double> t;

		for (int i = 0; i < N; i++) {
			double time = static_cast<double>(i) / fs;
			t.push_back(time);
		}
		return t;
	}
	static vector<double> translateToDecibels(vector<double> x) {
		vector<double> result;
		for (int i = 0; i < x.size(); i++) {
			result.push_back(10 * log10(x[i]));
		}
		return result;
	}
	static void adaptToComplex(vector<double> x, complex<double>* a) {
		for (int i = 0; i < N; i++) {
			a[i] = complex<double>(x[i], 0);
		}
	}
};
 class FFT {
 public:
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
	 void transform(Iter_T a, Iter_T b, int log2n)
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
 };
 class spectrum {
 public:
	 complex<double>* a;
	 complex<double>* b;
	 FFT fft;
	 spectrum() {
		 a = new complex<double>[N];
		 b = new complex<double>[N];
	 }
	 void drawSpectrum(vector<double> x, string title) {
		 global::adaptToComplex(x, a);
		 fft.transform(a, b, lg2n);
		 vector<double> freq;
		 vector<double> spectrum;
		 for (int i = 0; i < N_freq; i++) {
			 freq.push_back(i * fs / N);
			 spectrum.push_back(abs(b[i]) / N);
		 }
		 spectrum = global::translateToDecibels(spectrum);
		 plt::plot(freq, spectrum);
		 plt::grid();

		 plt::savefig(title);
		 plt::clf();
	 }
 };

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
}
int main() {
	complex<double>* a = new complex<double>[N];
	complex<double>* b = new complex<double>[N];
	vector<double> x;
	vector<double> t = global::time();
	x = ASK(t);
	spectrum s;
	s.drawSpectrum(x, "spectrum.png");
	

}