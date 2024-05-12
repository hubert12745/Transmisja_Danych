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
	magnitude = translateToDecibels(magnitude, value);
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

vector<int> bitsVec() {
	vector<int> s(N);
	int z = (int)Tbp;

	for (int i = 0; i < B.size(); i++) {
		int start = i * z;
		int end = start + z;
		for (int j = start; j < end; j++) {
			if (B[i] == 0) {
				s[j] = 0;
			}
			else {
				s[j] = 1;
			}
		}
	}

	return s;
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

vector<double> multiplyBySine(vector<double> signal, int mode = 0) {
	vector<double> time = timeVector();
	double fn1 = (double)(W + 1) / Tb;
	double fn2 = (double)(W + 2) / Tb;
	int z = (int)Tbp;
	int size = signal.size();
	vector<double> s(size);
	if (mode == 0) {
		for (int i = 0; i < size; i++) {
			double a = sin(2 * M_PI * fn * time[i]);
			s[i] = a;
		}
	}
	else if (mode == 1) {
		for (int i = 0; i < size; i++) {
			double a = sin(2 * M_PI * fn1 * time[i]);
			s[i] = a;
		}
	}
	else if (mode == 2) {
		for (int i = 0; i < size; i++) {
			double a = sin(2 * M_PI * fn2 * time[i]);
			s[i] = a;
		}
	}
	for (int i = 0; i < size; i++) {
		double a = signal[i] * s[i];
		s[i] = a;
	}
	return s;
}
void multiplyPlot(vector<vector<double>>(*func)()) {
	vector<vector<double>> result = func();
	result[1] = multiplyBySine(result[1]);
	plt::plot(result[0], result[1]);
	plt::show();
}

//calculaetSum zwraca wektor wektorow, gdzie pierwszy wektor to czas, a drugi to sygnal dla funkcji podanej jako argument
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
		for (int j = start; j < end; j++) {
			double a = s[j];
			sum += a;
			sumSignal[j] = sum;
		}
	}
	data[1] = sumSignal;
	return data;

}
//calculateSum zwraca wektor wektorow, gdzie pierwszy wektor to czas, a drugi to sygnal dla funkcji FSK, jeœli parametrem jest tryb, a nie funkcja. Tryb 
vector<vector<double>> calculateSum(int mode) {
	vector<vector<double>> data = FSK();
	vector<double> s = multiplyBySine(data[1], mode);
	vector<double> sumSignal(s.size());
	double sum = 0;
	int z = (int)Tbp;
	for (int i = 0; i < B.size(); i++) {
		int start = i * z;
		int end = start + z;
		sum = 0;
		for (int j = start; j < end; j++) {
			double a = s[j];
			sum += a;
			sumSignal[j] = sum;
		}
	}
	data[1] = sumSignal;
	return data;

}
void drawPlots(bool show = false) {
	vector<vector<double>> ask = ASK();
	plt::plot(ask[0], ask[1]);
	plt::grid(true);
	plt::title("ASK_z");
	if (show) plt::show();
	plt::savefig("plots/zad5/ask_z.png");
	plt::clf();
	vector<vector<double>> psk = PSK();
	plt::plot(psk[0], psk[1]);
	plt::grid(true);
	plt::title("PSK_z");
	if (show) plt::show();
	plt::savefig("plots/zad5/psk_z.png");
	plt::clf();
	vector<vector<double>> fsk = FSK();
	plt::plot(psk[0], fsk[1]);
	plt::grid(true);
	plt::title("FSK_z");
	if (show) plt::show();
	plt::savefig("plots/zad5/fsk_z.png");
	plt::clf();
}
void drawSpectrums(int val, bool save = 0) {
	double oldFn = fn;
	fn = 250;
	vector<vector<double>> ask = ASK();
	if (val != 0) cout << "ASK: ";
	drawSpectrum(ask[1], "za_widmo.png", val, save);
	vector<vector<double>> psk = PSK();
	if (val != 0) cout << "PSK: ";
	drawSpectrum(psk[1], "zp_widmo.png", val, save);
	vector<vector<double>> fsk = FSK();
	if (val != 0) cout << "FSK: ";
	drawSpectrum(fsk[1], "zf_widmo.png", val, save);
	fn = oldFn;
}

//mode 0 - wartoœæ progowa 0, p(t) > 0(FSK),		mode 1 - wartoœæ progowa 0, p(t) < 0(PSK),		mode 2 - wartoœæ progowa h, p(t) > h(ASK)
vector<int> binarizeData(vector<double> x, int mode) {
	vector<int> result(x.size());
	double trigger;

	switch (mode) {
	case 0: {
		trigger = 0;
		for (int i = 0; i < x.size(); i++)
			result[i] = (x[i] < trigger) ? 1 : 0;
		break;
	}
	case 1: {
		trigger = 0;
		for (int i = 0; i < x.size(); i++)
			result[i] = (x[i] > trigger) ? 1 : 0;
		break;
	}

	case 2: {
		trigger = (*max_element(x.begin(), x.end()) / 2)-1;
		for (int i = 0; i < x.size(); i++)
			result[i] = (x[i] > trigger) ? 1 : 0;
		break;
	}
	}
	return result;


}
vector<vector<double>> combineFSK() {
	vector<vector<double>> fsk1 = calculateSum(1);
	vector<vector<double>> fsk2 = calculateSum(2);
	vector<double> fsk(fsk2[1].size());
	for (int i = 0; i < fsk2[1].size(); i++)
	{
		double sum = -fsk1[1][i] + fsk2[1][i];
		fsk[i] = sum;
	}
	vector<vector<double>> result;
	result.push_back(fsk2[0]);
	result.push_back(fsk);
	return result;
}
vector<int> fitToBits(vector<double> x, int mode) {
	vector<int> signal = binarizeData(x, mode);
	vector<int> result(signal.size());
	int z = (int)Tbp;
	int half = (z / 2);
	for (int i = 0; i < B.size(); i++) {
		int start = i * z;
		int end = start + z;
		int sum = 0;
		for (int j = start; j < end; j++) {
			sum += signal[j];

		}

		if (sum > half) {
			for (int j = start; j < end; j++) {
				result[j] = 1;
			}
		}
		else {
			for (int j = start; j < end; j++) {
				result[j] = 0;
			}
		}
	}
	return result;
}
void z() {
	drawPlots();
}
void x() {
	vector<vector<double>> ask = calculateSum(ASK);
	plt::plot(ask[0], ask[1]);
	plt::title("ASK_x");
	plt::grid(true);
	plt::savefig("plots/zad5/ask_x.png");
	plt::clf();

	vector<vector<double>> psk = calculateSum(PSK);
	plt::plot(psk[0], psk[1]);
	plt::title("PSK_x");
	plt::grid(true);
	plt::savefig("plots/zad5/PSK_x.png");
	plt::clf();

	vector<vector<double>> fsk1 = calculateSum(1);
	plt::plot(fsk1[0], fsk1[1]);
	plt::title("FSK_x1");
	plt::grid(true);
	plt::savefig("plots/zad5/FSK_x1.png");
	plt::clf();

	vector<vector<double>> fsk2 = calculateSum(2);
	plt::plot(fsk2[0], fsk2[1]);
	plt::title("FSK_x2");
	plt::grid(true);
	plt::savefig("plots/zad5/FSK_x2.png");
	plt::clf();
}

void p() {
	vector<vector<double>> ask = calculateSum(ASK);
	vector<int> askB = binarizeData(ask[1], 2);
	plt::plot(ask[0], askB);
	plt::title("ASK_p");
	plt::grid(true);
	plt::savefig("plots/zad5/ask_p.png");
	plt::clf();

	vector<vector<double>> psk = calculateSum(PSK);
	vector<int> pskB = binarizeData(psk[1], 1);
	plt::plot(psk[0], pskB);
	plt::title("PSK_p");
	plt::grid(true);
	plt::savefig("plots/zad5/PSK_p.png");
	plt::clf();

	vector<vector<double>> fsk1 = calculateSum(1);
	vector<int> fsk1B = binarizeData(fsk1[1], 0);
	plt::plot(fsk1[0], fsk1B);
	plt::title("FSK_p1");
	plt::grid(true);
	plt::savefig("plots/zad5/FSK_p1.png");
	plt::clf();

	vector<vector<double>> fsk2 = calculateSum(2);
	vector<int> fsk2B = binarizeData(fsk2[1], 0);
	plt::plot(fsk2[0], fsk2B);
	plt::title("FSK_p2");
	plt::grid(true);
	plt::savefig("plots/zad5/FSK_p2.png");
	plt::clf();

	vector<vector<double>> fsk = combineFSK();
	vector<int> fskB = binarizeData(fsk[1], 0);
	plt::plot(fsk[0], fskB);
	plt::title("FSK_p");
	plt::grid(true);
	plt::savefig("plots/zad5/FSK_p.png");
	plt::clf();

}
void c() {
		vector<vector<double>> ask = calculateSum(ASK);

	vector<int> askF = fitToBits(ask[1], 2);
	plt::plot(ask[0], askF);
	plt::title("ASK_c");
	plt::grid(true);
	plt::savefig("plots/zad5/ask_c.png");
	plt::clf();

	vector<vector<double>> psk = calculateSum(PSK);

	vector<int> pskF = fitToBits(psk[1], 1);
	plt::plot(psk[0], pskF);
	plt::title("PSK_c");
	plt::grid(true);
	plt::savefig("plots/zad5/PSK_c.png");
	plt::clf();

	vector<vector<double>> fsk = combineFSK();
	vector<int> fskF = fitToBits(fsk[1], 0);
	plt::plot(fsk[0], fskF);
	plt::title("FSK_c");
	plt::grid(true);
	plt::savefig("plots/zad5/FSK_c.png");
	plt::clf();	
}

void cCompared() {
	vector<vector<double>> ask = calculateSum(ASK);

	vector<int> askF = fitToBits(ask[1], 2);
	plt::plot(ask[0], askF,"green");
	vector<int> bits = bitsVec();
	plt::plot(ask[0], bits, { {"c", "red"}, {"ls", ":"} });
	plt::title("porownanie ASK");
	plt::grid(true);
	plt::savefig("plots/zad5/porownanie_ask.png");
	plt::clf();

	vector<vector<double>> psk = calculateSum(PSK);
	vector<int> pskF = fitToBits(psk[1], 1);
	plt::plot(psk[0], pskF, "green");
	plt::plot(psk[0], bits, { {"c", "red"}, {"ls", ":"} });
	plt::title("porownanie PSK");
	plt::grid(true);
	plt::savefig("plots/zad5/porownanie_psk.png");
	plt::clf();

	vector<vector<double>> fsk = combineFSK();
	vector<int> fskF = fitToBits(fsk[1], 0);
	plt::plot(fsk[0], fskF, "green");
	plt::plot(fsk[0], bits, { {"c", "red"}, {"ls", ":"} });
	plt::title("porownanie FSK");
	plt::grid(true);
	plt::savefig("plots/zad5/porownanie_fsk.png");
	plt::clf();

}


int main() {

	z();
	x();
	p();
	c();
	cCompared();
	
	return 0;
}

