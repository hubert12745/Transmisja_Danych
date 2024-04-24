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
int N = fs * Tc;
int n = 1;
int W = 2;
int M = n * 8;
double Tb = (double)Tc / M;
double Tbp = (double)Tb * fs;
double fn = W * 1 / Tb;
vector<int> B = { 1,0,1,1,0,1,1,0 };
double A1 = 1 / 2;
double A2 = 2;


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
	cout << z << endl;
	cout<<"_____________________"<<endl;

	for (int i = 0; i < B.size(); i++) {
		int start = i * z;
		cout << start << endl;
		int end = start + z;
		cout << end << endl;
		if (B[i] == 0) {
			for (int j = start; j < end; j++) {
				s.push_back(A1 * sin(2 * M_PI * fn * time[j]));
			}
		}
		else {
			for (int j = start; j < end; j++) {
				s.push_back(A2 * sin(2 * M_PI * fn * time[j]));
			}
		}
	}
	return s;
}

int main() {
	vector<double> t = time();
	vector<double> ask = ASK(t);
	plt::plot(t, ask);
	plt::show();
	return 0;
}

