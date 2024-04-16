#include <iostream>
#include <vector>
#include <cmath>
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

#define M_PI 3.14159265358979323846
float t = 0;//s
float fs = 10000;//Hz
float Tc = 2;//s
int N = (int)floor(Tc * fs);
float T = 1 / fs;//
float f = 1 / T;
std::vector<float> x;
std::vector<float> y;

void func()
{
	for (int n = 0; n < N; n++) {
		t = n / fs;
		x.push_back(t);
		float x = cos(2 * M_PI * f * t) * cos(2.5 * pow(t, 0.2) * M_PI);
		y.push_back(x);
	}
}

int main()
{
	func();
	plt::plot(x, y);
	plt::savefig("x.png");
	plt::show();
}