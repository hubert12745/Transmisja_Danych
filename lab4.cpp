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
int M = 8;
int n = (int)floor(Tc * fs);
int lg2n = (int)log2(n);
const int N = pow(2, lg2n);
const int N_freq = (int)N / 2;
float Tb = Tc / M;
float Tbp = N / M;
float fn = 2 * (1/Tb);