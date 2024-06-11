#include <vector>
#include <cmath>
#include <iostream>
#include <complex>
#include <algorithm>
#include <iterator>
#include "matplotlibcpp.h"
#include "Eigen/Dense"


using namespace std;
using namespace Eigen;

namespace plt = matplotlibcpp;
#define M_PI 3.14159265358979323846


class Hamming {
public:
	static const int k = 11, n = 15;
	static const int m = n - k;
	int x1, x2, x4, x8;


	MatrixXi I = Matrix<int, k, k>::Identity();
	MatrixXi P;
	Matrix<int, k, n> G;
	MatrixXi H;

	void removeRow(MatrixXi& matrix, unsigned int rowToRemove) {
		unsigned int numRows = matrix.rows() - 1;
		unsigned int numCols = matrix.cols();

		if (rowToRemove < numRows)
			matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);

		matrix.conservativeResize(numRows, numCols);
	}
	void generateH() {
		MatrixXi i = Matrix<int, n - k, n - k>::Identity();
		MatrixXi pt = P.transpose();
		H.resize(m, n);
		H << i, pt;
	}
	void generateG() {
		G.resize(k, n);
		G << P, I;

	}


	void fillP() {
		P = Matrix<int, n + 1, m>::Zero(n + 1, m);
		for (int i = 0; i < n + 1; i++) {
			for (int j = 0; j < m; j++) {
				P(i, j) = (i >> j) & 1;
			}
		}
		removeRow(P, 0);
		for (int i = 1, delCount = 1; i < n + 1; i *= 2) {
			int idx = i - delCount;
			removeRow(P, i - delCount);
			delCount++;
		}
	}


	Hamming() {
		fillP();
		generateG();
		generateH();

	}

	void printI() {
		cout << I << endl;
	}
	void printP() {
		cout << P << endl;
	}
	void printG() {
		cout << G << endl;
	}
	void printH() {
		cout << H << endl;
	}
	MatrixXi moduloMatrix(MatrixXi matrix, int num) {
		MatrixXi result = MatrixXi::Zero(matrix.rows(), matrix.cols());
		result = (matrix.array() - (num * (matrix.array() / num))).matrix();
		return result;
	}

	MatrixXi encode(MatrixXi data) {
		Matrix<int, 1, n> result = (data * G);
		Matrix<int, 1, n> modArr = moduloMatrix(result, 2);
	/*	cout << "Data: " << data << endl;
		cout << "Encoded data: " << modArr << endl;*/
		return modArr;
	}
	MatrixXi decode(MatrixXi encodedMessage) {

		MatrixXi syndrome = moduloMatrix(encodedMessage * H.transpose(), 2);
		MatrixXi message = encodedMessage;

		int errorPos = 0;
		for (int j = 0; j < syndrome.cols(); j++) {
			errorPos += syndrome(0, j) << j;
		}

		/*cout << "Syndrome as binary: ";
		for (int i = 0; i < syndrome.cols(); ++i) {
			cout << syndrome(0, i) << " ";
		}
		cout << "\nCalculated error position : " << errorPos << endl;*/
		


		if (errorPos > 0) {
			message(0, errorPos - 1) ^= 1;
			//cout << "Error detected and corrected at position " << errorPos << endl;
		}
		else {
			//cout << "No errors detected." << endl;
		}

		MatrixXi originalData = message.block(0, m, 1, k);
		//cout << "Decoded data: ";
		//for (int i = 0; i < originalData.cols(); ++i) {
		//	cout << originalData(0, i) << " ";
		//}
		//cout << endl;

		return originalData;
	}
};
class Modulator {
private:
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
	int M;
	double Tb;
	double Tbp;
	double fn;
	double A1 = 1.0 / 2.0;
	double A2 = 2.0;
	vector<int> B;
	complex<double>* a = new complex<double>[N];
	complex<double>* b = new complex<double>[N];
	public:
		Modulator(vector<int> B) {
			this->B = B;
			M = w * B.size();
			Tb = (double)Tc / M;
			Tbp = (double)Tb * fs;
			fn = W * 1 / Tb;
		}
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
			vector<int> tempB = B;
			int z = (int)Tbp;
			int bSize = tempB.size();


			for (int i = 0; i < bSize; i++) {
				int start = i * z;
				int end = start + z;
				for (int j = start; j < end; j++) {
					if (tempB[i] == 0) {

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
		//mode = 0 dla ASK i PSK, mode = 1 dla FSK fn1, mode = 2 dla FSK fn2
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
		//calculaetSum zwraca wektor wektorow, gdzie pierwszy wektor to czas, a drugi to sygnal dla funkcji podanej jako argument
		vector<vector<double>> calculateSum(vector<vector<double>> data) {
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
		//calculateSum zwraca wektor wektorow, gdzie pierwszy wektor to czas, a drugi to sygnal dla funkcji FSK, je�li parametrem jest tryb, a nie funkcja. Tryb 1 dla fn1 , tryb 2 dla fn2
		vector<vector<double>> calculateSum(vector<vector<double>> data, int mode) {
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
		//mode 0 - wartosc progowa 0, p(t) > 0(FSK),		mode 1 - warto�� progowa 0, p(t) < 0(PSK),		mode 2 - warto�� progowa h, p(t) > h(ASK)
		vector<int> binarizeData(vector<double> x, int mode) {
			vector<int> result(x.size());
			double trigger;

			switch (mode) {
			case 0: {
				trigger = 0;
				for (int i = 0; i < x.size(); i++)
					result[i] = (x[i] > trigger) ? 1 : 0;
				break;
			}
			case 1: {
				trigger = 0;
				for (int i = 0; i < x.size(); i++)
					result[i] = (x[i] < trigger) ? 1 : 0;
				break;
			}

			case 2: {
				trigger = (*max_element(x.begin(), x.end()) / 2) - 1;
				for (int i = 0; i < x.size(); i++)
					result[i] = (x[i] > trigger) ? 1 : 0;
				break;
			}
			}
			return result;


		}
		vector<vector<double>> combineFSK(vector<vector<double>> data) {
			vector<vector<double>> fsk1 = calculateSum(data, 1);
			vector<vector<double>> fsk2 = calculateSum(data, 2);
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
			vector<int> result(B.size());
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
					result[i] = 1;
				}
				else {
					result[i] = 0;
				}
			}
			return result;
		}
		//mode = 0 dla ASK, mode = 1 dla PSK, mode = 2 dla FSK
		vector<vector<double>> modulateSignal(int mode) {
			vector<vector<double>> result;
			switch (mode) {
			case 0: {
					result = ASK();
					break;
				}
			case 1: {
					result = PSK();
					break;
				}
			case 2: {
					result = FSK();
					break;
				}
			}
			return result;
		}

		//mode = 1 dla FSK fn1, mode = 2 dla FSK fn2
		vector<int> demodulateFSKSignal(vector<vector<double>> data) {
				vector<vector<double>> fsk = combineFSK(data);
				vector<int> result = fitToBits(fsk[1], 0);
				return result;
		}
		//mode 1 - wartosc progowa 0, p(t) < 0(PSK),	mode 2 - wartosc progowa h, p(t) > h(ASK)
		vector<int> demodulateSignal(vector<vector<double>> data, int mode) {
			data = calculateSum(data);
			vector<int> result = fitToBits(data[1], mode);
			return result;
		}
};


void compareSignals(vector<vector<int>> x, vector<int> y) {
	int errors = 0;
	vector<int> temp;
	for (int i = 0; i < x.size(); i++) {
		for (int j : x[i]) {
			temp.push_back(j);
		}
	}
	for (int i = 0; i < temp.size(); i++) {
		if (temp[i] != y[i]) errors++;
	}
	cout << "Errors: " << errors << endl;
}
vector<int> encodeFrames(vector<vector<int>> frames) {
	vector<int> result;
	for (int i = 0; i < frames.size(); i++) {
		MatrixXi data(1, 11);
		for (int j = 0; j < frames[i].size(); j++) {
			data(0, j) = frames[i][j];
		}
		Hamming hamming;
		MatrixXi encoded = hamming.encode(data);
		for (int j = 0; j < encoded.cols(); j++) {
			result.push_back(encoded(0, j));
		}
	}
	return result;
}


int main() {
	vector<int> message = { 1,1,0,0,0,0,1,0,1,0,1,0,0,0,0,1,1,1,1,0,0,0,0,1,0,0,0,1,1,1,1,1,0,0,1,0,1,1,1,1,1,0,0,0 };
	if (message.size() % 11 != 0) { cout << "B.size() nie jest podzielne przez 11" << endl; return -1; }
	int frameSize = 11;
	vector<vector<int>> frames;
	for (int i = 0; i < message.size(); i += frameSize) {
		vector<int> frame;
		for (int j = i; j < i + frameSize; j++) {
			frame.push_back(message[j]);
		}
		frames.push_back(frame);
	}
	vector<vector<int>> result;
	for (int i = 0; i < frames.size(); i++) {
		vector<int> encodedVec;
		MatrixXi data(1, 11);
		for (int j = 0; j < frames[i].size(); j++) {
			data(0, j) = frames[i][j];
		}
		Hamming hamming;
		MatrixXi encoded = hamming.encode(data);
		for (int j = 0; j < encoded.cols(); j++) {
			encodedVec.push_back(encoded(0, j));
		}

		Modulator mod(encodedVec);
		//test dla ASK
		vector<vector<double>> signal = mod.modulateSignal(0);
		vector<int> demodulated = mod.demodulateSignal(signal, 2);
		//cout << "demodulated: " << endl;
		//for (int j = 0; j < demodulated.size(); j++) {
		//	cout << demodulated[j] << " ";
		//}
		MatrixXi demodulatedMatrix(1, 15);
		for (int j = 0; j < demodulated.size(); j++) {
			demodulatedMatrix(0, j) = demodulated[j];
		}
		vector<int> demodulatedVec;
		MatrixXi decoded = hamming.decode(demodulatedMatrix);
		for (int j = 0; j < decoded.cols(); j++) {
			demodulatedVec.push_back(decoded(0, j));
		}
		result.push_back(demodulatedVec);
	}
	compareSignals(result, message);
	//vector<int> encoded = encodeFrames(frames);
	//cout << "Encoded Size: " << encoded.size() << endl;
	//for (int i = 0; i < encoded.size(); i++) {
	//	cout << encoded[i] << " ";
	//}


	return 0;
}

