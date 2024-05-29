#include <vector>
#include <iostream>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;
vector<bool> encoder(vector<bool> input) {
	if (input.size() != 4) {
		throw std::invalid_argument("received invalid sized vector");
	}
	vector<bool> output(7);
	bool x1, x2, x4;
	output[2] = input[0];
	output[4] = input[1];
	output[5] = input[2];
	output[6] = input[3];
	x1 = output[2] ^ output[4] ^ output[6];
	x2 = output[2] ^ output[5] ^ output[6];
	x4 = output[4] ^ output[5] ^ output[6];
	output[0] = x1;
	output[1] = x2;
	output[3] = x4;

	for (auto b : output) {
		cout << b << " ";
	}
	return output;
}

vector<bool> decoder(vector<bool> input) {
	if (input.size() != 7) {
		throw std::invalid_argument("received invalid sized vector");
	}
	vector<bool> temp = input;
	bool x1p, x2p, x4p, x1s, x2s, x4s;
	x1p = input[2] ^ input[4] ^ input[6];
	x2p = input[2] ^ input[5] ^ input[6];
	x4p = input[4] ^ input[5] ^ input[6];
	x1s = input[0] ^ x1p;
	x2s = input[1] ^ x2p;
	x4s = input[3] ^ x4p;
	int s = x1s + x2s * 2 + x4s * 4;
	if (s == 0) {
		cout << "Output:" << endl;
		for (auto b : temp) {
			cout << b << " ";
		}
		cout << "No errors detected" << endl;
	}
	else {
		cout << "Output:" << endl;
		for (auto b : temp) {
			cout << b << " ";
		}

		cout << endl << "Error detected at position " << s << endl;
		cout << "Corrected output: ";
		temp[s - 1] = !temp[s - 1];
		for (auto b : temp) {
			cout << b << " ";
		}
		cout << endl;
	}
	return temp;
}
class Hamming {
public:
	static const int k = 11, n = 15;
	static const int m = n - k;

	MatrixXi I = Matrix<int, k, k>::Identity();
	Matrix<int, k, m> P;
	Matrix<int, k, n> G;
	MatrixXi data = MatrixXi::Zero(1,k);
	void generateG() {
		G << I, P;
	}
	void fillP() {
		P << 0, 0, 1, 1,
			0, 1, 0, 1,
			0, 1, 1, 0,
			0, 1, 1, 1,
			1, 0, 0, 1,
			1, 0, 1, 0,
			1, 0, 1, 1,
			1, 1, 0, 0,
			1, 1, 0, 1,
			1, 1, 1, 0,
			1, 1, 1, 1;
	}

	Hamming(MatrixXi data) {
		fillP();
		generateG();
		this->data = data;
	}
	//mgm::mat<k, k, int> I;
	//mgm::mat<k, m, int> P;
	//mgm::mat<11, 15, int> G;


	void printI() {
		cout << I << endl;
	}
	void printP() {
		cout << P << endl;
	}
	void printG() {
		cout << G << endl;
	}


	MatrixXi encode() {
		MatrixXi result = MatrixXi::Zero(1, n);
		cout << "Data: " << data << endl;
		cout << "G: " << G << endl;
		result = data * G;
		for (int i = 0; i < n; i++) {
			result(0, i) = result(0, i) % 2;
		}
		cout << "Encoded data: " << result << endl;
		return result;
	}
};


int main() {
	vector<bool> input = { 1, 1, 0, 1 };
	vector<bool> output = encoder(input);
	cout << endl;
	output[5] = !output[5];
	decoder(output);

	MatrixXi data(1,11);
	data << 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1;
	Hamming h(data);
	//h.printI();
	//h.fillP();
	//h.printP();
	//h.generateG();
	//h.printG();
	h.encode();
	//cout << endl;
	//h.eye();
	//h.print(h.I);
	return 0;
}