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
	int x1, x2, x4, x8;


	MatrixXi I = Matrix<int, k, k>::Identity();
	MatrixXi P;
	Matrix<int, k, n> G;
	MatrixXi data = MatrixXi::Zero(1,k);
	MatrixXi H;

	void removeRow(MatrixXi& matrix, unsigned int rowToRemove) {
		unsigned int numRows = matrix.rows() - 1;
		unsigned int numCols = matrix.cols();

		if (rowToRemove < numRows)
			matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);

		matrix.conservativeResize(numRows, numCols);
	}
	void generateH() {
		MatrixXi i = Matrix<int, n-k,n-k>::Identity();
		MatrixXi pt = P.transpose();
		H.resize(m, n);
		H << i,pt;
	}
	void generateG() {
		G.resize(k, n);
		G << P, I;

	}


	void fillP() {
		P = Matrix<int,n+1,m>::Zero(n + 1, m);
		for (int i = 0; i < n + 1; i++) {
			for (int j = 0; j < m; j++) {
				P(i, j) = (i >> j) & 1;
			}
		}
		removeRow(P, 0);
		for (int i = 1,  delCount = 1; i < n + 1; i *= 2) {
			int idx = i - delCount;
			removeRow(P, i - delCount);
			delCount++;
		}
	}


	Hamming(MatrixXi data) {
		fillP();
		generateG();
		generateH();
		this->data = data;
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

	MatrixXi encode() {
		Matrix<int, 1, n> result = (data * G);
		Matrix<int, 1, n> modArr = moduloMatrix(result, 2);
		cout << "Data: " << data << endl;
		cout << "Encoded data: " << modArr << endl;
		return modArr;
	}
	MatrixXi decode(MatrixXi& encodedMessage) {
	
		MatrixXi syndrome = moduloMatrix(encodedMessage * H.transpose(), 2);
		MatrixXi message = encodedMessage;

		int errorPos = 0;
		for (int j = 0; j < syndrome.cols(); j++) {
			errorPos += syndrome(0, j) << j;
		}

		cout << "Syndrome as binary: ";
		for (int i = 0; i < syndrome.cols(); ++i) {
			cout << syndrome(0, i) << " ";
		}
		cout << "\nCalculated error position : " << errorPos << endl;


		if (errorPos > 0) {
			message(0, errorPos - 1) ^= 1;
			cout << "Error detected and corrected at position " << errorPos << endl;
		}
		else {
			cout << "No errors detected." << endl;
		}

		MatrixXi originalData = message.block(0, m, 1, k);
		cout << "Decoded data: ";
		for (int i = 0; i < originalData.cols(); ++i) {
			cout << originalData(0, i) << " ";
		}
		cout << endl;

		return originalData;
	}
};


int main() {
	vector<bool> input = { 1, 1, 0, 1 };
	vector<bool> output = encoder(input);
	cout << endl;
	output[5] = !output[5];
	decoder(output);

	MatrixXi data(1,11);
	data << 1,1,0,0,0,0,1,0,0,1,0;
	Hamming h(data);
	//h.printI();
	//h.fillP();
	//h.printP();
	//h.generateG();
	//h.printG();
	MatrixXi result(1, 15);
	cout << endl;
	result = h.encode();

	//cout<< "result: "<<endl;
	//cout << result << endl;
	result(0, 3) ^= 1;
	h.decode(result);
	//cout << "Simulating error: " << endl;
	//for (int i = 0; i < h.n;) {
	//	cout <<endl << "Error at index " << i << endl;
	//	cout << "Result without error: " << endl;
	//	cout << result << endl;
	//	result(0, i) ^= 1;
	//	cout << "Result with error: " << endl;
	//	cout << result << endl;
	//	h.decode(result);
	//	result(0, i) ^= 1;
	//	i++;
	//}
	//h.printP();
	//h.printG();
	//h.printH();
	//cout << "Result: " << endl;
	/*cout << result << endl;
	MatrixXi decoded = h.decode(result);
	cout << "Decoded data: " << endl;
	cout << decoded << endl;
	cout << "Data: " << endl;
	 cout << data << endl;*/
	//cout << endl;
	//h.eye();
	//h.print(h.I);
	return 0;
}