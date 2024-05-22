#include <vector>
#include <iostream>
using namespace std;

vector<bool> coder(vector<bool> input) {
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
	bool x1p, x2p, x4p,x1s,x2s,x4s;
	x1p = input[2] ^ input[4] ^ input[6];
	x2p = input[2] ^ input[5] ^ input[6];
	x4p = input[4] ^ input[5] ^ input[6];
	x1s = input[0] ^ x1p;
	x2s = input[1] ^ x2p;
	x4s = input[3] ^ x4p;
	int s = x1s + x2s*2 + x4s*4;
	if (s == 0) {
		cout << "Output:" << endl;
		for (auto b : temp) {
			cout << b << " ";
		}
		cout << "No errors detected" << endl;
	}
	else {
		cout<<"Output:"<<endl;
		for (auto b : temp) {
			cout << b << " ";
		}

		cout <<endl<< "Error detected at position " << s << endl;
		cout << "Corrected output: ";
		temp[s-1] = !temp[s-1];
		for (auto b : temp) {
			cout << b << " ";
		}
	}
	return temp;
}

int main() {
	vector<bool> input = {1, 1, 0, 1};
	vector<bool> output = coder(input);
	cout << endl;
	output[5] = !output[5];
	decoder(output);
	return 0;
}