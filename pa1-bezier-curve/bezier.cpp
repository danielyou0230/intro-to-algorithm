#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

typedef struct 
{
	int ini;
	int end;
} parameters;

double *sample_x;
double *sample_y;
int n_pts;
double step;
double *input_x;
double *input_y;
double **binCoeff;

int main(int argc, char const *argv[])
{
	ifstream f_in;
	f_in.open(argv[1], ifstream::in);
	// Read number of control points in the curve
	// int n_pts;
	f_in >> n_pts;
	cout << "Total number of points: " << n_pts << endl;
	// Binomial Coefficient tree (Dynamic Programming)
	// double **binCoeff = new double*[n_pts];
	binCoeff = new double*[n_pts];
	binCoeff[0] = new double; 
	binCoeff[0][0] = 1; // Root
	for (int i = 1; i < n_pts; ++i) {
		binCoeff[i] = new double[i + 1]; // bottom-up
		for (int j = 0; j < i + 1; ++j) {
			if (j == 0 || j == i) // Rightmost and Leftmost nodes in Binomial Coefficient tree
				binCoeff[i][j] = 1;
			else // Sum of two parents 
				binCoeff[i][j] = binCoeff[i - 1][j - 1] + binCoeff[i - 1][j];
		}
	}
	// Read control points from file
	input_x = new double[n_pts];
	input_y = new double[n_pts];
	for (int i = 0; i < n_pts; ++i) {
		f_in >> input_x[i] >> input_y[i];
	}
	// Number of sample points
	int n_sample;
	f_in >> n_sample;
	cout << "Sample points: " << n_sample << endl;
	f_in.close();
	// Show 
	/* cout << "Binomial Coefficient: " << endl;
	for (int i = 0; i < n_pts; ++i) {
		for (int j = 0; j < i + 1; ++j)
			cout << binCoeff[i][j] << " ";
		cout << endl;
	} */
	// Allocate memory for output points
	sample_x = new double[n_sample]();
	sample_y = new double[n_sample]();
	// Sample resolution (time)
	// double step = 1.0 / (n_sample - 1);
	step = 1.0 / (n_sample - 1);
	// Calculate coordinates for each sample points
	for (int i = 0; i < n_sample; ++i) {
		// Substitute numerical values for (1 - t) and t
		for (int k = 0; k < n_pts; ++k) {
			// Sum of power of (1 - t) and t equals to (n - 1)
			// Polynomial = (1 - t) ^ (n - 1 - k) * t ^ (k)
			// e.g. n_pts = 5, input nodes = {a, b, c, d, e}
			//      Let x = (1 - t), y = t
			//      Curve = Binomial coefficient * input node * numerical t
			//            = 1 a (x ^ 4) (y ^ 0)
			//            + 4 b (x ^ 3) (y ^ 1) 
			//            + 6 c (x ^ 2) (y ^ 2) 
			//            + 4 d (x ^ 1) (y ^ 3) 
			//            + 1 e (x ^ 0) (y ^ 4) 
			double tmp = binCoeff[n_pts - 1][k] * pow(1 - step * i, n_pts - 1 - k) * pow(step * i, k);
			// Add up all control points base on it weights (numerical value of t)
			sample_x[i] += input_x[k] * tmp;
			sample_y[i] += input_y[k] * tmp;
		}
	}
	// Output result
	ofstream f_out;
	f_out.open (argv[2], ios::out);
	// cout << "Algorithm Elapsed Time: " << (double)clock() / CLOCKS_PER_SEC << "s" << endl;
	for (int i = 0; i < n_sample; ++i) {
		f_out << setprecision(2) << fixed;
		f_out << sample_x[i] << "\t" << sample_y[i] << endl;
	}
	f_out.close();
	// cout << "Execution Time        : " << (double)clock() / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}
