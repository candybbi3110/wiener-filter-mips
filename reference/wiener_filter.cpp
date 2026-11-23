// C++ reference model for the MIPS Wiener filter.
// Mirrors the assembly step by step (same rounding points) so outputs can be
// compared value by value.
//
// Usage: wiener_ref [input.txt] [desired.txt] [output.txt]

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace std;

const int N = 10;  // signal length
const int M = 10;  // filter length

double input_signal[N], desired_signal[N], output_signal[N] = {0};
double R[M][M] = {0};
double optimize_coefficient[M] = {0};
double GXX[M] = {0}, GDX[M] = {0};

// Reads up to N values; size returns how many numbers the file really holds.
bool loadInputFile(const string& filename, double arr[], int& size) {
    ifstream fin(filename);
    if (!fin) return false;

    size = 0;
    double value;
    while (fin >> value) {
        if (size < N) arr[size] = value;
        size++;
    }
    return true;
}

// Biased finite-sample estimates of autocorrelation and cross-correlation.
void build_Correlation(int m) {
    double sum_gxx = 0.0, sum_gdx = 0.0;
    for (int n = m; n < N; n++) {
        sum_gxx += input_signal[n] * input_signal[n - m];
        sum_gdx += desired_signal[n] * input_signal[n - m];
    }
    GXX[m] = sum_gxx / N;
    GDX[m] = sum_gdx / N;
}

// Toeplitz autocorrelation matrix: R[i][j] = GXX[|i - j|].
void build_R_matrix() {
    for (int i = 0; i < M; i++)
        for (int j = 0; j < M; j++)
            R[i][j] = GXX[abs(i - j)];
}

// Solves R h = GDX with Gaussian elimination (no pivoting) + back substitution.
void calculate_H_matrix() {
    for (int k = 0; k < M; k++) {
        for (int i = k + 1; i < M; i++) {
            double factor = R[i][k] / R[k][k];
            for (int j = k; j < M; j++) R[i][j] -= factor * R[k][j];
            GDX[i] -= factor * GDX[k];
        }
    }
    for (int i = M - 1; i >= 0; i--) {
        double sum = GDX[i];
        for (int j = i + 1; j < M; j++) sum -= R[i][j] * optimize_coefficient[j];
        optimize_coefficient[i] = sum / R[i][i];
    }
}

double round_1_digit(double value) { return round(value * 10.0) / 10.0; }

// y[n] = sum_k h[k] * x[n - k], rounded after every accumulation (as in MIPS).
void calculate_output() {
    for (int n = 0; n < N; n++) {
        for (int k = 0; k <= n && k < M; k++) {
            output_signal[n] += optimize_coefficient[k] * input_signal[n - k];
            output_signal[n] = round_1_digit(output_signal[n]);
        }
    }
}

double calculate_MMSE() {
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        double err = round_1_digit(output_signal[i] - desired_signal[i]);
        sum = round_1_digit(sum + round_1_digit(err * err));
    }
    return round_1_digit(sum / N);
}

int main(int argc, char* argv[]) {
    string inputPath   = argc > 1 ? argv[1] : "input.txt";
    string desiredPath = argc > 2 ? argv[2] : "desired.txt";
    string outputPath  = argc > 3 ? argv[3] : "output.txt";

    int sizeInput = 0, sizeDesired = 0;
    if (!loadInputFile(inputPath, input_signal, sizeInput) ||
        !loadInputFile(desiredPath, desired_signal, sizeDesired)) {
        cout << "Error : Cannot Open File !\n";
        return 1;
    }
    if (sizeInput != N || sizeDesired != N) {
        cout << "Error : size not match \n";
        ofstream(outputPath) << "Error : size not match \n";
        return 1;
    }

    for (int m = 0; m < M; m++) build_Correlation(m);
    build_R_matrix();
    calculate_H_matrix();
    calculate_output();
    double mmse = calculate_MMSE();

    ofstream fout(outputPath);
    if (!fout) {
        cout << "Error : Cannot create output file!\n";
        return 1;
    }
    for (ostream* os : {static_cast<ostream*>(&cout), static_cast<ostream*>(&fout)}) {
        *os << fixed << setprecision(1) << "Filtered output: ";
        for (int n = 0; n < N; n++) *os << output_signal[n] << " ";
        *os << "\nMMSE: " << mmse << "\n";
    }
    return 0;
}
