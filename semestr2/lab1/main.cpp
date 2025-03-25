#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#define M_PI 3.141592653

using namespace std;

//Гладкая функция
double f1(double x) { return exp(-pow(x, 2)); }
double df1(double x) { return -2 * x * exp(-pow(x, 2)); }

//Разрывная функция 
double f2(double x) { return fabs(x); }
double df2(double x) { return (x > 0) ? 1 : ((x < 0) ? -1 : 0); }

// равномерная сетка
vector<double> uniform_grid(int k, double a, double b) {
    vector<double> x_h(k);
    double h = (b - a) / (k - 1);
    for (int i = 0; i < k; i++) x_h[i] = a + i * h;
    return x_h;
}

// Чебышевская сетка
vector<double> chebyshev_grid(int k, double a, double b) {
    vector<double> x_h(k);
    for (int i = 0; i < k; i++)
        x_h[i] = 0.5 * (a + b) + 0.5 * (b - a) * cos((2 * i + 1) * M_PI / (2 * k));
    return x_h;
}

// вычисление коэффициентов полинома 
vector<vector<double>> hermite_coefficients(const vector<double>& x_h, double (*f)(double), double (*df)(double)) {
    int n = 2 * x_h.size();
    vector<vector<double>> table(n, vector<double>(n, 0));
    for (int i = 0; i < x_h.size(); i++) {
        int j = 2 * i;
        table[j][0] = f(x_h[i]);
        table[j + 1][0] = f(x_h[i]);
        table[j + 1][1] = df(x_h[i]);
        if (j > 0) table[j][1] = (table[j][0] - table[j - 1][0]) / (x_h[i] - x_h[i - 1]);
    }
    for (int col = 2; col < n; col++)
        for (int row = col; row < n; row++)
            table[row][col] = (table[row][col - 1] - table[row - 1][col - 1]) / (x_h[row / 2] - x_h[(row - col) / 2]);
    return table;
}

// вычисление полинома в точке x
double hermite_polynomial(double x, const vector<double>& x_h, const vector<vector<double>>& table) {
    double result = 0, term = 1;
    int n = 2 * x_h.size();
    for (int i = 0; i < n; i++) {
        result += table[i][i] * term;
        term *= (x - x_h[i / 2]);
    }
    return result;
}

double compute_error(const vector<double>& x_h, const vector<vector<double>>& table, double (*f)(double), double a, double b, int num_points) {
    double max_error = 0, step = (b - a) / (num_points - 1);
    for (int i = 0; i < num_points; i++) {
        double x = a + i * step;
        double error = abs(f(x) - hermite_polynomial(x, x_h, table));
        max_error = max(max_error, error);
    }
    return max_error;
}

int main() {
    double a = -2, b = 2;
    vector<int> k_values = { 5, 7, 9, 11 };  
    ofstream file("results.txt"); 

    if (!file) {
        cerr << "err" << endl;
        return 1;
    }

    file << "k сетка функция ошибка\n";

    for (int k : k_values) {
        vector<double> x_uniform_f1 = uniform_grid(k, a, b);
        vector<vector<double>> coeffs_uniform_f1 = hermite_coefficients(x_uniform_f1, f1, df1);
        double error_uniform_f1 = compute_error(x_uniform_f1, coeffs_uniform_f1, f1, a, b, 100);

        vector<double> x_chebyshev_f1 = chebyshev_grid(k, a, b);
        vector<vector<double>> coeffs_chebyshev_f1 = hermite_coefficients(x_chebyshev_f1, f1, df1);
        double error_chebyshev_f1 = compute_error(x_chebyshev_f1, coeffs_chebyshev_f1, f1, a, b, 100);

        file << k << " uniform exp(-x^2) " << error_uniform_f1 << "\n";
        file << k << " chebyshev exp(-x^2) " << error_chebyshev_f1 << "\n";

        vector<double> x_uniform_f2 = uniform_grid(k, a, b);
        vector<vector<double>> coeffs_uniform_f2 = hermite_coefficients(x_uniform_f2, f2, df2);
        double error_uniform_f2 = compute_error(x_uniform_f2, coeffs_uniform_f2, f2, a, b, 100);

        vector<double> x_chebyshev_f2 = chebyshev_grid(k, a, b);
        vector<vector<double>> coeffs_chebyshev_f2 = hermite_coefficients(x_chebyshev_f2, f2, df2);
        double error_chebyshev_f2 = compute_error(x_chebyshev_f2, coeffs_chebyshev_f2, f2, a, b, 100);

        file << k << " uniform |x| " << error_uniform_f2 << "\n";
        file << k << " chebyshev |x| " << error_chebyshev_f2 << "\n";
    }

    file.close();



    
    int k_demo = 5; 
    ofstream interp_file("interpolation_data.txt");
    if (!interp_file) {
        cerr << "err" << endl;
        return 1;
    }

    interp_file << "x y_real y_hermite сетка функция\n";

    vector<double> x_uniform_f1 = uniform_grid(k_demo, a, b);
    vector<vector<double>> coeffs_uniform_f1 = hermite_coefficients(x_uniform_f1, f1, df1);

    vector<double> x_chebyshev_f1 = chebyshev_grid(k_demo, a, b);
    vector<vector<double>> coeffs_chebyshev_f1 = hermite_coefficients(x_chebyshev_f1, f1, df1);

    vector<double> x_uniform_f2 = uniform_grid(k_demo, a, b);
    vector<vector<double>> coeffs_uniform_f2 = hermite_coefficients(x_uniform_f2, f2, df2);

    vector<double> x_chebyshev_f2 = chebyshev_grid(k_demo, a, b);
    vector<vector<double>> coeffs_chebyshev_f2 = hermite_coefficients(x_chebyshev_f2, f2, df2);

    int num_points = 100;
    double step = (b - a) / (num_points - 1);

    for (int i = 0; i < num_points; i++) {
        double x = a + i * step;

        double y_real_f1 = f1(x);
        double y_hermite_uniform_f1 = hermite_polynomial(x, x_uniform_f1, coeffs_uniform_f1);
        double y_hermite_chebyshev_f1 = hermite_polynomial(x, x_chebyshev_f1, coeffs_chebyshev_f1);

        interp_file << x << " " << y_real_f1 << " " << y_hermite_uniform_f1 << " uniform exp(-x^2)\n";
        interp_file << x << " " << y_real_f1 << " " << y_hermite_chebyshev_f1 << " chebyshev exp(-x^2)\n";

        double y_real_f2 = f2(x);
        double y_hermite_uniform_f2 = hermite_polynomial(x, x_uniform_f2, coeffs_uniform_f2);
        double y_hermite_chebyshev_f2 = hermite_polynomial(x, x_chebyshev_f2, coeffs_chebyshev_f2);

        interp_file << x << " " << y_real_f2 << " " << y_hermite_uniform_f2 << " uniform |x|\n";
        interp_file << x << " " << y_real_f2 << " " << y_hermite_chebyshev_f2 << " chebyshev |x|\n";
    }

    interp_file.close();
    return 0;
}
