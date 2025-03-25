#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

double f(double x) {
    return 2 * x * sin(x);
}

double right_rectangle_integral(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.0;

    for (int i = 1; i <= n; ++i) {
        double x_i = a + i * h;
        sum += f(x_i);
    }

    return sum * h;
}

int main() {
    double a = 0, b = 5;
    std::ofstream out_file("results.txt");

    out_file << std::fixed << std::setprecision(12);
    out_file << "n\tresult\n";

    for (int n = 100; n <= 10'000'000; n *= 3) {
        double result = right_rectangle_integral(a, b, n);
        out_file << n << "\t" << result << "\n";
    }
    std::cout << std::setprecision(10) << -(2 * sin(0) - 2 * 0 * cos(0)) + (2 * sin(5) - 2 * 5 * cos(5))<< std::endl;

    out_file.close();
    return 0;
}