#include "header.h"

/* трёхчленная формула метода Ричардсона с чебышевскими параметрами
 A - симметричная положительно определённая матрица, 0 < alpha < beta - границы спектра, delta - требуемая норма вектора невязки
 b - вектор, x - начальный вектор*/

void richardson() {
    unsigned int n, i;

    double alpha = 0, beta = 0, delta = 0, d = 10;

    ifstream fin;
    ofstream fout;
    fin.open("input.txt", fstream::in);
    fout.open("output.txt", fstream::out);

    fout.precision(8);
    fixed(fout);

    fin >> n >> alpha >> beta >> delta;

    vector<double> b(n, 0), c(n, 0), omega(2, 0), temp1(n, 0), temp2(n, 0);
    vector<vector<double> > a(n, c), x(2, c);

    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            fin >> a[k][j];
        }
    }

    for (int k = 0; k < n; ++k) {
        fin >> b[k];
    }

    for (int k = 0; k < n; ++k) {
        fin >> x[0][k];
    }

    vector<double> b1 = b;

    matrix_vector_mult(x[1], a, x[0]);
    gauss_scalar_vector_mult(-1, b);
    gauss_vector_sum(x[1], b);
    gauss_scalar_vector_mult(-2 / (beta + alpha), x[1]);
    gauss_vector_sum(x[1], x[0]);

    omega[1] = -1 * ((beta - alpha) / (beta + alpha));

    for (i = 1; d > delta; ++i) {
        x.push_back(c);
        omega.push_back(0);

        omega[i + 1] = 1 / ((2 / omega[1]) - omega[i]);

        temp1 = x[i];
        gauss_scalar_vector_mult(-1, x[i - 1]);
        gauss_vector_sum(temp1, x[i - 1]);
        gauss_scalar_vector_mult(omega[i] * omega[i + 1], temp1);

        temp2 = x[i];
        matrix_vector_mult(temp2, a, temp2);
        gauss_vector_sum(temp2, b);
        gauss_scalar_vector_mult(((-2 * (1 + omega[i] * omega[i + 1])) / (beta + alpha)), temp2);

        gauss_vector_sum(x[i + 1], x[i]);
        gauss_vector_sum(x[i + 1], temp1);
        gauss_vector_sum(x[i + 1], temp2);

        d = norm_residual(a, x[i], b1);
    }

    fout << --i << endl;
    display_vector(x[i], fout);
}