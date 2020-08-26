#include "header.h"

/*метод прогонки для трёхдиагональных матриц*/

void tridiagonal_matrix_algorithm() {
    unsigned int n;

    ifstream fin;
    ofstream fout;
    fin.open("input.txt", fstream::in);
    fout.open("output.txt", fstream::out);

    fin >> n;
    vector<double> c(n + 1, 0), x(n, 0), p(n, 0), q(n, 0);
    vector<vector<double> > a(n, c);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= n; ++j) {
            fin >> a[i][j];
        }
    }

    p[0] = a[0][1] / a[0][0];
    q[0] = a[0][n] / a[0][0];

    for (int i = 1; i < n - 1; ++i) {
        p[i] = a[i][i + 1] / (a[i][i] - p[i - 1] * a[i][i - 1]);
    }

    for (int i = 1; i < n; ++i) {
        q[i] = (a[i][n] - q[i - 1] * a[i][i - 1]) / (a[i][i] - p[i - 1] * a[i][i - 1]);
    }

    x[n - 1] = q[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = q[i] - p[i] * x[i + 1];
    }

    display_vector(x, fout);
}