#include "header.h"

/*метод Гаусса: решает систему Ax=b, находит детерминант A и обратную матрицу к A*/

void gauss2() {
    unsigned int n;
    double det = 1;

    ifstream fin;
    ofstream fout;
    fin.open("input.txt", fstream::in);
    fout.open("output.txt", fstream::out);

    fin >> n;
    vector<double> c(2 * n + 1, 0);
    vector<vector<double> > a(n, c);
    double div;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= n; ++j) {
            fin >> a[i][j];
        }
    }

    for (int i = 0; i < n; ++i) {
        for (int j = n + 1; j <= 2 * n; ++j) {
            if (j - n - 1 == i) {
                a[i][j] = 1;
            } else {
                a[i][j] = 0;
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        det *= a[i][i];
        div = 1 / a[i][i];
        gauss_scalar_vector_mult(div, a[i]);
        for (int j = i + 1; j < n; ++j) {
            c = a[i];
            gauss_scalar_vector_mult(-1 * a[j][i], c);
            gauss_vector_sum(a[j], c);
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        for (int j = 0; j < i; ++j) {
            c = a[i];
            gauss_scalar_vector_mult(-1 * a[j][i], c);
            gauss_vector_sum(a[j], c);
        }
    }

    fout.precision(8);
    fixed(fout);

    for (int i = 0; i < n; ++i) {
        fout << a[i][n] << ' ';
    }

    fout << endl << det << endl;

    for (int i = 0; i < n; ++i) {
        for (int j = n + 1; j <= 2 * n; ++j) {
            fout << a[i][j] << ' ';
        }
        fout << endl;
    }
}
