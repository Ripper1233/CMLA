#include "header.h"

/*метод Гаусса: решает систему Ax=b*/

void gauss1() {
    unsigned int n;

    ifstream fin;
    ofstream fout;
    fin.open("input.txt", fstream::in);
    fout.open("output.txt", fstream::out);

    fin >> n;
    vector<double> c(n + 1, 0);
    vector<vector<double> > a(n, c);
    double div;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= n; ++j) {
            fin >> a[i][j];
        }
    }

    for (int i = 0; i < n; ++i) {
        div = 1 / a[i][i];
        gauss_scalar_vector_mult(div, a[i]);
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                c = a[i];
                gauss_scalar_vector_mult(-1 * a[j][i], c);
                gauss_vector_sum(a[j], c);
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        fout << a[i][n] << ' ';
    }
}
