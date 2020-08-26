#include "header.h"

/*метод Гаусса: решает систему Ax=b, находит детерминант A*/

void gauss3() {
    int n;
    double sum, det = 1;

    ifstream fin;
    ofstream fout;
    fin.open("input.txt", fstream::in);
    fout.open("output.txt", fstream::out);

    fin >> n;
    vector<double> c(n + 1, 0), x(n, 0);
    vector<vector<double> > a(n, c);
    double div;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= n; ++j) {
            fin >> a[i][j];
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
        sum = 0;

        for (int j = i + 1; j < n; ++j) {
            sum += a[i][j] * x[j];
        }

        x[i] = a[i][n] - sum;
    }

    fout.precision(8);
    fixed(fout);

    display_vector(x, fout);

    fout << det << endl;
}