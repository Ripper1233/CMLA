#include "header.h"

/*qr-разложение с использованием отражений Хаусхолдера*/

void householder() {
    int n;

    ifstream fin;
    ofstream fout;
    fin.open("input.txt", fstream::in);
    fout.open("output.txt", fstream::out);

    fin >> n;
    vector<double> b(n, 0), x(n, 0), w(n, 0), s(n, 0);
    vector<vector<double> > a(n, b), r(n, b);
    vector<vector<vector<double> > > h(n, a);
    double norm_s, sum;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            fin >> a[i][j];
        }
    }

    for (int i = 0; i < n; ++i) {
        fin >> b[i];
    }

    for (int i = 0; i < n - 1; ++i) {
        clear(s);

        for (int j = 0; j < n; ++j) {
            if (j < i) {
                s[j] = 0;
            } else {
                s[j] = a[j][i];
            }
        }

        norm_s = norm_2(s);

        if (norm_s == 0) {
            clear(w);
            householder_matrix(h[i], w);
        } else {
            for (int j = 0; j < n; ++j) {
                if (i == j) {
                    w[j] = s[j] + sgn_h(a[i][i]) * norm_s;
                } else {
                    w[j] = s[j];
                }
            }

            vector_norm(w);

            householder_matrix(h[i], w);
        }

        square_matrix_mult(a, h[i], a);
    }

    for (int i = 0; i < n - 1; ++i) {
        matrix_vector_mult(b, h[i], b);
    }

    for (int i = n - 1; i >= 0; --i) {
        sum = 0;

        for (int j = i + 1; j < n; ++j) {
            sum += a[i][j] * x[j];
        }

        x[i] = (b[i] - sum) / a[i][i];
    }

    fout.precision(8);
    fixed(fout);

    display_vector(x, fout);
}