#include "header.h"

/*метод итераций Рэлея*/

void householder_rayleigh(vector<vector<double> > a, vector<double> &x, vector<double> b) {
    int s = a.size();

    vector<double> w(s, 0), c(s, 0);
    vector<vector<double> > h_final(s, c);
    vector<vector<vector<double> > > h_cur(s, a);

    double norm_s = 0;

    for (int i = 0; i < s - 1; ++i) {
        clear(c);

        for (int j = 0; j < s; ++j) {
            if (j < i) {
                c[j] = 0;
            } else {
                c[j] = a[j][i];
            }
        }

        norm_s = norm_2(c);

        if (norm_s == 0) {
            clear(w);
            householder_matrix(h_cur[i], w);
        } else {
            for (int j = 0; j < s; ++j) {
                if (i == j) {
                    w[j] = c[j] + sgn_h(a[i][i]) * norm_s;
                } else {
                    w[j] = c[j];
                }
            }

            vector_norm(w);

            householder_matrix(h_cur[i], w);
        }

        square_matrix_mult(a, h_cur[i], a);
    }

    for (int i = 0; i < s - 1; ++i) {
        matrix_vector_mult(b, h_cur[i], b);
    }

    for (int i = s - 1; i >= 0; --i) {
        double sum = 0;

        for (int j = i + 1; j < s; ++j) {
            sum += a[i][j] * x[j];
        }

        x[i] = (b[i] - sum) / a[i][i];
    }
}

void inv_iter_rayleigh() {
    int n, i;

    double delta = 0, d = 10, lambda;

    ifstream fin;
    ofstream fout;
    fin.open("input.txt", fstream::in);
    fout.open("output.txt", fstream::out);

    fout.precision(8);
    fixed(fout);

    fin >> n >> delta >> lambda;

    vector<double> c(0, 0), cl(n, 0);
    vector<vector<double> > a(n, cl), x(1, cl);

    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            fin >> a[k][j];
        }
    }

    for (int k = 0; k < n; ++k) {
        fin >> x[0][k];
    }

    vector_norm(x[0]);

    vector<vector<double> > a_cur(a);

    for (i = 0; d > delta; ++i) {
        for (int k = 0; k < n; ++k) {
            a_cur[k][k] = a[k][k] - lambda;
        }

        x.push_back(cl);
        householder_rayleigh(a_cur, x[i + 1], x[i]);

        c.push_back(scalar_product(x[i + 1], x[i]) / scalar_product(x[i], x[i]));

        vector_norm(x[i + 1]);

        if (i > 0)
            d = abs(c[i] - c[i - 1]);
    }
    --i;
    fout << i << ' ' << (1 / c[i]) + lambda;
}