#include "header.h"

/*qr-алгоритм для нахождения спектра матрицы*/

void householder(vector<vector<double> > a, vector<vector<double> > &h, vector<vector<double> > &r) {
    int s = a.size();

    vector<double> w(s, 0), c(s, 0);
    vector<vector<double> > h_final(s, c);
    vector<vector<vector<double> > > h_cur(s, a);

    double norm_s = 0;

    for (int i = 0; i < s; ++i) {
        for (int j = 0; j < s; ++j) {
            h[i][j] = int(i == j);
        }
    }

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

        square_matrix_mult(h, h, h_cur[i]);
    }

    r = a;
}

double error(vector<vector<double> > a) {
    double sum = 0;
    int s = a.size();

    for (int i = 1; i < s; ++i) {
        for (int j = 0; j < i; ++j) {
            sum += a[i][j] * a[i][j];
        }
    }

    return sqrt(sum);
}

void qr_algorithm() {
    unsigned int n, i;

    double delta = 0, err;

    ifstream fin;
    ofstream fout;
    fin.open("input.txt", fstream::in);
    fout.open("output.txt", fstream::out);

    fout.precision(8);
    fixed(fout);

    fin >> n >> delta;

    vector<double> c(n, 0);
    vector<vector<double> > lambda(n, c);
    vector<vector<vector<double> > > a(1, lambda), h(1, lambda), r(1, lambda);

    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            fin >> a[0][k][j];
        }
    }

    householder(a[0], h[0], r[0]);

    err = error(a[0]);

    for (i = 1; err > delta; ++i) {
        a.push_back(lambda);
        h.push_back(lambda);
        r.push_back(lambda);

        square_matrix_mult(a[i], r[i - 1], h[i - 1]);

        householder(a[i], h[i], r[i]);

        err = error(a[i]);
    }

    fout << --i << endl << endl;

    for (int j = 0; j < n; ++j) {
        fout << a[i][j][j] << ' ';
    }

    fout << endl << endl;

    display_matrix(a[i], fout);
}