#include "header.h"

/*метод сопряжённых градиентов для решения СЛАУ
 A - симметричная положительно определённая матрица, 0 < alpha < beta - границы спектра, delta - требуемая норма вектора невязки
 b - вектор, x - начальный вектор*/

double calc_gamma(vector<vector<double> > a, vector<double> r, vector<double> g) {
    unsigned int s = r.size();

    vector<double> temp(s, 0);

    matrix_vector_mult(temp, a, g);

    return (scalar_product(r, temp)) / (scalar_product(temp, g));
}

double calc_alpha(vector<vector<double> > a, vector<double> r, vector<double> g) {
    unsigned int s = r.size();

    vector<double> temp(s, 0);

    matrix_vector_mult(temp, a, g);

    return (scalar_product(r, g)) / (scalar_product(temp, g));
}

void calc_g(vector<double> &g_result, vector<double> r, double gamma, vector<double> g) {
    gauss_scalar_vector_mult(-1 * gamma, g);

    gauss_vector_sum(r, g);

    g_result = r;
}

void calc_x(vector<double> &x_result, vector<double> x, double alpha, vector<double> g) {
    gauss_scalar_vector_mult(-1 * alpha, g);

    gauss_vector_sum(x, g);

    x_result = x;
}

void calc_r(vector<double> &r_result, vector<double> r, double alpha, vector<vector<double> > a, vector<double> g) {
    matrix_vector_mult(g, a, g);

    gauss_scalar_vector_mult(-1 * alpha, g);

    gauss_vector_sum(r, g);

    r_result = r;
}

void residual(vector<double> &r, vector<vector<double> > a, vector<double> x, vector<double> b) {
    matrix_vector_mult(r, a, x);

    gauss_scalar_vector_mult(-1, b);

    gauss_vector_sum(r, b);
}

void conjugate_gradient_method() {
    unsigned int n, i = 2;

    double delta = 0, d;

    ifstream fin;
    ofstream fout;
    fin.open("input.txt", fstream::in);
    fout.open("output.txt", fstream::out);

    fout.precision(8);
    fixed(fout);

    fin >> n >> delta;

    vector<double> b(n, 0), c(n, 0), alpha(2, 0), gamma(1, 0), temp(n, 0);
    vector<vector<double> > a(n, c), x(2, c), g(2, c), r(2, c);

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

    residual(r[0], a, x[0], b);

    g[1] = r[0];

    alpha[1] = calc_alpha(a, r[0], g[1]);

    calc_x(x[1], x[0], alpha[1], g[1]);

    calc_r(r[1], r[0], alpha[1], a, g[1]);

    d = norm_2(r[1]);

    if (d > delta) {
        for (i = 2; d > delta; ++i) {
            x.push_back(c);
            r.push_back(c);
            g.push_back(c);

            gamma.push_back(calc_gamma(a, r[i - 1], g[i - 1]));
            calc_g(g[i], r[i - 1], gamma[i - 1], g[i - 1]);
            alpha.push_back(calc_alpha(a, r[i - 1], g[i]));
            calc_x(x[i], x[i - 1], alpha[i], g[i]);
            calc_r(r[i], r[i - 1], alpha[i], a, g[i]);

            d = norm_2(r[i]);
        }
    }
    fout << --i << endl;
    display_vector(x[i], fout);
}