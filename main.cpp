#include "header.h"

void gauss_scalar_vector_mult(double x, vector<double> &y) {
    int s = y.size();

    for (int i = 0; i < s; ++i) {
        y[i] *= x;
    }
}

void gauss_vector_sum(vector<double> &x, vector<double> y) {
    int s = x.size();

    for (int i = 0; i < s; ++i) {
        x[i] += y[i];
    }
}

void square_matrix_mult(vector<vector<double> > &z, vector<vector<double> > x, vector<vector<double> > y) {
    int s = x.size();
    double t;

    for (int i = 0; i < s; ++i) {
        for (int j = 0; j < s; ++j) {
            t = 0;

            for (int k = 0; k < s; ++k) {
                t += x[i][k] * y[k][j];
            }

            z[i][j] = t;

        }
    }
}

void matrix_vector_mult(vector<double> &z, vector<vector<double> > x, vector<double> y) {
    int s = x.size();
    double t;

    for (int i = 0; i < s; ++i) {
        t = 0;

        for (int j = 0; j < s; ++j) {
            t += x[i][j] * y[j];
        }

        z[i] = t;
    }
}

void matrix_sum(vector<vector<double> > &z, vector<vector<double> > x, vector<vector<double> > y) {
    int s = x.size();

    for (int i = 0; i < s; ++i) {
        for (int j = 0; j < s; ++j) {
            z[i][j] = x[i][j] + y[i][j];
        }
    }
}

double norm_2(vector<double> x) {
    int s = x.size();

    double c = 0;

    for (int i = 0; i < s; ++i) {
        c += x[i] * x[i];
    }

    return sqrt(c);
}

double sgn_h(double x) {
    if (x >= 0) {
        return 1;
    } else {
        return -1;
    }
}

void clear(vector<double> &x) {
    int s = x.size();

    for (int i = 0; i < s; ++i) {
        x[i] = 0;
    }
}

void householder_matrix(vector<vector<double> > &h, vector<double> w) {
    int s = w.size();

    for (int i = 0; i < s; ++i) {
        for (int j = 0; j < s; ++j) {
            if (i == j) {
                h[i][j] = 1 + (-2) * w[i] * w[j];
            } else {
                h[i][j] = -2 * w[i] * w[j];
            }
        }
    }
}

void display_vector(vector<double> x, ofstream &y) {
    int s = x.size();

    for (int i = 0; i < s; ++i) {
        y << x[i] << ' ';
    }

    y << endl << endl;
}

void display_matrix(vector<vector<double> > x, ofstream &y) {
    int s = x.size();

    for (int i = 0; i < s; ++i) {
        for (int j = 0; j < s; ++j) {
            y << x[i][j] << ' ';
        }
        y << endl;
    }
    y << endl;
}

void transpose(vector<vector<double> > &x) {
    int s = x.size();
    vector<vector<double> > y(x);

    for (int i = 0; i < s; ++i) {
        for (int j = 0; j < s; ++j) {
            y[j][i] = x[i][j];
        }
    }
    x = y;
}

void vector_norm(vector<double> &x) {
    double n = norm_2(x);
    int s = x.size();

    for (int i = 0; i < s; ++i) {
        x[i] /= n;
    }
}

double difference_norm(vector<double> x, vector<double> y) {
    gauss_scalar_vector_mult(-1, y);
    gauss_vector_sum(x, y);

    return norm_2(x);
}

double norm_residual(vector<vector<double> > a, vector<double> x, vector<double> b) {
    int s = x.size();

    vector<double> y(s, 0);

    matrix_vector_mult(y, a, x);

    return difference_norm(y, b);
}

double scalar_product(vector<double> x, vector<double> y) {
    int s = x.size();
    double t = 0;

    for (int i = 0; i < s; ++i) {
        t += x[i] * y[i];
    }

    return t;
}

int main() {
//    gauss1();
//    gauss2();
//    gauss3();
//    tridiagonal_matrix_algorithm();
//    householder();
//    richardson();
//    conjugate_gradient_method();
//    qr_algorithm();
//    inv_iter();
//    inv_iter_rayleigh();

    return 0;
}