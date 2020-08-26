#ifndef CMLA_HEADER_H
#define CMLA_HEADER_H


#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

void gauss_scalar_vector_mult(double x, vector<double> &y);
void gauss_vector_sum(vector<double> &x, vector<double> y);
void square_matrix_mult(vector<vector<double> > &z, vector<vector<double> > x, vector<vector<double> > y);
void matrix_vector_mult(vector<double> &z, vector<vector<double> > x, vector<double> y);
void matrix_sum(vector<vector<double> > &z, vector<vector<double> > x, vector<vector<double> > y);
double norm_2(vector<double> x);
double sgn_h(double x);
void clear(vector<double> &x);
void householder_matrix(vector<vector<double> > &h, vector<double> w);
void display_vector(vector<double> x, ofstream &y);
void display_matrix(vector<vector<double> > x, ofstream &y);
void transpose(vector<vector<double> > &x);
void vector_norm(vector<double> &x);
double difference_norm(vector<double> x, vector<double> y);
double norm_residual(vector<vector<double> > a, vector<double> x, vector<double> b);
double scalar_product(vector<double> x, vector<double> y);

void gauss1();
void gauss2();
void gauss3();
void tridiagonal_matrix_algorithm();
void householder();
void richardson();
void conjugate_gradient_method();
void qr_algorithm();
void inv_iter();
void inv_iter_rayleigh();

#endif //CMLA_HEADER_H