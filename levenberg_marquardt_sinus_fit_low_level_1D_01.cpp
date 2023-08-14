// levenberg_marquardt_sinus_fit_low_level_1D_01.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <cmath>
#include <math.h>

# define M_PI 3.141592653589793238462643383279502884197

using namespace std;

int levmar_sinus(double* t_data_inp, double* y_data_inp, const int M_inp, double* x0_inp, double* x_fit_outp, int k_max_inp, double eps_1_inp, double eps_2_inp, double tau_inp) {
    // initial iteration variable
    int k = 0;
    int ni = 2;
    // initialize Jacobian 1D matrix
    double* J = new double[3*M_inp];
    // fill Jacobian matrix
    for (int index_0 = 0; index_0 < M_inp; index_0++) {
        J[0 + index_0 * 3] = (-1.0f) * std::sin(t_data_inp[index_0]+x0_inp[1]);
        J[1 + index_0 * 3] = (-1.0f) * x0_inp[0] * std::cos(t_data_inp[index_0] + x0_inp[1]);
        J[2 + index_0 * 3] = -1.0f;
    }
    // initialize transpose of Jacobian matrix in 1D form
    double* J_transpose = new double[3*M_inp];
    // fill transpose of Jacobian matrix
    for (int index_0 = 0; index_0 < M_inp; index_0++) {
        for (int index_1 = 0; index_1 < 3; index_1++) {
            J_transpose[index_0 + index_1 * M_inp] = J[index_1 + index_0 * 3];
        }
    }
    // calculate A matrix
    // initialize A matrix
    double* A = new double[9];
    // initialize A matrix to zero values
    for (int index_0 = 0; index_0 < 9; index_0++) {
        A[index_0] = 0.0f;
    }
    // calculate A matrix as J_transpose * J
    // use multiplication of 2D matrices
    for (int index_0 = 0; index_0 < 3; index_0++) {
        for (int index_1 = 0; index_1 < 3; index_1++) {
            for (int index_2 = 0; index_2 < M_inp; index_2++) {
                A[index_1 + index_0 * 3] = A[index_1 + index_0 * 3] + J_transpose[index_2 + index_0 * M_inp] * J[index_1 + index_2 * 3];
            }
        }
    }
    // calculate f function
    // initialize f function
    double* f = new double[M_inp];
    // fill f function
    for (int index_0 = 0; index_0 < M_inp; index_0++) {
        f[index_0] = y_data_inp[index_0] - x0_inp[0] * std::sin(t_data_inp[index_0] + x0_inp[1]) - x0_inp[2];
    }
    // calculate transpose of f
    // initialize transpose of f
    double* f_transpose = new double[M_inp];
    // fill transpose of f
    for (int index_0 = 0; index_0 < M_inp; index_0++) {
        f_transpose[index_0] = f[index_0];
    }
    // calculate g as J_transpose * f_transpose
    // initialize g
    double* g = new double[3];
    // initialize g to zero values
    for (int index_0 = 0; index_0 < 3; index_0++) {
        g[index_0] = 0.0f;
    }
    for (int index_0 = 0; index_0 < 3; index_0++) {
        for (int index_1 = 0; index_1 < 1; index_1++) {
            for (int index_2 = 0; index_2 < M_inp; index_2++) {
                g[index_1 + index_0 * 1] = g[index_1 + index_0 * 1] + J_transpose[index_2 + index_0 * M_inp] * f_transpose[index_1 + index_2 * 1];
            }
        }
    }
    // calculate norm of g
    // initialize norm of g
    double g_norm = 0.0f;
    for (int index_0 = 0; index_0 < 3; index_0++) {
        g_norm = g_norm + g[index_0] * g[index_0];
    }
    g_norm = std::sqrt(g_norm);
    // boolean variable
    bool found_bool = (g_norm <= eps_1_inp);
    // initialize mi
    double mi = 0.0f;
    double A_diag_max = std::max(A[0], A[4]);
    A_diag_max = std::max(A_diag_max, A[8]);
    // calculate mi
    mi = tau_inp * A_diag_max;

    // initialize x vector
    double* x = new double[3];
    // fill x vector
    for (int index_0 = 0; index_0 < 3; index_0++) {
        x[index_0] = x0_inp[index_0];
    }
    // initialize x_new vector
    double* x_new = new double[3];
    // initialize x_new vector to zero values
    for (int index_0 = 0; index_0 < 3; index_0++) {
        x_new[index_0] = 0.0f;
    }
    // initialize transpose of g
    double* g_transpose = new double[3];
    // initialize transpose of g to zero values
    for (int index_0 = 0; index_0 < 3; index_0++) {
        g_transpose[index_0] = 0.0f;
    }
    // initialize B matrix
    double* B = new double[9];
    // initialize B matrix to zero values
    for (int index_0 = 0; index_0 < 9; index_0++) {
        B[index_0] = 0.0f;
    }
    // initialize inversion matrix of B
    double* B_inv = new double[9];
    // initialize inversion matrix of B to zero values
    for (int index_0 = 0; index_0 < 9; index_0++) {
        B_inv[index_0] = 0.0f;
    }
    // initialize adjoint matrix of B
    double* B_adj = new double[9];
    // initialize adjoint matrix of B to zero values
    for (int index_0 = 0; index_0 < 9; index_0++) {
        B_adj[index_0] = 0.0f;
    }
    // initialize value of determinant of B
    double B_det = 0.0f;
    // initialize minors of matrix B
    double B_minor_11 = 0.0f;
    double B_minor_12 = 0.0f;
    double B_minor_13 = 0.0f;
    double B_minor_21 = 0.0f;
    double B_minor_22 = 0.0f;
    double B_minor_23 = 0.0f;
    double B_minor_31 = 0.0f;
    double B_minor_32 = 0.0f;
    double B_minor_33 = 0.0f;
    // initialize h_lm vector
    double* h_lm = new double[3];
    // initialize h_lm vector to zero values
    for (int index_0 = 0; index_0 < 3; index_0++) {
        h_lm[index_0] = 0.0f;
    }
    // initialize transpose of h_lm vector
    double* h_lm_transpose = new double[3];
    // initialize transpose of h_lm vector to zero values
    for (int index_0 = 0; index_0 < 3; index_0++) {
        h_lm_transpose[index_0] = 0.0f;
    }
    // initialize mi * h_lm_transpose - g
    double* mi_h_lm_trans_minus_g = new double[3];
    // initialize mi * h_lm_transpose - g to zero values
    for (int index_0 = 0; index_0 < 3; index_0++) {
        mi_h_lm_trans_minus_g[index_0] = 0.0f;
    }
    // initialize norm of h_lm vector
    double h_lm_norm = 0.0f;
    // initilize norm of x vector
    double x_norm = 0.0f;
    // initialize F_x value
    double F_x = 0.0f;
    // initialize F_x_new value
    double F_x_new = 0.0f;
    // initialize ro_denominator
    double ro_denominator = 0.0f;
    // intialize ro value - gain ratio
    double ro = 0.0f;

    while (!found_bool && (k < k_max_inp)) {
        // increase iteration by one
        k++;
        // calculate matrix B
        // insert matrix A to matrix B
        for (int index_0 = 0; index_0 < 9; index_0++) {
            B[index_0] = A[index_0];
        }
        B[0] = B[0] + mi;
        B[4] = B[4] + mi;
        B[8] = B[8] + mi;
        // calculate transpose of g
        for (int index_0 = 0; index_0 < 3; index_0++) {
            g_transpose[index_0] = g[index_0];
        }
        // calculate inversion of B
        // calculate minor values
        B_minor_11 = (+1.0f) * (B[4] * B[8] - B[5] * B[7]);
        B_minor_12 = (-1.0f) * (B[3] * B[8] - B[5] * B[6]);
        B_minor_13 = (+1.0f) * (B[3] * B[7] - B[4] * B[6]);
        B_minor_21 = (-1.0f) * (B[1] * B[8] - B[2] * B[7]);
        B_minor_22 = (+1.0f) * (B[0] * B[8] - B[2] * B[6]);
        B_minor_23 = (-1.0f) * (B[0] * B[7] - B[1] * B[6]);
        B_minor_31 = (+1.0f) * (B[1] * B[5] - B[2] * B[4]);
        B_minor_32 = (-1.0f) * (B[0] * B[5] - B[2] * B[3]);
        B_minor_33 = (+1.0f) * (B[0] * B[4] - B[1] * B[3]);
        // calculate adjoint matrix of B
        B_adj[0] = B_minor_11; 
        B_adj[1] = B_minor_21;
        B_adj[2] = B_minor_31;
        B_adj[3] = B_minor_12;
        B_adj[4] = B_minor_22;
        B_adj[5] = B_minor_32;
        B_adj[6] = B_minor_13;
        B_adj[7] = B_minor_23;
        B_adj[8] = B_minor_33;
        // calculate determinant value of matrix B
        B_det = B[0] * B_minor_11 + B[1] * B_minor_12 + B[2] * B_minor_13;
        // calculate inversion of matrix B
        for (int index_0 = 0; index_0 < 9; index_0++) {
            B_inv[index_0] = B_adj[index_0] / B_det;
        }
        // calculate h_lm vector
        // initialize h_lm vector to zero values
        for (int index_0 = 0; index_0 < 3; index_0++) {
            h_lm[index_0] = 0.0f;
        }
        for (int index_0 = 0; index_0 < 1; index_0++) {
            for (int index_1 = 0; index_1 < 3; index_1++) {
                for (int index_2 = 0; index_2 < 3; index_2++) {
                    h_lm[index_1 + index_0 * 3] = h_lm[index_1 + index_0 * 3] + (-1.0f) * g_transpose[index_2 + index_0 * 3] * B_inv[index_1 + index_2 * 3];
                }
            }
        }
        // calculate norm of h_lm vector
        // fill norm of h_lm vector with zero value
        h_lm_norm = 0.0f;
        for (int index_0 = 0; index_0 < 3; index_0++) {
            h_lm_norm = h_lm_norm + h_lm[index_0] * h_lm[index_0];
        }
        h_lm_norm = std::sqrt(h_lm_norm);
        // calculate norm of x vector
        // fill norm of x vector with zero value
        x_norm = 0.0f;
        for (int index_0 = 0; index_0 < 3; index_0++) {
            x_norm = x_norm + x[index_0] * x[index_0];
        }
        x_norm = std::sqrt(x_norm);
        // main condition
        if (h_lm_norm <= eps_2_inp * (x_norm + eps_2_inp)) {
            found_bool = true;
        }
        else {
            // calculate vector x_new
            for (int index_0 = 0; index_0 < 3; index_0++) {
                x_new[index_0] = x[index_0] + h_lm[index_0];
            }
            // print iteration result
            std::cout << k << " " << double(x_new[0]) << " " << double(x_new[1]) << " " << double(x_new[2]) << std::endl;
            // calculate F(x)
            // caclulate function f
            for (int index_0 = 0; index_0 < M_inp; index_0++) {
                f[index_0] = y_data_inp[index_0] - x[0] * std::sin(t_data_inp[index_0] + x[1]) - x[2];
            }
            // calculate F_x value
            // initialize F_x to zero value 
            F_x = 0.0f;
            for (int index_0 = 0; index_0 < M_inp; index_0++) {
                F_x = F_x + f[index_0] * f[index_0];
            }
            F_x = 0.5f * F_x;

            // calculate F(x_new)
            // calculate function f
            for (int index_0 = 0; index_0 < M_inp; index_0++) {
                f[index_0] = y_data_inp[index_0] - x_new[0] * std::sin(t_data_inp[index_0] + x_new[1]) - x_new[2];
            }
            // calculate F_x_new value
            // initialize F_x_new to zero value
            F_x_new = 0.0f;
            for (int index_0 = 0; index_0 < M_inp; index_0++) {
                F_x_new = F_x_new + f[index_0] * f[index_0];
            }
            F_x_new = 0.5f * F_x_new;

            // calculate ro_denominator part
            // initialize transpose of h_lm vector
            for (int index_0 = 0; index_0 < 3; index_0++) {
                h_lm_transpose[index_0] = h_lm[index_0];
            }
            // calculate mi * h_lm_transpose - g
            // initialize mi * h_lm_transpose - g to zero values
            for (int index_0 = 0; index_0 < 3; index_0++) {
                mi_h_lm_trans_minus_g[index_0] = 0.0f;
            }
            for (int index_0 = 0; index_0 < 3; index_0++) {
                mi_h_lm_trans_minus_g[index_0] = mi * h_lm_transpose[index_0] - g[index_0];
            }
            // calculate ro_denominator
            ro_denominator = 0.0f;
            for (int index_0 = 0; index_0 < 1; index_0++) {
                for (int index_1 = 0; index_1 < 1; index_1++) {
                    for (int index_2 = 0; index_2 < 3; index_2++) {
                        ro_denominator = ro_denominator + h_lm[index_2 + index_0 * 3] * mi_h_lm_trans_minus_g[index_1 + index_2 * 1];
                    }
                }
            }
            ro_denominator = 0.5f * ro_denominator;
            // calculate ro value - gain ratio
            ro = (F_x - F_x_new) / ro_denominator;
            if (ro > 0.0f) {
                // insert vector x_new into the vector x
                for (int index_0 = 0; index_0 < 3; index_0++) {
                    x[index_0] = x_new[index_0];
                }
                // fill Jacobian matrix
                for (int index_0 = 0; index_0 < M_inp; index_0++) {
                    J[0 + index_0 * 3] = (-1.0f) * std::sin(t_data_inp[index_0] + x[1]);
                    J[1 + index_0 * 3] = (-1.0f) * x[0] * std::cos(t_data_inp[index_0] + x[1]);
                    J[2 + index_0 * 3] = -1.0f;
                }
                // fill transpose of Jacobian matrix
                for (int index_0 = 0; index_0 < M_inp; index_0++) {
                    for (int index_1 = 0; index_1 < 3; index_1++) {
                        J_transpose[index_0 + index_1 * M_inp] = J[index_1 + index_0 * 3];
                    }
                }
                // calculate A matrix
                // initialize A matrix to zero values
                for (int index_0 = 0; index_0 < 9; index_0++) {
                    A[index_0] = 0.0f;
                }
                // calculate A matrix as J_transpose * J
                // use multiplication of 2D matrices
                for (int index_0 = 0; index_0 < 3; index_0++) {
                    for (int index_1 = 0; index_1 < 3; index_1++) {
                        for (int index_2 = 0; index_2 < M_inp; index_2++) {
                            A[index_1 + index_0 * 3] = A[index_1 + index_0 * 3] + J_transpose[index_2 + index_0 * M_inp] * J[index_1 + index_2 * 3];
                        }
                    }
                }
                // calculate f function
                // fill f function
                for (int index_0 = 0; index_0 < M_inp; index_0++) {
                    f[index_0] = y_data_inp[index_0] - x[0] * std::sin(t_data_inp[index_0] + x[1]) - x[2];
                }
                // calculate transpose of f
                // fill transpose of f
                for (int index_0 = 0; index_0 < M_inp; index_0++) {
                    f_transpose[index_0] = f[index_0];
                }
                for (int index_0 = 0; index_0 < 3; index_0++) {
                    g[index_0] = 0.0f;
                }
                // calculate g as J_transpose * f_transpose
                for (int index_0 = 0; index_0 < 3; index_0++) {
                    for (int index_1 = 0; index_1 < 1; index_1++) {
                        for (int index_2 = 0; index_2 < M_inp; index_2++) {
                            g[index_1 + index_0 * 1] = g[index_1 + index_0 * 1] + J_transpose[index_2 + index_0 * M_inp] * f_transpose[index_1 + index_2 * 1];
                        }
                    }
                }
                // calculate norm of g
                // initialize norm of g
                g_norm = 0.0f;
                for (int index_0 = 0; index_0 < 3; index_0++) {
                    g_norm = g_norm + g[index_0] * g[index_0];
                }
                g_norm = std::sqrt(g_norm);
                // calculate boolean variable
                found_bool = (g_norm <= eps_1_inp);
                // calculate mi
                mi = mi * std::max(double(0.33333333333333333333f), double(1.0f - std::pow((2.0f * ro - 1), 3.0f)));
                // define ni
                ni = 2;
            }
            else {
                // calculate mi
                mi = mi * double(ni);
                // calculate ni
                ni = 2 * ni;
            }
        }
    }

    //std::cout << x_new << std::endl;
    // convert phase shift from fitting to interval (0, 2*pi)
    if (x_new[1] > 0.0f && x_new[0] > 0.0f) {
        x_new[1] = x_new[1] - (2.0f * M_PI) * int(x_new[1] / (2 * M_PI));
    }
    else if (x_new[1] < 0.0f && x_new[0] > 0.0f) {
        x_new[1] = x_new[1] - (2.0f * M_PI) * (int(x_new[1] / (2 * M_PI)) - 1);
    }
    else if (x_new[1] > 0.0f && x_new[0] < 0.0f) {
        x_new[0] = (-1.0f) * x_new[0];
        x_new[1] = x_new[1] + 1.0f * M_PI;
        x_new[1] = x_new[1] - (2.0f * M_PI) * int(x_new[1] / (2 * M_PI));
    }
    else if (x_new[1] < 0.0f && x_new[0] < 0.0f) {
        x_new[0] = (-1.0f) * x_new[0];
        x_new[1] = x_new[1] - 1.0f * M_PI;
        x_new[1] = x_new[1] - (2.0f * M_PI) * (int(x_new[1] / (2 * M_PI)) - 1);
    }
    else {
        x_new[1] = 0.0f;
    }
    //std::cout << x_new << std::endl;

    // store fitting results to output 1D double array
    x_fit_outp[0] = x_new[0];
    x_fit_outp[1] = x_new[1];
    x_fit_outp[2] = x_new[2];

    delete[] J;
    delete[] J_transpose;
    delete[] A;
    delete[] f;
    delete[] f_transpose;
    delete[] g;
    delete[] g_transpose;
    delete[] x;
    delete[] x_new;
    delete[] B;
    delete[] B_inv;
    delete[] B_adj;
    delete[] h_lm;
    delete[] h_lm_transpose;
    delete[] mi_h_lm_trans_minus_g;

    return 0;
}

int main()
{
    // define experimental sinus parameters
    // experimental amplitude
    double x1 = 306.523f;
    // experimental phase shift
    double x2 = double(M_PI / 3.0f);
    // experimental offset
    double x3 = 17.682f;
    // define experimental sinus parameter vector
    double* x = new double[3];
    // fill experimental parameters vector
    x[0] = x1;
    x[1] = x2;
    x[2] = x3;
    // define intial parameters of sinusoidal function
    // intial amplitude
    double x01 = 1.0f;
    // initial phase shift
    double x02 = double(M_PI / 2.0f);
    // initial offset
    double x03 = 0.0f;
    // define initial parameter vector
    double* x0 = new double[3];
    // fill initial parameters vector
    x0[0] = x01;
    x0[1] = x02;
    x0[2] = x03;
    // define number of steps in fringe scanning
    // define number of experimental points
    const int M = 5;
    // define input t_data 1D array, on the t axis 
    double* t_data = new double[M];
    for (int index_0 = 0; index_0 < M; index_0++) {
        t_data[index_0] = double(index_0) * ((2.0f * M_PI) / double(M));
    }
    // define input y_data 1D array
    double* y_data = new double[M];
    for (int index_0 = 0; index_0 < M; index_0++) {
        y_data[index_0] = x[0] * std::sin(t_data[index_0] + x[1]) + x[2];
    }
    // maximum number of iterations in fitting
    int k_max = 1000;
    // auxiliar fitting variable epsilon 1
    double eps_1 = 1.0E-8f;
    // auxiliar fitting variable epsilon 2
    double eps_2 = 1.0E-8f;
    // auxiliar fitting variable tau
    double tau = 1.0E-3f;
    // create output 1D array where fitting results will be stored
    double* x_fit = new double[3];

    // perform Levenberg-Marquardt sinusoidal fitting
    levmar_sinus(t_data, y_data, M, x0, x_fit, k_max, eps_1, eps_2, tau);

    // print fitting results on the screen
    std::cout << "The amplitude of sinus curve is: " << x_fit[0] << std::endl;
    std::cout << "The phase of sinus curve is: " << x_fit[1] << std::endl;
    std::cout << "The offset of sinus curve is: " << x_fit[2] << std::endl;

    // delete buffers
    delete[] x;
    delete[] x0;
    delete[] t_data;
    delete[] y_data;
    delete[] x_fit;

    return 0;
}
