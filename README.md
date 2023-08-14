# Levenberg-Marquardt-sinus-fitting-low-level-only-1D-arrays
I post here the Levenberg-Marquardt algorithm for fitting sinus function in C++ in low-level form, using only 1D arrays and for cycles. 

The sinusoidal function has form a*sin(x+b)+c, where a is the amplitude, b is the phase, and c is the offset. The algorithm uses Jacobian matrix J, in which the first column contains the derivation of the sinus function according to the amplitude, the second column contains the derivation of the sinus function according to the phase, and the third column contains the derivation of the sinus function according to the offset. The Jacobian matrix J is here defined as a 1D array. All matrices and vectors are here defined as 1D arrays. The algorithm uses the multiplication of the matrices in 1D array form. In addition, the Levenberg-Marquardt algorithm contains the calculation of the inversion matrix in 1D array form. The Levenberg-Marquardt algorithm uses only cycles. Therefore, it should be very fast, suitable for fast calculations. Enjoy this fast Levenberg-Marquardt algorithm in your applications where sinusoidal non-linear curve fit is necessary. 
