*********************************************************
Michal Lepek, wepas@wp.eu, michal.lepek.dokt@pw.edu.pl
02 Dec 2019
Read me file - How to run the code for plotting theory
in the exact combinatorial approach to coagulation
*********************************************************

Here, is the code in C++ for calculations of the theoretical predictions for the work "Combinatorial solutions to condensation, electrorheological and other aggregation kernels".

The only file that I wrote is "theory_kernels.cpp". To run this code you need two dependencies:

1. GNU Multiple Precision library for float numbers.
It is needed because Bell polynomials easily got extremely high values, exceeding the standard floating point memory. The instruction how to install this library is a separate document, "How to install GMP". This will show you the way of installing on Linux Ubuntu. I used the version GNU MPFR Library 4.0.1 from the webpage https://www.mpfr.org/ - it is attached in this repository ("mpfr-4.0.1.tar"). 

2. C++ interface for MPFR library, MPFRC++.
You can download it from http://www.holoborodko.com/pavel/mpfr/#download. Luckily, the interface is one header file, so no library building needed. This file ("mpreal.h") is also provided here in the repository as same as whole mpfrc++ package ("mpfrc++-3.6.2").

I worked on Linux Ubuntu 12.04, but the code should be cross-platform (if only you manage to install GNU MPFR on Windows or Mac - I did not try). To compile the code on Linux, I used (gcc 4.9.4):

g++ theory_kernels.cpp -o theory -lmpfr -lgmp -w -std=c++11 -O2

In the code, there are comments explaining what is going on there.

If you use this code for your research, please cite the main paper:

"Combinatorial solutions to coagulation kernel for linear chains"
by M. Lepek, A. Fronczak & P. Fronczak, Physica D: Nonlinear Phenomena, Volume 415, January 2021, 132756.
https://doi.org/10.1016/j.physd.2020.132756

Best regards!
M.L.