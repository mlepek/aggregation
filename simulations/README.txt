*********************************************************
Michal Lepek, wepas@wp.eu, michal.lepek.dokt@pw.edu.pl
02 Dec 2019
Read me file - How to run the code for simulations
of the coagulation process with different kernels
*********************************************************

Here, is the code in C++ for simulating the process of coagulation with various kernels which I used in the work "Combinatorial solutions to condensation, electrorheological and other aggregation kernels". This code is generic and allows you to simulate any kernel you want, if only you can write the equation for aggregation rate, K(i,j).

There are two source files:
simulation_kernels.cpp - main function and all simulations here
functions.cpp - few useful functions which I use in main

I have written this code and used it on Windows but it is self-contained and it should be cross-platform. I used Code::Blocks 17.12 environment for developing the project, so the ".cbp" project is available in the repository. I used only standard libraries. I compiled the project with -std=c++14 specifier. In the code, there are comments explaining what is going on there.

If you use this code for your research, please cite the main paper:

"Combinatorial solutions to coagulation kernel for linear chains"
by M. Lepek, A. Fronczak & P. Fronczak, Physica D: Nonlinear Phenomena, Volume 415, January 2021, 132756.
https://doi.org/10.1016/j.physd.2020.132756

Best regards!
M.L.