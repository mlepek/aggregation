/******************************************
The code for the numerical simulations of coagulating systems
Author: Michal Lepek
Date: 02 Dec 2019
If using this code for your research, please refer to the main work:
"Combinatorial solutions to condensation, electrorheological and other aggregation kernels"
by M. Lepek, A. Fronczak & P. Fronczak.
******************************************/

#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>


typedef std::vector<int> IVector;
typedef std::vector<double> DVector;

inline IVector histogram(IVector input, int maksymalny_rozmiar_klastra);
inline DVector histogramD(IVector input, int maksymalny_rozmiar_klastra);

/*
std::string to_string(int num)
{
    std::ostringstream str1;
    str1 << num;
    return str1.str();
}
*/

std::string to_string(float num)
{
    std::ostringstream str1;
    str1 << num;
    return str1.str();
}

int int_pow(int N, int power)
{
    int result = 1;
    for(int i=0; i<power; i++)
        result *= N;
    return result;
}

inline IVector histogram(IVector input, int maksymalny_rozmiar_klastra)
{
	IVector output(maksymalny_rozmiar_klastra);

	for (auto i : input)
		output.at(i - 1)++;

	return output;
}

inline DVector histogramD(IVector input, int maksymalny_rozmiar_klastra)
{
	DVector output(maksymalny_rozmiar_klastra);

	for (auto i : input)
		output.at(i - 1)++;

	return output;
}

#endif
