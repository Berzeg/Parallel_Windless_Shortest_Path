/*
Hashem Shawqi (hashemshawqi@cmail.carleton.ca)
School of Computer Science, Carleton University

This program receives an input of a sparse vector field, it performs a sequential kriging algorithm to establish a geostatistical correlation with regards to the distance between two points and the difference in magnitude of each.

The output of the program is a 2D array with original and interpolated vector values at each point.
*/
#include "../headers/WindVector.h"
#include "../headers/Semivariogram.h"

// Given a list of the initial data points this
// function calculates the semivariance between
// each pair and outputs it in linked-list format
void calculate_all_semivariances()
{
	for(int i = 0; i < size_of_data - 2; i++)
	{
		for(int j = i + 1; j < size_of_data - 1; j++)
		{

		}
	}
}

int main()
{

}