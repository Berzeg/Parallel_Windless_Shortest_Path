/*
Hashem Shawqi (hashemshawqi@cmail.carleton.ca)
School of Computer Science, Carleton University

This program receives an input of a sparse vector field, it performs a sequential kriging algorithm to establish a geostatistical correlation with regards to the distance between two points and the difference in magnitude of each.

The output of the program is a 2D array with original and interpolated vector values at each point.
*/
#include "../headers/WindVector.h"
#include "../headers/Semivariance.h"
#include <vector> /* vector */
#include <cmath> /* sqrt */

using namespace std;

// Given a list of the initial data points this
// function calculates the semivariance between
// each pair and outputs it in linked-list format
vector < Semivariance* >* calculate_all_semivariances( vector< WindVector* >& input_data )
{
	vector < Semivariance* > semivariances;

	for(int i = 0; i < input_data.size - 2; i++)
	{
		for(int j = i + 1; j < input_data.size - 1; j++)
		{
			// find the displacement of the two samples using their position coordinates
			float displacement = sqrt( pow( input_data[ i ]->x - input_data[ j ]->x, 2 ) + pow( input_data[ i ]->y - input_data[ j ]->y, 2 ) );
		
			float gamma[4];

			for( int k = 0; k < 4; k++ )
			{
				gamma[ k ] = 0.5 * ( pow( input_data[ i ]->velocity[ k ] - input_data[ j ]->velocity[ k ], 2 ) );
			}
		}
	}
}

int main()
{

}