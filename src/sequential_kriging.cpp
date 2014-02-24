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
vector < Semivariance* > calculate_all_semivariances( vector< WindVector* > &input_data )
{
	vector < Semivariance* > semivariances;

	for(int i = 0; i < input_data.size() - 1; i++)
	{
		for(int j = i + 1; j < input_data.size(); j++)
		{
			// variables for pythagoras
			float x1 = input_data[ i ]->x;
			float x2 = input_data[ j ]->x;
			float y1 = input_data[ i ]->y;
			float y2 = input_data[ j ]->y;

			// find the displacement of the two samples using their position coordinates
			// using pythagoras
			float displacement = sqrt( pow( x1 - x2, 2 ) + pow( y1 - y2, 2 ) );
		
			float gamma[2];

			for( int k = 0; k < 2; k++ )
			{
				// variables for the semivariance equation
				float v1 = input_data[ i ]->velocity[ k ];
				float v2 = input_data[ j ]->velocity[ k ];

				// semivariance equation for points i and j
				gamma[ k ] = 0.5 * pow( v1 - v2, 2 );
			}

			// Create the semivariance object and add it to the list
			Semivariance* newSemivariance = new Semivariance( displacement, gamma );

			cout << "The semivariance created has displacement: " << newSemivariance->displacement << "\n";
			semivariances.push_back( newSemivariance );
		}
	}

	cout << "size of semivariances (in function): " << semivariances.size() << "\n";
	return semivariances;
}

// This function takes in a list of semivariances buckets them
// according to their displacement. Then it calculates the average
// semivariance for each bucket and returns a list of averaged
// semivariances. 
// The max distance reflects the largest gap between any two 
// measurements in our data. 
// The step-size reflects the granularity of our buckets.
vector < Semivariance* > bucket_sort( int max_distance, int step_size, vector < Semivariance* > &input_semivariances )
{
	// The index_array keeps track of the index for the semivariances of each lag class
	// The count_array keeps count of the number of semivariances in each class
	// The bucketed_semivariances vector will hold the averaged value for all SVs of a
	// certain lag class.
	int num_buckets = max_distance / step_size + 1;
	int index_array[ num_buckets ];
	int count_array[ num_buckets ];
	vector < Semivariance* > bucketed_semivariances;

	// fill the index_array and count_array with a default value of 0
	std::fill_n( index_array, num_buckets, 0 );
	std::fill_n( count_array, num_buckets, 0 );

	// traverse the input SVs and accumulate their data into new SV objects and push them
	// into the SV bucket vector
	for (int i = 0; i < input_semivariances.size(); i++ )
	{
		Semivariance* sv = input_semivariances[i];
		int lag_class = sv->displacement / step_size;

		// get the index of our sv bucket
		int sv_insert_index = index_array[ lag_class ];

		// update the SV object in the bucket
		// if there is no object yet, then create one
		if ( count_array[ lag_class ] == 0 )
		{
			int new_sv_displacement = ( lag_class + 1 ) * step_size;
			Semivariance* new_sv = new Semivariance( new_sv_displacement, sv->semivariance );

			bucketed_semivariances.insert( bucketed_semivariances.begin() + sv_insert_index, new_sv );

			// update the index of all buckets after the current SV's lag class because we've
			// added a new bucket
			for ( int j = lag_class + 1; j < num_buckets; j++ )
			{
				index_array[ j ]++;
			}

		} else {

			// If the SV is already there, then add to the current SV components
			bucketed_semivariances[ sv_insert_index ]->semivariance[ X_COMPONENT ] += sv->semivariance[ X_COMPONENT ];
			bucketed_semivariances[ sv_insert_index ]->semivariance[ Y_COMPONENT ] += sv->semivariance[ Y_COMPONENT ];

		}

		// increment the count for bucket we're inserting the SV into
		count_array[ lag_class ]++;
	}


	// Now that we have summed all semivariances for each lag class we will divide the
	// SV values by the number of contributing SVs
	for ( int i = 0; i < num_buckets; i++ )
	{
		float count = float( count_array[ i ] );

		// if there were any SVs contributing to this bucket at all
		if ( count > 0 )
		{
			int index = index_array[ i ];
			Semivariance* bucket = bucketed_semivariances[ index ];

			// All average calculations have division, this is our division
			bucket->semivariance[ X_COMPONENT ] /= count;
			bucket->semivariance[ Y_COMPONENT ] /= count;
		}
	}

	return bucketed_semivariances;
}


// This is a function that takes a vector of semivariances, a pivot index,
// a pointer to left mean value for the x-component (XLMV) and another to 
// the right mean value for the x-component (XRMV), a YLMV, and a YRMV. It
// calculates the mean of all x-semivariances up to and including the pivot
// and stores it in the (X/Y)LMV, it does the same to the right side of the
// pivot and stores the result in (X/Y)RMV.
void calculate_left_and_right_means( vector< Semivariance* > semivariances, int pivot, float* XLMV, float* XRMV, float* YLMV, float* YRMV )
{
	// if this is true we'll end up dividing the right values by 0
	if ( semivariances.size() - pivot - 1 == 0 )
	{
		return;
	}

	// set the mean variables to 0
	*XLMV = 0;
	*YLMV = 0;
	*XRMV = 0;
	*YRMV = 0;

	for ( int i = 0; i <= pivot; i++ )
	{
		*XLMV += semivariances[ i ]->semivariance[ X_COMPONENT ];
		*YLMV += semivariances[ i ]->semivariance[ Y_COMPONENT ];
	}

	*XLMV /= ( pivot + 1 );
	*YLMV /= ( pivot + 1 );

	for ( int i = pivot + 1; i < semivariances.size(); i++ )
	{
		*XRMV += semivariances[ i ]->semivariance[ X_COMPONENT ];
		*YRMV += semivariances[ i ]->semivariance[ Y_COMPONENT ];
	}

	*XRMV /= ( semivariances.size() - pivot - 1 );
	*YRMV /= ( semivariances.size() - pivot - 1 );
}

// This function takes a vector of semivariances, a pivot index,
// the XLMV, YLMV, XRMV, YRMV from the calculate_left_and_right_means
// method. It also takes pointers to the x component left variance
// (xl_variance), yl_variance, xr_variance, and yr_variance; it
// calculates each variance and assigns the value to the parameters.
void calculate_left_and_right_variances( vector< Semivariance* > semivariances, int pivot, float XLMV, float XRMV, float YLMV, float YRMV, float* xl_variance, float* yl_variance, float* xr_variance, float* yr_variance)
{
	// if this is true we'll end up dividing the right values by 0
	if ( semivariances.size() - pivot - 1 == 0 )
	{
		return;
	}

	// set the variance variables to 0
	*xl_variance = 0;
	*yl_variance = 0;
	*xr_variance = 0;
	*yr_variance = 0;


	for ( int i = 0; i <= pivot; i++ )
	{
		*xl_variance += pow( XLMV - semivariances[ i ]->semivariance[ X_COMPONENT ], 2 );
		*yl_variance += pow( YLMV - semivariances[ i ]->semivariance[ Y_COMPONENT ], 2 );
	}

	*xl_variance /= ( pivot + 1 );
	*yl_variance /= ( pivot + 1 );

	for ( int i = pivot + 1; i < semivariances.size(); i++ )
	{
		*xr_variance += pow( XRMV - semivariances[ i ]->semivariance[ X_COMPONENT ], 2 );
		*yr_variance += pow( YRMV - semivariances[ i ]->semivariance[ Y_COMPONENT ], 2 );
	}

	*xr_variance /= ( semivariances.size() - pivot - 1 );
	*yr_variance /= ( semivariances.size() - pivot - 1 );
}

int main()
{

}