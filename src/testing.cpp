#include "../headers/WindVector.h"
#include "../headers/Semivariance.h"
#include <iostream>
#include <vector>
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
void calculate_left_and_right_means( vector< Semivariance* > semivariances, int pivot, int* XLMV, int* XRMV, int* YLMV, int* YRMV )
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

int main(int argc, char * argv[])
{
	/*float velocity[] = { 9.0, 1.2 };
	WindVector* wv = new WindVector( 1, 1, velocity );

	cout << "testing WindVector class:\n";
	cout << "x: " << wv->x << "\ny: " << wv-> y << "\n";
	for (int i = 0; i < 2; i++ )
	{
		cout << "velocity[" << i << "]: " << wv->velocity[i] << "\n";
	}	


	float gamma[] = { 1.2, 3.4 };
	Semivariance* sv = new Semivariance( 3.1, gamma );

	cout << "Testing semivariance class: \n";
	cout << "displacement: " << sv->displacement << "\n";
	for (int i = 0; i < 2; i++ )
	{
		cout << "semivariance[" << i << "]: " << sv->semivariance[i] << "\n";
	}


	// testing the semivariance method
	// create another wind vector in order to establish a semivariance with wv
	float velocity2[] = { 4.2, 3.3 };
	WindVector* wv2 = new WindVector( 5, 7, velocity2 );
	float velocity3[] = { 4.4, 9.1 };
	WindVector* wv3 = new WindVector( 3, 12, velocity3 );
	WindVector* wv4 = new WindVector( 12, 5, velocity3 );

	vector< WindVector* > test_input_vectors;
	test_input_vectors.push_back( wv );
	test_input_vectors.push_back( wv2 );
	test_input_vectors.push_back( wv3 );
	test_input_vectors.push_back( wv4 );

	vector< Semivariance* > semivariance_vector = calculate_all_semivariances( test_input_vectors );
	cout << "Testing the semivariance calculation method\n";
	cout << "The number of semivariances calculated: " << semivariance_vector.size();
	cout << "\nThe lag of the first semivariance: " << semivariance_vector[0]->displacement;
	cout << "\nThe value of the northern semivariance: " << semivariance_vector[0]->semivariance[NORTH] << "\n";
*/	

	// Testing the bucketing function with 6 semivariances to which I've precalculated the 
	// correct output as: 
	// x1 = 56.7, y1 = 42
	// x2 = 50.3, y2 = 58.3

	float component_sv1[] = { 12, 34 };
	float component_sv2[] = { 90, 13 };
	float component_sv3[] = { 68, 79 };
	float component_sv4[] = { 43, 65 };
	float component_sv5[] = { 17, 28 };
	float component_sv6[] = { 91, 82 };

	Semivariance* sv1 = new Semivariance( 2, component_sv1 );
	Semivariance* sv2 = new Semivariance( 3, component_sv2 );
	Semivariance* sv3 = new Semivariance( 4, component_sv3 );
	Semivariance* sv4 = new Semivariance( 5, component_sv4 );
	Semivariance* sv5 = new Semivariance( 6, component_sv5 );
	Semivariance* sv6 = new Semivariance( 9, component_sv6 );

	Semivariance* sv_array[] = { sv1, sv2, sv3, sv4, sv5, sv6 };
	vector< Semivariance* > sv_collection( sv_array, sv_array + sizeof( sv_array ) / sizeof( Semivariance* ) );

	vector <Semivariance*> average_SVs = bucket_sort( 10, 5, sv_collection );

	//cout << "The components of the first SV are:\n";
	//cout << "X_COMPONENT: " << average_SVs[0]->semivariance[ X_COMPONENT ] << ", Y_COMPONENT: " <<  average_SVs[0]->semivariance[ Y_COMPONENT ] <<  "\n";
	//cout << "The components of the second SV are:\n";
	//cout << "X-COMPONENT: " << average_SVs[1]->semivariance[ X_COMPONENT ] << ", Y-COMPONENT: " <<  average_SVs[1]->semivariance[ Y_COMPONENT ] << "\n";


	// Testing the graph fitting functions
	// The output should be the same as the above result
	int XLMV, YLMV, XRMV, YRMV;

	calculate_left_and_right_means( sv_collection, 2, &XLMV, &XRMV, &YLMV, &YRMV );

	cout << "After calculating the left and right means we've found:\n";
	cout << "XMLV = " << XLMV << ", YLMV = " << YLMV << ", XRMV = " << XRMV << ", YRMV = " << YRMV << "\n";
}