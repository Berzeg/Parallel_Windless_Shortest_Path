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
		
			float gamma[4];

			for( int k = 0; k < 4; k++ )
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

int main(int argc, char * argv[])
{
	float velocity[] = { 9.0, 1.2, 3.4, 5.6 };
	WindVector* wv = new WindVector( 1, 1, velocity );

	cout << "testing WindVector class:\n";
	cout << "x: " << wv->x << "\ny: " << wv-> y << "\n";
	for (int i = 0; i < 4; i++ )
	{
		cout << "velocity[" << i << "]: " << wv->velocity[i] << "\n";
	}	


	float gamma[] = { 1.2, 3.4, 5.6, 7.8 };
	Semivariance* sv = new Semivariance( 3.1, gamma );

	cout << "Testing semivariance class: \n";
	cout << "displacement: " << sv->displacement << "\n";
	for (int i = 0; i < 4; i++ )
	{
		cout << "semivariance[" << i << "]: " << sv->semivariance[i] << "\n";
	}


	// testing the semivariance method
	// create another wind vector in order to establish a semivariance with wv
	float velocity2[] = { 4.2, 3.3, 1.0, 2 };
	WindVector* wv2 = new WindVector( 5, 7, velocity2 );
	float velocity3[] = { 4.4, 9.1, 6.7, 2 };
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
}