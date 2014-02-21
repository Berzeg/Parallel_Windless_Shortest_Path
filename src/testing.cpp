#include "../headers/WindVector.h"
#include "../headers/Semivariance.h"
#include <iostream>

using namespace std;

int main(int argc, char * argv[])
{
	WindVector* wv = new WindVector( 1, 1, 5.0, 33 );

	cout << "testing WindVector class:\n";
	cout << "x: " << wv->x << "\ny: " << wv-> y << "\nmagnitude: " << wv->magnitude << "\nbearing: " << wv->bearing << "\n";  

	float gamma[] = { 1.2, 3.4, 5.6, 7.8 };
	float velocity[] = { 9.0, 1.2, 3.4, 5.6 };
	Semivariance* sv = new Semivariance( 3.1, gamma, velocity );

	cout << "Testing semivariance class: \n";
	cout << "displacement: " << sv->displacement << "\n";
	for (int i = 0; i < 4; i++ )
	{
		cout << "semivariance[" << i << "]: " << sv->semivariance[i] << "\n";
	}

	for (int i = 0; i < 4; i++ )
	{
		cout << "velocity[" << i << "]: " << sv->velocity[i] << "\n";
	}	
}