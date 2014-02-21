#include "../headers/WindVector.h"
#include "../headers/Semivariance.h"
#include <iostream>

using namespace std;

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
}