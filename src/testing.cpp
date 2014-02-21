#include "../headers/WindVector.h"
#include "../headers/Semivariance.h"
#include <iostream>

using namespace std;

int main(int argc, char * argv[])
{
	WindVector* wv = new WindVector( 1, 1, 5.0, 33 );

	cout << "testing WindVector class:\n";
	cout << "x: " << wv->x << "\ny: " << wv-> y << "\nmagnitude: " << wv->magnitude << "\nbearing: " << wv->bearing << "\n";  
}