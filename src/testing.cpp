#include "../headers/WindVector.h"
#include "../headers/Semivariogram.h"

int main(int argc, char * argv[])
{
	WindVector* wv = new WindVector( 1, 1, 5, 33 );

	cout << "testing WindVector class:\n"
	cout << "x: " << wv->x << "\ny: " << wv-> y << "\nmagnitude: " << wv->magnitude << "\nbearing: " << wv->bearing << "\n";  
}