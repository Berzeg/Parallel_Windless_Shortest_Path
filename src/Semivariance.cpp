#include "../headers/Semivariance.h"


Semivariance::Semivariance( double h, double gamma[] )
{
	displacement = h;
	
	for (int i = 0; i < 2; i++ )
	{
		semivariance[ i ] = gamma[ i ];
	}
}

Semivariance::Semivariance( double h, double gamma[] )
{
	displacement = double( h );
	
	for (int i = 0; i < 2; i++ )
	{
		semivariance[ i ] = gamma[ i ];
	}
}

Semivariance::~Semivariance() {}