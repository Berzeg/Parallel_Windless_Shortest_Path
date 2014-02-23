#include "../headers/Semivariance.h"


Semivariance::Semivariance( float h, float gamma[] )
{
	displacement = h;
	
	for (int i = 0; i < 2; i++ )
	{
		semivariance[ i ] = gamma[ i ];
	}
}

Semivariance::Semivariance( int h, float gamma[] )
{
	displacement = float( h );
	
	for (int i = 0; i < 2; i++ )
	{
		semivariance[ i ] = gamma[ i ];
	}
}

Semivariance::~Semivariance() {}