#include "../headers/Semivariance.h"

Semivariance::Semivariance( float h, float gamma[], float v[] )
{
	displacement = h;
	
	for (int i = 0; i < 4; i++ )
	{
		semivariance[ i ] = gamma[ i ];
		velocity[ i ] = v[ i ];
	}
}

Semivariance::~Semivariance() {}