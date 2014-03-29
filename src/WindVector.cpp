#include "../headers/WindVector.h"

WindVector::WindVector( int x, int y, double v[] )
{
	this->x = x;
	this->y = y;

	for( int i = 0; i < 2; i++ )
	{
		velocity[ i ] = v[ i ];
	}
}

WindVector::~WindVector() {}