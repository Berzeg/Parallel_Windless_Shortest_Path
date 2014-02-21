#include "../headers/WindVector.h"

WindVector::WindVector( int x, int y, float v[] )
{
	this->x = x;
	this->y = y;

	for( int i = 0; i < 4; i++ )
	{
		velocity[ i ] = v[ i ];
	}
}

WindVector::~WindVector() {}