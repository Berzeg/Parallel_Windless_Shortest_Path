#include "../headers/WindVector.h"

WindVector::WindVector( int x, int y, float magnitude, int bearing )
{
	this->x = x;
	this->y = y;
	this->magnitude = magnitude;
	this->bearing = bearing;
}

WindVector::~WindVector() {}