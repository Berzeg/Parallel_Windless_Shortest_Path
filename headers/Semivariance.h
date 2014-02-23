#ifndef SEMIVARIANCE
#define SEMIVARIANCE
#include <iostream>

/*
Hashem Shawqi (hashemshawqi@cmail.carleton.ca)
School of Computer Science, Carleton University

This class represents a semivariogram for each of the 4 main compass directions. However, if the velocity in one direction is negative then the velocity for that direction is set to 0.
*/

// increments with clockwise movement starting from north
enum direction { X_COMPONENT, Y_COMPONENT };

class Semivariance
{
public:
	float displacement;
	float semivariance[2];

	Semivariance( float, float[] );
	Semivariance( int, float[] );
	~Semivariance();
};

#endif