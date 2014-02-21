#ifndef SEMIVARIANCE
#define SEMIVARIANCE

/*
Hashem Shawqi (hashemshawqi@cmail.carleton.ca)
School of Computer Science, Carleton University

This class represents a semivariogram for each of the 4 main compass directions. However, if the velocity in one direction is negative then the velocity for that direction is set to 0.
*/

enum direction { NORTH = 0, EAST, SOUTH, WEST };

class Semivariance
{
public:
	float displacement;
	float semiVariance[4];
	float velocity[4];
};

#endif