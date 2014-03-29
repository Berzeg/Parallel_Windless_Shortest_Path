#ifndef WIND_VECTOR
#define WIND_VECTOR

/*
Hashem Shawqi (hashemshawqi@cmail.carleton.ca)
School of Computer Science, Carleton University

This class represents individual wind vectors that make up the vector field.
*/

class WindVector
{
public:
	int x;
	int y;
	double velocity[2];

	WindVector( int, int, double[] );
	~WindVector();
};

#endif