#include "../headers/WindVector.h"
#include "../headers/Semivariance.h"
#include <vector> /* vector */
#include <cmath> /* sqrt */
#include <string> /* string */
#include <fstream> /* ifstream, ofstream */
#include <sstream> /* istringstream, stringstream */



using namespace std;

// Given a list of the initial data points this
// function calculates the semivariance between
// each pair and outputs it in linked-list format
vector < Semivariance* > calculate_all_semivariances( vector< WindVector* > &input_data )
{
	vector < Semivariance* > semivariances;

	for(int i = 0; i < input_data.size() - 1; i++)
	{
		for(int j = i + 1; j < input_data.size(); j++)
		{
			// variables for pythagoras
			double x1 = input_data[ i ]->x;
			double x2 = input_data[ j ]->x;
			double y1 = input_data[ i ]->y;
			double y2 = input_data[ j ]->y;

			// find the displacement of the two samples using their position coordinates
			// using pythagoras
			double displacement = sqrt( pow( x1 - x2, 2 ) + pow( y1 - y2, 2 ) );
		
			double gamma[2];

			for( int k = 0; k < 2; k++ )
			{
				// variables for the semivariance equation
				double v1 = input_data[ i ]->velocity[ k ];
				double v2 = input_data[ j ]->velocity[ k ];

				// semivariance equation for points i and j
				gamma[ k ] = 0.5 * pow( v1 - v2, 2 );
			}

			// Create the semivariance object and add it to the list
			Semivariance* newSemivariance = new Semivariance( displacement, gamma );

			semivariances.push_back( newSemivariance );
		}
	}

	return semivariances;
}

// This function takes in a list of semivariances buckets them
// according to their displacement. Then it calculates the average
// semivariance for each bucket and returns a list of averaged
// semivariances. 
// The max distance reflects the largest gap between any two 
// measurements in our data. 
// The step-size reflects the granularity of our buckets.
vector < Semivariance* > bucket_sort( int max_distance, int step_size, vector < Semivariance* > &input_semivariances )
{
	// The index_array keeps track of the index for the semivariances of each lag class
	// The count_array keeps count of the number of semivariances in each class
	// The bucketed_semivariances vector will hold the averaged value for all SVs of a
	// certain lag class.
	int num_buckets = max_distance / step_size + 1;
	int index_array[ num_buckets ];
	int count_array[ num_buckets ];
	vector < Semivariance* > bucketed_semivariances;

	// fill the index_array and count_array with a default value of 0
	std::fill_n( index_array, num_buckets, 0 );
	std::fill_n( count_array, num_buckets, 0 );

	// traverse the input SVs and accumulate their data into new SV objects and push them
	// into the SV bucket vector
	for (int i = 0; i < input_semivariances.size(); i++ )
	{
		Semivariance* sv = input_semivariances[i];
		int lag_class = sv->displacement / step_size;

		// get the index of our sv bucket
		int sv_insert_index = index_array[ lag_class ];

		// update the SV object in the bucket
		// if there is no object yet, then create one
		if ( count_array[ lag_class ] == 0 )
		{
			int new_sv_displacement = ( lag_class + 1 ) * step_size;
			Semivariance* new_sv = new Semivariance( new_sv_displacement, sv->semivariance );

			bucketed_semivariances.insert( bucketed_semivariances.begin() + sv_insert_index, new_sv );

			// update the index of all buckets after the current SV's lag class because we've
			// added a new bucket
			for ( int j = lag_class + 1; j < num_buckets; j++ )
			{
				index_array[ j ]++;
			}

		} else {

			// If the SV is already there, then add to the current SV components
			bucketed_semivariances[ sv_insert_index ]->semivariance[ X_COMPONENT ] += sv->semivariance[ X_COMPONENT ];
			bucketed_semivariances[ sv_insert_index ]->semivariance[ Y_COMPONENT ] += sv->semivariance[ Y_COMPONENT ];

		}

		// increment the count for bucket we're inserting the SV into
		count_array[ lag_class ]++;
	}


	// Now that we have summed all semivariances for each lag class we will divide the
	// SV values by the number of contributing SVs
	for ( int i = 0; i < num_buckets; i++ )
	{
		double count = double( count_array[ i ] );

		// if there were any SVs contributing to this bucket at all
		if ( count > 0 )
		{
			int index = index_array[ i ];
			Semivariance* bucket = bucketed_semivariances[ index ];

			// All average calculations have division, this is our division
			bucket->semivariance[ X_COMPONENT ] /= count;
			bucket->semivariance[ Y_COMPONENT ] /= count;
		}
	}

	return bucketed_semivariances;
}


// This is a function that takes a vector of semivariances, a pivot index,
// a pointer to left mean value for the x-component (XLMV) and another to 
// the right mean value for the x-component (XRMV), a YLMV, and a YRMV. It
// calculates the mean of all x-semivariances up to and including the pivot
// and stores it in the (X/Y)LMV, it does the same to the right side of the
// pivot and stores the result in (X/Y)RMV.
void calculate_left_and_right_means( vector< Semivariance* > semivariances, int pivot, double* XLMV, double* XRMV, double* YLMV, double* YRMV )
{
	// if this is true we'll end up dividing the right values by 0
	if ( semivariances.size() - pivot - 1 == 0 )
	{
		return;
	}

	// set the mean variables to 0
	*XLMV = 0;
	*YLMV = 0;
	*XRMV = 0;
	*YRMV = 0;

	for ( int i = 0; i <= pivot; i++ )
	{
		*XLMV += semivariances[ i ]->semivariance[ X_COMPONENT ];
		*YLMV += semivariances[ i ]->semivariance[ Y_COMPONENT ];
	}

	*XLMV /= ( pivot + 1 );
	*YLMV /= ( pivot + 1 );

	for ( int i = pivot + 1; i < semivariances.size(); i++ )
	{
		*XRMV += semivariances[ i ]->semivariance[ X_COMPONENT ];
		*YRMV += semivariances[ i ]->semivariance[ Y_COMPONENT ];
	}

	*XRMV /= ( semivariances.size() - pivot - 1 );
	*YRMV /= ( semivariances.size() - pivot - 1 );
}

// This function takes a vector of semivariances, a pivot index,
// the XLMV, YLMV, XRMV, YRMV from the calculate_left_and_right_means
// method. It also takes pointers to the x component left variance
// (xl_variance), yl_variance, xr_variance, and yr_variance; it
// calculates each variance and assigns the value to the parameters.
void calculate_left_and_right_variances( vector< Semivariance* > semivariances, int pivot, double XLMV, double XRMV, double YLMV, double YRMV, double* xl_variance, double* yl_variance, double* xr_variance, double* yr_variance)
{
	// if this is true we'll end up dividing the right values by 0
	if ( semivariances.size() - pivot - 1 == 0 )
	{
		return;
	}

	// set the variance variables to 0
	*xl_variance = 0;
	*yl_variance = 0;
	*xr_variance = 0;
	*yr_variance = 0;


	for ( int i = 0; i <= pivot; i++ )
	{
		*xl_variance += pow( XLMV - semivariances[ i ]->semivariance[ X_COMPONENT ], 2 );
		*yl_variance += pow( YLMV - semivariances[ i ]->semivariance[ Y_COMPONENT ], 2 );
	}

	*xl_variance /= ( pivot + 1 );
	*yl_variance /= ( pivot + 1 );

	for ( int i = pivot + 1; i < semivariances.size(); i++ )
	{
		*xr_variance += pow( XRMV - semivariances[ i ]->semivariance[ X_COMPONENT ], 2 );
		*yr_variance += pow( YRMV - semivariances[ i ]->semivariance[ Y_COMPONENT ], 2 );
	}

	*xr_variance /= ( semivariances.size() - pivot - 1 );
	*yr_variance /= ( semivariances.size() - pivot - 1 );
}

// This calculates the variance-to-mean ratio (VMR) for all semivariances
// to the left and to the right of the pivot value, and for all possible
// pivot values. The function calculate the difference between VMRleft and
// VMRright, and returns the pivot index i at which the maximum difference
// is observed ( this is our SV model's range, the semi-variance with that
// index is the model's sill )
// Source: Chen, Qi, and Peng Gong. "Automatic variogram parameter extraction for textural classification of the panchromatic IKONOS imagery." Geoscience and Remote Sensing, IEEE Transactions on 42.5 (2004): 1106-1115.
void find_range_and_sill( vector < Semivariance* > semivariances, double* x_range, double* x_sill, double* x_nugget, double* y_range, double* y_sill, double* y_nugget )
{
	int x_max_index = 0;
	int y_max_index = 0;
	double x_max_value = 0;
	double y_max_value = 0;

	double XLMV, YLMV, XRMV, YRMV, xl_variance, yl_variance, xr_variance, yr_variance;

	// the nuggets are initially set to the first semivariance. We traverse the list
	// of SVs to make sure they're they represent the smallest SV in the lsit.
	*x_nugget = semivariances[ 0 ]->semivariance[ X_COMPONENT ];
	*y_nugget = semivariances[ 0 ]->semivariance[ Y_COMPONENT ];

	for( int pivot = 0; pivot < semivariances.size() - 2; pivot++ )
	{
		calculate_left_and_right_means( semivariances, pivot, &XLMV, &XRMV, &YLMV, &YRMV );
		calculate_left_and_right_variances( semivariances, pivot, XLMV, XRMV, YLMV, YRMV, &xl_variance, &yl_variance, &xr_variance, &yr_variance);

		// calculate the ratio of the variance to the mean for the left and right, and for the x and y components
		double XVMRleft = xl_variance / XLMV;
		double XVMRright = xr_variance / XRMV;

		double YVMRLeft = yl_variance / YLMV;
		double YVMRright = yr_variance / YRMV;

		// calculate the difference between the VMR to the left and right of the pivot for x and y
		double XDV = XVMRleft - XVMRright;
		double YDV = YVMRLeft - YVMRright;

		// if we have a new max then update the recorded extrema
		if ( XDV > x_max_value )
		{
			x_max_value = XDV;
			x_max_index = pivot;
		}

		if ( YDV > y_max_value )
		{
			y_max_value = YDV;
			y_max_index = pivot;
		}

		// find the smallest semivariance, it will be our nugget
		double x_SV = semivariances[ pivot ]->semivariance[ X_COMPONENT ];
		double y_SV = semivariances[ pivot ]->semivariance[ Y_COMPONENT ];

		*x_nugget = x_SV < *x_nugget ? x_SV : *x_nugget;
		*y_nugget = y_SV < *y_nugget ? y_SV : *y_nugget;
	}

	*x_range = semivariances[ x_max_index ]->displacement;
	*y_range = semivariances[ y_max_index ]->displacement;

	*x_sill = semivariances[ x_max_index ]->semivariance[ X_COMPONENT ];
	*y_sill = semivariances[ y_max_index ]->semivariance[ Y_COMPONENT ];

}

// This function takes the parameters of the exponential semivariance model.
// The function is also provided with the lag for which we will be estimating
// the semi-variance. The returned value is equal to the estimated semi-variance
// for that lag.
// The nugget is assumed to equal our smallest bucketed reading.
double interpolate_exponential( double displacement, double nugget, double sill, double range )
{
	if( displacement >= 0 && range != 0 )
	{
		double e = exp( - displacement / range );

		return nugget + sill * ( 1 - e );
	}
}

// Matrix Multiplication
void matrix_multiply( double** A, double** B, double** C, int n, int m )
{
	int rowA, colA, rowB, colB;

	for ( rowA = 0; rowA < n; rowA++ )
	{
		for ( colB = 0; colB < m; colB++ )
		{
			C[rowA][colB] = 0;

			for ( colA = 0; colA < n; colA++ )
			{
				C[ rowA ][ colB ] += A[ rowA ][ colA ] * B[ colA ][ colB ];
			}
		}
	}
}

// Print matrix as string
void matrixPrint( double** M, int n, int m )
{
	int i, j;

	for ( i = 0; i < n; i++ )
	{
		for ( j = 0; j < m; j++ )
		{
			cout << M[i][j] << ", ";
		}
		cout << "\n";
	}
}


// Matrix inversion
// This is the work of Paul Bourke
// taken from: http://www.cs.rochester.edu/~brown/Crypto/assts/projects/adj.html

/*
   Recursive definition of determinant using expansion by minors.
*/
double determinant(double **a,int n)
{
   int i,j,j1,j2;
   double det = 0;
   double **m = NULL;

   if (n < 1) { /* Error */

   } else if (n == 1) { /* Shouldn't get used */
      det = a[0][0];
   } else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
      	m = new double* [ n - 1 ];
         for (i=0;i<n-1;i++)
        	m[ i ] = new double [ n - 1 ];
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         det += pow(-1.0,j1+2.0) * a[0][j1] * determinant(m,n-1);
         for (i=0;i<n-1;i++)
            delete [] m[ i ];
         delete [] m;
      }
   }
   return(det);
}

/*
   Find the cofactor matrix of a square matrix
*/
void coFactor(double **a,int n,double **b)
{
   int i,j,ii,jj,i1,j1;
   double det;
   double **c;

   c = new double*[ n ];
   for (i=0;i<n;i++)
	 c[i] = new double[ n ];

   for (j=0;j<n;j++) {
      for (i=0;i<n;i++) {

         /* Form the adjoint a_ij */
         i1 = 0;
         for (ii=0;ii<n;ii++) {
            if (ii == i)
               continue;
            j1 = 0;
            for (jj=0;jj<n;jj++) {
               if (jj == j)
                  continue;
               c[i1][j1] = a[ii][jj];
               j1++;
            }
            i1++;
         }

         /* Calculate the determinate */
         det = determinant(c,n-1);

         /* Fill in the elements of the cofactor */
         b[i][j] = pow(-1.0,i+j+2.0) * det;
      }
   }
   for (i=0;i<n;i++)
      delete [] c[ i ];
   delete [] c;
}

/*
   Transpose of a square matrix, do it in place
*/
void transpose(double **a,int n)
{
   int i,j;
   double tmp;

   for (i=1;i<n;i++) {
      for (j=0;j<i;j++) {
         tmp = a[i][j];
         a[i][j] = a[j][i];
         a[j][i] = tmp;
      }
   }
}

// A matrix invert function that applies the above 3 algorithms in order
void matrixInverse(double **matrix, double **inverse, int n)
{
	int i, j;
	double det;	
	
	det = determinant( matrix, n );

	coFactor( matrix, n, inverse );

	// we don't need to transpose since the kriging matrix is symmetric
	for ( i = 0; i < n; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			inverse[i][j] = inverse[i][j] / det;
		}
	}
}

// This method takes the parameters for the wind model, and the input wind data.
// It also takes the coordinates of a point and returns the wind vector component
// for that point.
double find_wind_component_at_point( int wind_component, double x, double y, double range, double sill, double nugget, vector< WindVector* > wind_data )
{
	// Create the matrices that are involved in the calculation of the wind data
	int wind_data_size = wind_data.size();
	double** inter_sv_matrix = new double* [ wind_data_size ];
	double** inter_sv_inverse = new double* [ wind_data_size ];
	double** zero_sv_matrix = new double* [ wind_data_size ];
	double** weight_matrix = new double* [ wind_data_size ];

	// fill the matrices with the semivariances for the relevant points
	for ( int i = 0; i < wind_data_size; i++ )
	{
		inter_sv_matrix[ i ] = new double [ wind_data_size ];
		inter_sv_inverse[ i ] = new double [ wind_data_size ];
		zero_sv_matrix[ i ] = new double [ 1 ];
		weight_matrix[ i ] = new double [ 1 ];

		double x1 = wind_data[ i ]->x;
		double y1 = wind_data[ i ]->y;

		for ( int j = 0; j < wind_data_size; j++ )
		{
			// calculate semivariance between points i and j (i.e. sv_i,j)
			double x2 = wind_data[ j ]->x;
			double y2 = wind_data[ j ]->y;

			double displacement = sqrt( pow( x1 - x2, 2 ) + pow( y1 - y2, 2 ) );

			double temp_sv = interpolate_exponential( displacement, nugget, sill, range );

			inter_sv_matrix[ i ][ j ] = temp_sv;
		}

		double displacement = sqrt( pow( x1 - x, 2 ) + pow( y1 - y, 2 ) );

		double temp_sv = interpolate_exponential( displacement, nugget, sill, range );

		zero_sv_matrix[ i ][ 0 ] = temp_sv;
	}

	// now we want to find the weights. Invert the inter_sv_matrix and multiply it with
	// the zero_sv_matrix from the left
	matrixInverse( inter_sv_matrix, inter_sv_inverse, wind_data_size );

	matrix_multiply( inter_sv_inverse, zero_sv_matrix, weight_matrix, wind_data_size, 1 );

	double wind_reading = 0;

	for ( int i = 0; i < wind_data_size; i++ )
	{
		wind_reading += ( weight_matrix[ i ][ 0 ] * wind_data[ i ]->velocity[ wind_component ] );
	}

	// free memory of matrices involved
	for (int i = 0; i < wind_data_size; i++ )
	{
		delete [] inter_sv_matrix[ i ];
		delete [] inter_sv_inverse[ i ];
		delete [] zero_sv_matrix[ i ];
		delete [] weight_matrix[ i ];
	}

	delete [] inter_sv_inverse;
	delete [] inter_sv_matrix;
	delete [] zero_sv_matrix;
	delete [] weight_matrix;

	return wind_reading;
}

double*** calculate_wind_at_all_points( int width, int height, double x_range, double x_sill, double x_nugget, double y_range, double y_sill, double y_nugget, vector< WindVector* > wind_data )
{
	vector< WindVector* > complete_wind_data;
	double*** windMaps = new double** [ 2 ];
	windMaps[ 0 ] = new double*[ height ];
	windMaps[ 1 ] = new double*[ height ];

	for ( int row = 0; row < height; row++ )
	{
		windMaps[ 0 ][ row ] = new double[ width ];
		windMaps[ 1 ][ row ] = new double[ width ];

		for ( int col = 0; col < width; col++ )
		{
			windMaps[0][row][col] = find_wind_component_at_point( X_COMPONENT, col, row, x_range, x_sill, x_nugget, wind_data );
			windMaps[1][row][col] = find_wind_component_at_point( Y_COMPONENT, col, row, y_range, y_sill, y_nugget, wind_data );
		}
	}

	return windMaps;
}


// Create a width x height x 8 matrix that holds the weights for each cell and its 8 neighbours
double*** create_cost_matrix( int width, int height, double*** wind_model )
{
	double*** cost_matrix = new double** [ height ];

	double max_wind = 0;

	// find the maximum wind reading to avoid having negative weights
	for ( int i = 0; i < height; i++ )
	{
		for ( int j = 0; j < width; j++ )
		{
			double current_wind = sqrt( pow( wind_model[0][i][j], 2 ) + pow( wind_model[1][i][j], 2 ) );
		
			max_wind = max_wind > current_wind ? max_wind : current_wind;	
		}
	}

	// set the weights for each direction for each node
	for ( int i = 0; i < height; i++ )
	{
		cost_matrix[ i ] = new double* [ width ];

		for ( int j = 0; j < width; j++ )
		{
			cost_matrix[ i ][ j ] = new double [ 9 ];

			for ( int k = -1; k < 2; k++ )
			{
				for ( int l = -1; l < 2; l++ )
				{
					int new_row = i + k;
					int new_col = j + l;
					int direction_index = (k + 1) * 3 + (l + 1);

					if ( !( k==0 && l==0 ) && ( new_row >= 0 && new_row < height ) && ( new_col >= 0 && new_col < width ))
					{
						double current_wind = sqrt( pow( wind_model[0][i][j], 2 ) + pow( wind_model[1][i][j], 2 ) );
						double displacement = direction_index % 2 == 0 ? sqrt( 2 ) : 1;

						// find the projection of the wind onto the path of movement
						double wind_projection_along_path = l * wind_model[0][i][j] + k * wind_model[1][i][j];
						int sign = wind_projection_along_path >= 0 ? 1 : -1;

						// find the angle between the wind and the path of movement
						double theta = acos( wind_projection_along_path / ( displacement * current_wind ));

						// find the projection of the wind perpendicular to the path of movement
						double wind_projection_perpendicular_to_path = current_wind * displacement * cos( 180 - theta );

						cost_matrix[ i ][ j ][ direction_index ] = 1 + sign * pow( wind_projection_along_path / ( displacement * max_wind ), 2 ) + 1 + pow( wind_projection_perpendicular_to_path / ( displacement * max_wind ), 2);
					}
				}
			}
		}
	}

	return cost_matrix;
}

// given the weight matrix find the shortest path linking the source and destination points
// using Dijkstra's algorithm
// I followed (and modified) the wikipedia pseudocode for Dijkstra: http://en.wikipedia.org/wiki/Dijkstra%27s_algorithm
vector< vector< int > > find_shortest_path( int width, int height, int source_x, int source_y, int dest_x, int dest_y, double*** weight_matrix )
{
	int num_properties = 3;
	enum property { DISTANCE, PREV_X, PREV_Y };
	double*** properties_matrix = new double** [ num_properties ];
	vector<int> Q;

	for ( int property = 0; property < num_properties; property++ )
	{
		properties_matrix[ property ] = new double* [ height ];
	}

	// set up the default starting values for each node
	for ( int i = 0; i < height; i++ )
	{
		for ( int property = 0; property < num_properties; property++ )
		{
			properties_matrix[ property ][ i ] = new double [ width ];
		}

		for ( int j = 0; j < width; j++ )
		{
			properties_matrix[ DISTANCE ][ i ][ j ] = INFINITY;
			properties_matrix[ PREV_X ][ i ][ j ] = -1;
			properties_matrix[ PREV_Y ][ i ][ j ] = -1;

			Q.push_back( i * height + j );
		}
	}

	properties_matrix[ DISTANCE ][ source_y ][ source_x ] = 0;


	// implement Dijkstra
	while ( Q.size() != 0 )
	{
		int min_row = Q[ 0 ] / height;
		int min_col = Q[ 0 ] % width;
		int min_index = 0;

		double min_dist = properties_matrix[ DISTANCE ][ min_row ][ min_col ];

		for ( int i = 0; i < Q.size(); i++ )
		{
			int row = Q[ i ] / height;
			int col = Q[ i ] % width;
			double current_dist = properties_matrix[ DISTANCE ][ row ][ col ];


			if ( min_dist > current_dist )
			{
				min_row = row;
				min_col = col;
				min_index = i;
				min_dist = current_dist;
			}
		}

		Q.erase( Q.begin() + min_index );

		if ( min_dist == INFINITY )
		{
			break;
		}

		for ( int k = -1; k < 2; k++ )
		{
			for ( int l = -1; l < 2; l++ )
			{
				int new_row = min_row + k;
				int new_col = min_col + l;
				int direction_index = (k + 1) * 3 + (l + 1);

				if ( !( k==0 && l==0 ) && ( new_row >= 0 && new_row < height ) && ( new_col >= 0 && new_col < width ))
				{
					double alt_dist = min_dist + weight_matrix[ new_row ][ new_col ][ direction_index ];
					double v_dist = properties_matrix[ DISTANCE ][ new_row ][ new_col ];

					if ( alt_dist < v_dist )
					{
						properties_matrix[ DISTANCE ][ new_row ][ new_col ] = alt_dist;
						properties_matrix[ PREV_X ][ new_row ][ new_col ] = min_col;
						properties_matrix[ PREV_Y ][ new_row ][ new_col ] = min_row;
					}
				}
			}
		}
	}

	vector< int > x_path;
	vector< int > y_path;

	x_path.push_back( dest_x );
	y_path.push_back( dest_y );

	int prev_x = properties_matrix[ PREV_X ][ dest_y ][ dest_x ];
	int prev_y = properties_matrix[ PREV_Y ][ dest_y ][ dest_x ];

	while ( !( prev_x == source_x ) && !( prev_y == source_y ))
	{
		x_path.insert( x_path.begin(), prev_x );
		y_path.insert( y_path.begin(), prev_y );

		prev_x = properties_matrix[ PREV_X ][ prev_y ][ prev_x ];
		prev_y = properties_matrix[ PREV_Y ][ prev_y ][ prev_x ];
	}

	x_path.insert( x_path.begin(), source_x );
	y_path.insert( y_path.begin(), source_y );

	vector< vector< int > > path;
	path.push_back( x_path );
	path.push_back( y_path );

	return path;
}


int main(int argc, char * argv[])
{
	// double velocity[] = { 9.0, 1.2 };
	// WindVector* wv = new WindVector( 1, 1, velocity );

	// cout << "testing WindVector class:\n";
	// cout << "x: " << wv->x << "\ny: " << wv-> y << "\n";
	// for (int i = 0; i < 2; i++ )
	// {
	// 	cout << "velocity[" << i << "]: " << wv->velocity[i] << "\n";
	// }	


	// double gamma[] = { 1.2, 3.4 };
	// Semivariance* sv = new Semivariance( 3.1, gamma );

	// cout << "Testing semivariance class: \n";
	// cout << "displacement: " << sv->displacement << "\n";
	// for (int i = 0; i < 2; i++ )
	// {
	// 	cout << "semivariance[" << i << "]: " << sv->semivariance[i] << "\n";
	// }


	// // testing the semivariance method
	// // create another wind vector in order to establish a semivariance with wv
	// double velocity2[] = { 4.2, 3.3 };
	// WindVector* wv2 = new WindVector( 5, 7, velocity2 );
	// double velocity3[] = { 4.4, 9.1 };
	// WindVector* wv3 = new WindVector( 3, 12, velocity3 );
	// WindVector* wv4 = new WindVector( 12, 5, velocity3 );

	// vector< WindVector* > test_input_vectors;
	// test_input_vectors.push_back( wv );
	// test_input_vectors.push_back( wv2 );
	// test_input_vectors.push_back( wv3 );
	// test_input_vectors.push_back( wv4 );

	// vector< Semivariance* > semivariance_vector = calculate_all_semivariances( test_input_vectors );
	// cout << "Testing the semivariance calculation method\n";
	// cout << "The number of semivariances calculated: " << semivariance_vector.size();
	// cout << "\nThe lag of the first semivariance: " << semivariance_vector[0]->displacement;
	// cout << "\nThe value of the northern semivariance: " << semivariance_vector[0]->semivariance[NORTH] << "\n";
	

	// Testing the bucketing function with 6 semivariances to which I've precalculated the 
	// correct output as: 
	// x1 = 56.7, y1 = 42
	// x2 = 50.3, y2 = 58.3

	// double component_sv1[] = { 12, 34 };
	// double component_sv2[] = { 90, 13 };
	// double component_sv3[] = { 68, 79 };
	// double component_sv4[] = { 43, 65 };
	// double component_sv5[] = { 17, 28 };
	// double component_sv6[] = { 91, 82 };

	// Semivariance* sv1 = new Semivariance( 2, component_sv1 );
	// Semivariance* sv2 = new Semivariance( 3, component_sv2 );
	// Semivariance* sv3 = new Semivariance( 4, component_sv3 );
	// Semivariance* sv4 = new Semivariance( 5, component_sv4 );
	// Semivariance* sv5 = new Semivariance( 6, component_sv5 );
	// Semivariance* sv6 = new Semivariance( 9, component_sv6 );

	// Semivariance* sv_array[] = { sv1, sv2, sv3, sv4, sv5, sv6 };
	// vector< Semivariance* > sv_collection( sv_array, sv_array + sizeof( sv_array ) / sizeof( Semivariance* ) );

	// vector <Semivariance*> average_SVs = bucket_sort( 10, 5, sv_collection );

	// cout << "The components of the first SV are:\n";
	// cout << "X_COMPONENT: " << average_SVs[0]->semivariance[ X_COMPONENT ] << ", Y_COMPONENT: " <<  average_SVs[0]->semivariance[ Y_COMPONENT ] <<  "\n";
	// cout << "The components of the second SV are:\n";
	// cout << "X-COMPONENT: " << average_SVs[1]->semivariance[ X_COMPONENT ] << ", Y-COMPONENT: " <<  average_SVs[1]->semivariance[ Y_COMPONENT ] << "\n";


	// // Testing the graph fitting functions
	// // The output should be the same as the above result
	// double XLMV, YLMV, XRMV, YRMV;

	// calculate_left_and_right_means( sv_collection, 2, &XLMV, &XRMV, &YLMV, &YRMV );

	// cout << "After calculating the left and right means we've found:\n";
	// cout << "XMLV = " << XLMV << ", YLMV = " << YLMV << ", XRMV = " << XRMV << ", YRMV = " << YRMV << "\n";

	// double xl_variance, yl_variance, xr_variance, yr_variance;

	// calculate_left_and_right_variances( sv_collection, 2, XLMV, XRMV, YLMV, YRMV, &xl_variance, &yl_variance, &xr_variance, &yr_variance );

	// cout << "After calculating the left and right variances we've found:\n";
	// cout << "xl_variance = " << xl_variance << ", yl_variance = " << yl_variance << ", xr_variance = " << xr_variance << ", yr_variance = " << yr_variance << "\n";

	// double x_range, y_range, x_sill, y_sill;

	// find_range_and_sill( sv_collection, &x_range, &x_sill, &y_range, &y_sill );

	// cout << "After calling the find_range_and_sill function on the sv_collection we got:\n";
	// cout << "x_range = " << x_range << ", y_range = " << y_range << ", x_sill = " << x_sill << ", y_sill = " << y_sill << "\n"; 


	// // testing matrix multiplication
	// int i, j, n;
	// double **A, **B, **C, **Ci;

	// n = 3;

	// A = new double*[ n ];
	// B = new double*[ n ];
	// C = new double*[ n ];

	// for ( i = 0; i < n; i++ )
	// {
	// 	A[ i ] = new double[ n ];
	// 	B[ i ] = new double[ n ];
	// 	C[ i ] = new double[ n ];
	// }

	// A[ 0 ][ 0 ] = 1;
	// A[ 0 ][ 1 ] = 4;
	// A[ 0 ][ 2 ] = 1;

	// A[ 1 ][ 0 ] = 4;
	// A[ 1 ][ 1 ] = 2;
	// A[ 1 ][ 2 ] = 0;

	// A[ 2 ][ 0 ] = 1;
	// A[ 2 ][ 1 ] = 0;
	// A[ 2 ][ 2 ] = -5;

	// matrix_multiply( A, A, C, n, n );	

	// cout << "Testing matrix multiplication method\nMultiplying A with itself. A:\n";
	// matrixPrint( A, n );
	// cout << "result of multiplication:\n";
	// matrixPrint( C, n );

	// cout << "Testing matrix inverse method. The inverse matrix is:\n";
	// cout << "Remember I'm not transposing the matrix because we only have symmetric matrices in our application\n";

	// matrixInverse( A, B, n );
	// matrixPrint( B, n );
	// matrix_multiply( A, B, C, n, n );

	// cout << "This should be the identity matrix after multiplying the matrix with its inverse:\n";
	// matrixPrint( C, n );


	// Testing complete functionality on a 10x12 map
	// first create the wind vectors
	int size = 6;
	const double Xs[] = { 1, 3, 3, 7, 9, 9 };
	const double Ys[] = { 7, 11, 1, 9, 4, 1 };
	const double windX[] = { 6, -3, 1, -3, 12, 7 };
	const double windY[] = { 1, 1, 4, -10, 8, 2 };

	vector<WindVector*> wind_vectors;

	cout << "About to start filling wind_vector\n";

	for ( int i = 0; i < size; i++ )
	{
		double* temp_wind = new double [2];
		temp_wind[ 0 ] = windX[ i ];
		temp_wind[ 1 ] = windY[ i ];

		WindVector* wv = new WindVector( Xs[ i ], Ys[ i ], temp_wind );
		wind_vectors.push_back( wv );
	}

	cout << "created wind_vectors\n";

	vector< Semivariance* > semivariance_list = calculate_all_semivariances( wind_vectors );
	vector< Semivariance* > bucketed_semivariances;
	int bucket_size = 1; // default bucket
	double max_displacement = 0;
	double x_range, x_sill, x_nugget, y_range, y_sill, y_nugget;

	for( int i=0; i < semivariance_list.size(); i++ )
	{
		double current_displacement = semivariance_list[ i ]->displacement;
		max_displacement = current_displacement > max_displacement ? current_displacement : max_displacement;
	}

	bucketed_semivariances = bucket_sort( int( max_displacement ), bucket_size, semivariance_list );

	cout << "Bucketed semivariances\n";

	find_range_and_sill( bucketed_semivariances, &x_range, &x_sill, &x_nugget, &y_range, &y_sill, &y_nugget );

	cout << "Modelled semivariances\n";

	double*** complete_wind_model = calculate_wind_at_all_points( 10, 12, x_range, x_sill, x_nugget, y_range, y_sill, y_nugget, wind_vectors );

	cout << "Calculated whole wind model\n";

	//matrixPrint( complete_wind_model[0], 12, 10 );
	//matrixPrint( complete_wind_model[1], 12, 10 );

	// Testing out the weighting for the wind_model
	double*** weights_matrix = create_cost_matrix( 10, 12, complete_wind_model );

	cout << "The weights for the wind_model are:\n\n";

	// for ( int i = 0; i < 12; i++ )
	// {
	// 	for ( int j = 0; j < 10; j++ )
	// 	{
	// 		cout << "looking at cell " << i << "x" << j << ": \n"; 

	// 		for ( int k = -1; k < 2; k++ )
	// 		{
	// 			for ( int l = -1; l < 2; l++ )
	// 			{
	// 				int new_row = i + k;
	// 				int new_col = j + l;
	// 				int direction_index = (k + 1) * 3 + (l + 1);

	// 				if ( !( k==0 && l==0 ) && ( new_row >= 0 && new_row < 12 ) && ( new_col >= 0 && new_col < 10 ))
	// 				{
	// 					cout << weights_matrix[ i ][ j ][ direction_index ] << "; ";
	// 				}
	// 			}

	// 			cout << "\n";
	// 		}

	// 		cout << "-----------\n";
	// 	}
	// }

	// testing the dijkstra function
	cout << "This function returns the dijkstra shortest path from point (2, 1) to (7, 8)\n";

	vector< vector< int > > path =  find_shortest_path( 10, 12, 2, 1, 7, 8, weights_matrix );

	for ( int i = 0; i < path[ 0 ].size(); i++ )
	{
		cout << "The " << i << "th node is: (" << path[ 0 ][ i ] << ", " << path[ 1 ][ i ] << ")\n";
	}
}