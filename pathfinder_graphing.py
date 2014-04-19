import numpy as np
import matplotlib.pyplot as plt

path = []
wind = []

# open the wind and path files and read their content unto a programmable datastructure
wind_file = open( "output_wind.txt", 'r' ).readlines()
path_file = open( "output_path.txt", 'r' ).readlines()

for line in wind_file[1:]:
	row_of_wind_vectors = []

	for wind_reading in line.split(' '):
		if ( wind_reading != '\n' ):
			x, y = wind_reading.split(',')

			row_of_wind_vectors.append( ( float(x) , float(y) ) )

	wind.append( row_of_wind_vectors )

for line in path_file:
	x, y = line.split(',')
	path.append( ( int(x), int(y) ) )

