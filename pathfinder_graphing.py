import numpy as np
import matplotlib.pyplot as plt
import pylab as P

path = []
wind = []
recordings = []

# open the wind and path files and read their content unto a programmable datastructure
wind_file = open( "output_wind.txt", 'r' ).readlines()
path_file = open( "output_path.txt", 'r' ).readlines()
input_file = open( "input.txt", 'r' ).readlines()

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

for line in input_file[3:]:
	x, y = line.split(' ')[0:2]

	recordings.append( ( int(x), int(y)) )


# find the maximum wind speed
max_wind = 0.0

for row in range(len(wind)):
	for col in range(len(wind[0])):
		current_wind = ( wind[ row ][ col ][ 0 ]**2 + wind[ row ][ col ][ 1 ]**2 )**0.5
		max_wind = current_wind if current_wind > max_wind else max_wind


# draw the shortest path on the graph
# FIRST: draw red circles on each point
# for row in range(len(wind)):
# 	for col in range(len(wind[0])):
# 		plt.Circle((col, row), radius=0.1, color='r')

# SECOND: draw green circles on the points included in the path
# and link those points with lines
for i in range(len(path)):
	point = path[ i ]
	plt.Circle((point[0], point[1]), radius=0.1, color='g')

	if (i + 1 < len(path)):
		next_point = path[ i + 1 ]
		actual_dx = next_point[ 0 ] - point[ 0 ]
		actual_dy = next_point[ 1 ] - point[ 1 ]

		dx = actual_dx * 3 / 4.0
		dy = actual_dy * 3 / 4.0

		head_length = ( ( actual_dx / 4.0 )**2 + ( actual_dy / 4.0 )**2 )**0.5
		head_width = head_length / 2.0

		kwargs = { 'color' : '#000000' }

		P.arrow( point[0], point[1], dx, dy, head_width=head_width, head_length=head_length, alpha=0.25, **kwargs )


	# draw the arrows across the graph
for row in range(len(wind)):
	for col in range(len(wind[0])):
		dx = wind[row][col][0] / max_wind
		dy = wind[row][col][1] / max_wind

		head_length = ( dx**2 + dy**2 )**0.5 / 3.0
		head_width = head_length / 1.5

		kwargs = { 'color' : '#0FD67D' }

		if (col, row) in recordings:
			kwargs = { 'color' : '#D60F58' }

		P.arrow( col, row, dx * (2 / 3.0), dy * (2 / 3.0), head_width=head_width, head_length=head_length, **kwargs )

# set plot parameters
P.xlim( -1.5, len(wind[0]) + 0.5 )
P.ylim( -1.5, len(wind) + 0.5 )

# P.draw()
P.show()
# plt.savefig('testgraph.pdf')