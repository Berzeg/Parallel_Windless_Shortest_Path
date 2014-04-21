from __future__ import print_function
import random
import os
import subprocess
import re
# import numpy as np
# import matplotlib as plt

parameters = []
times = {}

dimensions = [(10,10), (10,20), (15,20), (20, 20), (20,25)]
input_readings = [ 3, 4, 5, 6, 7, 8 ]

# compile the cpp file once
subprocess.call( "make all", shell=True )

# dim is the size of the map (x, y)
# source is the coordinate of the starting point (x, y)
# dest is the coordinate of the destination point (x, y)
# input_size is an integer representing the # of input recorded values
def create_input_file( dim, source, dest, input_size ):
	input_file = open("input.txt", 'w')
	print( "{0} {1}".format( dim[0], dim[1] ), file=input_file )
	print( "{0} {1}".format( source[0], source[1] ), file=input_file )
	print( "{0} {1}".format( dest[0], dest[1] ), file=input_file )

	reading_coordinates = []

	# find |sim_num| unique coordinates
	for i in range( input_size ) :
		x = random.randint( 0, dim[0] - 1 )
		y = random.randint( 0, dim[1] - 1 )

		while ( (x, y) in reading_coordinates ):
			x = random.randint( 0, dim[0] )
			y = random.randint( 0, dim[1] )

		reading_coordinates.append( (x, y) )

		# add a random wind reading for each of the pre-specified coordinates
		x_wind = random.uniform( -40, 40 )
		y_wind = random.uniform( -40, 40 )

		print( "{0} {1} {2} {3}".format( x, y, x_wind, y_wind ), file=input_file )

	input_file.close()


# compile and run the kriging/pathfinder cpp file
def run_pathfinder():
	subprocess.call( "./sequential_kriging", shell=True )

def read_times():
	times_file = open( "output_times.txt", 'r' ).readlines()

	for line in times_file:
		name, time = line.split(',')

		if name not in times:
			times[ name ] = []

		times[ name ].append( int(time) )

# def draw_barchart():
# 	# Plot bars and create text labels for the table
# 	cell_text = []
# 	for i in range(len(times[times.keys()[0]])):
# 		y_offset = 0

# 		for method in times:
# 		    plt.bar(0, method[ i ], 0.5, bottom=y_offset )#color=colors[row])
# 		    y_offset = y_offset + method[ i ]
# 		    cell_text.append( '%1.1f' % x for x in method[ i ])
# 	# Reverse colors and text labels to display the last value at the top.
# 	# colors = colors[::-1]
# 	cell_text.reverse()

# 	the_table = plt.table(cellText=cell_text,
# 	                      # rowLabels=rows,
# 	                      # rowColours=colors,
# 	                      # colLabels=columns,
# 	                      loc='bottom')

# 	plt.subplots_adjust(left=0.2, bottom=0.2)

# 	plt.ylabel("Time (\mu s)")
# 	# plt.yticks(values * value_increment, ['%d' % val for val in values])
# 	# plt.xticks([])
# 	plt.title("Time Elapsed to Compute the Optimum Path\n with a 2D graph containing 500 nodes")

# 	plt.show()

# for dim in dimensions:
# 	for input_size in input_readings:	
# 		create_input_file( dim, (0, 0), ( dim[0] - 1, dim[1] - 1 ), input_size )
# 		run_pathfinder()
# 		read_times()

# for input_size in input_readings:	
# 	create_input_file( (20,25), (0, 0), ( 19, 24 ), input_size )
# 	run_pathfinder()
# 	read_times()


# print( "times: " + str(times) )

# # remove any leftover clutter
# os.system( "make clean" )