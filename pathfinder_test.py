from __future__ import print_function
import random
import os
import subprocess
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-dim", type=str, help="-dim specifies the dimensions of the map in question using any discrete integral values. e.g. \"100, 100\"",
		default="100, 100")
parser.add_argument("-source", type=str, help="-source specifies the coordinates of the path's starting point within the grid defined in -dim. e.g. \"20, 80\"",
		default="0,0")
parser.add_argument("-dest", type=str, help="-dest specifies the coordinates of the path's finishing point within the grid defined in -dim. e.g. \"20, 80\"",
		default="0,0")
parser.add_argument("-sim", type=int, help="-sim accepts truth values of 0 or 1 (1 being true). It state whether you want the program to create random input wind readings to a file 'input.txt'",
		choices=[0,1], default=0)
parser.add_argument("--sim-num", type=int, help="If -sim has a truth value of 1 then --sim-num determines the number of random input readings to write to 'input.txt'",
		default=10)

args = parser.parse_args()

dim = args.dim.split(',')
map_width = int( dim[ 0 ] )
map_height = int( dim[ 1 ] )

source = args.source.split(',')
source_x = int( source[ 0 ] )
source_y = int( source[ 1 ] )

dest = args.dest.split(',')
dest_x = int( dest[ 0 ] )
dest_y = int( dest[ 1 ] )

# if the user wants to simulate the input readings then create the input file
if ( args.sim == 1 ):
	input_file = open("input.txt", 'w')
	print( "{0} {1}".format( map_width, map_height ), file=input_file )
	print( "{0} {1}".format( source_x, source_y ), file=input_file )
	print( "{0} {1}".format( dest_x, dest_y ), file=input_file )

	reading_coordinates = []

	# find |sim_num| unique coordinates
	for i in range( args.sim_num ) :
		x = random.randint( 0, map_width - 1 )
		y = random.randint( 0, map_height - 1 )

		while ( (x, y) in reading_coordinates ):
			x = random.randint( 0, map_width )
			y = random.randint( 0, map_height )

		reading_coordinates.append( (x, y) )

		# add a random wind reading for each of the pre-specified coordinates
		x_wind = random.uniform( -40, 40 )
		y_wind = random.uniform( -40, 40 )

		print( "{0} {1} {2} {3}".format( x, y, x_wind, y_wind ), file=input_file )

	input_file.close()

# compile the kriging/pathfinder cpp file
subprocess.call( "make all", shell=True )
subprocess.call( "./sequential_kriging", shell=True )
os.system( "rm ./sequential_kriging" )