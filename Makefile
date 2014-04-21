CC := g++
TARGET := sequential_kriging

SOURCE := sequential_kriging.cpp WindVector.cpp Semivariance.cpp
FIXSOURCE := $(addprefix src/, $(SOURCE))

VPATH := src : include

vpath %.h include
vpath %.cpp src

## Clean Rule
clean:
	@rm -f $(TARGET) $(OBJECTS)
 
all: $(SOURCE)
	# $(CC) $(FIXSOURCE) -o $(TARGET) 
	$(CC) $(FIXSOURCE) -o $(TARGET) -std=gnu++0x