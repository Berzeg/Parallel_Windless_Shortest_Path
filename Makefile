CC := g++
TARGET := test

SOURCE := testing.cpp WindVector.cpp
FIXSOURCE := $(addprefix src/, $(SOURCE))

VPATH := src : include

vpath %.h include
vpath %.cpp src

## Clean Rule
clean:
	@rm -f $(TARGET) $(OBJECTS)
 
all: $(SOURCE)
	$(CC) $(FIXSOURCE) -o $(TARGET)