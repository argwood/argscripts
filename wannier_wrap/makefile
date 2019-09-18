OBJECTS = main.o read_xyz.o wannier.o
HEADERS = read_xyz.h wannier.h
CFLAGS = -Wall -g -std=c++11 
LDFLAGS = -lm
CC = g++

default: main

main: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)

%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

clean: 
	rm -f $(OBJECTS) main
