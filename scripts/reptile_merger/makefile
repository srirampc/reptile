# Makefile for reptile_merger, 15 Sept 2010

SRCDIR = .
CXX = g++
CXXFLAGS = -Wall -I$(SRCDIR)
LIBS = -lz

.PHONY: all
all : reptile_merger

.PHONY: clean
clean: 
	rm -f reptile_merger

reptile_merger:
	$(CXX) $(CXXFLAGS) reptile_merger.cpp -o $@ $(LIBS)

