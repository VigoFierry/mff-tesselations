# on Linux OS open terminal in this directory and type "make", the following will make the binary

# set variables:
CXX=g++
CFLAGS = -std=c++11
.PHONY: all

all: exampleLinux clean
# link binaries into one file:
exampleLinux: exampleLinux.o Linux/Release/Tessellation.so
	$(CXX) $(CFLAGS) -o $@ $^
# create exampleLinux binary:
exampleLinux.o: exampleLinux.cpp
	$(CXX) $(CFLAGS) -o $@ -c $^
# delete support binaries:
clean:
	rm -f exampleLinux.o
