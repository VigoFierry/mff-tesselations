Author: Milan Pultar
Email: milan.pultar@gmail.com
Date: 2017-07-11

Provided is the library to create and modify power tessellation. You may override any class to suit your needs.
Inheritance is a nice basic c++ feature that enables this. Voronoi tessellation is a special case of power tessellation
and can be created by adding all grains with fixed radius. Other distance functions may be used to construct the corresponding
tessellation if you override the (Grain)DistanceTo and (Plane)ComputeD methods. Adding a grain takes approx. 24 ms on
average laptop while removing a grain is somewhat slower.

You may use the method (Grain)to_connected() to get the connected representation (no duplicates) of the Grain.
It is safe to insert/remove a grain to the tessellation after you transformed some of its grains to the connected form (presented in the example).

Included is one example as well as both Linux and Windows builds (debug and release).
Linux build compiler: g++ version 6.3.0, -std=c++11
Windows build: MS Visual Studio 2017.

MatlabView folder contains a script to view the output of the program. Try it on viewMe.txt.

example23bitRelease.exe is a built exampleWin.cpp code on Windows32bit that demonstrates the functions of the library.
If you try to compile the exampleWin.cpp code yourself, be aware that you have to link the library correctly.

exampleLinux is a build of exampleLinux.cpp.
