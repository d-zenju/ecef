#!/bin/sh

g++ -std=c++11 -o ecef/ecef.o -c ecef/ecef.cpp
g++ -std=c++11 ecef/ecef.o main.cpp 
