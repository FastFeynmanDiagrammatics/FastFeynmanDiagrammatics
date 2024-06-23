#!/usr/bin/env bash

g++ -std=c++17 -Wall -O0 -I /Users/rrossi/Dropbox\ \(Simons\ Foundation\)/projects/18_fermi/ test.cpp -o test.x
g++ -std=c++17 -Wall -Ofast -DNDEBUG -I /Users/rrossi/Dropbox\ \(Simons\ Foundation\)/projects/18_fermi/ test.cpp -o FastTest.x
./test.x
