#!/usr/bin/env bash

g++ -std=c++17 -Ofast -Wpedantic -Wall -DNDEBUG -g -I /Users/rrossi/Dropbox\ \(Simons\ Foundation\)/projects/18_fermi/ HubbardAtom.cpp -o Atom.x
clang++ -std=c++17 -Ofast -Wpedantic -Wall -DNDEBUG -g -I /Users/rrossi/Dropbox\ \(Simons\ Foundation\)/projects/18_fermi/ HubbardAtom.cpp -o CAtom.x
./Atom.x $1 $2

