#!/usr/bin/env bash

g++ -std=c++17 -Ofast test.cpp -o test.x
./test.x $1
