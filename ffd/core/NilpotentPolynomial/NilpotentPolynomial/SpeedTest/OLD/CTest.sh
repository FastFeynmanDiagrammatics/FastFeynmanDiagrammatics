#!/usr/bin/env bash

clang++ -std=c++17 -Ofast test.cpp -o ctest.x
./ctest.x
