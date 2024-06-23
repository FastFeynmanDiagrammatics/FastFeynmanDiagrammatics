#!/usr/bin/env bash

g++ -std=c++17 -O0 -Wall -g -DFFD_CHEBYSHEVPOLYNOMIAL_UNIT_TEST_FLAG -I ../../../../ UnitTest.cpp -o UnitTest.x
./UnitTest.x
