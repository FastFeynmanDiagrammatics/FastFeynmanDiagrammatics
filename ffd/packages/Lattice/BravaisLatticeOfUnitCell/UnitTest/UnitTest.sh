#!/usr/bin/env bash

g++ -std=c++17 -O0 -Wall -DFFD_LATTICE_UNIT_TEST_FLAG -I $FFD_LIB_PATH UnitTest.cpp -o UnitTest.x
./UnitTest.x
