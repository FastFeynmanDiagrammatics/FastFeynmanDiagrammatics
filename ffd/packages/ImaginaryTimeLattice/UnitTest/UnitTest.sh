#!/usr/bin/env bash

g++ -std=c++17 -O0  -DFFD_UNIT_TEST_FLAG -I ../../../../ UnitTest.cpp -o UnitTest.x
./UnitTest.x
