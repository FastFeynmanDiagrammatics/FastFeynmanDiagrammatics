#!/usr/bin/env bash

g++ -std=c++17 -O0 -Wall -Wfatal-errors -Wno-unused-parameter -DFFD_UNIT_TEST_FLAG -I ../../../../ UnitTest.cpp -o UnitTest.x
./UnitTest.x
