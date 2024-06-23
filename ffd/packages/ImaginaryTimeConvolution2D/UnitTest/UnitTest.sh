#!/usr/bin/env bash

unit_test_flag="-DFFD_IMAGINARYTIMECONVOLUTION2D_UNIT_TEST_FLAG"

g++ -std=c++17 -O0 -Wall -Wfatal-errors $unit_test_flag -I ../../../../ UnitTest.cpp -o .UnitTest.x
g++ -std=c++17 -Ofast -Wall -Wfatal-errors $unit_test_flag -I ../../../../ UnitTest.cpp -o .FastUnitTest.x

./.UnitTest.x
