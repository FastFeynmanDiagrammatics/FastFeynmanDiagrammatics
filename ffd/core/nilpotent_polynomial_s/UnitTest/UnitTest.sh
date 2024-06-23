#!/usr/bin/env bash

g++-10 -std=c++17 -O0 -Wall -Wextra -Wpedantic -DFFD_NILPOTENTPOLYNOMIAL_S_UNIT_TEST_FLAG -I ../../../../ UnitTest.cpp -o .UnitTest.x
./.UnitTest.x
