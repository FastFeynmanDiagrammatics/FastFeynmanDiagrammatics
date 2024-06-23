#!/usr/bin/env bash

g++ -std=c++17 -O0 -Wall -Wfatal-errors -DFFD_UNIT_TEST_FLAG -I ../../../../ UnitTest.cpp -o .UnitTest.x
g++ -std=c++17 -Ofast -Wall -Wfatal-errors -DFFD_UNIT_TEST_FLAG -I ../../../../ UnitTest.cpp -o .FastUnitTest.x

./.UnitTest.x
