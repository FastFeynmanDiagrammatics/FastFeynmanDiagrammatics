#!/usr/bin/env bash

g++ -std=c++17 -O0 -Wfatal-errors -Wall -DFFD_UNIT_TEST_FLAG -I ~/prj/19/ UnitTest.cpp -o UnitTest.x
./UnitTest.x
