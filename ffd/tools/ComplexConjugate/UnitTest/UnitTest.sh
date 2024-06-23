#!/usr/bin/env bash

g++ -std=c++17 -O0 -Wall -Wfatal-errors -DFFD_UNIT_TEST_FLAG -I ../../../../ UnitTest.cpp -o .UnitTest.x

echo slow version compiled.

g++ -std=c++17 -O0 -Wall -Wfatal-errors -DFFD_UNIT_TEST_FLAG -I ../../../../ UnitTest.cpp -o .FastUnitTest.x

echo fast version compiled. Add argument for fast execution.

if [ $# -eq 1 ]; then
	 echo fast execution...
	 ./.FastUnitTest.x
else
	 echo slow execution...
	 ./.UnitTest.x
fi
