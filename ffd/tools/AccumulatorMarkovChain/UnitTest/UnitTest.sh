#!/usr/bin/env bash

g++ -std=c++17 -O0 -Wall -Wfatal-errors -DFFD_UNIT_TEST_FLAG -I ../../../../ UnitTest.cpp -o .UnitTest.x

echo slow version compiled.

g++ -std=c++17 -Ofast -Wall -Wfatal-errors -DFFD_UNIT_TEST_FLAG -I ../../../../ UnitTest.cpp -o .FastUnitTest.x

echo faster version compiled. Add argument for fast execution.

if [ $# -eq 1 ]; then
	 echo faster execution. still not fastest because assertions are there
	 echo add -DNDEBUG to command to remove assertions from code
	 echo UnitTest evidently wont work with -DNDEBUG
	 time ./.FastUnitTest.x
else
	 echo slow execution...
	 time ./.UnitTest.x
fi