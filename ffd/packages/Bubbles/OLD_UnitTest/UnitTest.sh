#!/usr/bin/env bash

g++ -std=c++17 -O0 -Wall -Wfatal-errors -DFFD_UNIT_TEST_FLAG -I ../../../../ UnitTest.cpp -o .UnitTest.x

g++ -std=c++17 -Ofast -Wall -Wfatal-errors -DFFD_UNIT_TEST_FLAG -I ../../../../ UnitTest.cpp -o .FastUnitTest.x

echo slow and fast version compiled. Add argument for fast execution

if [ $# -eq 1 ]
then
    echo fast execution
    ./.FastUnitTest.x
else
    echo slow execution
    ./.UnitTest.x
fi

    
   
