#!/usr/bin/env bash

order0=10 
orderfast=100

g++ -std=c++17 -O0 -Wall -DFFD_UNIT_TEST_FLAG -I /Users/rrossi/prj/19/ main.cpp -o integ_test_O0.x
./integ_test_O0.x $order0
g++ -std=c++17 -Ofast -DNDEBUG -Wall -I /Users/rrossi/prj/19/ main.cpp -o integ_test_fast.x
./integ_test_fast.x $orderfast
