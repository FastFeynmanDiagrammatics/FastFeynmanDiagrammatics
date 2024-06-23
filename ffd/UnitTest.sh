#!/usr/bin/env bash

cpp_comp="g++-14"
cpp_std="-std=c++20"
warnings_flags="-Wpedantic -Wextra -Wall -Werror"
unit_test_flag="-DFFD_UNIT_TEST_FLAG"
debug_flag="-g"
slow_flag="-O0"
fast_flag="-O3"
ffd_path="../"
cpp_file="UnitTest.cpp"
slow_file=".UnitTest.x"
fast_file=".FastUnitTest.x"

echo "using "$cpp_comp" compiler"
#time $cpp_comp $cpp_std $debug_flag $warning_flags $unit_test_flag $slow_flag -I $ffd_path $cpp_file -o $slow_file
echo slow version compiled.
time $cpp_comp $cpp_std $warning_flags $unit_test_flag $fast_flag -I $ffd_path $cpp_file -o $fast_file
echo fast version compiled and running. It is OK if no output
./$fast_file
