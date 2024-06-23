#!/usr/bin/env bash

cpp_comp="g++-8"
cpp_std="-std=c++17"
warnings_flags="-Wpedantic -Wextra -Wall -Werror"
unit_test_flag="-DFFD_FFT_UNIT_TEST_FLAG"
max_errors="3"
debug_flag="-g"
slow_flag="-O0"
fast_flag="-O0 -g"
ffd_path="../../../../"
cpp_file="UnitTest.cpp"
slow_file=".UnitTest.x"
fast_file=".FastUnitTest.x"

# $cpp_comp $cpp_std $error_limit_cmd=$max_errors $slow_flag $warning_flags $debug_flag $unit_test_flag -I $ffd_path $cpp_file -o $slow_file

# echo slow version compiled.

$cpp_comp $cpp_std $fast_flag $warning_flags $unit_test_flag -I $ffd_path $cpp_file -o $fast_file

echo faster version compiled.

time ./$fast_file