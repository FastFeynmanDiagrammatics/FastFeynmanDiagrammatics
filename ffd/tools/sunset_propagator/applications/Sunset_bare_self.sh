#!/usr/bin/env bash

cpp_comp="g++"
cpp_std="-std=c++17"
warnings_flags="-Wpedantic -Wextra -Wall -Werror"
unit_test_flag="-DFFD_SUNSET_PROPAGATOR_UNIT_TEST_FLAG"
max_errors="1"
debug_flag="-g"
slow_flag="-O0"
fast_flag="-Ofast -march=native"
ffd_path="../../../../"
cpp_file="sunset_bare_self.cpp"
fast_file=".sunset_self.x"
$cpp_comp $cpp_std $fast_flag $warning_flags  $unit_test_flag -I $ffd_path $cpp_file -o $fast_file

echo faster version compiled

time ./$fast_file
