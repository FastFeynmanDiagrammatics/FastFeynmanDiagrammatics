#!/usr/bin/env bash

cpp_comp="g++"
cpp_std="-std=c++17"
warnings_flags="-Wpedantic -Wextra -Wall -Werror"
unit_test_flag="-DFFD_YRZ_PROPAGATOR_UNIT_TEST_FLAG"
error_limit_cmd="-ferror-limit"
max_errors="1"
debug_flag="-g"
slow_flag="-O0"
fast_flag="-Ofast -march=native -pipe"
ffd_path="../../../../"
cpp_file="main.cpp"
fast_file=".main.x"


$cpp_comp -DNDEBUG $cpp_std $fast_flag $warning_flags $unit_test_flag -I $ffd_path $cpp_file -o $fast_file

time ./$fast_file
