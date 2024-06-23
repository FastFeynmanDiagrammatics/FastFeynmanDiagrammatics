#!/usr/bin/env bash

cpp_comp="g++"
cpp_std="-std=c++17"
warnings_flags="-Wpedantic -Wextra -Wall -Werror"
unit_test_flag="-DFFD_SIGMADET_HUBBARD_UNIT_TEST_FLAG"
error_limit_cmd="-ferror-limit"
max_errors="1"
debug_flag="-g"
slow_flag="-O0"
fast_flag="-Ofast"
ffd_path="../../../../../"
cpp_file="compute_propagator_atom2.cpp"
fast_file=".propagator_atom2.x"

$cpp_comp $cpp_std $fast_flag $warning_flags $error_limit_cmd=$max_errors $unit_test_flag -I $ffd_path $cpp_file -o $fast_file

echo faster version compiled. 

time ./$fast_file
