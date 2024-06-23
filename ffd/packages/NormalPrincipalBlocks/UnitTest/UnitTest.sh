#!/usr/bin/env bash

cpp_comp="g++-9"
cpp_std="-std=c++17"
warnings_flags="-Wpedantic -Wextra -Wall -Werror"
unit_test_flag="-DFFD_NORMALPRINCIPALBLOCKS_UNIT_TEST_FLAG"
error_limit_cmd="-ferror-limit"
max_errors="1"
debug_flag="-g"
slow_flag="-O0"
fast_flag="-Ofast"
ffd_path=$FFD_LIB_PATH
cpp_file="UnitTest.cpp"
slow_file=".UnitTest.x"
fast_file=".FastUnitTest.x"

$cpp_comp $cpp_std $slow_flag $warning_flags $debug_flag  $unit_test_flag -I $ffd_path $cpp_file -o $slow_file

echo slow version compiled.

$cpp_comp $cpp_std $fast_flag $warning_flags  $unit_test_flag -I $ffd_path $cpp_file -o $fast_file

echo faster version compiled. Add argument for fast execution.

if [ $# -eq 1 ]; then
	 echo faster execution. still not fastest because assertions are there
	 echo add -DNDEBUG to command to remove assertions from code
	 echo UnitTest evidently wont work with -DNDEBUG
	 time ./$fast_file
else
	 echo slow execution...
	 time ./$slow_file
fi
