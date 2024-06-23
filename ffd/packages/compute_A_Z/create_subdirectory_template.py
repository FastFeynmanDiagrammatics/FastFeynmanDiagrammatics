#!/usr/bin/env python3

import os, sys
import subprocess

if len(sys.argv) > 1:
	 subdir_name = sys.argv[1]
if len(sys.argv) > 2:
	 print("second argument has been discarded!!")

include_prj = open("include.hpp", "a")
include_prj.write("#include\""+subdir_name+".hpp\"\n")
subdir_prj = open(subdir_name+".hpp", "x")
subdir_prj.write("#include\""+subdir_name+"/include.hpp\"\n")
os.mkdir(subdir_name)

include_sub_prj = open(subdir_name+"/include.hpp", "x")
subprocess.run(["cp", "create_header_template.py", subdir_name])

