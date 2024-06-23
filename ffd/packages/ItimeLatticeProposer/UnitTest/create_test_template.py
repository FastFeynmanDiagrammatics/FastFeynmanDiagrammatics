#!/usr/bin/env python3

import os, sys
import subprocess

name_space = "itime_lattice_proposer::unit_test"
prj_name = "include"

if len(sys.argv) != 2:
	 print("Usage: enter file_name")
else:
	 file_name = sys.argv[1]

	 file = open(file_name+".hpp", "x")
	 file.write("namespace ffd::"+name_space+"{\n\n\n\n")
	 file.write("}//namespace\n")

	 include_prj = open(prj_name+".hpp", "a")
	 include_prj.write("#include\""+file_name+".hpp\"\n")