#!/usr/bin/env python

import os, sys
import subprocess

name_space = "itime_proposer"
prj_name = "include"

if len(sys.argv) != 2:
	 print("Usage: enter file_name")
else:
	 file_name = sys.argv[1]

	 file = open(file_name+".hpp", "w")
	 file.write("namespace ffd::"+name_space+"{\n\n\n\n")
	 file.write("}//namespace\n")

	 include_prj = open(prj_name+".hpp", "a")
	 include_prj.write("#include\""+file_name+".hpp\"\n")
