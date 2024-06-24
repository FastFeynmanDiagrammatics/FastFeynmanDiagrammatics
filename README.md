Authors: Riccardo Rossi and Fedor Šimkovic IV

# FastFeynmanDiagrammatics: self-energy and susceptibility codes
This code has been used to obtain the numerical results of the paper [https://arxiv.org/abs/2209.09237](https://arxiv.org/abs/2209.09237).

## Location of the codes and important files

The code for the self-energy is at
```
ffd/examples/SquareTriangular1PIC
```
while the code for the susceptibility is at
```
ffd/examples/SquareTriangularSSC
```
In these two folders, the important files are 
```
main.cpp
``` 
which contains the C++ code,
```
parameters.hpp
```
which contains the parameters, and
```
README.md
```
which contains the instructions on how to run and interpret the results.

The postprocessing code to merge results of different processes is at
```
ffd/postprocessing/mergeruns_fix_density.cpp
```

The self-energy and the susceptibility codes have similar usage and conventions.

## System requirements
- A C++ compiler with a standard of at least C++20