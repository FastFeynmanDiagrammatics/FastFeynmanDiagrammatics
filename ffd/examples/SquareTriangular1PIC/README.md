# FastFeynmanDiagrammatics: self-energy code

*Authors: Riccardo Rossi and Fedor Šimkovic IV*

## What does the code compute?

The codes outputs the self-energy on a certain momentum-frequency grid.

### The model

We consider, for the purpose of this example, the Hubbard model
$$
\hat{H}=\hat{H}_0+ U \sum_i n_{i\uparrow}n_{i\downarrow}
$$
where 
$$
\hat{H}_0=[-t\sum_{<i,j>}c_i^\dagger c_j  -t'\sum_{<<i,j>>}c_i^\dagger c_j + h.c.] -\mu \sum_{i,\sigma } n_{i\sigma}
$$


### A rough sketch of the output

The code outputs stochastic estimates of the self-energy double-expansion coefficients
$$
\Sigma_{j,l}(\mathbf{k},i\omega_n)
$$
where $$\Sigma_{j,l}(\mathbf{k},i\omega_n)$$ is the $$j$$-th insertion of $$U$$ and $$l$$-th insertion of a chemical potential shift (which this code assumes of order $$U$$). This data can be used to compute quantities at fixed density by inserting the right chemical potential shift in the post-processing. 

## Important files
```
main.cpp
```
contains the c++ code, and
```
parameters.hpp
```
contains the parameters.

The postprocessing code is at 
```
ffd/postprocessing/mergeruns_fix_density.cpp
```


## Compilation

### System requirements

- A C++ compiler with a standard of at least C++20

### Running the code

```
$mkdir build; cd build
```
```
$cmake ..
```
```
$cmake --build .
```
```
$./main.x
```
The output is saved in the `out` folder.


## Parameters of the codes

### Location

The parameters are in the file `parameters.hpp`

### Physical parameters

For both codes
$$
\texttt{Beta} = \beta= \text{"Inverse temperature"}
$$ 
$$
\texttt{Lx} = \text{"number of lattice sites in the x direction"}
$$
$$
\texttt{Ly} = \text{"number of lattice sites in the y direction"}
$$
$$
\texttt{geometry} = \text{"Lattice geometry"}
$$
`geometry` can be equal to `square` (which includes rectangular geometry) or `triangular`
$$
\texttt{n0} = \text{"average number of particles per site"}
$$
$$
\texttt{order} = \text{"maximal expansion order of the calculation"}
$$
$$
\texttt{omega\_begin},\texttt{omega\_end} 
$$
the code computes Matsubara frequencies $n\in${`omega_begin`,`omega_begin+1`,...,`omega_end-1`}
$$
\texttt{t\_hopping} = \text{"t parameter of the Hubbard model"}
$$
$$
\texttt{tprime\_hopping} = \text{"t' parameter of the Hubbard model"}
$$

### Important numerical parameters

$$
\texttt{k\_slice} = \text{type of momentum accumulation}
$$
if `k_slice` is equal to `BZ`, a uniform grid in the Brillouin zone is used. If `k_slice` is equal to other choices, a path in momentum space is chosen.
$$
\texttt{n\_k\_points\_max}=\text{"Maximal number of k points to save"}
$$
if `n_k_points` is big enough, all the points in the BZ are chosen (if big, the number of BZ points to be chosen will be roughly a eight of the BZ for square-lattice symmetry).
$$
\texttt{verbose\_file} = \text{how much metadata I write to files}
$$
$$
\texttt{verbose\_stdout} = \text{how much metadata I write to screen}
$$
`verbose_file` and `verbose_stdout` can be put to `0` to not output anything not-strictly-useful.
$$
\texttt{print\_interval} = \text{Time in seconds between a print to file}
$$
The code writes stochastically to file, with an (average) period `print_interval` that might be set high (like 3600 for one hour) for clusters that are slow to write. 
$$
\texttt{pid\_counter} = \text{"Process ID counter"}
$$
A process is defined as a set of identical parallel runs (not only on the same CPU). `pid_counter` is used to merge different runs belonging to exactly the same set of numerical (including cpu time) and physical data. At each run, `pid_counter` must be increased by one in order to avoid merging data with different running times.

## Running it on a cluster with slurm

### Checking the output

The output is saved in the out directory. If the main file is run in the folder `ffd/examples/SquareTriangular1PIC`, then the output is

```
ffd/examples/SquareTriangular1PIC/out
```
Inside, directories corresponding to different physical parameters are created. Inside these directories, there are output data and metadata about the calculation. 

We can look at an example: a 12x12 square-lattice Beta 5 estimation of the real part (it is the `r` at the end) of the self energy at  $k_x=1.5708$, $k_y=0.523599$ and at zero Matsubara frequency
```ls ffd/examples/SquareTriangular1PIC/out/square_12x12_B5_n0-0.9456_t1_tp-0.3/kx1.5708_ky0.523599_om0_r```

```
34_OHgte51DNp.ds  34_iivkddQqmK.ds  39_3yjK281XTX.ds	39_KmphnoXnBH.ds   parameters
```
34 and 39 are PID: it means that we have done two runs, each run of two threads (the random string after 34 and 39 is the thread ID of each run).
We can take a look at what they look like
```
$cat 34_OHgte51DNp.ds
-1 -1 0.01708555356 0.0001328640543
0 0 0 0
1 0 0 0
2 0 -5.8872743e-05 6.911e-07
2 1 0 0
3 0 9.161273e-06 9.837e-08
3 1 -0.00029400867 4.069e-06
3 2 -6.7070159e-20 2.021e-19
4 0 -1.2512831e-05 1.348e-07
4 1 5.9150667e-05 1.522e-06
4 2 -0.00040060298 1.61e-05
5 0 1.1599469e-06 2.838e-08
5 1 -2.6178687e-05 1.026e-06
6 0 -3.1160658e-07 2.495e-08
```
so the format is order-in-U   order-in-chemical-potential-shift   value    error. -1 means the normalization sector.

The `.ds` files are raw files: they need to be merged together , **at postprocessing**, to create a result for each PID (the two 34 and the two 39 will be merged together), and then the result from each process will be merged together to create a single file with values and errorbars for each coefficients.

### Changing parameters

1. Always increase the PID counter `pid_counter` in `parameters.hpp`
2. Change the other parameters as needed
3. Recompile the code

### Minimal SLURM script

```
#!/usr/bin/env bash
local -i nnodes=1
local -i ncores=10
#SBATCH -N$nnodes --exclusive --ntasks-per-node=$ncores
#SBATCH --time=run_time
srun main.x
```
where `nnodes` are the number of nodes and `ncores` is the number of threads per node (can be higher than the number of cores by a factor of 2 in some cases). Each thread will write files independently: **no MPI needed**.


## Gathering data and analysis

### Data gathering

A script is provided, `ffd/postprocessing/mergeruns_fix_density.cpp`, to perform the postprocessing data gathering in the `out` folder.
After compiling `mergeruns_fix_density.cpp`, calling the executable `mergeruns_fix_density.x`, it can prepare recursively any folder by 
```
./mergeruns_fix_density.x -v --fix-density target_folder
```
or it can also be launched recursively on everything starting from current folder (might be slow though)
```
./mergeruns_fix_density.x -v --fix-density *
```

### Resummation

The series coefficients can also be futher resummed (code not provided here).
