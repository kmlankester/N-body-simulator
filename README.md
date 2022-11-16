Implementation of an N-body simulation for a Parallel Scientific Computing module. Details of the assignment can be seen below.



# Coursework: PSCI
## Overview

In this assignment, you will write a very simple N-body solver which runs efficiently on a single node of a supercomputer.
It should simulate N objects (bodies) in space -- think of suns or molecules -- which all move and fuse when they hit each other.

The assignment is broken down into steps of increasing difficulty. Each step’s result is to be handed in separately, so they can be marked independently. The assignment steps mirror content from the lecture. If you get stuck, perhaps revise the lecture slides and notes. The challenges within the assignment cover both topics from Parallel Scientific Computing: numerical algorithms and parallelisation.

- Please add your CIS ID and name to this file. If we cannot associate your submission with your name we cannot mark it.
- This assignment is to be completed and handed in via GitHub Classroom. That means your changes must be committed and pushed to your repository before the submission date.
- For deadlines, please consult the relevant level rubric on Blackboard.

## Tasks

You can find the instructions and the marking scheme for the assignment in the file [Assignment2021.pdf](Assignment2021.pdf) which is included in the repository.

## Compiling and executing the code

You will be modifying the file `assignment-code.cpp` in this repository. You can compile the C++ source code with the command
```bash
make
```
which generates two programmes, `assignment-gcc` (using the GNU Compiler Collection) and `assignment-icpc` (using the Intel Compiler). If you are working in an environment where only one of the two compilers is available, you can use `make assignment-gcc` or `make assignment-icpc` to only produce the binary file of interest, but please make sure that your final submission generates both files without errors, as both will be used for marking. Note that without command-line parameters, the two programs only print a usage message.

Running the code produces a number of output files for visualization with ParaView (extensions `.vtp` and `.pvd`). These can be removed with
```bash
make clean_paraview
```
or with
```bash
make cleanall
```
which will also remove the files `assignment-gcc` and `assignment-icpc`.
