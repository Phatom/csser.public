# CSSER
C++ Stochastic Series Expansion Runner
======================================
This is a free software; you may use/distribute the code freely but all work must acknowledge Ushnish Ray as an author.
There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

This library uses the SSE algorithm to compute observables of sign-free model Hamiltonians. Capable of handling very large 
systems, for instance see: Physical Review A 91 (2), 023625, where a system of 200,000 Bosonic particles were simulated to map an
experimental system. The API is flexible to handle a large class of model Hamiltonians, including spin systems, 
disordered systems etc.

Minimum Requirements
--------------------
g++ 4.8 with c++0x support

Gnu Scientfic Library

OpenMP

OpenMPI 1.6.5

Main Library Installation Instructions
--------------------------------------
Create the Makefile using
python paramGenerator/makefileCreate.py MakeFileDirList Makefile

This will create the makefile. You will need to add the linker locations for gsl and gslcblas in the Makefile. The default Makefile assumes that the libraries are loaded under standard LD_LIBRARY_PATH locations. Once the paths are setup correctly run make to create the library. 

Model Installation
------------------
The different models that use the main library are located under Models. Each model requires a separate installation and creation of executables. These are much simpler so the Makefile has already been provided. Run make to create the executables. The parameter file etc. will be model specific so refer to the documentation in the model folders for further instructions. 
