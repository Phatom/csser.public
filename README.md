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
OpenMP 1.
OpenMPI 1.6.5
