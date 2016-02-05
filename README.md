# CSSER
C++ Stochastic Series Expansion Runner
======================================
Author: Ushnish Ray
email: ushnish.qmc@gmail.com

This is a free software; you may use/distribute the code freely but all work must acknowledge Ushnish Ray as an author.
There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

This library uses the SSE algorithm [1]-[3] to compute observables of sign-free model Hamiltonians. Capable of handling very large systems, for instance see [4] where a system of 200,000 Bosonic particles were simulated to map an
experimental system. The API is flexible to handle a large class of model Hamiltonians, including spin systems, 
disordered systems etc.

* [1] A. W. Sandvik, Phys. Rev. B 59, R14157 (1999). 
* [2] A. Dorneich and M. Troyer, Phys. Rev. E 64, 066701 (2001). 
* [3] O. F. Syljuasen and A. W. Sandvik, Phys. Rev. E 66, 046701 (2002).
* [4] McKay, D., Ray, U., Natu, S., Russ, P., Ceperley, D., & DeMarco, B. Phys. Rev. A, 91(2), 023625 (2015).

Minimum Requirements
--------------------
*g++ 4.8 with c++0x support
*Gnu Scientfic Library 
*OpenMP
*OpenMPI 1.6.5

Main Library Installation Instructions
--------------------------------------
Create a folder objLibrary using: mkdir objLibrary
Create the Makefile using
python paramGenerator/makefileCreate.py MakeFileDirList Makefile

This will create the makefile and the target location for the object files that are used to create the static linker libarary libsse.a used by all models. You will need to add the linker locations for gsl and gslcblas in the Makefile. The default Makefile assumes that the libraries are loaded under standard LD_LIBRARY_PATH locations. Once the paths are setup correctly run make. 

Model Installation
------------------
The different models that use the main library are located under Models. Each model requires a separate installation and creation of executables. These are much simpler so the Makefile has already been provided. Run make to create the executables. The parameter file etc. will be model specific so refer to the documentation in the model folders for further instructions. 

For instructions on the how the different elemets of the code fits and the API email me directly. I have the UML diagrams as well as sequence diagrams that will enable users to understand the mechanisms of the code. 
