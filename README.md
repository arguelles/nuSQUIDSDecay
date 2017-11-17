# nuSQUIDSDecay

A specialization of the nuSQUIDS class that includes neutrino decay.
For physics details, see Appendix A of arXiv:1711.05921

The class is defined completely in its header file, include/nusquids_decay.h
The examples reside in the examples/ directory.
//----------------------------------Dependencies------------------------------------//

-The SQUIDS library (https://github.com/jsalvado/SQuIDS)
-The nuSQUIDS library (https://github.com/arguelles/nuSQuIDS)
-The HDF5 library (https://support.hdfgroup.org/HDF5/)
	(note: some packaged versions of HDF5 do not include support for the
	c++ API. The API is detailed at https://support.hdfgroup.org/HDF5/doc/cpplus_RM/
	and compiling the library from source w/ proper c++ options should 
	resolve any HDF5/c++ issues should they arise).
-The Doxygen documentation generator (only for compiling docs),
	see http://www.stack.nl/~dimitri/doxygen/ for details.

-Note: Library paths are managed with the pkg-config utility. See the Makefile for 
	its usage. Depending on the OS and the method used to install the libraries
	listed above, you may need to write pkg-config files which point to these 
	libraries. See https://people.freedesktop.org/~dbn/pkg-config-guide.html
	for details.

//-------------------------------Compilation----------------------------------//

Simply running "make" should compile both the coupling and partial rate examples.
These examples are explained in their respective source files, and should, together
with the class documentation, provide a workable understanding of how to use 
the nuSQUIDSDecay class. To compile the documentation, run "doxygen Doxyfile".
Use your favorite browser to open doc/html/index.html. The class documentation
is then under Namespaces/nusquids/nuSQUIDSDecay.

//--------------------------------Execution-----------------------------------//

To run the examples, change directories to examples/ and execute the examples
there. This is done to satisfy the relative paths pointing to the fluxes/
and output/ directories. Both examples will read input fluxes from the fluxes/
directory, and write both initial and final fluxes from kaon and pion channels
to text files in output/. 
The format of each line of the output file is:
cos(zenith angle)   neutrino energy(eV)   nu_mu flux   nu_mu_bar flux

More flux flavors can be output simply by modifying the WriteFlux() function
in the example source files appropriately.

//----------------------------------------------------------------------------//

For more information please email:

Alexander (Zander) Moss (zander@caltech.edu)
Marjon Moulai (marjon@mit.edu)
Carlos Arguelles (caad@mit.edu)
Janet Conrad (conrad@mit.edu)
