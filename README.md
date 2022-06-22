
                         AMG4PSBLAS  
 Algebraic Multigrid Package  based on PSBLAS (Parallel Sparse BLAS version 3.8)
    
Salvatore Filippone    (University of Rome Tor Vergata and IAC-CNR)
Pasqua D'Ambra         (IAC-CNR, Naples, IT)
Fabio Durastante       (IAC-CNR, Naples, IT)

---------------------------------------------------------------------

AMG4PSBLAS is a package of Algebraic MultiGrid (AMG)
preconditioners for  the iterative solution of large and sparse linear systems.

It is an evolution of MLD2P4 (see LICENSE.MLD2P4), but it has been
thoroughly reworked, and it is sufficiently different to warrant a new
project name.


MAIN REFERENCES:

     

P. D'Ambra, D. di Serafino, S. Filippone,
MLD2P4: a Package of Parallel Algebraic Multilevel Domain Decomposition
Preconditioners in Fortran 95,
ACM Transactions on Mathematical Software, 37 (3), 2010, art. 30,
doi: 10.1145/1824801.1824808.


TO COMPILE

0. Unpack the tar file in a directory of your choice (preferrably
   outside the main PSBLAS directory).
1. run configure --with-psblas=<ABSOLUTE path of the PSBLAS install directory>
   adding the options for MUMPS, SuperLU, SuperLU_Dist, UMFPACK as desired.
   See MLD2P4 User's and Reference Guide (Section 3) for details.
2. Tweak Make.inc if you are not satisfied.
3. make; 
4. Go into the test subdirectory and build the examples of your choice.
5. (if desired): make install 


NOTES

- The single precision version is supported only by MUMPS and SuperLU;
  thus, even if you specify at configure time to use UMFPACK or SuperLU_Dist, 
  the corresponding preconditioner options will be available only from
  the double precision version.


The AMG4PSBLAS team. 
---------------
Salvatore Filippone
Pasqua     D'Ambra
Fabio     Durastante
