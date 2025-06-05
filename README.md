# AMG4PSBLAS  v1.2
Algebraic Multigrid Package  based on [PSBLAS](https://github.com/sfilippone/psblas3) (Parallel Sparse BLAS version 3.9)

AMG4PSBLAS is a package of parallel algebraic multilevel preconditioners included in the PSCToolkit (Parallel Sparse Computation Toolkit) software framework.

It is a progress of a software development project started in 2007, named MLD2P4, which originally implemented a multilevel version of some domain decomposition preconditioners of additive-Schwarz type and was based on a parallel decoupled version of the well known smoothed aggregation method to generate the multilevel hierarchy of coarser matrices.

In the last years the package was extended for including new algorithms and functionalities for the setup and application new AMG preconditioners with the final aims of improving efficiency and scalability when tens of thousands cores are used and of boosting reliability in dealing with general symmetric positive definite linear systems.

It is an evolution of MLD2P4 (see [LICENSE.MLD2P4](LICENSE.MLD2P4)), but due to the significant number of changes and the increase in scope, we decided to rename the package as AMG4PSBLAS.

AMG4PSBLAS has been designed to provide scalable and easy-to-use preconditioners in the context of the PSBLAS (Parallel Sparse Basic Linear Algebra Subprograms) computational framework and can be used in conjuction with the Krylov solvers available in this framework. Our package is based on a completely algebraic approach; therefore users level interfaces assume that the system matrix and preconditioners are represented as PSBLAS distributed sparse matrices.

AMG4PSBLAS enables the user to easily specify different features of an algebraic multilevel preconditioner, thus allowing to experiment with different preconditioners for the problem and parallel computers at hand.

The package employs object-oriented design techniques in Fortran 2008, with interfaces to additional third party libraries such as MUMPS, UMFPACK, SuperLU, and SuperLU_Dist, which can be exploited in building multilevel preconditioners. The parallel implementation is based on a Single Program Multiple Data (SPMD) paradigm; the inter-process communication is based on MPI and is managed mainly through PSBLAS.

## Main Refrerences:

The main reference for this project is
> D'Ambra, P., Durastante, F., & Filippone, S. (2021). AMG preconditioners for linear solvers towards extreme scale. SIAM Journal on Scientific Computing, 43(5), S679-S703.

AMG4PSBLAS is the suite of preconditioners for the Parallel Sparse Computation Toolkit ([PSCToolkit](https://psctoolkit.github.io/)) suite of libraries. See the paper:
> D’Ambra, P., Durastante, F., & Filippone, S. (2023). Parallel Sparse Computation Toolkit. Software Impacts, 15, 100463.

The main reference for features inherited from MLD2P4 is
> P. D'Ambra, D. di Serafino, S. Filippone,
> MLD2P4: a Package of Parallel Algebraic Multilevel Domain Decomposition
> Preconditioners in Fortran 95,
> ACM Transactions on Mathematical Software, 37 (3), 2010, art. 30,
> doi: 10.1145/1824801.1824808.

## Installing

Installation requires having a working version of the [PSBLAS](https://github.com/sfilippone/psblas3) library installed.
AMG4PSBLAS has several interfaces to third-party libraries that can be used in the construction and application phases of preconditioners. 
In particular, it is possible to link AMG4PSBLAS with the libraries: MUMPS, SuperLU, SuperLU_Dist, UMFPACK. This is _not mandatory_ and the library can run 
in isolation and without these features.

0. Unpack the tar file in a directory of your choice (preferrably
   outside the main PSBLAS directory).
1. run configure `--with-psblas=<ABSOLUTE path of the PSBLAS install directory>`
   adding the options for MUMPS, SuperLU, SuperLU_Dist, UMFPACK as desired.
   See [AMG4PSBLAS User's and Reference Guide](docs/amg4psblas_1.0-guide.pdf) (Section 3) for details.
2. Tweak `Make.inc` if you are not satisfied.
3. run `make`; 
4. Go into the test subdirectory and build the examples of your choice.
5. (if desired): `make install` 

>[!CAUTION]
>The single precision version is supported only by MUMPS and SuperLU;
>thus, even if you specify at configure time to use UMFPACK or SuperLU_Dist,
>the corresponding preconditioner options will be available only from
>the double precision version.

### CUDA, OpeMP, OpenACC

CUDA, OpenMP and OpenACC features are transparently inherited by PSBLAS installation. If PSBLAS has been configured (and installed) with these supports then AMG4PSBLAS will transparently inherit them. It will then be possible to move the computation to GPU accelerator simply by selecting the appropriate variable types. If these have not been activated or installed for PSBLAS then they will not be available for AMG4PSBLAS either and the operation will be purely on CPU/MPI. See also the samples/cuda folder.

### EoCoE - Software as service portal

In the European project “Energy oriented Center of Excellence: toward exascale for energy” we made available a software as service portal: [https://eocoe.psnc.pl/](https://eocoe.psnc.pl/). This permits to test several cutting-edge computational methods for accelerating the transition to the production, storage and management of clean, decarbonized energy. Among them you have the possibility of running PSBLAS+AMG4PSBLAS on some test problems to become familiar with using the software.

## TODO and bugs
 
- [X] Fix all reamining bugs. Bugs? We dont' have any ! 🤓

> [!NOTE]
> To report bugs 🐛 or issues ❓ please use the [GitHub issue system](https://github.com/sfilippone/amg4psblas/issues).

## The AMG4PSBLAS team. 

- Pasqua D'Ambra         (IAC-CNR, Naples, IT)
- Fabio Durastante       (University of Pisa and IAC-CNR, IT)
- Salvatore Filippone    (University of Rome Tor Vergata and IAC-CNR, IT)

**Contributors** (_roughly reverse cronological order_):

- Luca       Pepè Sciarria
- Andea      Di Iorio
- Ambra	     Abdullahi Hassan
- Alfredo    Buttari
