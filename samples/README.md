# Samples

This folder contains several example for the AMG4PSBLAS library. After having compiled (and if needed installed) the library
these example can be compiled and run to become familiar with the use of the library and the preconditioners implemented in it.

## Simple

To compile the examples it is sufficient to enter the folder and run `make`. The executables will be moved to the corresponding `run` subdirectory.

### pdegen

This folder contains two main examples
- `amg_[s/d]example_ml.f90`
- `amg_[s/d]example_1lev.f90`

The difference between the `s` and the `d` variants is the use of the single or double precision.

#### Example `amg_dexample_ml.f90`

This sample program solves a linear system obtained by discretizing a PDE with Dirichlet boundary conditions. The solver used is **Flexible Conjugate Gradient (FCG)**, coupled with one of the following multi-level preconditioners, as explained in **Section 4.1** of the *AMG4PSBLAS User's and Reference Guide*:

##### Available Preconditioner Choices:

- **Choice = 1** (Default multi-level preconditioner):  
  V-cycle with decoupled smoothed aggregation, **1 hybrid forward/backward Gauss-Seidel** sweep as pre/post-smoother, and **UMFPACK** as the coarsest-level solver.  
  *(See Section 4.1, Listing 1)*  

- **Choice = 2**:  
  V-cycle preconditioner with **1 block-Jacobi sweep** (using ILU(0) on the blocks) as pre/post-smoother, and **8 block-Jacobi sweeps** (with ILU(0) on the blocks) as the coarsest-level solver.  
  *(See Section 4.1, Listing 2)*  

- **Choice = 3**:  
  W-cycle preconditioner based on coupled aggregation relying on **matching**, with:
  - Maximum aggregate size of **8**  
  - Smoothed prolongators  
  - **2 hybrid forward/backward Gauss-Seidel sweeps** as pre/post-smoother  
  - A **distributed coarsest matrix**  
  - Preconditioned **Flexible Conjugate Gradient** as the coarsest-level solver  
  *(See Section 4.1, Listing 3)*  

##### Input Data  
The matrix and the right-hand side (RHS) are read from files. If an RHS is not available, a **unit RHS** is set.

#### The PDE Formulation  
The PDE is a general second-order equation in **3D**:
$$- \left( a_1 \frac{d^2 u}{dx^2} + a_2 \frac{d^2 u}{dy^2} + a_3 \frac{d^2 u}{dz^2} \right) + \left( b_1 \frac{du}{dx} + b_2 \frac{du}{dy} + b_3 \frac{du}{dz} \right) + c u = f$$

with **Dirichlet boundary conditions**: $$u = g$$ on the unit cube: $$0 \leq x,y,z \leq 1$$

##### Special Case: Laplace Equation  
If $b_1 = b_2 = b_3 = c = 0$, the PDE reduces to the **Laplace equation**.

##### Computational Domain and Data Distribution  
In this sample program:
1. The index space of the **discretized computational domain** is numbered sequentially in a standard way.
2. The corresponding vector is then distributed according to a **BLOCK data distribution**.

#### Example `amg_dexample_1lev.f90`

This sample program solves a linear system obtained by discretizing a PDE with Dirichlet boundary conditions. The solver used is **BiCGStab**, preconditioned by **Restricted Additive Schwarz (RAS)** with overlap **2** and **ILU(0)** on the local blocks, as explained in **Section 4.1** of the *AMG4PSBLAS User's and Reference Guide*.

##### The PDE Formulation  
The PDE is a general second-order equation in **3D**:

$$- \left( a_1 \frac{d^2 u}{dx^2} + a_2 \frac{d^2 u}{dy^2} + a_3 \frac{d^2 u}{dz^2} \right) + \left( b_1 \frac{du}{dx} + b_2 \frac{du}{dy} + b_3 \frac{du}{dz} \right)  + c u = f$$

with **Dirichlet boundary conditions**: $$ u = g $$ on the unit cube: $$0 \leq x,y,z \leq 1$$

##### Special Case: Laplace Equation  
If $b_1 = b_2 = b_3 = c = 0$, the PDE reduces to the **Laplace equation**.

### fileread

This sample program `amg_[s/d/c/z]example_1lev.f90` solves a linear system using **BiCGStab**, preconditioned by **Restricted Additive Schwarz (RAS)** with overlap **2** and **ILU(0)** on the local blocks, as explained in **Section 4.1** of the *AMG4PSBLAS User's and Reference Guide*.

##### Input Data  
The matrix and the right-hand side (RHS) are read from files. If an RHS is not available, a **unit RHS** is set.

## newlsv
This folder contains a simple program to demonstrate how to define a new **solver** object. The actual code is simply a copy of the **ILU(0)** solver, but it demonstrates the integration process, which can be achieved even at the level of the user program without touching the main library. The program solves a simple discretization of the Poisson equation with Dirichlet boundary conditions

## cuda
This folder contains a simple program to demonstrate how to integrate CUDA-enabled data structures in your application, if available. The program will compile and run even if the main PSBLAS library has been compiled without CUDA support; it builds the same problem as in the **newslv** folder.

## advanced

This folder contains more complicated examples where you can choose, by setting them from the input files, most of the options available inside the AMG4PSBLAS library, it is a good starting point to test the different combinations on a finite difference discretization of a simple differential equation or on matrices read from files.

### pdegen

This folder contains four examples:

- `amg_[s/d]_pde[2/3]d.f90` 

##### The 3D Case

This sample program solves a linear system obtained by discretizing a PDE with **Dirichlet boundary conditions**.

##### The PDE Formulation  
The PDE is a general second-order equation in **3D**:  
$$- \left( a_1 \frac{d^2 u}{dx^2} + a_2 \frac{d^2 u}{dy^2} + a_3 \frac{d^2 u}{dz^2} \right) + \left( b_1 \frac{du}{dx} + b_2 \frac{du}{dy} + b_3 \frac{du}{dz} \right) + c u = f$$  

with **Dirichlet boundary conditions**:  
$$u = g$$  

on the unit cube:  
$$0 \leq x,y,z \leq 1$$  

##### Special Case: Laplace Equation  
If $$b_1 = b_2 = b_3 = c = 0$$, the PDE reduces to the **Laplace equation**.

##### Data Distribution Choices  
There are three available choices for data distribution:

1. **Simple BLOCK distribution**  
2. **Arbitrary index assignment** (typically from a graph partitioner)  
3. **3D distribution** where the unit cube is partitioned into subcubes, each assigned to a process.

##### The 2D Case 

This sample program solves a linear system obtained by discretizing a PDE with **Dirichlet boundary conditions**.

##### The PDE Formulation  
The PDE is a general second-order equation in **2D**:  
$$- \left( a_1 \frac{d^2 u}{dx^2} + a_2 \frac{d^2 u}{dy^2} \right) + \left( b_1 \frac{du}{dx} + b_2 \frac{du}{dy} \right) + c u = f$$  

with **Dirichlet boundary conditions**:  
$$u = g$$  

on the unit square:  
$$0 \leq x,y \leq 1$$  

##### Special Case: Laplace Equation  
If $$b_1 = b_2 = c = 0$$, the PDE reduces to the **Laplace equation**.

##### Data Distribution Choices  
There are three available choices for data distribution:

1. **Simple BLOCK distribution**  
2. **Arbitrary index assignment** (typically from a graph partitioner)  
3. **2D distribution** where the unit square is partitioned into rectangles, each assigned to a process.  

### Using pdegen for scalability studies
You can use the programs in pdegen to perform scalability studies, but you need to consider the following aspects:
1. **Sizing your test case for strong scalability**  In a strong scalability study you solve the exact same problem varying the number or processors (and matching processes). The main issue is the size of the problem: with modern hardware you need to have a substantial size, say 1 million equations per node: therefore it is advisable to start with the largest possibile problem size that would fit onto a single node/processor, and proceed from that; the output from the programs will give an indication of the amount of memory occupied. Once the local size N/NP drops below a certain threshold, performance will also start to flat out or drop. The exact value of the threshold depends on the ratio between the speed of the computing cores and the speed of the network: the faster the network, the smaller the threshold.
2. **Sizing your test case for weak scalability** In a weak scalability study you solve a problem whose size scales with the number of processes/processing cores: if you have a problem of size N=1M on 1 process, then you run with N=2M on 2 processes, and so on. The relationship between N and the IDIM parameter in the input file depends on whether you are running a 2D problem (N=IDIM^2) or a 3D problem *N=IDIM^3). At the time of this writing, a size of 1 or a few million equations per process is a reasonable starting point.
3. **What to measure in a scalability study** The test programs print out a detailed description of the time to set up the problem, the time to set up the preconditioner and the time to solve the problem.     The latter (time to solution) is usually the first parameter to be studied: it is further split into number of iterations and time per iteration. Ideally the number of iteration should be constant    across different configurations (number of processes and/or problem size), whereas the time per iteration should
    1. Scale linearly with the number of processes in a strong scaling experiment;
    2. Remain constant in a weak scaling experiment.
       
  Of course these ideal conditions will not be met exactly in practice, but it is important to have a fine level split of the two factors contributing.

### fileread

The Fortran source code in `amg_[s/d/c/z]f_sample.f90` demonstrates how to read a sparse matrix and its right-hand side (RHS) from files, set up an algebraic multigrid (AMG) preconditioner, and solve a linear system using an iterative solver.

1. **Initialization and Setup**  
   The program initializes the MPI environment and sets up the AMG4PSBLAS parameters. It processes input options that configure solver and preconditioner settings.

2. **File Reading for Matrix and RHS**  
   The code reads the matrix and RHS from files. If the RHS file is missing, it defaults to a unit RHS (i.e. a vector with all entries equal to 1). This enables the formulation of the linear system:  
   $$Ax = b$$  
   where $$A$$ is the matrix and $$b$$ is the right-hand side vector.

3. **AMG Preconditioner Construction**  
   After reading the matrix, the program sets up an AMG preconditioner. This preconditioner creates a hierarchy of coarser grids that improves the convergence of the iterative solver when applied to large, sparse systems.

4. **Iterative Solver Execution**  
   With the preconditioner in place, the code employs an iterative Krylov subspace method (such as Conjugate Gradient or BiCGStab) to solve the system. The AMG preconditioner is used within the iterative loop to accelerate convergence.

5. **Output and Finalization**  
   Upon convergence, the program outputs key information such as the number of iterations and the residual norm. Finally, it finalizes the MPI environment and properly terminates the execution.

Overall, the sample program serves as a practical demonstration of using AMG4PSBLAS in a parallel computing environment. It guides the user through initializing the computation, reading and distributing problem data, configuring the AMG preconditioner, executing an iterative solver, and finalizing the computation—all of which are crucial steps for efficiently solving large sparse linear systems.
