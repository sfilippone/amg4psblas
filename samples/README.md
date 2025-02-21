# Samples

This folder contains several example for the AMG4PSBLAS library. After having compiled (and if needed installed) the library
these example can be compiled and run to become familiar with the use of the library and the preconditioners implemented in it.

## Simple

### pdegen

This folder contains three main examples
- `amg_dexample_ml.f90`
- ``
- ``

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

\[
- \left( a_1 \frac{d^2 u}{dx^2} + a_2 \frac{d^2 u}{dy^2} + a_3 \frac{d^2 u}{dz^2} \right) 
+ \left( b_1 \frac{du}{dx} + b_2 \frac{du}{dy} + b_3 \frac{du}{dz} \right) 
+ c u = f
\]

with **Dirichlet boundary conditions**:

\[
u = g
\]

on the unit cube:  
\[
0 \leq x,y,z \leq 1
\]

##### Special Case: Laplace Equation  
If \( b_1 = b_2 = b_3 = c = 0 \), the PDE reduces to the **Laplace equation**.

##### Computational Domain and Data Distribution  
In this sample program:
1. The index space of the **discretized computational domain** is numbered sequentially in a standard way.
2. The corresponding vector is then distributed according to a **BLOCK data distribution**.


