/*
 * 
 *                            MLD2P4  version 2.1
 *   MultiLevel Domain Decomposition Parallel Preconditioners Package
 *              based on PSBLAS (Parallel Sparse BLAS version 3.5)
 *   
 *   (C) Copyright 2008-2018
 *
 *      Salvatore Filippone                              
 *      Pasqua D'Ambra                                  
 *      Daniela di Serafino   
 * 
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *    1. Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *    2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions, and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *    3. The name of the MLD2P4 group or the names of its contributors may
 *       not be used to endorse or promote products derived from this
 *       software without specific written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
 *  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 * 
 *
 * File: amg_slu_interface.c
 *
 * Functions: amg_dslu_fact, amg_dslu_solve, amg_dslu_free.
 *
 * This file is an interface to the SuperLU routines for sparse factorization and
 * solve. It was obtained by modifying the c_fortran_dgssv.c file from the SuperLU
 * source distribution; original copyright terms are reproduced below.
 * 
 */ 


/*		=====================

Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

(1) Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer. 
(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution. 
(3) Neither the name of Lawrence Berkeley National Laboratory, U.S. Dept. of
Energy nor the names of its contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
  
*/

/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */

#ifdef Have_SLU_
#include "slu_ddefs.h"

#define HANDLE_SIZE  8

typedef struct {
    SuperMatrix *L;
    SuperMatrix *U;
    int *perm_c;
    int *perm_r;
} factors_t;


#else

#include <stdio.h>

#endif



int amg_dslu_fact(int n, int nnz, double *values,
		  int *colptr, int *rowind, void **f_factors)
{
/* 
 * This routine can be called from Fortran.
 *  performs LU decomposition.
 *
 * f_factors (input/output)  
 *      On  output contains the pointer pointing to
 *       the structure of the factored matrices.
 *
 */
 
#ifdef Have_SLU_
    SuperMatrix A, AC;
    SuperMatrix *L, *U;
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    int *etree;  /* column elimination tree */
    SCformat *Lstore;
    NCformat *Ustore;
    int      i, panel_size, permc_spec, relax;
    trans_t  trans;
    mem_usage_t   mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    factors_t *LUfactors;
    GlobalLU_t Glu;   /* Not needed on return. */
    int info;

    trans = NOTRANS;

    
    /* Set the default input options. */
    set_default_options(&options);
    
    /* Initialize the statistics variables. */
    StatInit(&stat);
    
    dCreate_CompCol_Matrix(&A, n, n, nnz, values, rowind, colptr,
			   SLU_NC, SLU_D, SLU_GE);
    L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
    U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
    if ( !(perm_r = intMalloc(n)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
    
    /*
	 * Get column permutation vector perm_c[], according to permc_spec:
	 *   permc_spec = 0: natural ordering 
	 *   permc_spec = 1: minimum degree on structure of A'*A
	 *   permc_spec = 2: minimum degree on structure of A'+A
	 *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
	 */    	
    options.ColPerm=2;
    permc_spec = options.ColPerm;
    get_perm_c(permc_spec, &A, perm_c);
    
    sp_preorder(&options, &A, perm_c, etree, &AC);
    
    panel_size = sp_ienv(1);
    relax = sp_ienv(2);
#if defined(SLU_VERSION_5)
    dgstrf(&options, &AC, relax, panel_size, etree,
	   NULL, 0, perm_c, perm_r, L, U, &Glu, &stat, &info);
#elif defined(SLU_VERSION_4)
    dgstrf(&options, &AC, relax, panel_size, etree,
	   NULL, 0, perm_c, perm_r, L, U, &stat, &info);
#else
    choke_on_me;
#endif
    
    if ( info == 0 ) {
      Lstore = (SCformat *) L->Store;
      Ustore = (NCformat *) U->Store;
      dQuerySpace(L, U, &mem_usage);
#if 0
      printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
      printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
      printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz);
      printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
	     mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
#endif
    } else {
      printf("dgstrf() error returns INFO= %d\n", info);
      if ( info <= n ) { /* factorization completes */
	dQuerySpace(L, U, &mem_usage);
	printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
	       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
      }
    }
    
    /* Save the LU factors in the factors handle */
    LUfactors = (factors_t*) SUPERLU_MALLOC(sizeof(factors_t));
    LUfactors->L = L;
    LUfactors->U = U;
    LUfactors->perm_c = perm_c;
    LUfactors->perm_r = perm_r;
    *f_factors = (void *) LUfactors;
    
    /* Free un-wanted storage */
    SUPERLU_FREE(etree);
    Destroy_SuperMatrix_Store(&A);
    Destroy_CompCol_Permuted(&AC);
    StatFree(&stat);
    return(info);
#else
    fprintf(stderr," SLU Not Configured, fix make.inc and recompile\n");
    return(-1);
#endif
}


int amg_dslu_solve(int itrans, int n, int nrhs, double *b, int ldb,
		   void *f_factors)
{
  /* 
   * This routine can be called from Fortran.
   *      performs triangular solve
   *
   */
  int info;
#ifdef Have_SLU_ 
  SuperMatrix B;
  SuperMatrix *L, *U;
  int *perm_r; /* row permutations from partial pivoting */
  int *perm_c; /* column permutation vector */
  int *etree;  /* column elimination tree */
  SCformat *Lstore;
  NCformat *Ustore;
  int      i, panel_size, permc_spec, relax;
  trans_t  trans;
  SuperLUStat_t stat;
  factors_t *LUfactors;
  
  if (itrans == 0) {
    trans = NOTRANS;
  } else if (itrans ==1) {
    trans = TRANS;
  } else if (itrans ==2) {
    trans = CONJ;
  } else {
    trans = NOTRANS;
  }
  /* Initialize the statistics variables. */
  StatInit(&stat);
  
  /* Extract the LU factors in the factors handle */
  LUfactors = (factors_t*) f_factors;
  L      = LUfactors->L;
  U      = LUfactors->U;
  perm_c = LUfactors->perm_c;
  perm_r = LUfactors->perm_r;
  
  dCreate_Dense_Matrix(&B, n, nrhs, b, ldb, SLU_DN, SLU_D, SLU_GE);
  /* Solve the system A*X=B, overwriting B with X. */
  dgstrs(trans, L, U, perm_c, perm_r, &B, &stat, &info);
  
  Destroy_SuperMatrix_Store(&B);
  StatFree(&stat);
#else
  fprintf(stderr," SLU Not Configured, fix make.inc and recompile\n");
  info=-1;
#endif
  return(info);
}


int amg_dslu_free(void *f_factors)
{
  /* 
   * This routine can be called from Fortran.
   *
   *      free all storage in the end
   *
   */
#ifdef Have_SLU_ 
  factors_t *LUfactors; 
  
  /* Free the LU factors in the factors handle */
  LUfactors = (factors_t*) f_factors;
  if (LUfactors != NULL) {
    SUPERLU_FREE (LUfactors->perm_r);
    SUPERLU_FREE (LUfactors->perm_c);
    Destroy_SuperNode_Matrix(LUfactors->L);
    Destroy_CompCol_Matrix(LUfactors->U);
    SUPERLU_FREE (LUfactors->L);
    SUPERLU_FREE (LUfactors->U);
    SUPERLU_FREE (LUfactors);
  }
  return(0);
#else
  fprintf(stderr," SLU Not Configured, fix make.inc and recompile\n");
  return(-1);
#endif
}


