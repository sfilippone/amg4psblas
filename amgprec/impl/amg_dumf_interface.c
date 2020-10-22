/*
 * 
 *                            AMG4PSBLAS version 1.0
 *   Algebraic Multigrid Package
 *              based on PSBLAS (Parallel Sparse BLAS version 3.5)
 *   
 *   (C) Copyright 2020
 *
 *      Salvatore Filippone                              
 *      Pasqua D'Ambra                                  
 *      Fabio Durastante        
 * 
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *    1. Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *    2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions, and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *    3. The name of the AMG4PSBLAS group or the names of its contributors may
 *       not be used to endorse or promote products derived from this
 *       software without specific written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AMG4PSBLAS GROUP OR ITS CONTRIBUTORS
 *  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 * 
 *
 * File: amg_umf_interface.c
 *
 * Functions: amg_dumf_fact_, amg_dumf_solve_, amg_umf_free_.
 *
 * This file is an interface to the UMFPACK routines for sparse factorization and
 * solve. It was obtained by adapting umfpack_di_demo under the original UMFPACK
 * copyright terms reproduced below.
 * 
 */

/*		=====================
UMFPACK Version 4.4 (Jan. 28, 2005), Copyright (c) 2005 by Timothy A.
Davis.  All Rights Reserved.

UMFPACK License:

    Your use or distribution of UMFPACK or any modified version of
    UMFPACK implies that you agree to this License.

    THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
    EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

    Permission is hereby granted to use or copy this program, provided
    that the Copyright, this License, and the Availability of the original
    version is retained on all copies.  User documentation of any code that
    uses UMFPACK or any modified version of UMFPACK code must cite the
    Copyright, this License, the Availability note, and "Used by permission."
    Permission to modify the code and to distribute modified code is granted,
    provided the Copyright, this License, and the Availability note are
    retained, and a notice that the code was modified is included.  This
    software was developed with support from the National Science Foundation,
    and is provided to you free of charge.

Availability:

    http://www.cise.ufl.edu/research/sparse/umfpack

*/


#include <stdio.h>
#ifdef Have_UMF_		 
#include "umfpack.h"
#endif

int amg_dumf_fact(int n, int nnz,
		  double *values, int *rowind, int *colptr,
		  void **symptr, void **numptr,
		  long long int *ssize,
		  long long int *nsize)

{
 
#ifdef Have_UMF_
  double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
  void *Symbolic, *Numeric ;
  int i, info;
  
  
  umfpack_di_defaults(Control);
  
  info = umfpack_di_symbolic (n, n, colptr, rowind, values, &Symbolic,
				Control, Info);
  
    
  if ( info == UMFPACK_OK ) {
    info = 0;
  } else {
    printf("umfpack_di_symbolic() error returns INFO= %d\n", info);
    umfpack_di_report_status(Control,info);
    *symptr = (void *) NULL; 
    *numptr = (void *) NULL; 
    return -11;
  }
    
  *symptr = Symbolic; 
  *ssize  = Info[UMFPACK_SYMBOLIC_SIZE]; 
  *ssize *= Info[UMFPACK_SIZE_OF_UNIT]; 

  info = umfpack_di_numeric (colptr, rowind, values, Symbolic, &Numeric,
				Control, Info) ;
  
    
  if ( info == UMFPACK_OK ) {
    info = 0;
    *numptr =  Numeric; 
    *nsize  = Info[UMFPACK_NUMERIC_SIZE]; 
    *nsize *= Info[UMFPACK_SIZE_OF_UNIT]; 

  } else {
    printf("umfpack_di_numeric() error returns INFO= %d\n", info);
    umfpack_di_report_status(Control,info);
    info = -12;
    *numptr =  NULL; 
  }


  return info;
  
#else
    fprintf(stderr," UMF Not Configured, fix make.inc and recompile\n");
    return -1;
#endif    
}


int amg_dumf_solve(int itrans, int n,  
                 double *x,  double *b, int ldb,
		 void *numptr)

{
#ifdef Have_UMF_ 
  double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
  void *Symbolic, *Numeric ;
  int i,trans, info;
  
  
  umfpack_di_defaults(Control);
  Control[UMFPACK_IRSTEP]=0;


  if (itrans == 0) {
    trans = UMFPACK_A;
  } else if (itrans ==1) {
    trans = UMFPACK_At;
  } else {
    trans = UMFPACK_A;
  }

  info = umfpack_di_solve(trans,NULL,NULL,NULL,
			  x,b,numptr,Control,Info);
  return info;
#else
  fprintf(stderr," UMF Not Configured, fix make.inc and recompile\n");
  return -1;
#endif
    
}


int amg_dumf_free(void *symptr,	 void *numptr)

{
#ifdef Have_UMF_ 
  void *Symbolic, *Numeric ;
  Symbolic = symptr;
  Numeric  = numptr;
  
  if (numptr != NULL) umfpack_di_free_numeric(&Numeric);
  if (symptr != NULL) umfpack_di_free_symbolic(&Symbolic);
  return 0;
#else
  fprintf(stderr," UMF Not Configured, fix make.inc and recompile\n");
  return -1;
#endif
}


