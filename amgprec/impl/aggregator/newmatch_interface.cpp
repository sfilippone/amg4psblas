
#include <string.h>
#include <stdio.h>

#include "bcm.h"

bcm_CSRMatrix bootCMatch(bcm_CSRMatrix *C, int *match_algorithm, int *n_sweeps, int *max_nlevels, int *max_csize, bcm_Vector *w);
bcm_CSRMatrix bootCMatch(bcm_CSRMatrix *C, int *match_algorithm, int *n_sweeps, int *max_nlevels, int *max_csize, bcm_Vector *w){
   bcm_Vector *w_temp;
   int info;
   //double *w_inp;
   //w_inp=bcm_VectorData(w); 

   bcm_CSRMatrix *P;
   bcm_CSRMatrix *Ac;
   int ftcoarse=1;
   int cr_it=0, cr_relax_type=0;
   double cr_relax_weight=0.0;
   // Here I am building Ac but I won't use it. 
   Ac=bcm_CSRMatchingAgg(C, &w, &P, *match_algorithm, *n_sweeps,*max_csize,  *max_nlevels,  &ftcoarse,
			 cr_it, cr_relax_type, cr_relax_weight);
   //w_inp=bcm_VectorData(w);
   bcm_CSRMatrixDestroy(Ac);
 return *P;
}

int mld_bootCMatch_if(bcm_CSRMatrix *C, int match_algorithm, int n_sweeps,
		      int max_nlevels, int max_csize, bcm_Vector *w,
		      int isz, int ilaggr[], double valaggr[], int *num_cols){
   bcm_Vector *w_temp;
   int info;
   //double *w_inp;
   //w_inp=bcm_VectorData(w); 

   bcm_CSRMatrix *P;
   bcm_CSRMatrix *Ac;
   int *irp, *ja, nr, nz, nc,i,j;
   double *val;
   int ftcoarse=1;
   int cr_it=0, cr_relax_type=0;
   double cr_relax_weight=0.0;

   // Sanity checks
   nr = bcm_CSRMatrixNumRows(C);
   nc = bcm_VectorSize(w);
//   fprintf(stderr,"Sanity check:  %d   %d \n",nr,nc);

   // Here I am building Ac but I won't use it. 
   Ac=bcm_CSRMatchingAgg(C, &w, &P, match_algorithm, n_sweeps, max_csize, max_nlevels, &ftcoarse,
			 cr_it, cr_relax_type, cr_relax_weight);
   irp = bcm_CSRMatrixI(P);
   ja  = bcm_CSRMatrixJ(P);
   val = bcm_CSRMatrixData(P);
   nr  = bcm_CSRMatrixNumRows(P);
   nc  = bcm_CSRMatrixNumCols(P);
   nz  = bcm_CSRMatrixNumNonzeros(P);

   if (isz < nr) return(-1);
   if (nz != nr) return(-2);
   /* loop here only makes sense when nr==nz */
   for (i=0; i< nr; i++) {
     for (j=irp[i]; j<irp[i+1]; j++) {
       ilaggr[i]  = ja[j] + 1;
       valaggr[i] = val[j];
     }
   }
   *num_cols = nc;
   //w_inp=bcm_VectorData(w);
   bcm_CSRMatrixDestroy(Ac);
   bcm_CSRMatrixDestroy(P);
   return(0);
}
