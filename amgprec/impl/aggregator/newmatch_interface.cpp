#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "psb_base_cbind.h"
#include "MatchingAlgorithms.h"

#ifdef __cplusplus
extern "C" {
#endif

psb_i_t dnew_Match_If(psb_i_t ipar, psb_i_t matching, psb_d_t lambda,
		      psb_i_t nr, psb_i_t irp[], psb_i_t ja[],
		      psb_d_t val[], psb_d_t diag[],
		      psb_d_t w[], psb_i_t mate[]);

#ifdef __cplusplus
  }
#endif
  
psb_i_t dnew_Match_If(psb_i_t ipar, psb_i_t matching, psb_d_t lambda,
		      psb_i_t nr, psb_i_t irp[], psb_i_t ja[],
		      psb_d_t val[], psb_d_t diag[], psb_d_t w[],
		      psb_i_t mate[])
{
   psb_i_t info;
   psb_i_t i,j,k;
   psb_i_t ftcoarse=1;
   psb_i_t cr_it=0, cr_relax_type=0;
   psb_d_t cr_relax_weight=0.0;

   vector<NODE_T> s;
   vector<NODE_T> t;
   vector<VAL_T> weights;
   vector<NODE_T> mateNode;
   NODE_T u,v;
   VAL_T weight;
   psb_i_t preprocess = matching; // 0 no greedy   1 greedy
   psb_i_t romaInput  = ipar; // 1 sequential  2 parallel 
   // VAL_T   lambda     = 2; // positive real value
   psb_d_t aii, ajj, aij, wii, wjj, tmp1, tmp2, minabs, edgnrm;
   psb_i_t nt;   // number of threads, got with 1 for testing purposes.
   psb_d_t timeDiff;
   MatchStat pstat;
   double eps=1e-16;
   double minweight,maxweight;
   char *numthreadsenv;

   numthreadsenv=getenv("OMP_NUM_THREADS");
   if (numthreadsenv) {
     sscanf(numthreadsenv,"%d",&nt);
   } else {
     nt = 1;
   }

   minabs = 1e300;
   //   fprintf(stderr,"Sanity check:  %d   %d \n",nr,nc);
   k=0;
   for (i=1; i<nr; i++) {
     for (j=irp[i-1]; j<irp[i]; j++)  {
       v = i-1;        // I 
       u = ja[j-1] - 1;    // J 
       if (v>u) {
	 // Define Ahat entry
	 aij = val[j-1];
	 aii = diag[v];
	 ajj = diag[u];
	 wii = w[v];
	 wjj = w[u];
	 edgnrm = aii*(wii*wii) + ajj*(wjj*wjj);
	 if (edgnrm > eps) {
	   weight = abs(1.0 - (2*1.0*aij*wii*wjj)/(aii*(wii*wii) + ajj*(wjj*wjj)));
	 } else {
	   weight = eps;
	 }
	 //
	 s.push_back(u);
	 t.push_back(v);
	 weights.push_back(weight);
	 k = k + 1 ;
	 if (weight<minabs) minabs=weight;
	 
       }
     }
   }
   
   maxweight = eps;
   minweight = 1e300;
   //fprintf(stderr,"minabs %g\n",minabs);
   for (i=0; i<k; i++) {
     weights[i] = log(weights[i]/(0.999*minabs));
     if (weights[i]>maxweight) maxweight=weights[i];
     if (weights[i]<minweight) minweight=weights[i];
   }

   if (lambda<0.0){
     lambda = maxweight-2.0*minweight+eps;
     if (lambda<0.0) lambda=eps;
   } else if (lambda >= 0 && lambda <= 1.0){
     lambda = lambda*eps + (1.0-lambda)*(fmax(maxweight-2.0*minweight,0.0) );
   }
   fprintf(stderr,"Calling matching:  pre %d nt %d  lambda %g   %g  %g\n",
	   preprocess,nt,lambda,maxweight,minweight);
   
   runRomaWrapper(s,t,weights, nr, mateNode,preprocess,romaInput,lambda ,nt, pstat, timeDiff);
   /* loop here only makes sense when nr==nz */
   for (i=0; i< nr; i++) {
     //fprintf(stderr,"From runRomaWrapper: %d   %d\n",i,mateNode[i]);
     if (mateNode[i]>=0) {
       mate[i] = mateNode[i]+1;
     } else {
       mate[i] = mateNode[i];
       //fprintf(stderr,"From runRomaWrapper: %d   %d\n",i,mateNode[i]);
     }
   }
   return(0);
}
