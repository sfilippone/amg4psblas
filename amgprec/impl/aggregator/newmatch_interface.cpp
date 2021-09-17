#include <string.h>
#include <stdio.h>
#include "psb_base_cbind.h"
#include "MatchingAlgorithms.h"

#ifdef __cplusplus
extern "C" {
#endif

psb_i_t dnew_Match_If(psb_i_t nr, psb_i_t irp[], psb_i_t ja[],
		      psb_d_t val[], psb_d_t diag[],
		      psb_d_t w[], psb_i_t mate[])

#ifdef __cplusplus
  }
#endif
  
psb_i_t dnew_Match_If(psb_i_t nr, psb_i_t irp[], psb_i_t ja[],
		      psb_d_t val[], psb_d_t diag[], psb_d_t w[],
		      psb_i_t mate[])
{
   psb_i_t info;
   psb_i_t i,j;
   psb_i_t ftcoarse=1;
   psb_i_t cr_it=0, cr_relax_type=0;
   psb_d_t cr_relax_weight=0.0;

   vector<NODE_T> s;
   vector<NODE_T> t;
   vector<VAL_T> weights;
   vector<NODE_T> mateNode;
   NODE_T u,v;
   VAL_T weight;
   psb_i_t preprocess = atoi(argv[2]);
   psb_i_t romaInput = atoi(argv[3]);
   VAL_T lambda = atof(argv[4]);  
   psb_d_t aii, ajj, aij, wii, wjj, tmp1, tmp2, minabs, edgnrm;
   psb_i_t nt;
   psb_d_t timeDiff;
   MatchStat pstat;


   //   fprintf(stderr,"Sanity check:  %d   %d \n",nr,nc);
   for (i=1; i<nr; i++) {
     for (j=irp[i-1]; j<irp[i]; j++)  {
       v = i-1;        // I 
       u = ja[j-1];    // J 
       if (v>u) {
	 // Define Ahat entry
	 aij = val[j-1];
	 aii = diag[v];
	 ajj = diag[u];
	 wii = w[v];
	 wjj = w[u];
	 edgnrm = aii*(wii*wii) + ajj*(wjj*wjj);
	 if (edgnrm > eps) {
	   weight = 1.0 - (2*1.0*aij*wii*wjj)/(aii*(wii*wii) + ajj*(wjj*wjj));
	 } else {
	   weight = 1e-16;
	 }
	 //
	 s.push_back(u);
	 t.push_back(v);
	 weights.push_back(weight); 	 
       }
     }
   }
   runRomaWrapper(s,t,weights, nr, mateNode,preprocess,romaInput,lambda ,nt, pstat, timeDiff);
   if (isz < nr) return(-1);
   if (nz != nr) return(-2);
   /* loop here only makes sense when nr==nz */
   for (i=0; i< nr; i++) {
     mate[i] = mateNode[i]+1;
   }
   return(0);
}
