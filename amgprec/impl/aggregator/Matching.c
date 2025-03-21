#include "amg_config.h"
#if defined(PSB_SERIAL_MPI)
#include <stdlib.h>
#include "psb_types.h"

psb_l_t d_trymatch(psb_l_t rowindex, psb_l_t colindex,
		     psb_l_t nrows_W, psb_l_t *W_i, psb_l_t *W_j, psb_d_t *W_data,
		     psb_l_t *jrowindex, psb_l_t ljrowindex,
		     psb_l_t *jcolindex, psb_l_t ljcolindex, psb_l_t *rmatch)
{

  psb_l_t tryrowmatch, trycolmatch;
  psb_l_t i, j, k, nzrow_W, startj, kindex;
  psb_d_t cweight, nweight; 

  //  psb_l_t     *W_i = bcm_CSRMatrixI(W);
  //  psb_l_t     *W_j = bcm_CSRMatrixJ(W);
  //  psb_l_t     nrows_W  =  bcm_CSRMatrixNumRows(W);
  //  psb_d_t  *W_data=bcm_CSRMatrixData(W);

  k=-1;
  i=0;
  while (i<ljrowindex && k ==-1)
    {
      if(jrowindex[i]==colindex) k=i;
      i++;
    }
  if(k >= 0)
    {
      ljrowindex=ljrowindex-1;
      for(i=k; i<ljrowindex; ++i) jrowindex[i]=jrowindex[i+1];
    }

  k=-1;
  i=0;
  while(i<ljcolindex && k==-1)
    {
      if(jcolindex[i]==rowindex) k=i;
      i++;
    }
  if(k >= 0)
    {
      ljcolindex=ljcolindex-1;
      for(i=k; i<ljcolindex; ++i) jcolindex[i]=jcolindex[i+1];
    }

  nzrow_W=W_i[rowindex+1]-W_i[rowindex];
  startj=W_i[rowindex];
  kindex=-1;
  k=0;
  while(k<nzrow_W && kindex==-1) 
    {
      if(W_j[startj+k]==colindex) kindex=startj+k;
      k++;
    }
  cweight=W_data[kindex];

  while((rmatch[rowindex]==-1 && rmatch[colindex]==-1)
	&& (ljrowindex !=0 || ljcolindex !=0))    
    {

      if (rmatch[rowindex] == -1 && ljrowindex !=0)
	{
	  tryrowmatch=jrowindex[0];
	  nzrow_W=W_i[rowindex+1]-W_i[rowindex];
	  startj=W_i[rowindex];
	  k=0;
	  kindex=-1;
	  while(k<nzrow_W && kindex==-1) 
	    {
	      if(W_j[startj+k]==tryrowmatch) kindex=startj+k;
	      k++;
	    }
	  nweight=W_data[kindex];

	  ljrowindex=ljrowindex-1;
	  for(i=0; i<ljrowindex; ++i) jrowindex[i]=jrowindex[i+1];
        
	  if(nweight > cweight && rmatch[tryrowmatch]==-1)
	    {
          
              nzrow_W=W_i[tryrowmatch+1]-W_i[tryrowmatch];
	      psb_l_t *trymatchindexrow;
	      trymatchindexrow= (psb_l_t *) calloc(nzrow_W, sizeof(psb_l_t));

	      startj=W_i[tryrowmatch];
	      for(k=0; k<nzrow_W; ++k) trymatchindexrow[k]=W_j[startj+k]; 

	      d_trymatch(rowindex,tryrowmatch,
			   nrows_W, W_i, W_j, W_data,
			   jrowindex,ljrowindex,
			   trymatchindexrow,nzrow_W,rmatch);
	      free(trymatchindexrow);
	    }
	}

      if(rmatch[colindex]==-1 && ljcolindex!=0)
	{
	  trycolmatch=jcolindex[0];
	  nzrow_W=W_i[colindex+1]-W_i[colindex];
	  startj=W_i[colindex];
	  k=0;
	  kindex=-1;
	  while(k<nzrow_W && kindex==-1) 
	    {
	      if(W_j[startj+k]==trycolmatch) kindex=startj+k;
	      k++;
	    }
	  nweight=W_data[kindex];

	  ljcolindex=ljcolindex-1;
	  for(i=0; i<ljcolindex; ++i) jcolindex[i]=jcolindex[i+1];
        
	  if(nweight > cweight && rmatch[trycolmatch]==-1)
	    {
	      nzrow_W=W_i[trycolmatch+1]-W_i[trycolmatch];
	      psb_l_t *trymatchindexcol;
	      trymatchindexcol= (psb_l_t *) calloc(nzrow_W, sizeof(psb_l_t));

	      startj=W_i[trycolmatch];
	      for(k=0; k<nzrow_W; ++k) trymatchindexcol[k]=W_j[startj+k];

	      d_trymatch(colindex,trycolmatch,
			   nrows_W, W_i, W_j, W_data,			   
			   jcolindex,ljcolindex,
			   trymatchindexcol,nzrow_W,rmatch);
	      free(trymatchindexcol);
	    }
	}
    }
        
  if(rmatch[rowindex]==-1 & rmatch[colindex]==-1)
    {
      rmatch[rowindex]=colindex;
      rmatch[colindex]=rowindex;
    }

  return 0;
}




psb_l_t *d_CSRMatrixHMatch( psb_l_t  nrows_B,  psb_l_t  ncols_B,
			       psb_l_t  nnz_B,  psb_l_t     *B_i,
			       psb_l_t     *B_j, psb_d_t  *B_data)
{

  psb_l_t i, j, k, *rmatch;
  psb_d_t  *c, alpha;
  psb_l_t jbp, nzrows_B;
  psb_d_t tmp=0.0;
  psb_l_t rno, cno, nzrows_cno, startj, ljrowindex, ljjrowindex;
  psb_l_t *jrowindex, *jjrowindex;

#if 0
  psb_l_t     *B_i = bcm_CSRMatrixI(B);
  psb_l_t     *B_j = bcm_CSRMatrixJ(B);
  psb_d_t  *B_data = bcm_CSRMatrixData(B);
  psb_l_t  nrows_B  =  bcm_CSRMatrixNumRows(B);
  psb_l_t  ncols_B  =  bcm_CSRMatrixNumCols(B);
  psb_l_t  nnz_B  =  bcm_CSRMatrixNumNonzeros(B);
#endif
  //  assert(nrows_B==ncols_B);

  rmatch = (psb_l_t *) calloc(nrows_B, sizeof(psb_l_t));

  for(i=0; i<nrows_B; ++i)  rmatch[i]=-1;

  jrowindex = (psb_l_t *) calloc(nrows_B, sizeof(psb_l_t));
  jjrowindex = (psb_l_t *) calloc(nrows_B, sizeof(psb_l_t));

  jbp=0;
  for(i=0; i< nrows_B; ++i)
    { 
      nzrows_B=B_i[i+1]-B_i[i];

      for(j=0; j<nzrows_B; ++j) jrowindex[j]=B_j[jbp+j];
      for(j=0; j<nzrows_B; ++j)
	{
	  rno=i;
	  cno=B_j[jbp+j];

	  startj=B_i[cno];
	  nzrows_cno=B_i[cno+1]-startj;

	  for(k=0; k<nzrows_cno; ++k) jjrowindex[k]=B_j[startj+k];

	  if(rmatch[rno] == -1 && rmatch[cno] == -1)
	    d_trymatch(rno,cno,nrows_B,B_i,B_j,B_data,
			 jrowindex,nzrows_B,jjrowindex,nzrows_cno,rmatch);
        
	}

      jbp=jbp+nzrows_B;
    }
   
  free(jrowindex);
  free(jjrowindex);
  
  return rmatch;
}


void dMatching(psb_l_t NLVer, psb_l_t NLEdge,
	       psb_l_t *verLocPtr, psb_l_t *verLocInd, psb_d_t *edgeLocWeight,
	       psb_l_t *verDistance, psb_l_t *Mate)
{
  psb_l_t *lmate;
  lmate = d_CSRMatrixHMatch( NLVer, NLVer,   NLEdge,
			       verLocPtr, verLocInd,edgeLocWeight);
  for (psb_l_t i=0; i<NLVer; i++){
    Mate[i] = lmate[i];
  }
  free(lmate);
  return;
}

void sMatching(psb_l_t NLVer, psb_l_t NLEdge,
	       psb_l_t *verLocPtr, psb_l_t *verLocInd, psb_s_t *edgeLocWeight,
	       psb_l_t *verDistance, psb_l_t *Mate)
{

  return;
}

#endif
