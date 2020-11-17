#include <stdlib.h>
#include "amg_c_dprec.h"

amg_c_dprec* amg_c_dprec_new()
{
  amg_c_dprec* temp;
  
  temp=(amg_c_dprec *) malloc(sizeof(amg_c_dprec));
  temp->dprec=NULL;
  return(temp);
}


int amg_c_dprec_delete(amg_c_dprec* p)
{
  int iret;
  iret=amg_c_dprecfree(p);
  if (iret ==0) free(p);
  return(iret);
}

