#include <stdlib.h>
#include "mld_c_dprec.h"

mld_c_dprec* mld_c_dprec_new()
{
  mld_c_dprec* temp;
  
  temp=(mld_c_dprec *) malloc(sizeof(mld_c_dprec));
  temp->dprec=NULL;
  return(temp);
}


int mld_c_dprec_delete(mld_c_dprec* p)
{
  int iret;
  iret=mld_c_dprecfree(p);
  if (iret ==0) free(p);
  return(iret);
}

