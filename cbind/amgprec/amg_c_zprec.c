#include <stdlib.h>
#include "amg_c_zprec.h"

amg_c_zprec* amg_c_new_zprec()
{
  amg_c_zprec* temp;
  
  temp=(amg_c_zprec *) malloc(sizeof(amg_c_zprec));
  temp->dprec=NULL;
  return(temp);
}


psb_c_i_t amg_c_delete_zprec(amg_c_zprec* p)
{
  int iret;
  iret=amg_c_zprecfree(p);
  if (iret ==0) free(p);
  return(iret);
}

