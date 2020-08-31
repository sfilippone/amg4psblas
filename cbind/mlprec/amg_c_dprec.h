#ifndef MLD_C_DPREC_
#define MLD_C_DPREC_

#include "mld_const.h"
#include "psb_base_cbind.h"
#include "psb_prec_cbind.h"
#include "psb_krylov_cbind.h"

/* Object handle related routines */
/* Note:  mld_get_XXX_handle returns:  <= 0  unsuccessful */
/*                                     >0    valid handle */
#ifdef __cplusplus
extern "C" {
#endif
  typedef struct MLD_C_DPREC {
    void *dprec;
  } mld_c_dprec; 
  
  mld_c_dprec* mld_c_dprec_new();
  psb_i_t mld_c_dprec_delete(mld_c_dprec* p);
 
  psb_i_t mld_c_dprecinit(psb_i_t ictxt, mld_c_dprec *ph, const char *ptype);
  psb_i_t mld_c_dprecseti(mld_c_dprec *ph, const char *what, psb_i_t val);
  psb_i_t mld_c_dprecsetc(mld_c_dprec *ph, const char *what, const char *val);
  psb_i_t mld_c_dprecsetr(mld_c_dprec *ph, const char *what, double val);
  psb_i_t mld_c_dprecbld(psb_c_dspmat *ah, psb_c_descriptor *cdh, mld_c_dprec *ph);
  psb_i_t mld_c_dhierarchy_build(psb_c_dspmat *ah, psb_c_descriptor *cdh, mld_c_dprec *ph);
  psb_i_t mld_c_dsmoothers_build(psb_c_dspmat *ah, psb_c_descriptor *cdh, mld_c_dprec *ph);
  psb_i_t mld_c_dprecfree(mld_c_dprec *ph);
  psb_i_t mld_c_dprecbld_opt(psb_c_dspmat *ah, psb_c_descriptor *cdh, 
			  mld_c_dprec *ph, const char *afmt);
  psb_i_t mld_c_ddescr(mld_c_dprec *ph);

  psb_i_t mld_c_dkrylov(const char *method, psb_c_dspmat *ah, mld_c_dprec *ph, 
		  psb_c_dvector *bh, psb_c_dvector *xh,
		  psb_c_descriptor *cdh, psb_c_SolverOptions *opt);


#ifdef __cplusplus
}
#endif

#endif
