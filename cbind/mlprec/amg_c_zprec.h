#ifndef MLD_C_ZPREC_
#define MLD_C_ZPREC_

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
  typedef struct MLD_C_ZPREC {
    void *dprec;
  } mld_c_zprec; 
  
  mld_c_zprec* mld_c_zprec_new();
  psb_i_t mld_c_zprec_delete(mld_c_zprec* p);
 
  psb_i_t mld_c_zprecinit(psb_i_t ictxt, mld_c_zprec *ph, const char *ptype);
  psb_i_t mld_c_zprecseti(mld_c_zprec *ph, const char *what, psb_i_t val);
  psb_i_t mld_c_zprecsetc(mld_c_zprec *ph, const char *what, const char *val);
  psb_i_t mld_c_zprecsetr(mld_c_zprec *ph, const char *what, double val);
  psb_i_t mld_c_zprecbld(psb_c_dspmat *ah, psb_c_descriptor *cdh, mld_c_zprec *ph);
  psb_i_t mld_c_zhierarchy_build(psb_c_dspmat *ah, psb_c_descriptor *cdh, mld_c_zprec *ph);
  psb_i_t mld_c_zsmoothers_build(psb_c_dspmat *ah, psb_c_descriptor *cdh, mld_c_zprec *ph);
  psb_i_t mld_c_zprecfree(mld_c_zprec *ph);
  psb_i_t mld_c_zprecbld_opt(psb_c_zspmat *ah, psb_c_descriptor *cdh, 
			  mld_c_zprec *ph, const char *afmt);

  psb_i_t mld_c_zdescr(mld_c_zprec *ph);

  psb_i_t mld_c_zkrylov(const char *method, psb_c_zspmat *ah, mld_c_zprec *ph, 
		  psb_c_zvector *bh, psb_c_zvector *xh,
		  psb_c_descriptor *cdh, psb_c_SolverOptions *opt);


#ifdef __cplusplus
}
#endif

#endif
