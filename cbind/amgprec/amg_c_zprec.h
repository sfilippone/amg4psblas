#ifndef AMG_C_ZPREC_
#define AMG_C_ZPREC_

#include "amg_const.h"
#include "psb_base_cbind.h"
#include "psb_prec_cbind.h"
#include "psb_krylov_cbind.h"

/* Object handle related routines */
/* Note:  amg_get_XXX_handle returns:  <= 0  unsuccessful */
/*                                     >0    valid handle */
#ifdef __cplusplus
extern "C" {
#endif
  typedef struct AMG_C_ZPREC {
    void *dprec;
  } amg_c_zprec; 
  
  amg_c_zprec* amg_c_zprec_new();
  psb_i_t amg_c_zprec_delete(amg_c_zprec* p);
 
  psb_i_t amg_c_zprecinit(psb_c_ctxt cctxt, amg_c_zprec *ph, const char *ptype);
  psb_i_t amg_c_zprecseti(amg_c_zprec *ph, const char *what, psb_i_t val);
  psb_i_t amg_c_zprecsetc(amg_c_zprec *ph, const char *what, const char *val);
  psb_i_t amg_c_zprecsetr(amg_c_zprec *ph, const char *what, double val);
  psb_i_t amg_c_zprecbld(psb_c_dspmat *ah, psb_c_descriptor *cdh, amg_c_zprec *ph);
  psb_i_t amg_c_zhierarchy_build(psb_c_dspmat *ah, psb_c_descriptor *cdh, amg_c_zprec *ph);
  psb_i_t amg_c_zsmoothers_build(psb_c_dspmat *ah, psb_c_descriptor *cdh, amg_c_zprec *ph);
  psb_i_t amg_c_zprecfree(amg_c_zprec *ph);
  psb_i_t amg_c_zprecbld_opt(psb_c_zspmat *ah, psb_c_descriptor *cdh, 
			  amg_c_zprec *ph, const char *afmt);

  psb_i_t amg_c_zdescr(amg_c_zprec *ph);

  psb_i_t amg_c_zkrylov(const char *method, psb_c_zspmat *ah, amg_c_zprec *ph, 
		  psb_c_zvector *bh, psb_c_zvector *xh,
		  psb_c_descriptor *cdh, psb_c_SolverOptions *opt);


#ifdef __cplusplus
}
#endif

#endif
