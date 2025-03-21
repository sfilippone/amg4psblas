#include "amg_config.h"

#if defined(PSB_SERIAL_MPI)
#include "psb_types.h"
void dMatching(psb_l_t NLVer, psb_l_t NLEdge,
	       psb_l_t *verLocPtr, psb_l_t *verLocInd, psb_d_t *edgeLocWeight,
	       psb_l_t *verDistance, psb_l_t *Mate);
void sMatching(psb_l_t NLVer, psb_l_t NLEdge,
	       psb_l_t *verLocPtr, psb_l_t *verLocInd, psb_s_t *edgeLocWeight,
	       psb_l_t *verDistance, psb_l_t *Mate);
#endif
