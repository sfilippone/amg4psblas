/*
                             AMG4PSBLAS version 1.2
    Algebraic Multigrid Package
               based on PSBLAS (Parallel Sparse BLAS version 3.9)

    (C) Copyright 2021

        Salvatore Filippone
        Pasqua D'Ambra
        Fabio Durastante

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:
      1. Redistributions of source code must retain the above copyright
         notice, this list of conditions and the following disclaimer.
      2. Redistributions in binary form must reproduce the above copyright
         notice, this list of conditions, and the following disclaimer in the
         documentation and/or other materials provided with the distribution.
      3. The name of the AMG4PSBLAS group or the names of its contributors may
         not be used to endorse or promote products derived from this
         software without specific written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AMG4PSBLAS GROUP OR ITS CONTRIBUTORS
    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.

    Includes material from BootCMatch, see copyright below.
*/

/*
                BootCMatch
     Bootstrap AMG based on Compatible weighted Matching, version 0.9
    (C) Copyright 2017
                       Pasqua D'Ambra         IAC-CNR, IT
                       Panayot S. Vassilevski Portland State University, OR USA

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
    1. Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above copyright
       notice, this list of conditions, and the following disclaimer in the
       documentation and/or other materials provided with the distribution.
    3. The name of the BootCMatch group or the names of its contributors may
       not be used to endorse or promote products derived from this
       software without specific written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE BootCMatch GROUP OR ITS CONTRIBUTORS
  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.

*/
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
