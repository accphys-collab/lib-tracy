/* Tracy-3

   J. Bengtsson, BNL 2007

   NO   1   link to the linear TPSA (nv_tps = 1)
          >1   link to Berz' TPSA

*/

// C standard library
#include <stdio.h>
#include <stddef.h>
#include <setjmp.h>
#include <time.h>
#include <memory.h>
#include <malloc.h>
//#include <execinfo.h>


// C++ standard library
#include <cstdlib>
#include <cfloat>
#include <cctype>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
using namespace std;

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_multimin.h>
#include "gslport.h"

// Tracy-3
#include "field.h"
#include "mathlib.h"

#if NO == 1
  // linear TPSA
  #include "tpsa_lin.h"
  #include "tpsa_lin_pm.h"
#else
  // interface to M. Berz' TPSA
  #include "tpsa_for.h"
  #include "tpsa_for_pm.h"
#endif

#include "tracy.h"
#include "tracy_global.h"
#include "ety.h"

#include "num_rec.h"

#include "radia2tracy.h"
#include "pascalio.h"

#include "t2elem.h"
#include "t2cell.h"
#include "t2lat.h"
#include "t2ring.h"

#include "fft.h"

#include "physlib.h"
#include "nsls-ii_lib.h"

#include "lsoc.h"

#include "modnaff.h"

#include "naffutils.h"
#include "complexeheader_naff.h"

#include "soleillib.h"

#include "rdmfile.h"
#include "prtmfile.h"


// Truncated Power Series Algebra (TPSA)
extern const int  nv_tps, nd_tps, iref_tps;
extern int        no_tps, ndpt_tps;
extern double     eps_tps;

extern ElemFamType ElemFam[];

extern CellType Cell[];

extern globvalrec globval;
