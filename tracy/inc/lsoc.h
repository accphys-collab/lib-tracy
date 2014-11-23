/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/

void prt_gcmat(const int plane);

void gcmat(const int plane);

void gcmat(const int n_bpm, const long int bpms[],
	   const int n_corr, const long int corrs[], const int plane,
	   const bool svd);

void gcmat(const int bpm, const int corr, const int plane);

void gcmat1(const int bpm, const int corr, const int plane);

void lsoc(const int niter, const int plane);

void lsoc(const int niter, const int bpm, const int corr, const int plane);

void lsoc1(const int niter, const int bpm, const int corr, const int plane);
