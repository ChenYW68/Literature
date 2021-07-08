#include <cmath>
#include <Rcpp.h>
#include <string>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#include <limits>


#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

void zeros(double *a, int n){
  for(int i = 0; i < n; i++)
    a[i] = 0.0;
}
// [[Rcpp::export]]
SEXP theta2_fun_C(SEXP m_r, SEXP n_r, 
                  SEXP theta2_r, SEXP ds_r,
                  SEXP D_r, SEXP bandKernel_r,
                  SEXP Q_r, SEXP S00_r, SEXP S01_r, 
				  SEXP theta1_r, SEXP sigTheta1_r, 
				  SEXP nThreads_r){
	int m = INTEGER(m_r)[0];
	int n = INTEGER(n_r)[0];	
	double *theta2 = REAL(theta2_r);			  
	double ds = REAL(ds_r)[0];			  
	double *D = REAL(D_r);	
	double *bandKernel = REAL(bandKernel_r);	
	double *Q = REAL(Q_r);	
	double *S00 = REAL(S00_r);	
	double *S01 = REAL(S01_r);
	double theta1 = REAL(theta1_r)[0];
    double sigTheta1 = REAL(sigTheta1_r)[0];
    int nThreads = INTEGER(nThreads_r)[0];	

   int threadID = 0;


 #ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif

    const int inc = 1;
    int i, j, k, l, nProtect=0;
	double *M = (double *) R_alloc(m*m*nThreads, sizeof(double)); //zeros(M, m*m*nThreads);
	double *Qt = (double *) R_alloc(m*m*nThreads, sizeof(double)); //zeros(Qt, m*m*nThreads);
	double *Qm = (double *) R_alloc(m*m*nThreads, sizeof(double)); //zeros(Qm, m*m*nThreads);
	double *M_Q_M = (double *) R_alloc(m*m*nThreads, sizeof(double)); //zeros(M_Q_M, m*m*nThreads);
	double *Q_S0 = (double *) R_alloc(m*m*nThreads, sizeof(double)); //zeros(Q_S0, m*m*nThreads);
	double *Q_S1 = (double *) R_alloc(m*m*nThreads, sizeof(double)); //zeros(Q_S1, m*m*nThreads);
	
	SEXP logLik_r;
	PROTECT(logLik_r = Rf_allocVector(REALSXP, n)); nProtect++;
	
	//double *f = (double *) R_alloc(n*nThreads, sizeof(double)); zeros(f, n*nThreads);
	char const *ntran = "N";
	const double one = 1.0;
	const double zero = 0.0;
	char Trans = 'T';
	double f1, f2;
	int mm = m*m;
	double f;
	
#ifdef _OPENMP
#pragma omp parallel for private(i, j, l, f, f1, f2, threadID)
#endif 
	for(k = 0; k < n; k++){
#ifdef _OPENMP
	threadID = omp_get_thread_num();
#endif 
		for(i = 0; i < m; i++){
			 for(j = i; j < m; j++){
			   M[m*m*threadID + i*m + j] = ds*exp(-pow(D[i*m + j], 2)/theta2[k])*bandKernel[i*m + j];
			   M[m*m*threadID + j*m + i] = M[m*m*threadID + i*m + j];
			}
		}
		//F77_NAME(dcopy)(&mm, Q, &inc, &Qt[m*m*threadID], &inc);
		//zeros(&Qm[m*m*threadID], m*m);
		
		F77_NAME(dgemm)(ntran, ntran, &m, &m, &m, &one, Q, &m, &M[m*m*threadID], &m, &zero, &Qm[m*m*threadID], &m);
		//Rprintf("f1 =  %3.10f; f2 =  %3.10f;  \n", Q_M[m*m*threadID], Q_M[m*m*threadID + 15]);
		F77_NAME(dgemm)(ntran, ntran, &m, &m, &m, &one, &Qm[m*m*threadID], &m, S01, &m, &zero, &Q_S0[m*m*threadID], &m);
		// zeros(M_Q_M, m*m*nThreads);
		F77_NAME(dgemm)(&Trans, ntran, &m, &m, &m, &one, &M[m*m*threadID], &m, &Qm[m*m*threadID], &m, &zero, &M_Q_M[m*m*threadID], &m);
		//zeros(Q_S1, m*m*nThreads);
		
		F77_NAME(dgemm)(ntran, ntran, &m, &m, &m, &one, &M_Q_M[m*m*threadID], &m, S00, &m, &zero, &Q_S1[m*m*threadID], &m);
		
		f1 = 0.0;
		f2 = 0.0;
		for(l = 0; l < m; l++){
		 // for(j = i; j < i + 1; j++){
		    f1 += theta1*Q_S0[m*m*threadID + l*m + l];
			f2 -= (theta1 + pow(sigTheta1, 2))*Q_S1[m*m*threadID + l*m + l]/2;
		 //}
		}
			//Rprintf("f1 =  %3.5f; f2 =  %3.5f;  \n", f1, f2);	
		//f[k*threadID] = f1 + f2;
		f =  f1 + f2;
		F77_NAME(dcopy)(&inc, &f, &inc, &REAL(logLik_r)[k], &inc);  
			
	}
		
		
	SEXP result_r, resultName_r;
    int nResultListObjs = 1;

    PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    SET_VECTOR_ELT(result_r, 0, logLik_r);
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("logLik"));
    Rf_namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
}


















