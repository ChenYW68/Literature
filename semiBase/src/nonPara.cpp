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
void error(std::string error_type)
{
  std::cout << error_type << '\n';
}
void zeros(double *a, int n){
  for(int i = 0; i < n; i++)
    a[i] = 0.0;
}

//which index of b equals a, where b is of length n
int which(int a, int *b, int n){
  int i;
  for(i = 0; i < n; i++){
    if(a == b[i]){
      return(i);
    }
  }
  
  error("c++ error: which failed");
  return -9999;
}

void mvrnorm(double *des, double *mu, double *cholCov, int dim){
  
  int i;
  int inc = 1;
  double one = 1.0;
  
  for(i = 0; i < dim; i++){
    des[i] = Rf_rnorm(0, 1);
  }
  
  F77_NAME(dtrmv)("L", "N", "N", &dim, cholCov, &dim, des, &inc);
  F77_NAME(daxpy)(&dim, &one, mu, &inc, des, &inc);
}

double logit(double theta, double a, double b){
  return log((theta-a)/(b-theta));
}

double logitInv(double z, double a, double b){
  return b-(b-a)/(1+exp(z));
}

double dist2(double &a1, double &a2, double &b1, double &b2){
  return(sqrt(pow(a1-b1,2)+pow(a2-b2,2)));
}

std::string getCorName(int i){
  
  if(i == 0){
    return "exponential";
  }else if(i == 1){
    return "spherical";
  }else if(i == 2){
    return "matern";
  }else if(i == 3){
    return "gaussian";
  }else{
    error("c++ error: cov.model is not correctly specified");
  }
  
}
void getNNIndx(int i, int m, int &iNNIndx, int &iNN){
  
  if(i == 0){
    iNNIndx = 0;//this should never be accessed
    iNN = 0;
    return;
  }else if(i < m){
    iNNIndx = static_cast<int>(static_cast<double>(i)/2*(i-1));
    iNN = i;
    return;
  }else{
    iNNIndx = static_cast<int>(static_cast<double>(m)/2*(m-1)+(i-m)*m);
    iNN = m;
    return;
  } 
}

void mkUIndx0(int n, int m, int* nnIndx, int* uIndx, int* uIndxLU){ 
  
  int iNNIndx, iNN, i, j, k, l, h;
  
  for(i = 0, l = 0; i < n; i++){    
    uIndxLU[i] = l; 
    for(j = 0, h = 0; j < n; j++){   
      getNNIndx(j, m, iNNIndx, iNN);  
      for(k = 0; k < iNN; k++){      	
        if(nnIndx[iNNIndx+k] == i){
          uIndx[l+h] = j;
          h++;
        }    
      }
    }
    l += h;
    uIndxLU[n+i] = h;
    R_CheckUserInterrupt();
  }
}

void mkUIndx1(int n, int m, int* nnIndx, int* uIndx, int* uIndxLU){ 
  
  int iNNIndx, iNN, i, j, k, l, h;
  
  for(i = 0, l = 0; i < n; i++){    
    uIndxLU[i] = l; 
    for(j = n-1, h = 0; j > i; j--){   
      getNNIndx(j, m, iNNIndx, iNN);  
      for(k = 0; k < iNN; k++){      	
        if(nnIndx[iNNIndx+k] == i){
          uIndx[l+h] = j;
          h++;
        }    
      }
    }
    l += h;
    uIndxLU[n+i] = h;
    R_CheckUserInterrupt();
  }
}
void crs_csc(int n, int *i_A, int *j_A, int *i_B, int *j_B){
  
  int i, j, col, cumsum, temp, row, dest, last;
  
  int nnz = i_A[n];
  
  for(i = 0; i < n; i++){
    i_B[i] = 0;
  }
  
  for(i = 0; i < nnz; i++){            
    i_B[j_A[i]]++;
  }
  
  // cumsum the nnz per column to get i_B[]
  for(col = 0, cumsum = 0; col < n; col++){     
    temp  = i_B[col];
    i_B[col] = cumsum;
    cumsum += temp;
  }
  i_B[n] = nnz; 
  
  for(row = 0; row < n; row++){
    for(j = i_A[row]; j < i_A[row+1]; j++){
      col  = j_A[j];
      dest = i_B[col];
      
      j_B[dest] = row;
      i_B[col]++;
    }
  }  
  
  for(col = 0, last = 0; col <= n; col++){
    temp  = i_B[col];
    i_B[col] = last;
    last = temp;
  }
} 


void mkUIndx2(int n, int m, int* nnIndx, int *nnIndxLU, int* uIndx, int* uIndxLU){ 
  
  int i, k;
  int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
  
  //int *j_A = new int[nIndx]; is nnIndx
  int *i_nnIndx = new int[n+1];
  //int *j_A_csc = new int[nIndx];//uIndx
  int *i_A_csc = new int[n+1];
  
  for(i = 0, k = 0; i < n; i++){
    if(nnIndxLU[n+i] == 0){//excludes rows with no elements, i.e., the first row because it is zero by design A[0,0] = 0
      i_nnIndx[0] = 0;
    }else{
      i_nnIndx[k] = i_nnIndx[k-1]+nnIndxLU[n+i-1];
    }
    k++;
  }
  i_nnIndx[n] = i_nnIndx[0]+nIndx;
  
  crs_csc(n, i_nnIndx, nnIndx, i_A_csc, uIndx);
  
  for(i = 0; i < n; i++){
    uIndxLU[i] = i_A_csc[i];
    uIndxLU[i+n] = i_A_csc[i+1]-i_A_csc[i];
  }
  
  delete[] i_nnIndx;
  delete[] i_A_csc;
  
}

///////////////////////////////////////////////////////////////////
//u index 
///////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
SEXP mkUIndx(SEXP n_r, SEXP m_r, SEXP nnIndx_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r, SEXP nnIndxLU_r, SEXP searchType_r){
  
  int n = INTEGER(n_r)[0];
  int m = INTEGER(m_r)[0];
  int *nnIndx = INTEGER(nnIndx_r);
  int *uIndx = INTEGER(uIndx_r);
  int *uIndxLU = INTEGER(uIndxLU_r);
  int *uiIndx = INTEGER(uiIndx_r);
  int *nnIndxLU = INTEGER(nnIndxLU_r);
  int searchType = INTEGER(searchType_r)[0];
  int i, j, k;
  //int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
  
  if(searchType == 0){
    mkUIndx0(n, m, nnIndx, uIndx, uIndxLU);
  }else if(searchType == 1){
    mkUIndx1(n, m, nnIndx, uIndx, uIndxLU);
  }else{
    mkUIndx2(n, m, nnIndx, nnIndxLU, uIndx, uIndxLU);
  }
  
  //u lists those locations that have the i-th location as a neighbor
  //then for each of those locations that have i as a neighbor, we need to know the index of i in each of their B vectors (i.e. where does i fall in their neighbor set)
  for(i = 0; i < n; i++){//for each i
    for(j = 0; j < uIndxLU[n+i]; j++){//for each location that has i as a neighbor
      k = uIndx[uIndxLU[i]+j];//index of a location that has i as a neighbor
      uiIndx[uIndxLU[i]+j] = which(i, &nnIndx[nnIndxLU[k]], nnIndxLU[n+k]);
    }
  }
  
  return R_NilValue;
}


///////////////////////////////////////////////////////////////////
//Brute force 
///////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
SEXP mkNNIndx(SEXP n_r, SEXP m_r, SEXP coords_r, SEXP nnIndx_r, SEXP nnDist_r, SEXP nnIndxLU_r, SEXP nThreads_r){
  
  int i, j, iNNIndx, iNN;
  double d;
  
  int n = INTEGER(n_r)[0];
  int m = INTEGER(m_r)[0];
  double *coords = REAL(coords_r);
  int *nnIndx = INTEGER(nnIndx_r);
  double *nnDist = REAL(nnDist_r);
  int *nnIndxLU = INTEGER(nnIndxLU_r);
  int nThreads = INTEGER(nThreads_r)[0];
  
#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif
  
  int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
  
  for(i = 0; i < nIndx; i++){
    nnDist[i] = std::numeric_limits<double>::infinity();
  }
  
#ifdef _OPENMP
#pragma omp parallel for private(j, iNNIndx, iNN, d)
#endif
  for(i = 0; i < n; i++){
    getNNIndx(i, m, iNNIndx, iNN);
    nnIndxLU[i] = iNNIndx;
    nnIndxLU[n+i] = iNN;   
    if(i != 0){  
      for(j = 0; j < i; j++){	
        d = dist2(coords[i], coords[n+i], coords[j], coords[n+j]);	
        if(d < nnDist[iNNIndx+iNN-1]){	  
          nnDist[iNNIndx+iNN-1] = d;
          nnIndx[iNNIndx+iNN-1] = j;
          rsort_with_index(&nnDist[iNNIndx], &nnIndx[iNNIndx], iNN);
        }	
      }
    }
  }
  
  return R_NilValue;
}


///////////////////////////////////////////////////////////////////
//code book
///////////////////////////////////////////////////////////////////

//Description: using the fast mean-distance-ordered nn search by Ra and Kim 1993
//Input:
//ui = is the index for which we need the m nearest neighbors
//m = number of nearest neighbors
//n = number of observations, i.e., length of u
//sIndx = the NNGP ordering index of length n that is pre-sorted by u
//u = x+y vector of coordinates assumed sorted on input
//rSIndx = vector or pointer to a vector to store the resulting nn sIndx (this is at most length m for ui >= m)
//rNNDist = vector or point to a vector to store the resulting nn Euclidean distance (this is at most length m for ui >= m) 


double dmi(double *x, double *c, int inc){
  return pow(x[0]+x[inc]-c[0]-c[inc], 2);
}

double dei(double *x, double *c, int inc){
  return pow(x[0]-c[0],2)+pow(x[inc]-c[inc],2);
}
void fastNN(int m, int n, double *coords, int ui, double *u, int *sIndx, int *rSIndx, double *rSNNDist){
  
  int i,j;
  bool up, down;
  double dm, de;
  
  //rSNNDist will hold de (i.e., squared Euclidean distance) initially.
  for(i = 0; i < m; i++){
    rSNNDist[i] = std::numeric_limits<double>::infinity();
  }
  
  i = j = ui;
  
  up = down = true;
  
  while(up || down){
    
    if(i == 0){
      down = false;
    }
    
    if(j == (n-1)){
      up = false;
    }
    
    if(down){
      
      i--;
      
      dm = dmi(&coords[sIndx[ui]], &coords[sIndx[i]], n);
      
      if(dm > 2*rSNNDist[m-1]){
        down = false;
        
      }else{
        de = dei(&coords[sIndx[ui]], &coords[sIndx[i]], n);
        
        if(de < rSNNDist[m-1] && sIndx[i] < sIndx[ui]){
          rSNNDist[m-1] = de;
          rSIndx[m-1] = sIndx[i];
          rsort_with_index(rSNNDist, rSIndx, m);
        }
        
      }
    }//end down
    
    if(up){
      
      j++;
      
      dm = dmi(&coords[sIndx[ui]], &coords[sIndx[j]], n);
      
      if(dm > 2*rSNNDist[m-1]){
        up = false;
        
      }else{
        de = dei(&coords[sIndx[ui]], &coords[sIndx[j]], n);
        
        if(de < rSNNDist[m-1] && sIndx[j] < sIndx[ui]){
          rSNNDist[m-1] = de;
          rSIndx[m-1] = sIndx[j];
          rsort_with_index(rSNNDist, rSIndx, m);
        }
        
      }
      
    }//end up
    
  }
  
  for(i = 0; i < m; i++){
    rSNNDist[i] = sqrt(rSNNDist[i]);
  }
  
  return;
}
// [[Rcpp::export]]
SEXP mkNNIndxCB(SEXP n_r, SEXP m_r, SEXP coords_r, SEXP nnIndx_r, SEXP nnDist_r, SEXP nnIndxLU_r, SEXP nThreads_r){
  
  int n = INTEGER(n_r)[0];
  int m = INTEGER(m_r)[0];
  double *coords = REAL(coords_r);
  int *nnIndx = INTEGER(nnIndx_r);
  double *nnDist = REAL(nnDist_r);
  int *nnIndxLU = INTEGER(nnIndxLU_r);
  int nThreads = INTEGER(nThreads_r)[0];
  
#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif
  
  int i, iNNIndx, iNN;
  
  int *sIndx = new int[n];
  double *u = new double[n];
  
  for(i = 0; i < n; i++){
    sIndx[i] = i;
    u[i] = coords[i]+coords[n+i];
  }
  
  rsort_with_index(u, sIndx, n); 
  
  //make nnIndxLU and fill nnIndx and d
#ifdef _OPENMP
#pragma omp parallel for private(iNNIndx, iNN)
#endif  
  for(i = 0; i < n; i++){ //note this i indexes the u vector
    getNNIndx(sIndx[i], m, iNNIndx, iNN);
    nnIndxLU[sIndx[i]] = iNNIndx;
    nnIndxLU[n+sIndx[i]] = iNN;   
    fastNN(iNN, n, coords, i, u, sIndx, &nnIndx[iNNIndx], &nnDist[iNNIndx]);
  } 
  
  return R_NilValue;
}

double spCor(double &D, double &phi, double &nu, int &covModel, double *bk){
  
  //0 exponential
  //1 spherical
  //2 matern
  //3 gaussian
  D = abs(D);
  if(covModel == 0){//exponential
    
    return exp(-D/phi);//;//0.3*pow(D/phi,1)
    
  }else if(covModel == 1){//spherical
    
    if(D > 0 && D <= phi){
      return 1.0 - 1.5*D/phi + 0.5*pow(D/phi,3);
    }else if(D >= 1.0*phi){
      return 0.0;
    }else{
      return 1.0;
    }
  }else if(covModel == 2){//matern
    
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*pi/2*(besselI(d*phi,-nu)-besselI(d*phi, nu))/sin(nu*pi), or
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)
    
    if(D > 0.0){
      return pow(D/phi, nu)/(pow(2, nu-1)*Rf_gammafn(nu))*Rf_bessel_k_ex(D/phi, nu, 1.0, bk);//thread safe bessel
    }else{
      return 1.0;
    } 
  }else if(covModel == 3){//gaussian
    
    return exp(-1.0*(pow(D, 2)/phi));
    
  }else{
    error("c++ error: cov.model is not correctly specified");
  }
}

////which index of b equals a, where b is of length n
// int which(int a, int *b, int n){
// int i;
// for(i = 0; i < n; i++){
// if(a == b[i]){
// return(i);
// }
// }

// error("c++ error: which failed");
// return -9999;
// }

//Description: computes the quadratic term.
double Q(double *B, double *F, double *u, double *v, int n, int *nnIndx, int *nnIndxLU){
  
  double a, b, q = 0;
  int i, j;
  
#ifdef _OPENMP
#pragma omp parallel for private(a, b, j) reduction(+:q)
#endif  
  for(i = 0; i < n; i++){
    a = 0;
    b = 0;
    for(j = 0; j < nnIndxLU[n+i]; j++){
      a += B[nnIndxLU[i]+j]*u[nnIndx[nnIndxLU[i]+j]];
      b += B[nnIndxLU[i]+j]*v[nnIndx[nnIndxLU[i]+j]];
    }
    q += (u[i] - a)*(v[i] - b)/F[i];
  }
  return(q);
}


// SEXP Qsq(double *B, double *F, double *y, double *X, int n, int p, int *nnIndx, int *nnIndxLU){
//   
//   double a, b, q = 0, nProtect=0;
//   int i, j, k1, k2;
//   SEXP uv_r;
//   PROTECT(uv_r = Rf_allocVector(REALSXP, p*p)); nProtect++;
//   double *uv = REAL(uv_r);
//   
// for(k1 = 0; k1 < p; k1++){ 
//   Z[k1] = Q(B, F, X[k1],y)
//   for(k2 = 0; k2 < p; k2++){ 
//     
//   }
// }
// SEXP result_r, resultName_r;
// int nResultListObjs = 1;
// 
// PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
// PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
// 
// SET_VECTOR_ELT(result_r, 0, uv_r);
// SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("uv")); 
// Rf_namesgets(result_r, resultName_r);
// 
// //unprotect
// UNPROTECT(nProtect);
// return(result_r);
// }


// void printMtrx(double *m, int nRow, int nCol){
// int i, j;
// for(i = 0; i < nRow; i++){
// Rprintf("\t");
// for(j = 0; j < nCol; j++){
// Rprintf("%.10f\t", m[j*nRow+i]);
// }
// Rprintf("\n");
// }
// }


// void printMtrxInt(int *m, int nRow, int nCol){
// int i, j;
// for(i = 0; i < nRow; i++){
// Rprintf("\t");
// for(j = 0; j < nCol; j++){
// Rprintf("%i\t", m[j*nRow+i]);
// }
// Rprintf("\n");
// }

// }


//Description: update B and F.
void updateBF(double *B, double *F, double *c, double *C, 
              double *coords, int *nnIndx, int *nnIndxLU, 
              int n, int m, double *sigmaSq, double phi, double nu, 
              int covModel, double *bk, double nuUnifb){
  
  int i, k, l;
  int info = 0;
  int inc = 1;
  double one = 1.0;
  double zero = 0.0;
  char lower = 'L';
  
  //bk must be 1+(int)floor(alpha) * nthread
  int nb = 1+static_cast<int>(floor(nuUnifb));
  int threadID = 0;
  double e;
  int mm = m*m;
  
#ifdef _OPENMP
#pragma omp parallel for private(k, l, info, threadID, e)
#endif
  for(i = 0; i < n; i++){
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif
    if(i > 0){
      for(k = 0; k < nnIndxLU[n+i]; k++){
        e = dist2(coords[i], coords[n+i], coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]]);
        c[m*threadID+k] = sigmaSq[i]*sigmaSq[nnIndx[nnIndxLU[i]+k]]*spCor(e, phi, nu, covModel, &bk[threadID*nb]);
        for(l = 0; l <= k; l++){
          e = dist2(coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]], coords[nnIndx[nnIndxLU[i]+l]], coords[n+nnIndx[nnIndxLU[i]+l]]); 
          C[mm*threadID+l*nnIndxLU[n+i]+k] = sigmaSq[nnIndx[nnIndxLU[i]+l]]*sigmaSq[nnIndx[nnIndxLU[i]+k]]*spCor(e, phi, nu, covModel, &bk[threadID*nb]); 
          // Rprintf("i: %i; k: %i; l: %i; C: %3.2f; phi = %3.2f, sigmaSq[i]: %3.2f\n", i, k, l, 
          //          C[mm*threadID+l*nnIndxLU[n+i]+k], phi, sigmaSq[i]);
          
          //sigmaSq[nnIndx[nnIndxLU[i]+l]]*sigmaSq[nnIndx[nnIndxLU[i]+k]]*
          //sigmaSq[i]*sigmaSq[nnIndx[nnIndxLU[i]+k]]*
          
        }
      }
      F77_NAME(dpotrf)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info); if(info != 0){error("11c++ error: dpotrf failed\n");}
      F77_NAME(dpotri)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info); if(info != 0){error("c++ error: dpotri failed\n");}
      F77_NAME(dsymv)(&lower, &nnIndxLU[n+i], &one, &C[mm*threadID], &nnIndxLU[n+i], &c[m*threadID], &inc, &zero, &B[nnIndxLU[i]], &inc);
      F[i] = sigmaSq[i]*sigmaSq[i] - F77_NAME(ddot)(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, &c[m*threadID], &inc); 
    }else{
      //B[i] = 0;
      F[i] = sigmaSq[i]*sigmaSq[i];
    }
    //Rprintf("-----------------------------------\n");
  }
}
//Description: update B and F.
void updateBF1(double *B, double *F, double *c, double *C, 
               double *coords, int *nnIndx, int *nnIndxLU, 
               int n, int m, double sigmaSq, double phi, double nu, 
               int covModel, double *bk, double nuUnifb){
  
  int i, k, l;
  int info = 0;
  int inc = 1;
  double one = 1.0;
  double zero = 0.0;
  char lower = 'L';
  
  //bk must be 1+(int)floor(alpha) * nthread
  int nb = 1+static_cast<int>(floor(nuUnifb));
  int threadID = 0;
  double e;
  int mm = m*m;
  
#ifdef _OPENMP
#pragma omp parallel for private(k, l, info, threadID, e)
#endif
  for(i = 0; i < n; i++){
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif
    if(i > 0){
      for(k = 0; k < nnIndxLU[n+i]; k++){
        e = dist2(coords[i], coords[n+i], coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]]);
        c[m*threadID+k] = sigmaSq*spCor(e, phi, nu, covModel, &bk[threadID*nb]);
        for(l = 0; l <= k; l++){
          e = dist2(coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]], coords[nnIndx[nnIndxLU[i]+l]], coords[n+nnIndx[nnIndxLU[i]+l]]); 
          C[mm*threadID+l*nnIndxLU[n+i]+k] = sigmaSq*spCor(e, phi, nu, covModel, &bk[threadID*nb]); 
          // Rprintf("i: %i; k: %i; l: %i; C: %3.2f; phi = %3.2f, sigmaSq[i]: %3.2f\n", i, k, l, 
          //          C[mm*threadID+l*nnIndxLU[n+i]+k], phi, sigmaSq[i]);
          
        }
      }
      F77_NAME(dpotrf)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info); if(info != 0){error("11c++ error: dpotrf failed\n");}
      F77_NAME(dpotri)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info); if(info != 0){error("c++ error: dpotri failed\n");}
      F77_NAME(dsymv)(&lower, &nnIndxLU[n+i], &one, &C[mm*threadID], &nnIndxLU[n+i], &c[m*threadID], &inc, &zero, &B[nnIndxLU[i]], &inc);
      F[i] = sigmaSq - F77_NAME(ddot)(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, &c[m*threadID], &inc); 
    }else{
      //B[i] = 0;
      F[i] = sigmaSq;
    }
    //Rprintf("-----------------------------------\n");
  }
}
// [[Rcpp::export]]
SEXP RdistC(SEXP coords_r, SEXP larCoords_r, SEXP n_r, SEXP m_r, 
			SEXP covModel_r, SEXP phi_r,  SEXP nu_r, SEXP nuUnifb_r,
			SEXP nThreads_r){
  int i, j, nProtect=0;
  int n = INTEGER(n_r)[0];
  int m = INTEGER(m_r)[0];
  int nThreads = INTEGER(nThreads_r)[0];
  int nuUnifb = INTEGER(nuUnifb_r)[0];
  double phi = REAL(phi_r)[0];
  double nu = REAL(nu_r)[0];
  double *coords = REAL(coords_r);
  double *larCoords = REAL(larCoords_r);
  int covModel = INTEGER(covModel_r)[0];
  //double *BasisIndx = (double *) R_alloc(n*p, sizeof(double));nProtect++; 
  //double *BasisIndxLU = (double *) R_alloc(n, sizeof(double));nProtect++;
  //double *H = (double *) R_alloc(n*p, sizeof(double));nProtect++;
#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif	

   //int nuUnifb = 0;
   double *bk = (double *) R_alloc((1.0+static_cast<int>(floor(nuUnifb))), sizeof(double));
   int nb = 1+static_cast<int>(floor(nuUnifb));
	  
  SEXP Dist_r, Corr_r;
  PROTECT(Dist_r = Rf_allocVector(REALSXP, n*m)); nProtect++;
  PROTECT(Corr_r = Rf_allocVector(REALSXP, n*m)); nProtect++;
  double *Dist = REAL(Dist_r);
  double *Corr = REAL(Corr_r);
  
#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
  for(i = 0; i < n; i++){
    for(j = 0; j < m; j++){
      Dist[m*i + j] = dist2(larCoords[j], larCoords[m+j], coords[i], coords[n + i]);///theta
	  
	  Corr[m*i + j] = spCor(Dist[m*i + j], phi, nu, covModel, &bk[nb]);
    }
  }
  
 
  //make return object
  SEXP result_r, resultName_r;
  int nResultListObjs = 2;
  
  PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  
  SET_VECTOR_ELT(result_r, 0, Dist_r);
  SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("Dist"));
  
  SET_VECTOR_ELT(result_r, 1, Corr_r);
  SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("Corr"));
  Rf_namesgets(result_r, resultName_r);
  
  //unprotect
  UNPROTECT(nProtect);
  return(result_r);
}

// [[Rcpp::export]]
SEXP ConExp(SEXP X_r, SEXP covZ_r, SEXP coords_r, 
            SEXP Q_r,
            SEXP n_r, SEXP p_r, SEXP Kernel_r, 
            SEXP h_r, SEXP GeomVariable_r, 
            SEXP nThreads_r){// SEXP lowB_r, SEXP invF_r,
  const int inc = 1;
  const double one = 1.0;
  const double zero = 0.0;
  // char const *lower = "L";
  // char const *ntran = "N";
  // char Trans = 'T';
  
  
  int n = INTEGER(n_r)[0];
  int GeomVariable = INTEGER(GeomVariable_r)[0];
  int s, i, j, info;
  double h = REAL(h_r)[0];
  int p = INTEGER(p_r)[0];
  int nThreads = INTEGER(nThreads_r)[0];
  int threadID = 0;
  int nProtect = 0;
  
  
  int Kernel = INTEGER(Kernel_r)[0];
  double *X = REAL(X_r);
  double *Z = REAL(covZ_r);
  double *coords = REAL(coords_r);
  double *Q = REAL(Q_r);
  // double *lowB = REAL(lowB_r);
  // double *invF = REAL(invF_r);
  
  double *Dist = (double *) R_alloc(n*nThreads, sizeof(double));
  double *K = (double *) R_alloc(n*nThreads, sizeof(double));
  double *Weight = (double *) R_alloc(n, sizeof(double));
  
#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif
  
  // char Trans = 'T';
  // char TransN = 'N';
  // double *Q_temp = (double *) R_alloc(n*n, sizeof(double));
  
  // SEXP Q_r;
  // PROTECT(Q_r = Rf_allocMatrix(REALSXP, n, n)); nProtect++;
  // double *Q = REAL(Q_r);zeros(Q, n*n);
 //// double *Q = (double *) R_alloc(n*n, sizeof(double));
  
  // F77_NAME(dgemm)(&Trans, &TransN, &n, &n, &n, &one, lowB, &n, invF, &n, &zero, Q_temp, &n);
  // F77_NAME(dgemm)(&TransN, &TransN, &n, &n, &n, &one, Q_temp, &n, lowB, &n, &zero, Q, &n);
  
  // weight
  for(i = 0; i < n; i++){
    for(j = i; j < i + 1; j++){
      Weight[i] = Q[j + n*i];
    }
  }
  
  SEXP XconZ_r;
  PROTECT(XconZ_r = Rf_allocMatrix(REALSXP, n, p)); nProtect++;
  double *XconZ = REAL(XconZ_r);
  //double *XconZ = (double *) R_alloc(n*p, sizeof(double));
  double WsumX_temp1, WsumX_temp2;
  
#ifdef _OPENMP
#pragma omp parallel for private(i, j, threadID, WsumX_temp1, WsumX_temp2)
#endif
  for(s = 0; s < n; s++){
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif
    for(j = 0; j < p; j++){
      WsumX_temp1 = 0.0;
      WsumX_temp2 = 0.0;
      for(i = 0; i < n; i++){
        if(GeomVariable == 1){
          Dist[n*threadID + i] = (Z[i]  - Z[s])/h;// pow(Z[i]  - TestZ[s], 2);
        }else{
          Dist[n*threadID + i] = dist2(coords[s], coords[n + s], coords[i], coords[n + i]);///theta  
        }
        if(Kernel == 1){
          K[n*threadID + i] = exp(-abs(Dist[n*threadID + i]))/h;  
        }else{
          K[n*threadID + i] = exp(-pow(Dist[n*threadID + i], 2)*h)/h; 
        }
        WsumX_temp1 += K[n*threadID + i]*Weight[i]*X[n*j + i];
        WsumX_temp2 += K[n*threadID + i]*Weight[i];
      }
      XconZ[n*j + s] = WsumX_temp1/WsumX_temp2; 
    }
  }
  
  // SEXP XconQ_r;
  // PROTECT(XconQ_r = Rf_allocMatrix(REALSXP, p, n)); nProtect++;
  // double *XconQ = REAL(XconQ_r);
  // 
  // double *S_temp = (double *) R_alloc(n*p, sizeof(double));
  // for(i = 0; i < n; i++){
  //   for(j = 0; j < p; j++){
  //     S_temp[n*j + i] = X[n*j + i] - XconZ[n*j + i];
  //   }
  // }
  // 
  // F77_NAME(dgemm)(&Trans, &TransN, &p, &n, &n, &one, S_temp, &n, Q, &n, &zero, XconQ, &p);
  // 
  
  
  //make return object
  SEXP result_r, resultName_r;
  int nResultListObjs = 1;
  PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  
  // SET_VECTOR_ELT(result_r, 0, XconQ_r);
  // SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("XconQ"));
  
  // SET_VECTOR_ELT(result_r, 0, Q_r);
  // SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("Q"));
  
  SET_VECTOR_ELT(result_r, 0, XconZ_r);
  SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("XconZ"));
  
  
  Rf_namesgets(result_r, resultName_r);
  
  //unprotect
  UNPROTECT(nProtect);
  return(result_r);
  
}

// [[Rcpp::export]]
SEXP local_kernel_est(SEXP y_r, SEXP covZ_r, 
                      SEXP coords_r, SEXP n_r,
                      SEXP covModel_r, SEXP h_r,
					  SEXP nu_r, SEXP nuUnifb_r, 
                      SEXP nThreads_r){
  
  const int inc = 1;
  const double one = 1.0;
  const double zero = 0.0;
  char const *lower = "L";
  char const *ntran = "N";
  char Trans = 'T';
  int nuUnifb = INTEGER(nuUnifb_r)[0];
  int n = INTEGER(n_r)[0];
  int s, i, j, info, p = 2;
  double h = REAL(h_r)[0];
  int np = n*p;
  int pp = 4;
  
  int nThreads = INTEGER(nThreads_r)[0];
  int threadID = 0;
  
  double nu = REAL(nu_r)[0];
	 int covModel = INTEGER(covModel_r)[0];
	//  int nuUnifb = 0;
	  double *bk = (double *) R_alloc(nThreads*(1.0+static_cast<int>(floor(nuUnifb))), sizeof(double));
	  int nb = 1+static_cast<int>(floor(nuUnifb));
	  
  //int Kernel = INTEGER(Kernel_r)[0];
  double *y = REAL(y_r);
  double *Z = REAL(covZ_r);
  double *coords = REAL(coords_r);
  
 double *Dist = (double *) R_alloc(n*nThreads, sizeof(double));
 double *K = (double *) R_alloc(n*nThreads, sizeof(double));
  double *X = (double *) R_alloc(p*n*nThreads, sizeof(double));
  double *KerX = (double *) R_alloc(p*n*nThreads, sizeof(double));
   double *tKerX = (double *) R_alloc(p*n*nThreads, sizeof(double));
  double *XtX = (double *) R_alloc(p*p*nThreads, sizeof(double));
   double *XtX_temp = (double *) R_alloc(p*p*nThreads, sizeof(double));
  
  //	double *KerY = (double *) R_alloc(n, sizeof(double));
  
  //int Kernel = 2;
  
  double *tmp_p = (double *) R_alloc(p*nThreads, sizeof(double));
  double d1, t;
  
  int nProtect=0;
  SEXP alpha_r, S_r;
  PROTECT(alpha_r = Rf_allocMatrix(REALSXP, p, n)); nProtect++; //Rf_allocVector
  //double *alpha = REAL(alpha_r);
  PROTECT(S_r = Rf_allocMatrix(REALSXP, p, n*n)); nProtect++; 
  
  
  //double *c =(double *) R_alloc(m[0]*nThreads, sizeof(double));
  double *alpha = (double *) R_alloc(p*nThreads, sizeof(double));
 // double *S_temp = (double *) R_alloc(n*p*n*nThreads, sizeof(double));
 
 
#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif
  
  //private(i, j, Dist, K, X, KerX, XtX, alpha_temp, tmp_p)
  
#ifdef _OPENMP
#pragma omp parallel for private(i, j, d1, t, threadID)
#endif
  for(s = 0; s < n; s++){
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif

    for(i = 0; i < n; i++){
        Dist[n*threadID + i] = (Z[i]  - Z[s])/h;// pow(Z[i]  - TestZ[s], 2);
		d1 = (Z[i]  - Z[s]);
		 // Rprintf("d1: %3.2f\n", d1);
      // if(Kernel == 1){
        // K[n*threadID + i] = exp(-abs(Dist[n*threadID + i]))/h;  
      // }else{
       K[n*threadID + i] = spCor(d1, h, nu, covModel, &bk[threadID*nb])/h;//exp(-pow(Dist[n*threadID + i], 2)*h)/h; 
		  // Rprintf("C++: %3.2f\n", spCor(d1, h, nu, covModel, &bk[threadID*nb])/h);
		  // Rprintf("K: %3.2f\n", exp(-pow(Dist[n*threadID + i], 2)*h)/h);
      // }
    }
    
    for(i = 0; i < n; i++){
      for(j = 0; j < p; j++){
        if(j==0){
          X[p*n*threadID + n*j + i] = 1;
          KerX[p*n*threadID + n*j + i] = K[n*threadID + i];
        }else{
          X[p*n*threadID + n*j + i] = Dist[n*threadID + i];  
          KerX[p*n*threadID + n*j + i] = K[n*threadID + i]*Dist[n*threadID + i];
        }
      }
      
    }
	
   
     for(j = 0; j < p; j++){
	   for(i = 0; i < n; i++){
		   if(j == 0){
			 tKerX[p*n*threadID + p*i + j] = K[n*threadID + i];  
		   }else{
			 tKerX[p*n*threadID + p*i + j] = K[n*threadID + i]*Dist[n*threadID + i];  
		   }		  
	    }
	   }
	
	 // for(i = 0; i < n; i++){
		 // d1 = (Z[i]  - Z[s]);
		 // t = spCor(d1, h, nu, covModel, &bk[threadID*nb])/h;
      // for(j = 0; j < p; j++){
        // if(j==0){
          // X[p*n*threadID + n*j + i] = 1;
          // KerX[p*n*threadID + n*j + i] = t;
        // }else{
          // X[p*n*threadID + n*j + i] = d1/h;  
          // KerX[p*n*threadID + n*j + i] = d1*t/h;
        // }
      // }
      
    // }
	
	   // for(j = 0; j < p; j++){
	   // for(i = 0; i < n; i++){
		   // d1 = (Z[i]  - Z[s]);
		   // t = spCor(d1, h, nu, covModel, &bk[threadID*nb])/h;
		   // if(j == 0){
			 // tKerX[p*n*threadID + p*i + j] = t;  
		   // }else{
			 // tKerX[p*n*threadID + p*i + j] = d1*t/h;  
		   // }		  
	    // }
	   // }
	
	
    //F77_NAME(dgemv)(&Trans, &n, &p, &one, &KerX[p*n*threadID], &n, y, &inc, &zero, &tmp_p[p*threadID], &inc);
    F77_NAME(dgemm)(&Trans, ntran, &p, &p, &n, &one, &KerX[p*n*threadID], &n, &X[p*n*threadID], &n, &zero, &XtX[p*p*threadID], &p);
    
	F77_NAME(dcopy)(&pp, &XtX[p*p*threadID], &inc, &XtX_temp[p*p*threadID], &inc);
	
	F77_NAME(dposv)(lower, &p, &n, &XtX_temp[p*p*threadID], &p, &tKerX[p*n*threadID], &p, &info);//if(info != 0){error("c++ error: dpotrf failed\n");}
	F77_NAME(dcopy)(&np, &tKerX[p*n*threadID], &inc, &REAL(S_r)[s*p*n], &inc);
	F77_NAME(dgemv)(ntran, &p, &n, &one, &tKerX[p*n*threadID], &p, y, &inc, &zero, &alpha[p*threadID], &inc);	
	
	
    // F77_NAME(dpotrf)(lower, &p, &XtX[p*p*threadID], &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
    // F77_NAME(dpotri)(lower, &p, &XtX[p*p*threadID], &p, &info); if(info != 0){error("c++ error: dpotri failed\n");}//solve(A)
    // F77_NAME(dsymv)(lower, &p, &one, &XtX[p*p*threadID], &p, &tmp_p[p*threadID], &inc, &zero, &alpha[p*threadID], &inc);
    
	
    F77_NAME(dcopy)(&p, &alpha[p*threadID], &inc, &REAL(alpha_r)[s*p], &inc);
	
	
	
	
	// F77_NAME(dgetrf)(&p, &p, &XtX_temp[p*p*threadID], &p, &p, &info);
	// F77_NAME(dgetrs)(ntran, &p, &n, &XtX_temp[p*p*threadID], &p, &p, &tKerX[p*n*threadID], &n, &info);
	
	// F77_NAME(dpotrf)(lower, &p, &XtX_temp[p*p*threadID], &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
    // F77_NAME(dpotri)(lower, &p, &XtX_temp[p*p*threadID], &p, &info); if(info != 0){error("c++ error: dpotri failed\n");}
	
	 // F77_NAME(dpotrf)(lower, &p, &XtX_temp[p*p*threadID], &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
    // F77_NAME(dpotrs)(lower, &p, &n, &XtX_temp[p*p*threadID], &p, &tKerX[p*n*threadID], &p, &info); if(info != 0){error("c++ error: dpotri failed\n");}
	
	
	
    // alpha[s] = alpha_temp[n*p*threadID];
  }
  
  //make return object
  SEXP result_r, resultName_r;
  int nResultListObjs = 2;
  PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  
  SET_VECTOR_ELT(result_r, 0, alpha_r);
  SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("alpha"));
  
  SET_VECTOR_ELT(result_r, 1, S_r);
  SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("S"));
  Rf_namesgets(result_r, resultName_r);
  
  //unprotect
  UNPROTECT(nProtect);
  // SEXP out;
  // PROTECT(out = Rf_allocVector(REALSXP, nIndx));
  // B = REAL(out);
  
  return(result_r);
  
}

// [[Rcpp::export]]
SEXP local_kernel_pred(SEXP y_r, 
                       SEXP covZ_r, SEXP TestCovZ_r, 
                       SEXP coords_r, SEXP TestCoords_r,
                       SEXP n_r, SEXP nTest_r,
                       SEXP covModel_r, SEXP h_r,
					   SEXP nu_r, SEXP nuUnifb_r, 
                       SEXP nThreads_r){
  const int inc = 1;
  const double one = 1.0;
  const double zero = 0.0;
  char const *lower = "L";
  char const *ntran = "N";
  char Trans = 'T';
  int nuUnifb = INTEGER(nuUnifb_r)[0];
  int n = INTEGER(n_r)[0];
  int nTest = INTEGER(nTest_r)[0];
  
  // int GeomVariable = INTEGER(GeomVariable_r)[0];
  int s, i, j, info, p = 2;
  double h = REAL(h_r)[0];
  int np = n*p;
  int pp = p*p;
  
  int nThreads = INTEGER(nThreads_r)[0];
  int threadID = 0;
  
    double nu = REAL(nu_r)[0];
	 int covModel = INTEGER(covModel_r)[0];
	 // int nuUnifb = 0;
	  double *bk = (double *) R_alloc(nThreads*(1.0+static_cast<int>(floor(nuUnifb))), sizeof(double));
	  int nb = 1+static_cast<int>(floor(nuUnifb));
	  
  //int Kernel = INTEGER(Kernel_r)[0];
  double *y = REAL(y_r);
  double *Z = REAL(covZ_r);
  double *TestCovZ = REAL(TestCovZ_r);
  double *coords = REAL(coords_r);
  double *TestCoords = REAL(TestCoords_r);

  double *Dist = (double *) R_alloc(n*nThreads, sizeof(double));
  double *K = (double *) R_alloc(n*nThreads, sizeof(double));
  double *X = (double *) R_alloc(p*n*nThreads, sizeof(double));
  double *KerX = (double *) R_alloc(p*n*nThreads, sizeof(double));
  double *XtX = (double *) R_alloc(p*p*nThreads, sizeof(double));
  
  double *tKerX = (double *) R_alloc(p*n*nThreads, sizeof(double));
  double *XtX_temp = (double *) R_alloc(p*p*nThreads, sizeof(double));
  
  
  
  //	double *KerY = (double *) R_alloc(n, sizeof(double));
  //double *tmp_p = (double *) R_alloc(p*nThreads, sizeof(double));
  int nProtect=0;
  SEXP alpha_r, S_r;
  PROTECT(alpha_r = Rf_allocMatrix(REALSXP, p, nTest)); nProtect++;
  PROTECT(S_r = Rf_allocMatrix(REALSXP, p, n*nTest)); nProtect++; 
   
  
 
  //double *c =(double *) R_alloc(m[0]*nThreads, sizeof(double));
  double *alpha = (double *) R_alloc(p*nThreads, sizeof(double));
  double *S_temp = (double *) R_alloc(n*p*n*nThreads, sizeof(double));
  double d1, t;
  //int Kernel = 2;
  
  
#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif
  
  //private(i, j, Dist, K, X, KerX, XtX, alpha_temp, tmp_p)
  
#ifdef _OPENMP
#pragma omp parallel for private(i, j, d1, t, threadID)
#endif
  for(s = 0; s < nTest; s++){
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif
    for(i = 0; i < n; i++){
        Dist[n*threadID + i] = (Z[i]  - TestCovZ[s])/h;// pow(Z[i]  - TestZ[s], 2);
		d1 = (Z[i]  - TestCovZ[s]);
      // if(Kernel == 1){
        // K[n*threadID + i] = exp(-abs(Dist[n*threadID + i]))/h;  
      // }else{
        // K[n*threadID + i] = spCor(d1, h, nu, covModel, &bk[threadID*nb])/pow(h, 1);//exp(-pow(Dist[n*threadID + i], 2)*h)/h; 
      // }
	   K[n*threadID + i] = spCor(d1, h, nu, covModel, &bk[threadID*nb])/h;
    }
    
    for(i = 0; i < n; i++){
      for(j = 0; j < p; j++){
        if(j==0){
          X[p*n*threadID + n*j + i] = 1;
          KerX[p*n*threadID + n*j + i] = K[n*threadID + i];
        }else{
          X[p*n*threadID + n*j + i] = Dist[n*threadID + i];  
          KerX[p*n*threadID + n*j + i] = K[n*threadID + i]*Dist[n*threadID + i]; 
        }
      }
      
    }
    
       for(j = 0; j < p; j++){
	   for(i = 0; i < n; i++){
		   if(j == 0){
			 tKerX[p*n*threadID + p*i + j] = K[n*threadID + i];  
		   }else{
			 tKerX[p*n*threadID + p*i + j] = K[n*threadID + i]*Dist[n*threadID + i];  
		   }		  
	    }
	   }
	   
	   // for(i = 0; i < n; i++){
		    // d1 = (Z[i]  - TestCovZ[s])/h;
			// t = spCor(d1, h, nu, covModel, &bk[threadID*nb])/pow(h, 1);
      // for(j = 0; j < p; j++){
        // if(j==0){
          // X[p*n*threadID + n*j + i] = 1;
          // KerX[p*n*threadID + n*j + i] = t;
        // }else{
          // X[p*n*threadID + n*j + i] = d1;  
          // KerX[p*n*threadID + n*j + i] = d1*t; 
        // }
      // }
      
    // }  
	   
	     // for(j = 0; j < p; j++){
	   // for(i = 0; i < n; i++){
		   // d1 = (Z[i]  - TestCovZ[s])/h;
		   // t = spCor(d1, h, nu, covModel, &bk[threadID*nb])/pow(h, 1);
		   // if(j == 0){
			 // tKerX[p*n*threadID + p*i + j] = t;  
		   // }else{
			 // tKerX[p*n*threadID + p*i + j] = d1*t;  
		   // }		  
	    // }
	   // } 
	   
	   
	   
   // F77_NAME(dgemv)(&Trans, &n, &p, &one, &KerX[p*n*threadID], &n, y, &inc, &zero, &tmp_p[p*threadID], &inc);
    F77_NAME(dgemm)(&Trans, ntran, &p, &p, &n, &one, &KerX[p*n*threadID], &n, &X[p*n*threadID], &n, &zero, &XtX[p*p*threadID], &p);
	
    F77_NAME(dcopy)(&pp, &XtX[p*p*threadID], &inc, &XtX_temp[p*p*threadID], &inc);
	// F77_NAME(dgetrf)(&p, &p, &XtX_temp[p*p*threadID], &p, &p, &info);
	// F77_NAME(dgetrs)(ntran, &p, &n, &XtX_temp[p*p*threadID], &p, &p, &tKerX[p*n*threadID], &n, &info);
	F77_NAME(dposv)(lower, &p, &n, &XtX_temp[p*p*threadID], &p, &tKerX[p*n*threadID], &p, &info);//if(info != 0){error("c++ error: dpotrf failed\n");}
	F77_NAME(dcopy)(&np, &tKerX[p*n*threadID], &inc, &REAL(S_r)[s*p*n], &inc);
	
	F77_NAME(dgemv)(ntran, &p, &n, &one, &tKerX[p*n*threadID], &p, y, &inc, &zero, &alpha[p*threadID], &inc);	
	F77_NAME(dcopy)(&p, &alpha[p*threadID], &inc, &REAL(alpha_r)[s*p], &inc);
	
    // F77_NAME(dpotrf)(lower, &p, &XtX[p*p*threadID], &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
    // F77_NAME(dpotri)(lower, &p, &XtX[p*p*threadID], &p, &info); if(info != 0){error("c++ error: dpotri failed\n");}//solve(A)
    // F77_NAME(dsymv)(lower, &p, &one, &XtX[p*p*threadID], &p, &tmp_p[p*threadID], &inc, &zero, &alpha[p*threadID], &inc);
    
   
    //alpha[s] = alpha_temp[n*p*threadID];
	
	
  }
  
  //make return object
  SEXP result_r, resultName_r;
  int nResultListObjs = 2;
  PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  
  SET_VECTOR_ELT(result_r, 0, alpha_r);
  SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("alpha"));
    SET_VECTOR_ELT(result_r, 1, S_r);
  SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("S"));
  Rf_namesgets(result_r, resultName_r);
  
  //unprotect
  UNPROTECT(nProtect);
  // SEXP out;
  // PROTECT(out = Rf_allocVector(REALSXP, nIndx));
  // B = REAL(out);
  
  return(result_r);
  
}



// [[Rcpp::export]]
SEXP spatial_LLE(SEXP y_r, SEXP X_r, 
				  SEXP coord_r, SEXP fsDensity_r, SEXP n_r,
				  SEXP p_r, SEXP covModel_r, SEXP h_r, 
				  SEXP nu_r, SEXP nuUnifb_r,
				  SEXP adWidth_r, SEXP mm_r,
				  SEXP nThreads_r){
	double *y = REAL(y_r);
    double *X = REAL(X_r);
	double *coords = REAL(coord_r);
	double *fsDensity = REAL(fsDensity_r);
	int *adWidth = INTEGER(adWidth_r);
	int n = INTEGER(n_r)[0];
	int mm = INTEGER(mm_r)[0];
    int p = INTEGER(p_r)[0];	
	//int Kernel_index = INTEGER(Kernel_r)[0];
	double h = REAL(h_r)[0];
	int nThreads = INTEGER(nThreads_r)[0];
	//int Kernel = INTEGER(Kernel_r)[0];
	 int nuUnifb = INTEGER(nuUnifb_r)[0];
	double nu = REAL(nu_r)[0];
	 int covModel = INTEGER(covModel_r)[0];
	 // int nuUnifb = 0;
	  double *bk = (double *) R_alloc(nThreads*(1.0+static_cast<int>(floor(nuUnifb))), sizeof(double));
	  int nb = 1+static_cast<int>(floor(nuUnifb));
	
	
	int s, i, j, k, info, Ps = 3*p*n;
	int p3 = 3*p;
	int p3s = p3*p3;
	int nProtect=0;
	const int inc = 1;
    const double one = 1.0;
    const double zero = 0.0;
    char const *lower = "L";
	char const *ntran = "N";
	char Trans = 'T';
	int threadID = 0;
    SEXP S_r, alpha_r, Mat_r, KK_r;
    PROTECT(S_r = Rf_allocMatrix(REALSXP, p3, n*n)); nProtect++;
	//PROTECT(Xs_r = Rf_allocMatrix(REALSXP, n, Ps)); nProtect++;
	//double *Xs = (double *) R_alloc(3*n*p*nThreads, sizeof(double));zeros(Xs, 3*n*p*nThreads);
	PROTECT(alpha_r = Rf_allocMatrix(REALSXP, p3, n)); nProtect++;
	//double *Xs = REAL(Xs_r);
	 //PROTECT(Mat_r = Rf_allocMatrix(REALSXP, p3, p3*n)); nProtect++;
	
	
	//PROTECT(KK_r = Rf_allocMatrix(REALSXP, n, p3*n)); nProtect++;
	
	// SEXP info_r;
    // PROTECT(info_r = Rf_allocVector(INTSXP, n)); nProtect++;
     // int *info = INTEGER(info_r);
	 //int *info = (int *) R_alloc(n*nThreads, sizeof(int));
	
	
	 //double *Dist = (double *) R_alloc(n*nThreads, sizeof(double));zeros(Dist, n*nThreads);
    // double *K = (double *) R_alloc(n*nThreads, sizeof(double));zeros(K, n*nThreads);
    //double *KerY = (double *) R_alloc(n*nThreads, sizeof(double));
    double *DagKer = (double *) R_alloc(n*n*nThreads, sizeof(double)); zeros(DagKer, n*n*nThreads);
  
  
  double *alpha = (double *) R_alloc(p3*nThreads, sizeof(double));
  double *Xs = (double *) R_alloc(Ps*nThreads, sizeof(double));zeros(Xs, Ps*nThreads);
  double *Q_temp = (double *) R_alloc(Ps*nThreads, sizeof(double)); zeros(Q_temp, Ps*nThreads);
  double *S = (double *) R_alloc(Ps*nThreads, sizeof(double)); zeros(S, Ps*nThreads);
  double *Mat = (double *) R_alloc(p3*p3*nThreads, sizeof(double)); zeros(Mat, p3*p3*nThreads);
  double d0, d1, t; 
 
  
#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif
  
  
#ifdef _OPENMP
#pragma omp parallel for private(i, j, k, d0, d1, t, threadID)
#endif
  for(s = 0; s < n; s++){
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif
         if(mm > 0){
			d0 = h*dist2(coords[s], coords[n + s], coords[adWidth[(mm - 1)*n + s]], coords[n + adWidth[(mm - 1)*n + s]]); 
		 }else{
			d0 = h;
		 }
		 for(i = 0; i < n; i++){
			   for(j = 0; j < p; j++){
					Xs[3*n*p*threadID + j*3*n + i] = X[j*n + i];
					Xs[3*n*p*threadID + (j*3 + 1)*n + i] = X[j*n + i]*(coords[i] - coords[s])/d0;
					Xs[3*n*p*threadID + (j*3 + 2)*n + i] = X[j*n + i]*(coords[n + i] - coords[n + s])/d0;

			} 
		 }
		 if(mm > 0){
			for(k = 0; k < n; k++){
				// for(j = k; j < k + 1; j++){					
					              // DagKer[n*n*threadID + k + n*k] = 0.0;	
				                // }
								 DagKer[n*n*threadID + k + n*k] = 0.0;	
				             }
			    for(i = 0; i < mm; i++){
						  d1 = dist2(coords[s], coords[n + s], coords[adWidth[i*n + s]], coords[n + adWidth[i*n + s]])*fsDensity[adWidth[i*n + s]];
				           //t = fsDensity[adWidth[i*n + s]]*pow(d1/d0, 2);
						   // t = fsDensity[adWidth[i*n + s]]*exp(-pow(d1, 1)/d0);
						   t = fsDensity[adWidth[i*n + s]]*spCor(d1, d0, nu, covModel, &bk[threadID*nb]);
				         if(t<1){DagKer[n*n*threadID + adWidth[i*n + s] + n*adWidth[i*n + s]] = t/pow(d0, 2);}//(1-t)
                     // Rprintf("adWidth[i*n + s] = %i ; K: %3.10f \n", adWidth[i*n + s], DagKer[n*n*threadID + adWidth[i*n + s] + n*adWidth[i*n + s]]);						 
				        }	
				
		}else{
			for(i = 0; i < n; i++){
			  d1 = fsDensity[i]*dist2(coords[s], coords[n + s], coords[i], coords[n + i]);///theta 

              t = fsDensity[i]*spCor(d1, h, nu, covModel, &bk[threadID*nb])/pow(h, 2);			  
			  // if(Kernel == 1){
				// t = fsDensity[i]*exp(-d1/h)/pow(h, 2);
				
				
				// int d;
				// d = 1 - pow(Dist[n*threadID + i], 2);
				// if(d>0){K[n*threadID + i] = 0.75*d/pow(h, 2);}else{K[n*threadID + i] = 0;}
				
			  // }else{
				// t = fsDensity[i]*exp(-pow(d1, 2)/h)/pow(h, 2); 
			  // }
			  //kernel function
					
			   DagKer[n*n*threadID + i + n*i] = t;		
				
			}
	}

     F77_NAME(dgemm)(&Trans, ntran, &p3, &n, &n, &one, &Xs[Ps*threadID], &n, &DagKer[n*n*threadID], &n, &zero, &S[Ps*threadID],&p3);	
     F77_NAME(dgemm)(ntran, ntran, &p3, &p3, &n, &one, &S[Ps*threadID], &p3, &Xs[Ps*threadID], &n, &zero, &Mat[p3*p3*threadID],&p3);
	 
	 
	 F77_NAME(dposv)(lower, &p3, &n, &Mat[p3*p3*threadID], &p3, &S[Ps*threadID], &p3, &info); //if(info != 0){error("c++ error: dpotrf failed\n");}
	 F77_NAME(dcopy)(&Ps, &S[Ps*threadID], &inc, &REAL(S_r)[s*Ps], &inc);
	 F77_NAME(dgemv)(ntran, &p3, &n, &one, &S[Ps*threadID], &p3, y, &inc, &zero, &alpha[p3*threadID], &inc);	
	 F77_NAME(dcopy)(&p3, &alpha[p3*threadID], &inc, &REAL(alpha_r)[s*p3], &inc);
	 
	 
	 
	 
	  //F77_NAME(dcopy)(&p3s, &Mat[p3*p3*threadID], &inc, &REAL(Mat_r)[s*p3s], &inc);	
	 //if(info != 0){error("c++ error: dpotrf failed\n");}
	 //F77_NAME(dcopy)(&inc, &info[s], &inc, &INTEGER(info_r)[s], &inc);

	
	 //F77_NAME(dgemm)(ntran, ntran, &p3, &n, &n, &one, &S[Ps*threadID], &p3, &DagKer[n*n*threadID], &n, &zero, &S[Ps*threadID],&p3);	
     

	
		//F77_NAME(dcopy)(&Ps, &Xs[3*n*p*threadID], &inc, &REAL(Xs_r)[s*Ps], &inc);
	 }
	//make return object
  SEXP result_r, resultName_r;
  int nResultListObjs = 2;
  PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  
  SET_VECTOR_ELT(result_r, 0, alpha_r);
  SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("alpha"));
  
  SET_VECTOR_ELT(result_r, 1, S_r);
  SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("S"));
  
   // SET_VECTOR_ELT(result_r, 2, KK_r);
   // SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("Mat"));
  
  Rf_namesgets(result_r, resultName_r);
  
  //unprotect
  UNPROTECT(nProtect);
  // SEXP out;
  // PROTECT(out = Rf_allocVector(REALSXP, nIndx));
  // B = REAL(out);
  
  return(result_r);
						  
}


// [[Rcpp::export]]
SEXP spatial_LLE_Pred(SEXP y_r, SEXP X_r, 
					  SEXP coord_r, 
					  SEXP fsDensity_r,
					  SEXP n_r,
					  SEXP pred_coord_r, SEXP m_r,
					  SEXP p_r, SEXP covModel_r, 
					  SEXP h_r, SEXP nu_r, 
					  SEXP nuUnifb_r, 
					  SEXP adWidth_r, 
					  SEXP mm_r,
					  SEXP nThreads_r){
	double *y = REAL(y_r);
    double *X = REAL(X_r);
	double *coords = REAL(coord_r);
	double *fsDensity = REAL(fsDensity_r);
	double *pred_coord = REAL(pred_coord_r);
	int *adWidth = INTEGER(adWidth_r);
	int n = INTEGER(n_r)[0];
	int m = INTEGER(m_r)[0];
	int mm = INTEGER(mm_r)[0];
    int p = INTEGER(p_r)[0];	
	//int Kernel_index = INTEGER(Kernel_r)[0];
	double h = REAL(h_r)[0];
	int nThreads = INTEGER(nThreads_r)[0];
	//int Kernel = INTEGER(Kernel_r)[0];
	int nuUnifb = INTEGER(nuUnifb_r)[0];
	double nu = REAL(nu_r)[0];
	 int covModel = INTEGER(covModel_r)[0];
	//  int nuUnifb = 0;
	  double *bk = (double *) R_alloc(nThreads*(1.0+static_cast<int>(floor(nuUnifb))), sizeof(double));
	  int nb = 1+static_cast<int>(floor(nuUnifb));
	
	
	
	
	int s, i, j, k, info, Ps = 3*p*n;
	int p3 = 3*p;
	int p3s = p3*p3;
	int nProtect=0;
	const int inc = 1;
    const double one = 1.0;
    const double zero = 0.0;
    char const *lower = "L";
	char const *ntran = "N";
	char Trans = 'T';
	int threadID = 0;
    SEXP S_r, alpha_r;
    PROTECT(S_r = Rf_allocMatrix(REALSXP, p3, n*m)); nProtect++;
	//PROTECT(Xs_r = Rf_allocMatrix(REALSXP, n, Ps)); nProtect++;
	//double *Xs = (double *) R_alloc(3*n*p*nThreads, sizeof(double));zeros(Xs, 3*n*p*nThreads);
	PROTECT(alpha_r = Rf_allocMatrix(REALSXP, p3, m)); nProtect++;
	//double *Xs = REAL(Xs_r);
	 //PROTECT(Mat_r = Rf_allocMatrix(REALSXP, p3, p3*n)); nProtect++;
	
	
	
	// double *Dist = (double *) R_alloc(n*nThreads, sizeof(double));zeros(Dist, n*nThreads);
    // double *K = (double *) R_alloc(n*nThreads, sizeof(double));zeros(K, n*nThreads);
    //double *KerY = (double *) R_alloc(n*nThreads, sizeof(double));
    double *DagKer = (double *) R_alloc(n*n*nThreads, sizeof(double)); zeros(DagKer, n*n*nThreads);
  
  
  double *alpha = (double *) R_alloc(p3*nThreads, sizeof(double));
  double *Xs = (double *) R_alloc(Ps*nThreads, sizeof(double));zeros(Xs, Ps*nThreads);
  double *Q_temp = (double *) R_alloc(Ps*nThreads, sizeof(double)); zeros(Q_temp, Ps*nThreads);
  double *S = (double *) R_alloc(Ps*nThreads, sizeof(double)); zeros(S, Ps*nThreads);
  double *Mat = (double *) R_alloc(p3*p3*nThreads, sizeof(double)); zeros(Mat, p3*p3*nThreads);
  double d0, d1, t;
  
  
#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif
  
  
#ifdef _OPENMP
#pragma omp parallel for private(i, j, d0, d1, t, threadID)
#endif
  for(s = 0; s < m; s++){
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif

        if(mm > 0){
			d0 = h*dist2(pred_coord[s], pred_coord[m + s], coords[adWidth[(mm - 1)*m + s]], coords[n + adWidth[(mm - 1)*m + s]]);
		 }else{
			d0 = h;
		 }
		 for(i = 0; i < n; i++){
			   for(j = 0; j < p; j++){
					Xs[3*n*p*threadID + j*3*n + i] = X[j*n + i];
					Xs[3*n*p*threadID + (j*3 + 1)*n + i] = X[j*n + i]*(coords[i] - pred_coord[s])/d0;
					Xs[3*n*p*threadID + (j*3 + 2)*n + i] = X[j*n + i]*(coords[n + i] - pred_coord[m + s])/d0;

			} 
		 }
		  if(mm > 0){
		      for(k = 0; k < n; k++){
				// for(j = k; j < k + 1; j++){					
					              // DagKer[n*n*threadID + j + n*k] = 0.0;	
				                // }
								 DagKer[n*n*threadID + k + n*k] = 0.0;	
				             }
			    for(i = 0; i < mm; i++){
						   d1 = dist2(pred_coord[s], pred_coord[m + s], coords[adWidth[i*m + s]], coords[n + adWidth[i*m + s]])*fsDensity[adWidth[i*m + s]];
				          // t = (pow(d1/d0, 2))*fsDensity[adWidth[i*m + s]];
						  t = fsDensity[adWidth[i*m + s]]*spCor(d1, d0, nu, covModel, &bk[threadID*nb]);
				         if(t < 1){DagKer[n*n*threadID + adWidth[i*m + s] + n*adWidth[i*m + s]] = (t)/pow(d0, 2);//(1-t)
						 }
				
				        }
		   
		}else{
			for(i = 0; i < n; i++){
			  d1 = dist2(pred_coord[s], pred_coord[m + s], coords[i], coords[n + i]);///theta  
			  t = fsDensity[i]*spCor(d1, h, nu, covModel, &bk[threadID*nb])/pow(h, 2);	
			  // if(Kernel == 1){
				// t = fsDensity[i]*exp(-abs(d1)/h)/pow(h, 2); 
			  // }else{ 
				// t = fsDensity[i]*exp(-pow(d1, 2)/h)/pow(h, 2); 
			  // }
			    DagKer[n*n*threadID + i + n*i] = t;		
			}
		}

     F77_NAME(dgemm)(&Trans, ntran, &p3, &n, &n, &one, &Xs[Ps*threadID], &n, &DagKer[n*n*threadID], &n, &zero, &S[Ps*threadID],&p3);	
     F77_NAME(dgemm)(ntran, ntran, &p3, &p3, &n, &one, &S[Ps*threadID], &p3, &Xs[Ps*threadID], &n, &zero, &Mat[p3*p3*threadID],&p3);		
		
	 //F77_NAME(dcopy)(&p3s, &Mat[p3*p3*threadID], &inc, &REAL(Mat_r)[s*p3s], &inc);	
		
	 F77_NAME(dposv)(lower, &p3, &n, &Mat[p3*p3*threadID], &p3, &S[Ps*threadID], &p3, &info);//if(info != 0){error("c++ error: dpotrf failed\n");}
		
	 //F77_NAME(dgemm)(ntran, ntran, &p3, &n, &n, &one, &S[Ps*threadID], &p3, &DagKer[n*n*threadID], &n, &zero, &S[Ps*threadID],&p3);	
     F77_NAME(dcopy)(&Ps, &S[Ps*threadID], &inc, &REAL(S_r)[s*Ps], &inc);
	
	 F77_NAME(dgemv)(ntran, &p3, &n, &one, &S[Ps*threadID], &p3, y, &inc, &zero, &alpha[p3*threadID], &inc);	
	 F77_NAME(dcopy)(&p3, &alpha[p3*threadID], &inc, &REAL(alpha_r)[s*p3], &inc);

	
		//F77_NAME(dcopy)(&Ps, &Xs[3*n*p*threadID], &inc, &REAL(Xs_r)[s*Ps], &inc);
	 }
	//make return object
  SEXP result_r, resultName_r;
  int nResultListObjs = 2;
  PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  
  SET_VECTOR_ELT(result_r, 0, alpha_r);
  SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("alpha"));
  
  SET_VECTOR_ELT(result_r, 1, S_r);
  SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("S"));
  Rf_namesgets(result_r, resultName_r);
  
  //unprotect
  UNPROTECT(nProtect);
  // SEXP out;
  // PROTECT(out = Rf_allocVector(REALSXP, nIndx));
  // B = REAL(out);
  
  return(result_r);
						  
}

// [[Rcpp::export]]
SEXP semiQLME(SEXP y_r, SEXP n_r, SEXP m_r, SEXP N_r, SEXP coords_r, SEXP covModel_r,
              SEXP nnIndx_r, SEXP nnIndxLU_r, SEXP sigmaSq_r, SEXP phi_r, 
			  SEXP nu_r, SEXP nThreads_r){
			   // SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r, SEXP TestCoords_r,
  int i, j, k, l, s, info, nProtect=0;
  const int inc = 1;
  const double one = 1.0;
  const double negOne = -1.0;
  const double zero = 0.0;
 
  
  //get args
  int n = INTEGER(n_r)[0];
  int m = INTEGER(m_r)[0];
  int N = INTEGER(N_r)[0];

  double *coords = REAL(coords_r);
  double *y = REAL(y_r);
  
  int *nnIndx = INTEGER(nnIndx_r);
  int *nnIndxLU = INTEGER(nnIndxLU_r);
  // int *uIndx = INTEGER(uIndx_r);
  // int *uIndxLU = INTEGER(uIndxLU_r);
  // int *uiIndx = INTEGER(uiIndx_r);

  
  int covModel = INTEGER(covModel_r)[0];
  std::string corName = getCorName(covModel);
  

  int nThreads = INTEGER(nThreads_r)[0];

  char Trans = 'T';
  char TransN = 'N';
  char lower = 'L';
  
#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif
  
  
  //starting	

  double *sigmaSq = REAL(sigmaSq_r);
  double *phi = REAL(phi_r);

  //allocate for the U index vector that keep track of which locations have the i-th location as a neighbor
  int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);

  //other stuff
  int mm = m*m;
  
  int nn = n*n;
  // SEXP neiB_r, varF_r;
  // PROTECT(neiB_r = Rf_allocVector(REALSXP, nIndx)); nProtect++;
  // PROTECT(varF_r = Rf_allocVector(REALSXP, n)); nProtect++;
  // double *B = REAL(neiB_r);
  // double *F = REAL(varF_r);
  double *B = (double *) R_alloc(nIndx, sizeof(double));
  double *F = (double *) R_alloc(n, sizeof(double));

  double *c =(double *) R_alloc(m*nThreads, sizeof(double));
  double *C = (double *) R_alloc(mm*nThreads, sizeof(double));
  int nuUnifb = 0;
  double *bk = (double *) R_alloc(nThreads*(1.0+static_cast<int>(floor(nuUnifb))), sizeof(double));

  double nu = 0;
  if(corName == "matern"){nu = REAL(nu_r)[0];}
  
  double sq_err, qle;
  SEXP Q_r, qleSum_r;
  PROTECT(Q_r = Rf_allocMatrix(REALSXP, n, n)); nProtect++;
  PROTECT(qleSum_r = Rf_allocVector(REALSXP, N)); nProtect++;
  double *Q_temp2 = REAL(Q_r);zeros(Q_temp2, n*n);
  double *qleSum = REAL(qleSum_r);
  
  double *Q_temp = (double *) R_alloc(n*n, sizeof(double));
  for(int s = 0; s < N; s++){ 
     updateBF(B, F, c, C, coords, nnIndx, nnIndxLU, n, m, sigmaSq, phi[s], nu, covModel, bk, nuUnifb);
    sq_err = Q(B, F, y, y, n, nnIndx, nnIndxLU);

    double *lowB = (double *) R_alloc(n*n, sizeof(double));zeros(lowB, n*n);
	double *invF = (double *) R_alloc(n*n, sizeof(double));zeros(invF, n*n);

	for(i = 0; i < n; i++){
		if(nnIndxLU[n + i]>0){		
		  for(j = 0; j < nnIndxLU[n + i]; j++){
			lowB[i + nnIndx[nnIndxLU[i] + j]*n] = - B[nnIndxLU[i] + j]; 		
		  }
		}
		
		 for(j = i; j < i + 1; j++){
			invF[j + n*i] = 1/F[i]; 
			lowB[j + n*i] = 1;		
		  }
	}


qle = 0.0;

 F77_NAME(dgemm)(&Trans, &TransN, &n, &n, &n, &one, lowB, &n, invF, &n, &zero, Q_temp,&n);
 F77_NAME(dgemm)(&TransN, &TransN, &n, &n, &n, &one, Q_temp, &n, lowB, &n, &zero, Q_temp2,&n);
 F77_NAME(dpotrf)(&lower, &n, Q_temp2, &n, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
 
  for(i = 0; i < n; i++){
    for(j = i; j < i + 1; j++){
     qle  += log(Q_temp2[j + n*i]);
    }
  }
  qleSum[s] = qle - sq_err;
  } 
  //make return object
  SEXP result_r, resultName_r;
  int nResultListObjs = 1;
  PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  
  SET_VECTOR_ELT(result_r, 0, qleSum_r);
  SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("qleSum"));
  
  // SET_VECTOR_ELT(result_r, 1, Q_r);
  // SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("Q"));
  
  
  Rf_namesgets(result_r, resultName_r);
  
  //unprotect
  UNPROTECT(nProtect);
  // SEXP out;
  // PROTECT(out = Rf_allocVector(REALSXP, nIndx));
  // B = REAL(out);
  
  return(result_r);
}















void BetaEst(double *alpha, double *A, double *B, double *F, double *u, double *y, int *nnIndx, int *nnIndxLU, int n, int p){
  int info;
  const int inc = 1;
  const double one = 1.0;
  const double zero = 0.0;
  char const *lower = "L";
  // SEXP A_r, b_r, alpha_r;
  // PROTECT(A_r = Rf_allocVector(REALSXP, p*p)); nProtect++;
  // PROTECT(b_r = Rf_allocVector(REALSXP, p)); nProtect++;
  // PROTECT(alpha_r = Rf_allocVector(REALSXP, p)); nProtect++;
  // double *A = REAL(A_r); 
  // double *b = REAL(b_r);
  // double *alpha = REAL(alpha_r);
  int pp = p*p;
  
  double *A0 = (double *) R_alloc(p*p, sizeof(double));
  //  double *alpha = (double *) R_alloc(p, sizeof(double));
  double *b = (double *) R_alloc(p, sizeof(double));
  
  int k0 =0;
  double *Xi = (double *) R_alloc(n, sizeof(double));
  double *Xj = (double *) R_alloc(n, sizeof(double));
  for(int k1 = 0; k1 < p; k1++){
    for(int i = 0; i < n; i++){
      Xi[i] = u[k1*n + i];
    }
    b[k1] = Q(B, F, Xi, y, n, nnIndx, nnIndxLU);
    for(int k2 = 0; k2 < p; k2++){
      for(int j = 0; j < n; j++){
        Xj[j] = u[k2*n+ j];
      } 	   
      A0[k0]= Q(B, F, Xi, Xj, n, nnIndx, nnIndxLU);
      k0++;
    }
  }
  
  F77_NAME(dcopy)(&pp, A0, &inc, A, &inc);
  //double *alpha_temp = (double *) R_alloc(p, sizeof(double));
  
  F77_NAME(dpotrf)(lower, &p, A0, &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
  F77_NAME(dpotri)(lower, &p, A0, &p, &info); if(info != 0){error("c++ error: dpotri failed\n");}//solve(A)
  F77_NAME(dsymv)(lower, &p, &one, A0, &p, b, &inc, &zero, alpha, &inc);
  // Rprintf("alpha: %3.10f \n", alpha_temp[0]);
  // Rprintf("alpha: %3.10f \n", alpha_temp[1]);
  // for(int k1 = 0; k1 < p; k1++){
  // alpha[k1] = alpha_temp[k1];
  // }
  //return(alpha[0]);
}


// [[Rcpp::export]]
SEXP bivariate_local_kernel_est(SEXP y_r, SEXP covZ_r, 
                      SEXP lon_r, SEXP lat_r, SEXP n_r,
                      SEXP Kernel_r, SEXP h_r, 
                      SEXP nThreads_r){
  
  const int inc = 1;
  const double one = 1.0;
  const double zero = 0.0;
  char const *lower = "L";
  char const *ntran = "N";
  char Trans = 'T';
  
  int n = INTEGER(n_r)[0];
  
  int s, i, j, info, p = 3;
  double *h = REAL(h_r);
  
  int nThreads = INTEGER(nThreads_r)[0];
  int threadID = 0;
  
  int Kernel = INTEGER(Kernel_r)[0];
  double *y = REAL(y_r);
  double *Z = REAL(covZ_r);
  double *lon = REAL(lon_r);
  double *lat = REAL(lat_r);
  
  double *Dist = (double *) R_alloc(n*nThreads*2, sizeof(double));
  //double *Dist2 = (double *) R_alloc(n*nThreads, sizeof(double));
  double *K = (double *) R_alloc(n*nThreads, sizeof(double));
  // double *X = (double *) R_alloc(p*n*nThreads, sizeof(double));
  // double *KerX = (double *) R_alloc(p*n*nThreads, sizeof(double));
  double *XtX = (double *) R_alloc(p*p*nThreads, sizeof(double));
  
  //	double *KerY = (double *) R_alloc(n, sizeof(double));
  int nProtect=0;
  SEXP alpha_r;
  PROTECT(alpha_r = Rf_allocMatrix(REALSXP, p, n)); nProtect++; //Rf_allocVector
  //double *alpha = REAL(alpha_r);
   SEXP XX_r, XK_r;
  PROTECT(XX_r = Rf_allocVector(REALSXP, p*n*nThreads)); nProtect++;
  PROTECT(XK_r = Rf_allocVector(REALSXP, p*n*nThreads)); nProtect++;
  double *X = REAL(XX_r); 
  double *KerX = REAL(XK_r);
  
  double *tmp_p = (double *) R_alloc(p*nThreads, sizeof(double));
  
  
  //double *c =(double *) R_alloc(m[0]*nThreads, sizeof(double));
  double *alpha_temp = (double *) R_alloc(n*p*nThreads, sizeof(double));
  
#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif
  
  //private(i, j, Dist, K, X, KerX, XtX, alpha_temp, tmp_p)
  
#ifdef _OPENMP
#pragma omp parallel for private(i, j, threadID)
#endif
  for(s = 0; s < n; s++){
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif
    for(i = 0; i < n; i++){
        Dist[n*threadID*2 + i] = (lon[i]  - lon[s])/h[0];// pow(Z[i]  - TestZ[s], 2);
        Dist[n*threadID*2 + n + i] = (lat[i]  - lat[s])/h[1];
		
		//Dist[n*threadID + i] = dist2(lon[i], lat[i], lon[s], lat[s]);///theta 
      if(Kernel == 1){
		  K[n*threadID + i] =  exp(-dist2(lon[i], lat[i], lon[s], lat[s])/h[0]);
		  
        //K[n*threadID + i] = exp(-abs(Dist[n*threadID*2 + i]))*exp(-abs(Dist[n*threadID*2 + n + i]))/(h[0] * h[1]);
        //K[n*threadID*2 + n + i] = exp(-abs(Dist[n*threadID*2 + n + i]))/h[1];		
      }else{
		  K[n*threadID + i] =  exp(-pow(dist2(lon[i], lat[i], lon[s], lat[s]), 2)/h[0]);
        //K[n*threadID + i] = exp(-pow(Dist[n*threadID*2 + i], 2)*h[0])*exp(-pow(Dist[n*threadID*2 + n + i], 2)*h[1])/(h[0] * h[1]); 
		//K[n*threadID*2 + n + i] = exp(-pow(Dist[n*threadID*2 + n + i], 2))/h[1]; 
      }
    }
    
    for(i = 0; i < n; i++){
      for(j = 0; j < p; j++){
        if(j==0){
          X[p*n*threadID + n*j + i] = Z[i];
          KerX[p*n*threadID + n*j + i] = K[n*threadID + i]*Z[i];
        }else{
          X[p*n*threadID + n*j + i] = Dist[n*threadID*2 + n*(j - 1) + i]*Z[i];  
          KerX[p*n*threadID + n*j + i] = K[n*threadID + i]*Dist[n*threadID*2 + n*(j - 1) + i]*Z[i]; 
        }
      }
      //KerY[i] = pow(K[i], 0.5)*y[i];
      
    }
    //Rprintf("s = %i ; KerX[0]: %3.10f \n", s, KerX[0]);
    
    F77_NAME(dgemv)(&Trans, &n, &p, &one, &KerX[p*n*threadID], &n, y, &inc, &zero, &tmp_p[p*threadID], &inc);
    F77_NAME(dgemm)(&Trans, ntran, &p, &p, &n, &one, &KerX[p*n*threadID], &n, &X[p*n*threadID], &n, &zero, &XtX[p*p*threadID], &p);
    
    F77_NAME(dpotrf)(lower, &p, &XtX[p*p*threadID], &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
    F77_NAME(dpotri)(lower, &p, &XtX[p*p*threadID], &p, &info); if(info != 0){error("c++ error: dpotri failed\n");}//solve(A)
    F77_NAME(dsymv)(lower, &p, &one, &XtX[p*p*threadID], &p, &tmp_p[p*threadID], &inc, &zero, &alpha_temp[n*p*threadID], &inc);
    
    F77_NAME(dcopy)(&p, &alpha_temp[n*p*threadID], &inc, &REAL(alpha_r)[s*p], &inc);
    // alpha[s] = alpha_temp[n*p*threadID];
  }
  
  //make return object
  SEXP result_r, resultName_r;
  int nResultListObjs = 3;
  PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  
  SET_VECTOR_ELT(result_r, 0, alpha_r);
  SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("alpha"));
  SET_VECTOR_ELT(result_r, 1, XX_r);
  SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("X"));
  SET_VECTOR_ELT(result_r, 2, XK_r);
  SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("KerX"));
  
  Rf_namesgets(result_r, resultName_r);
  
  //unprotect
  UNPROTECT(nProtect);
  // SEXP out;
  // PROTECT(out = Rf_allocVector(REALSXP, nIndx));
  // B = REAL(out);
  
  return(result_r);
  
}












// [[Rcpp::export]]
SEXP OputBF(SEXP n_r, SEXP m_r, SEXP coords_r, SEXP covModel_r, SEXP nnIndx_r, SEXP nnIndxLU_r,
            SEXP sigmaSq_r, SEXP phi_r, SEXP nu_r, SEXP nuUnifb_r, SEXP nThreads_r){
  // SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r, SEXP TestCoords_r,
  int i, j, nProtect=0;
  const double one = 1.0;
  const double zero = 0.0;
  
  //get args
  int n = INTEGER(n_r)[0];
  int m = INTEGER(m_r)[0];
  int nuUnifb = INTEGER(nuUnifb_r)[0];
  
  double *coords = REAL(coords_r);
  
  int *nnIndx = INTEGER(nnIndx_r);
  int *nnIndxLU = INTEGER(nnIndxLU_r);
  // int *uIndx = INTEGER(uIndx_r);
  // int *uIndxLU = INTEGER(uIndxLU_r);
  // int *uiIndx = INTEGER(uiIndx_r);
  
  
  int covModel = INTEGER(covModel_r)[0];
  std::string corName = getCorName(covModel);
  
  
  int nThreads = INTEGER(nThreads_r)[0];
  
  
#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif
  
  //parameters
  int nTheta, phiIndx, nuIndx;
  
  if(corName != "matern"){
    nTheta = 1;//sigma^2, tau^2, phi
    // sigmaSqIndx = 0; tauSqIndx = 1; 
    phiIndx = 1;
  }else{
    nTheta = 2;//sigma^2, tau^2, phi, nu,
    //sigmaSqIndx = 0; tauSqIndx = 1; 
    phiIndx = 1; nuIndx = 2;
  }
  
  //starting	
  
  double *theta = (double *) R_alloc(nTheta, sizeof(double));
  double *sigmaSq = REAL(sigmaSq_r);
  // theta[sigmaSqIndx] = REAL(sigmaSq_r)[0];
  //theta[tauSqIndx] = REAL(tauSq_r)[0];
  theta[phiIndx] = REAL(phi_r)[0];
  if(corName == "matern"){
    theta[nuIndx] = REAL(nu_r)[0];
  }
  
  //allocate for the U index vector that keep track of which locations have the i-th location as a neighbor
  int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
  
  //other stuff
  int mm = m*m;
  
  
  SEXP neiB_r, varF_r;
  PROTECT(neiB_r = Rf_allocVector(REALSXP, nIndx)); nProtect++;
  PROTECT(varF_r = Rf_allocVector(REALSXP, n)); nProtect++;
  double *B = REAL(neiB_r);
  double *F = REAL(varF_r);
  //double *B = (double *) R_alloc(nIndx, sizeof(double));
  //double *F = (double *) R_alloc(n, sizeof(double));
  
  double *c =(double *) R_alloc(m*nThreads, sizeof(double));
  double *C = (double *) R_alloc(mm*nThreads, sizeof(double));
  //int nuUnifb = 0;
  double *bk = (double *) R_alloc(nThreads*(1.0+static_cast<int>(floor(nuUnifb))), sizeof(double));
  
  double nu = 0;
  if(corName == "matern"){nu = theta[nuIndx];}
  
  updateBF(B, F, c, C, coords, nnIndx, nnIndxLU, n, m, sigmaSq, theta[phiIndx], nu, covModel, bk, nuUnifb);
  
  SEXP lowB_r, invF_r;
  PROTECT(lowB_r = Rf_allocMatrix(REALSXP, n, n)); nProtect++;
  PROTECT(invF_r = Rf_allocMatrix(REALSXP, n, n)); nProtect++;
  
  double *lowB = REAL(lowB_r);zeros(lowB, n*n);
  double *invF = REAL(invF_r);zeros(invF, n*n);
  
  // double *lowB = (double *) R_alloc(n*n, sizeof(double));zeros(lowB, n*n);
  // double *Fmat = (double *) R_alloc(n*n, sizeof(double));zeros(lowB, n*n);
  
  for(i = 0; i < n; i++){
    if(nnIndxLU[n + i]>0){		
      for(j = 0; j < nnIndxLU[n + i]; j++){
        lowB[i + nnIndx[nnIndxLU[i] + j]*n] = - B[nnIndxLU[i] + j]; 		
      }
    }
    // if(nnIndxLU[n + 1]==1){
    // lowB[1 + nnIndx[nnIndxLU[1]]*n] = B[nnIndxLU[1]]; 
    // Rprintf("lowB: %i., %3.5f. \n", nnIndxLU[1], B[nnIndxLU[1]]);		
    // }
    for(j = i; j < i + 1; j++){
      invF[j + n*i] = 1/F[i]; 
      lowB[j + n*i] = 1;		
    }
  }
  //Rprintf("lowB: %i. \n", kkk);
  char Trans = 'T';
  char TransN = 'N';
  
  
  double *Q_temp = (double *) R_alloc(n*n, sizeof(double));
  //double *Q = (double *) R_alloc(n*n, sizeof(double));
  SEXP Q_r;
  PROTECT(Q_r = Rf_allocMatrix(REALSXP, n, n)); nProtect++;
  double *Q = REAL(Q_r);zeros(Q, n*n);
  
  F77_NAME(dgemm)(&Trans, &TransN, &n, &n, &n, &one, lowB, &n, invF, &n, &zero, Q_temp,&n);
  F77_NAME(dgemm)(&TransN, &TransN, &n, &n, &n, &one, Q_temp, &n, lowB, &n, &zero, Q,&n);
  
  // PutRNGstate();
  //make return object
  SEXP result_r, resultName_r;
  int nResultListObjs = 5;
  PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  
  SET_VECTOR_ELT(result_r, 0, neiB_r);
  SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("neiB"));
  SET_VECTOR_ELT(result_r, 1, invF_r);
  SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("invF"));
  SET_VECTOR_ELT(result_r, 2, varF_r);
  SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("varF"));
  
  SET_VECTOR_ELT(result_r, 3, lowB_r);
  SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("lowB"));
  SET_VECTOR_ELT(result_r, 4, Q_r);
  SET_VECTOR_ELT(resultName_r, 4, Rf_mkChar("Q"));
  
  Rf_namesgets(result_r, resultName_r);
  
  //unprotect
  UNPROTECT(nProtect);
  // SEXP out;
  // PROTECT(out = Rf_allocVector(REALSXP, nIndx));
  // B = REAL(out);
  
  return(result_r);
}

// [[Rcpp::export]]
SEXP SemiAlphaPro(SEXP y_r, SEXP Z_r, SEXP n_r,	   
				   SEXP coords_r, 
				   SEXP Q_r,
				   SEXP nnIndx_r, 
				   SEXP nnIndxLU_r,
				   SEXP Kernel_r,
				   SEXP h_r,	
				   SEXP GeomVariable_r,
				   SEXP nThreads_r){
  int i, j, info, nProtect=0;
  const int inc = 1;
  const double one = 1.0;
  const double zero = 0.0;
  //get args
  int n = INTEGER(n_r)[0];
  int Kernel = INTEGER(Kernel_r)[0];
  int p = 2;  int pp = p*p; int np = n*p;
  double h = REAL(h_r)[0];
    
  int GeomVariable = INTEGER(GeomVariable_r)[0];
  
  double *y = REAL(y_r);
  // double *x = REAL(x_r);
  double *Z = REAL(Z_r);
  double *coords = REAL(coords_r);
  double *Q = REAL(Q_r);
  int *nnIndx = INTEGER(nnIndx_r);
  int *nnIndxLU = INTEGER(nnIndxLU_r);
  
  int nThreads = INTEGER(nThreads_r)[0];
  char const *lower = "L";
  char const *ntran = "N";
  char Trans = 'T';
  
#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif
  
  double *KerZ = (double *) R_alloc(n*p*nThreads, sizeof(double));
  double *tKerZ = (double *) R_alloc(p*n*nThreads, sizeof(double));
  double *Dist = (double *) R_alloc(n*nThreads, sizeof(double));
  double *K = (double *) R_alloc(n*nThreads, sizeof(double));
  //double *KerY = (double *) R_alloc(n*nThreads, sizeof(double));
  double *DagKer = (double *) R_alloc(n*n*nThreads, sizeof(double)); zeros(DagKer, n*n*nThreads);

  SEXP alpha_r, S_r;
  PROTECT(alpha_r = Rf_allocMatrix(REALSXP, p, n)); nProtect++;
  PROTECT(S_r = Rf_allocMatrix(REALSXP, p, n*n)); nProtect++;
  //double *Alpha = REAL(alpha_r);
  double *alpha = (double *) R_alloc(p*nThreads, sizeof(double));
  
  double *Q_temp = (double *) R_alloc(n*p*nThreads, sizeof(double)); zeros(Q_temp, p*n*nThreads);
  double *S = (double *) R_alloc(p*n*nThreads, sizeof(double)); zeros(S, p*n*nThreads);
  double *MAT = (double *) R_alloc(p*p*nThreads, sizeof(double)); zeros(MAT, p*p*nThreads);


 // SEXP MAT_r;
  // PROTECT(MAT_r = Rf_allocMatrix(REALSXP, n, p)); nProtect++;
// double *KerZ = REAL(MAT_r);
int threadID = 0;


#ifdef _OPENMP
#pragma omp parallel for private(i, j, threadID)
#endif
  for(int s = 0; s < n; s++){
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif
    for(i = 0; i < n; i++){
      if(GeomVariable == 1){
        Dist[n*threadID + i] = (Z[i]  - Z[s])/h;// pow(Z[i]  - TestZ[s], 2);
      }else{
        Dist[n*threadID + i] = dist2(coords[s], coords[n + s], coords[i], coords[n + i]);///theta  
      }
      if(Kernel == 1){
        K[n*threadID + i] = exp(-abs(Dist[n*threadID + i]))/h;  
      }else{
        K[n*threadID + i] = exp(-pow(Dist[n*threadID + i], 2)*h)/h; 
      }
    }
    
    for(i = 0; i < n; i++){
      for(j = 0; j < p; j++){
        if(j==0){
          //X[n*j + i] = 1;
          KerZ[n*p*threadID + n*j + i] = pow(K[n*threadID + i], 0.5);
        }else{
          //X[n*j + i] = Dist[i]/phi;  
          KerZ[n*p*threadID + n*j + i] = pow(K[n*threadID + i], 0.5)*Dist[n*threadID + i];
        }
      }
    }

	for(i = 0; i < n; i++){
		for(j = i; j < i + 1; j++){ 
		  DagKer[n*n*threadID + j + n*i] = pow(K[n*threadID + i], 0.5);		
		}
	}
     F77_NAME(dgemm)(&Trans, ntran, &p, &n, &n, &one, &KerZ[n*p*threadID], &n, Q, &n, &zero, &Q_temp[n*p*threadID],&p);	
     F77_NAME(dgemm)(ntran, ntran, &p, &p, &n, &one, &Q_temp[n*p*threadID], &p, &KerZ[n*p*threadID], &n, &zero, &MAT[p*p*threadID],&p);		
		
	 F77_NAME(dposv)(lower, &p, &n, &MAT[p*p*threadID], &p, &Q_temp[n*p*threadID], &p, &info);if(info != 0){error("c++ error: dpotrf failed\n");}
		
	 F77_NAME(dgemm)(ntran, ntran, &p, &n, &n, &one, &Q_temp[n*p*threadID], &p, &DagKer[n*n*threadID], &n, &zero, &S[n*p*threadID],&p);	
     F77_NAME(dcopy)(&np, &S[n*p*threadID], &inc, &REAL(S_r)[s*np], &inc);
	
	 F77_NAME(dgemv)(ntran, &p, &n, &one, &S[n*p*threadID], &p, y, &inc, &zero, &alpha[p*threadID], &inc);	
	 F77_NAME(dcopy)(&p, &alpha[p*threadID], &inc, &REAL(alpha_r)[s*p], &inc);
	
  }
  
  
  // PutRNGstate();
  //make return object
  SEXP result_r, resultName_r;
  int nResultListObjs = 2;
  PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  
  SET_VECTOR_ELT(result_r, 0, alpha_r);
  SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("alpha"));
  
  SET_VECTOR_ELT(result_r, 1, S_r);
  SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("S"));
  Rf_namesgets(result_r, resultName_r);
  
  //unprotect
  UNPROTECT(nProtect);
  // SEXP out;
  // PROTECT(out = Rf_allocVector(REALSXP, nIndx));
  // B = REAL(out);
  
  return(result_r);
}

// [[Rcpp::export]]
SEXP semiProPred(SEXP y_r, 
              SEXP Z_r, 
              SEXP TestZ_r, 
              SEXP n_r, 
              SEXP nTest_r,
              SEXP coords_r, 
              SEXP TestCoords_r,
			  SEXP Q_r,
              SEXP nnIndx_r, 
              SEXP nnIndxLU_r,
              SEXP Kernel_r,
              SEXP h_r,	
              SEXP GeomVariable_r,				 
              SEXP nThreads_r){
  // SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r, SEXP TestCoords_r,
  const int inc = 1;
  const double one = 1.0;
  const double zero = 0.0;
  char const *lower = "L";
  char const *ntran = "N";
  char Trans = 'T';
  int i, j, info, nProtect=0;
  int p = 2; 
  int pp = p*p; 
  //get args
  int n = INTEGER(n_r)[0];
  int np = n*p;
  int Kernel = INTEGER(Kernel_r)[0];
  int nTest = INTEGER(nTest_r)[0];
  double h = REAL(h_r)[0];
  int GeomVariable = INTEGER(GeomVariable_r)[0];
  double *y = REAL(y_r);
  // double *x = REAL(x_r);
  
  double *Z = REAL(Z_r);
  double *TestZ = REAL(TestZ_r);
  
  double *coords = REAL(coords_r);
  double *TestCoords = REAL(TestCoords_r);
  
  
  double *Q = REAL(Q_r);
  int *nnIndx = INTEGER(nnIndx_r);
  int *nnIndxLU = INTEGER(nnIndxLU_r);
  int nThreads = INTEGER(nThreads_r)[0];
  
  
#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif
  
  double *KerZ = (double *) R_alloc(p*n*nThreads, sizeof(double));
  double *Dist = (double *) R_alloc(n*nThreads, sizeof(double));
  double *K = (double *) R_alloc(n*nThreads, sizeof(double));
  double *DagKer = (double *) R_alloc(n*n*nThreads, sizeof(double)); zeros(DagKer, n*n*nThreads);
  double *A = (double *) R_alloc(p*p, sizeof(double));
  
  //double  phi0 = 0.1;//theta[phiIndx];
  SEXP alpha_r, S_r;
  PROTECT(alpha_r = Rf_allocMatrix(REALSXP, p, nTest)); nProtect++;
  PROTECT(S_r = Rf_allocMatrix(REALSXP, p, n*nTest)); nProtect++;
  //double *Alpha = REAL(alpha_r);
  double *alpha = (double *) R_alloc(p*nThreads, sizeof(double));

  double *Q_temp = (double *) R_alloc(n*p*nThreads, sizeof(double)); zeros(Q_temp, p*n*nThreads);
  double *S = (double *) R_alloc(p*n*nThreads, sizeof(double)); zeros(S, p*n*nThreads);
  double *MAT = (double *) R_alloc(p*p*nThreads, sizeof(double)); zeros(MAT, p*p*nThreads);
  
  
int threadID = 0;


#ifdef _OPENMP
#pragma omp parallel for private(i, j, threadID)
#endif
  for(int s = 0; s < nTest; s++){
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif
    for(i = 0; i < n; i++){
      if(GeomVariable == 1){
        Dist[n*threadID + i] = (Z[i]  - TestZ[s])/h;// pow(Z[i]  - TestZ[s], 2);
      }else{
        Dist[n*threadID + i] = dist2(TestCoords[s], TestCoords[nTest + s], coords[i], coords[n + i]);///theta  
      }
      if(Kernel == 1){
        K[n*threadID + i] = exp(-abs(Dist[n*threadID + i]))/h;  
      }else{
        K[n*threadID + i] = exp(-pow(Dist[n*threadID + i], 2)*h)/h; 
      }
    }
    
    for(i = 0; i < n; i++){
      for(j = 0; j < p; j++){
        if(j==0){
          KerZ[p*n*threadID + n*j + i] = pow(K[n*threadID +i], 0.5);
        }else{
          KerZ[p*n*threadID + n*j + i] = pow(K[n*threadID +i], 0.5)*Dist[n*threadID + i];
        }
      }
    }
	
    for(i = 0; i < n; i++){
		for(j = i; j < i + 1; j++){ 
		  DagKer[n*n*threadID + j + n*i] = pow(K[n*threadID +i], 0.5);		
		}
	}
	 
	F77_NAME(dgemm)(&Trans, ntran, &p, &n, &n, &one, &KerZ[n*p*threadID], &n, Q, &n, &zero, &Q_temp[n*p*threadID],&p);	
     F77_NAME(dgemm)(ntran, ntran, &p, &p, &n, &one, &Q_temp[n*p*threadID], &p, &KerZ[n*p*threadID], &n, &zero, &MAT[p*p*threadID],&p);		
		
	 F77_NAME(dposv)(lower, &p, &n, &MAT[p*p*threadID], &p, &Q_temp[n*p*threadID], &p, &info);if(info != 0){error("c++ error: dpotrf failed\n");}
		
	 F77_NAME(dgemm)(ntran, ntran, &p, &n, &n, &one, &Q_temp[n*p*threadID], &p, &DagKer[n*n*threadID], &n, &zero, &S[n*p*threadID],&p);	
     F77_NAME(dcopy)(&np, &S[n*p*threadID], &inc, &REAL(S_r)[s*np], &inc);
	
	 F77_NAME(dgemv)(ntran, &p, &n, &one, &S[n*p*threadID], &p, y, &inc, &zero, &alpha[p*threadID], &inc);	
	 F77_NAME(dcopy)(&p, &alpha[p*threadID], &inc, &REAL(alpha_r)[s*p], &inc);
  }
  

  //make return object
  SEXP result_r, resultName_r;
  int nResultListObjs = 2;
  PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  
  SET_VECTOR_ELT(result_r, 0, alpha_r);
  SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("alpha"));

  SET_VECTOR_ELT(result_r, 1, S_r);
  SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("S"));
  Rf_namesgets(result_r, resultName_r);
  //unprotect
  UNPROTECT(nProtect);
  // SEXP out;
  // PROTECT(out = Rf_allocVector(REALSXP, nIndx));
  // B = REAL(out);
  
  return(result_r);
}
















// [[Rcpp::export]]
SEXP SemiAlpha(SEXP y_r, SEXP Z_r, SEXP n_r,	   
               SEXP coords_r, 
               SEXP B_r, SEXP F_r, SEXP Q_r,
               SEXP nnIndx_r, 
               SEXP nnIndxLU_r,
               SEXP Kernel_r,
               SEXP h_r,	
               SEXP GeomVariable_r,
               SEXP nThreads_r){
  int i, j, info, nProtect=0;
  const int inc = 1;
  const double one = 1.0;
  const double zero = 0.0;
  //get args
  int n = INTEGER(n_r)[0];
  int Kernel = INTEGER(Kernel_r)[0];
  int p = 2;  int pp = p*p; int np = n*p;
  double h = REAL(h_r)[0];
  
  int GeomVariable = INTEGER(GeomVariable_r)[0];
  
  double *y = REAL(y_r);
  // double *x = REAL(x_r);
  double *Z = REAL(Z_r);
  double *coords = REAL(coords_r);
  double *B = REAL(B_r);
  double *F = REAL(F_r);
  double *Q = REAL(Q_r);
  int *nnIndx = INTEGER(nnIndx_r);
  int *nnIndxLU = INTEGER(nnIndxLU_r);
  
  int nThreads = INTEGER(nThreads_r)[0];
  char const *lower = "L";
  char const *ntran = "N";
  char Trans = 'T';
  
#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif
  
  double *KerZ = (double *) R_alloc(p*n, sizeof(double));
  double *tKerZ = (double *) R_alloc(p*n, sizeof(double));
  double *Dist = (double *) R_alloc(n, sizeof(double));
  double *K = (double *) R_alloc(n, sizeof(double));
  double *KerY = (double *) R_alloc(n, sizeof(double));
  double *DagKer = (double *) R_alloc(n*n, sizeof(double)); zeros(DagKer, n*n);
  double *A = (double *) R_alloc(p*p, sizeof(double));
  
  double *Q_temp1 = (double *) R_alloc(n*n, sizeof(double));
  double *Q_temp2 = (double *) R_alloc(p*n, sizeof(double));
  SEXP alpha_r, S_r;
  PROTECT(alpha_r = Rf_allocMatrix(REALSXP, p, n)); nProtect++;
  PROTECT(S_r = Rf_allocMatrix(REALSXP, p, n*n)); nProtect++;
  //double *Alpha = REAL(alpha_r);
  double *alpha = (double *) R_alloc(p, sizeof(double));
  
  // #ifdef _OPENMP
  // #pragma omp parallel for private(i, j, alpha)
  // #endif   
  for(int s = 0; s < n; s++){
    //solve distance from testing location to basis table
    // #ifdef _OPENMP
    // #pragma omp parallel 
    // #endif
    for(i = 0; i < n; i++){
      if(GeomVariable == 1){
        Dist[i] = (Z[i]  - Z[s])/h;// pow(Z[i]  - TestZ[s], 2);
      }else{
        Dist[i] = dist2(coords[s], coords[n + s], coords[i], coords[n + i]);///theta  
      }
      if(Kernel == 1){
        K[i] = exp(-abs(Dist[i]))/h;  
      }else{
        K[i] = exp(-pow(Dist[i], 2)*h)/h; 
      }
    }
    
    for(i = 0; i < n; i++){
      for(j = 0; j < p; j++){
        if(j==0){
          //X[n*j + i] = 1;
          KerZ[n*j + i] = pow(K[i], 0.5);
        }else{
          //X[n*j + i] = Dist[i]/phi;  
          KerZ[n*j + i] = pow(K[i], 0.5)*Dist[i];
        }
      }
      KerY[i] = pow(K[i], 0.5)*y[i];
      // Rprintf("alpha[i]: %3.10f \n", KerZ[i]);
    }
    
   // for(j = 0; j < p; j++){
     // for(i = 0; i < n; i++){
	   // if(j == 0){
		 // tKerZ[p*i + j] = pow(K[i], 0.5);  
	   // }else{
		 // tKerZ[p*i + j] = pow(K[i], 0.5)*Dist[i];  
	   // }		  
	// }
   // }
	for(i = 0; i < n; i++){
		for(j = i; j < i + 1; j++){ 
		  DagKer[j + n*i] = pow(K[i], 0.5);		
		}
	}
	
	
    BetaEst(alpha, A, B, F, KerZ, KerY, nnIndx, nnIndxLU, n, p);
	F77_NAME(dcopy)(&p, alpha, &inc, &REAL(alpha_r)[s*p], &inc);
		
		
	F77_NAME(dgemm)(ntran, ntran, &n, &n, &n, &one, Q, &n, DagKer, &n, &zero, Q_temp1,&n);
	 F77_NAME(dgemm)(&Trans, ntran, &p, &n, &n, &one, KerZ, &n, Q_temp1, &n, &zero, Q_temp2,&p);
	 F77_NAME(dposv)(lower, &p, &n, A, &p, Q_temp2, &p, &info);if(info != 0){error("c++ error: dpotrf failed\n");}
	F77_NAME(dcopy)(&np, Q_temp2, &inc, &REAL(S_r)[s*p*n], &inc);
	
   // F77_NAME(dcopy)(&pp, Q_temp2, &inc, &REAL(S_r)[s*pp], &inc);

    // for(int k1 = 0; k1 < p; k1++){
    // Alpha[s*p + k1] = alpha[k1];
    // Rprintf("B: %3.10f \n", Alpha[s*p + k1]);
    // Rprintf("alpha[k1]: %3.10f \n", alpha[k1]);
    // }
  }
  
  
  // PutRNGstate();
  //make return object
  SEXP result_r, resultName_r;
  int nResultListObjs = 2;
  PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  
  SET_VECTOR_ELT(result_r, 0, alpha_r);
  SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("alpha"));
  
  SET_VECTOR_ELT(result_r, 1, S_r);
  SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("S"));
  Rf_namesgets(result_r, resultName_r);
  
  //unprotect
  UNPROTECT(nProtect);
  // SEXP out;
  // PROTECT(out = Rf_allocVector(REALSXP, nIndx));
  // B = REAL(out);
  
  return(result_r);
}

// // [[Rcpp::export]]
// SEXP SemiAlpha1(SEXP y_r, SEXP Z_r, SEXP n_r, SEXP h_r,			   
//                 SEXP coords_r, 
//                 SEXP B_r, SEXP F_r, 
//                 SEXP nnIndx_r, 
//                 SEXP nnIndxLU_r,
//                 SEXP GeomVariable_r,
//                 SEXP Kernel_r,
//                 SEXP nThreads_r){
//   int i, j, nProtect=0;
//   const int inc = 1;
//   
//   //get args
//   int n = INTEGER(n_r)[0];
//   int Kernel = INTEGER(Kernel_r)[0];
//   int p = 2;
//   double h = REAL(h_r)[0];
//   
//   int GeomVariable = INTEGER(GeomVariable_r)[0];
//   
//   double *y = REAL(y_r);
//   // double *x = REAL(x_r);
//   double *Z = REAL(Z_r);
//   double *coords = REAL(coords_r);
//   double *B = REAL(B_r);
//   double *F = REAL(F_r);
//   int *nnIndx = INTEGER(nnIndx_r);
//   int *nnIndxLU = INTEGER(nnIndxLU_r);
//   
//   int nThreads = INTEGER(nThreads_r)[0];
//   
//   
// #ifdef _OPENMP
//   omp_set_num_threads(nThreads);
// #else
//   if(nThreads > 1){
//     warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
//     nThreads = 1;
//   }
// #endif
//   
//   double *KerZ = (double *) R_alloc(p*n, sizeof(double));
//   double *Dist = (double *) R_alloc(n, sizeof(double));
//   double *K = (double *) R_alloc(n, sizeof(double));
//   double *KerY = (double *) R_alloc(n, sizeof(double));
//   SEXP alpha_r;
//   PROTECT(alpha_r = Rf_allocMatrix(REALSXP, p, n)); nProtect++;
//   //double *Alpha = REAL(alpha_r);
//   double *alpha = (double *) R_alloc(p, sizeof(double));
//   
//   // #ifdef _OPENMP
//   // #pragma omp parallel for private(i, j, alpha)
//   // #endif   
//   for(int s = 0; s < n; s++){
//     //solve distance from testing location to basis table
//     // #ifdef _OPENMP
//     // #pragma omp parallel 
//     // #endif
//     for(i = 0; i < n; i++){
//       if(GeomVariable == 1){
//         Dist[i] = (Z[i]  - Z[s])/h;// pow(Z[i]  - TestZ[s], 2);
//       }else{
//         Dist[i] = dist2(coords[s], coords[n + s], coords[i], coords[n + i]);///theta  
//       }
//       if(Kernel == 1){
//         K[i] = exp(-abs(Dist[i]))/h;  
//       }else{
//         K[i] = exp(-pow(Dist[i], 2))/h; 
//       }
//     }
//     
//     for(i = 0; i < n; i++){
//       for(j = 0; j < p; j++){
//         if(j==0){
//           //X[n*j + i] = 1;
//           KerZ[n*j + i] = pow(K[i], 0.5);
//         }else{
//           //X[n*j + i] = Dist[i]/phi;  
//           KerZ[n*j + i] = pow(K[i], 0.5)*Dist[i];
//         }
//       }
//       KerY[i] = pow(K[i], 0.5)*y[i];
//       // Rprintf("alpha[i]: %3.10f \n", KerZ[i]);
//     }
//     
//     BetaEst(alpha, B, F, KerZ, KerY, nnIndx, nnIndxLU, n, p);
//     
//     F77_NAME(dcopy)(&p, alpha, &inc, &REAL(alpha_r)[s*p], &inc);
//     // for(int k1 = 0; k1 < p; k1++){
//     // Alpha[s*p + k1] = alpha[k1];
//     // Rprintf("B: %3.10f \n", Alpha[s*p + k1]);
//     // Rprintf("alpha[k1]: %3.10f \n", alpha[k1]);
//     // }
//   }
//   
//   
//   // PutRNGstate();
//   //make return object
//   SEXP result_r, resultName_r;
//   int nResultListObjs = 1;
//   PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
//   PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
//   
//   SET_VECTOR_ELT(result_r, 0, alpha_r);
//   SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("alpha"));
//   
//   Rf_namesgets(result_r, resultName_r);
//   
//   //unprotect
//   UNPROTECT(nProtect);
//   // SEXP out;
//   // PROTECT(out = Rf_allocVector(REALSXP, nIndx));
//   // B = REAL(out);
//   
//   return(result_r);
// }


// [[Rcpp::export]]
SEXP semiPred(SEXP y_r, 
              SEXP Z_r, 
              SEXP TestZ_r, 
              SEXP n_r, 
              SEXP nTest_r,
              SEXP coords_r, 
              SEXP TestCoords_r, 
              SEXP B_r,
              SEXP F_r, SEXP Q_r,
              SEXP nnIndx_r, 
              SEXP nnIndxLU_r,
              SEXP Kernel_r,
              SEXP h_r,	
              SEXP GeomVariable_r,				 
              SEXP nThreads_r){
  // SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r, SEXP TestCoords_r,
  const int inc = 1;
  const double one = 1.0;
  const double zero = 0.0;
  char const *lower = "L";
  char const *ntran = "N";
  char Trans = 'T';
  int i, j, info, nProtect=0;
  int p = 2; 
  int pp = p*p; 
  
  //get args
  int n = INTEGER(n_r)[0];
  int np = n*p;
  int Kernel = INTEGER(Kernel_r)[0];
  int nTest = INTEGER(nTest_r)[0];
  double h = REAL(h_r)[0];
  int GeomVariable = INTEGER(GeomVariable_r)[0];
  double *y = REAL(y_r);
  // double *x = REAL(x_r);
  
  double *Z = REAL(Z_r);
  double *TestZ = REAL(TestZ_r);
  
  double *coords = REAL(coords_r);
  double *TestCoords = REAL(TestCoords_r);
  
  double *B = REAL(B_r);
  double *F = REAL(F_r);
  double *Q = REAL(Q_r);
  int *nnIndx = INTEGER(nnIndx_r);
  int *nnIndxLU = INTEGER(nnIndxLU_r);
  int nThreads = INTEGER(nThreads_r)[0];
  
  
#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif
  
  double *KerZ = (double *) R_alloc(p*n, sizeof(double));
  double *Dist = (double *) R_alloc(n, sizeof(double));
  double *K = (double *) R_alloc(n, sizeof(double));
  double *KerY = (double *) R_alloc(n, sizeof(double));
  double *DagKer = (double *) R_alloc(n*n, sizeof(double)); zeros(DagKer, n*n);
  double *A = (double *) R_alloc(p*p, sizeof(double));
  
  double *Q_temp1 = (double *) R_alloc(n*n, sizeof(double));
  double *Q_temp2 = (double *) R_alloc(p*n, sizeof(double));
  //double  phi0 = 0.1;//theta[phiIndx];
  SEXP alpha_r, S_r;
  PROTECT(alpha_r = Rf_allocMatrix(REALSXP, p, nTest)); nProtect++;
  PROTECT(S_r = Rf_allocMatrix(REALSXP, p, n*nTest)); nProtect++;
  //double *Alpha = REAL(alpha_r);
  double *alpha = (double *) R_alloc(p, sizeof(double));
  // #ifdef _OPENMP
  // #pragma omp parallel for private(i, j, alpha)
  // #endif 
  
  for(int s = 0; s < nTest; s++){
    //solve distance from testing location to basis table
    // #ifdef _OPENMP
    // #pragma omp parallel 
    // #endif
    for(i = 0; i < n; i++){
      if(GeomVariable == 1){
        Dist[i] = (Z[i]  - TestZ[s])/h;// pow(Z[i]  - TestZ[s], 2);
      }else{
        Dist[i] = dist2(TestCoords[s], TestCoords[nTest + s], coords[i], coords[n + i]);///theta  
      }
      if(Kernel == 1){
        K[i] = exp(-abs(Dist[i]))/h;  
      }else{
        K[i] = exp(-pow(Dist[i], 2)*h)/h; 
      }
    }
    
    for(i = 0; i < n; i++){
      for(j = 0; j < p; j++){
        if(j==0){
          KerZ[n*j + i] = pow(K[i], 0.5);
        }else{
          KerZ[n*j + i] = pow(K[i], 0.5)*Dist[i];
        }
      }
      KerY[i] = pow(K[i], 0.5)*y[i];
    }
	
    for(i = 0; i < n; i++){
		for(j = i; j < i + 1; j++){ 
		  DagKer[j + n*i] = pow(K[i], 0.5);		
		}
	}
	 BetaEst(alpha, A, B, F, KerZ, KerY, nnIndx, nnIndxLU, n, p);
	 F77_NAME(dcopy)(&p, alpha, &inc, &REAL(alpha_r)[s*p], &inc);
	 
	 F77_NAME(dgemm)(ntran, ntran, &n, &n, &n, &one, Q, &n, DagKer, &n, &zero, Q_temp1,&n);
	 F77_NAME(dgemm)(&Trans, ntran, &p, &n, &n, &one, KerZ, &n, Q_temp1, &n, &zero, Q_temp2,&p);
	 F77_NAME(dposv)(lower, &p, &n, A, &p, Q_temp2, &p, &info);if(info != 0){error("c++ error: dpotrf failed\n");}
	 F77_NAME(dcopy)(&np, Q_temp2, &inc, &REAL(S_r)[s*p*n], &inc);
   // BetaEst(alpha, B, F, KerZ, KerY, nnIndx, nnIndxLU, n, p);
    //F77_NAME(dcopy)(&p, alpha, &inc, &REAL(alpha_r)[s*p], &inc);
  }
  

  //make return object
  SEXP result_r, resultName_r;
  int nResultListObjs = 2;
  PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
  
  SET_VECTOR_ELT(result_r, 0, alpha_r);
  SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("alpha"));

  SET_VECTOR_ELT(result_r, 1, S_r);
  SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("S"));
  Rf_namesgets(result_r, resultName_r);
  //unprotect
  UNPROTECT(nProtect);
  // SEXP out;
  // PROTECT(out = Rf_allocVector(REALSXP, nIndx));
  // B = REAL(out);
  
  return(result_r);
}

// // [[Rcpp::export]]
// SEXP SemiBeta(SEXP y_r, SEXP Z_r, 
//               SEXP n_r, SEXP m_r, SEXP h_r,			   
//               SEXP coords_r, 
//               SEXP GeomVariable_r,
//               SEXP covModel_r, SEXP nnIndx_r, SEXP nnIndxLU_r,
//               SEXP sigmaSq_r, SEXP tauSq_r,
//               SEXP phi_r, SEXP nu_r, SEXP nThreads_r){
//   // SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r, SEXP TestCoords_r,
//   int i, j, nProtect=0;
//   const int inc = 1;
//   const double zero = 0.0;
//   const double one = 1.0;
//   //get args
//   int n = INTEGER(n_r)[0];
//   int m = INTEGER(m_r)[0];
//   int p = 2;
//   double h = REAL(h_r)[0];
//   
//   
//   int GeomVariable = INTEGER(GeomVariable_r)[0];
//   
//   
//   double *y = REAL(y_r);
//   // double *x = REAL(x_r);
//   
//   double *Z = REAL(Z_r);
//   
//   double *coords = REAL(coords_r);
//   
//   int *nnIndx = INTEGER(nnIndx_r);
//   int *nnIndxLU = INTEGER(nnIndxLU_r);
//   // int *uIndx = INTEGER(uIndx_r);
//   // int *uIndxLU = INTEGER(uIndxLU_r);
//   // int *uiIndx = INTEGER(uiIndx_r);
//   
//   
//   int covModel = INTEGER(covModel_r)[0];
//   std::string corName = getCorName(covModel);
//   
//   
//   int nThreads = INTEGER(nThreads_r)[0];
//   
//   
// #ifdef _OPENMP
//   omp_set_num_threads(nThreads);
// #else
//   if(nThreads > 1){
//     warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
//     nThreads = 1;
//   }
// #endif
//   
//   //parameters
//   int nTheta, sigmaSqIndx, tauSqIndx, phiIndx, nuIndx;
//   
//   if(corName != "matern"){
//     nTheta = 3;//sigma^2, tau^2, phi
//     sigmaSqIndx = 0; tauSqIndx = 1; phiIndx = 2;
//   }else{
//     nTheta = 4;//sigma^2, tau^2, phi, nu,
//     sigmaSqIndx = 0; tauSqIndx = 1; phiIndx = 2; nuIndx = 3;
//   }
//   
//   //starting	
//   
//   double *theta = (double *) R_alloc(nTheta, sizeof(double));
//   
//   theta[sigmaSqIndx] = REAL(sigmaSq_r)[0];
//   theta[tauSqIndx] = REAL(tauSq_r)[0];
//   theta[phiIndx] = REAL(phi_r)[0];
//   if(corName == "matern"){
//     theta[nuIndx] = REAL(nu_r)[0];
//   }
//   
//   //allocate for the U index vector that keep track of which locations have the i-th location as a neighbor
//   int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
//   
//   //other stuff
//   int mm = m*m;
//   
//   SEXP B_r, F_r;
//   PROTECT(B_r = Rf_allocVector(REALSXP, nIndx)); nProtect++;
//   PROTECT(F_r = Rf_allocVector(REALSXP, n)); nProtect++;
//   double *B = REAL(B_r);
//   double *F = REAL(F_r);
//   //double *B = (double *) R_alloc(nIndx, sizeof(double));
//   //double *F = (double *) R_alloc(n, sizeof(double));
//   
//   double *c =(double *) R_alloc(m*nThreads, sizeof(double));
//   double *C = (double *) R_alloc(mm*nThreads, sizeof(double));
//   int nuUnifb = 0;
//   double *bk = (double *) R_alloc(nThreads*(1.0+static_cast<int>(floor(nuUnifb))), sizeof(double));
//   
//   double nu = 0;
//   if(corName == "matern"){nu = theta[nuIndx];}
//   
//   updateBF1(B, F, c, C, coords, nnIndx, nnIndxLU, n, m, 
//             theta[sigmaSqIndx], theta[phiIndx], nu, 
//             covModel, bk, nuUnifb);
//   // GetRNGstate();
//   
//   // SEXP UV_r, H_r;
//   // PROTECT(UV_r = Rf_allocVector(REALSXP, 1)); nProtect++;
//   // double *UV = REAL(UV_r);
//   // UV[0] = Q(B, F, x, y, n, nnIndx, nnIndxLU);
//   
//   // PROTECT(H_r = Rf_allocVector(REALSXP, nIndx)); nProtect++;
//   // double *H = REAL(H_r); 
//   
//   // H = B;
//   // for(i = 0; i < n; i++){
//   // if(nnIndxLU[n + i]>0){		
//   // for(j = 0; j < nnIndxLU[n + i]; j++){
//   // kkk++;
//   // H[i] = 1.0;//*B[nnIndxLU[i] + j]
//   //Rprintf("H: %3.10f \n", H[kkk]); nnIndxLU[n - 1] + m
//   // }
//   // Rprintf("B: %3.10f \n", B[nnIndxLU[i]]);
//   
//   // }else{
//   // H[0] = 0.0;	 
//   // }
//   // }
//   
//   
//   double *KerX = (double *) R_alloc(p*n, sizeof(double));
//   double *Dist = (double *) R_alloc(n, sizeof(double));
//   double *K = (double *) R_alloc(n, sizeof(double));
//   double *KerY = (double *) R_alloc(n, sizeof(double));
//   
//   //double  phi0 = 0.1;//theta[phiIndx];
//   
//   SEXP alpha_r;
//   PROTECT(alpha_r = Rf_allocMatrix(REALSXP, p, n)); nProtect++;
//   //double *Alpha = REAL(alpha_r);
//   double *alpha = (double *) R_alloc(p, sizeof(double));
//   
//   // #ifdef _OPENMP
//   // #pragma omp parallel for private(i, j, alpha)
//   // #endif   
//   for(int s = 0; s < n; s++){
//     //solve distance from testing location to basis table
//     // #ifdef _OPENMP
//     // #pragma omp parallel 
//     // #endif
//     for(i = 0; i < n; i++){
//       if(GeomVariable == 1){
//         Dist[i] = (Z[i]  - Z[s])/h;// pow(Z[i]  - TestZ[s], 2);
//       }else{
//         Dist[i] = dist2(coords[s], coords[n + s], coords[i], coords[n + i]);///theta  
//       }
//       K[i] = exp(-abs(Dist[i]))/h;  
//     }
//     
//     for(i = 0; i < n; i++){
//       for(j = 0; j < p; j++){
//         if(j==0){
//           //X[n*j + i] = 1;
//           KerX[n*j + i] = pow(K[i], 0.5);
//         }else{
//           //X[n*j + i] = Dist[i]/phi;  
//           KerX[n*j + i] = pow(K[i], 0.5)*Dist[i];
//         }
//       }
//       KerY[i] = pow(K[i], 0.5)*y[i];
//       // Rprintf("alpha[i]: %3.10f \n", KerX[i]);
//     }
//     
//     BetaEst(alpha, B, F, KerX, KerY, nnIndx, nnIndxLU, n, p);
//     
//     F77_NAME(dcopy)(&p, alpha, &inc, &REAL(alpha_r)[s*p], &inc);
//     // for(int k1 = 0; k1 < p; k1++){
//     // Alpha[s*p + k1] = alpha[k1];
//     // Rprintf("B: %3.10f \n", Alpha[s*p + k1]);
//     // Rprintf("alpha[k1]: %3.10f \n", alpha[k1]);
//     // }
//   }
//   
//   SEXP Bmat_r, Fmat_r, Qmat_r;
//   PROTECT(Bmat_r = Rf_allocVector(REALSXP, n*n)); nProtect++;
//   PROTECT(Fmat_r = Rf_allocVector(REALSXP, n*n)); nProtect++;
//   PROTECT(Qmat_r = Rf_allocVector(REALSXP, n*n)); nProtect++;
//   double *Bmat = REAL(Bmat_r);zeros(Bmat, n*n);
//   double *Fmat = REAL(Fmat_r);zeros(Fmat, n*n);
//   double *Qmat = REAL(Qmat_r);zeros(Qmat, n*n);
//   // double *Bmat = (double *) R_alloc(n*n, sizeof(double));zeros(Bmat, n*n);
//   // double *Fmat = (double *) R_alloc(n*n, sizeof(double));zeros(Bmat, n*n);
//   
//   for(i = 0; i < n; i++){
//     if(nnIndxLU[n + i]>0){		
//       for(j = 0; j < nnIndxLU[n + i]; j++){
//         
//         Bmat[i + nnIndx[nnIndxLU[i] + j]*n] = - B[nnIndxLU[i] + j]; 		
//       }
//     }
//     // if(nnIndxLU[n + 1]==1){
//     // Bmat[1 + nnIndx[nnIndxLU[1]]*n] = B[nnIndxLU[1]]; 
//     // Rprintf("Bmat: %i., %3.5f. \n", nnIndxLU[1], B[nnIndxLU[1]]);		
//     // }
//     for(j = i; j < i + 1; j++){
//       Fmat[j + n*i] = 1/F[i]; 
//       Bmat[j + n*i] = 1;		
//     }
//   }
//   //Rprintf("Bmat: %i. \n", kkk);
//   char Trans = 'T';
//   char TransN = 'N';
//   double *Q_temp = (double *) R_alloc(n*n, sizeof(double));
//   F77_NAME(dgemm)(&Trans, &TransN, &n, &n, &n, &one, Bmat, &n, Fmat, &n, &zero, Q_temp,&n);
//   F77_NAME(dgemm)(&TransN, &TransN, &n, &n, &n, &one, Q_temp, &n, Bmat, &n, &zero, Qmat,&n);
//   //F77_NAME(dcopy)(&p, alpha, &inc, &REAL(b_r), &inc);
//   
//   // double *XCondExpZ = (double *) R_alloc(n*p, sizeof(double));zeros(XCondExpZ, n*p);
//   // double DiagKer_Sum, Cxz;
//   
//   // int Px;
//   // int Xp;
//   // double *X = (double *) R_alloc(n*p, sizeof(double));
//   ////double *tildX = (double *) R_alloc(n*p, sizeof(double));
//   // for(int s = 0; s < n; s++){
//   // DiagKer_Sum = 0.0
//   // for(i = 0; i < n; i++){
//   // Dist[i] = dist2(coords[s], coords[n + s], coords[i], coords[n + i]);///theta  
//   // K[i] = exp(-Dist[i]/phi)/phi;  
//   // }
//   // for(Xp = 0; Xp < Px; Xp++){
//   // Cxz = 0.0;
//   // for(i = 0; i < n; i++){
//   // for(j = i; j < i + 1; j++){
//   // DiagKer_Sum += Qmat[j + n*i]*K[i];
//   // Cxz += X[i + n*Xp]*Qmat[j + n*i]*K[i];
//   // }			
//   // }
//   // XCondExpZ[s + n*Xp] = Cxz;
//   // }
//   // for(i = 0; i < n; i++){
//   // for(Xp = 0; Xp < Px; Xp++){
//   // XCondExpZ[i + n*Xp] =  XCondExpZ[i + n*Xp]/DiagKer_Sum;	
//   // }
//   // }
//   
//   
//   // }
//   
//   
//   // PutRNGstate();
//   //make return object
//   SEXP result_r, resultName_r;
//   int nResultListObjs = 5;
//   PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
//   PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
//   
//   SET_VECTOR_ELT(result_r, 0, alpha_r);
//   SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("alpha"));
//   SET_VECTOR_ELT(result_r, 1, Bmat_r);
//   SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("B"));
//   SET_VECTOR_ELT(result_r, 2, Fmat_r);
//   SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("D"));
//   
//   SET_VECTOR_ELT(result_r, 3, Qmat_r);
//   SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("Q"));
//   SET_VECTOR_ELT(result_r, 4, B_r);
//   SET_VECTOR_ELT(resultName_r, 4, Rf_mkChar("Biag"));
//   SET_VECTOR_ELT(result_r, 4, F_r);
//   SET_VECTOR_ELT(resultName_r, 4, Rf_mkChar("F"));
//   
//   Rf_namesgets(result_r, resultName_r);
//   
//   //unprotect
//   UNPROTECT(nProtect);
//   // SEXP out;
//   // PROTECT(out = Rf_allocVector(REALSXP, nIndx));
//   // B = REAL(out);
//   
//   return(result_r);
// }
// 
// // [[Rcpp::export]]
// SEXP semiPredict1(SEXP y_r, 
//                   SEXP Z_r, 
//                   SEXP TestZ_r, 
//                   SEXP n_r, 
//                   SEXP nTest_r, 
//                   SEXP m_r,
//                   SEXP h_r,	
//                   SEXP coords_r, 
//                   SEXP TestCoords_r, 
//                   SEXP GeomVariable_r,
//                   SEXP covModel_r, 
//                   SEXP nnIndx_r, 
//                   SEXP nnIndxLU_r,
//                   SEXP sigmaSq_r, 
//                   SEXP tauSq_r,
//                   SEXP phi_r, 
//                   SEXP nu_r, 
//                   SEXP nThreads_r){
//   // SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r, SEXP TestCoords_r,
//   int i, j, nProtect=0;
//   int p = 2;  
//   const int inc = 1;
//   //const double zero = 0.0;
//   
//   //get args
//   int n = INTEGER(n_r)[0];
//   int nTest = INTEGER(nTest_r)[0];
//   int m = INTEGER(m_r)[0];
//   double h = REAL(h_r)[0];
//   int GeomVariable = INTEGER(GeomVariable_r)[0];
//   double *y = REAL(y_r);
//   // double *x = REAL(x_r);
//   
//   double *Z = REAL(Z_r);
//   double *TestZ = REAL(TestZ_r);
//   
//   double *coords = REAL(coords_r);
//   double *TestCoords = REAL(TestCoords_r);
//   
//   int *nnIndx = INTEGER(nnIndx_r);
//   int *nnIndxLU = INTEGER(nnIndxLU_r);
//   
//   
//   int covModel = INTEGER(covModel_r)[0];
//   std::string corName = getCorName(covModel);
//   
//   
//   int nThreads = INTEGER(nThreads_r)[0];
//   
//   
// #ifdef _OPENMP
//   omp_set_num_threads(nThreads);
// #else
//   if(nThreads > 1){
//     warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
//     nThreads = 1;
//   }
// #endif
//   
//   //parameters
//   int nTheta, sigmaSqIndx, tauSqIndx, phiIndx, nuIndx;
//   
//   if(corName != "matern"){
//     nTheta = 3;//sigma^2, tau^2, phi
//     sigmaSqIndx = 0; tauSqIndx = 1; phiIndx = 2;
//   }else{
//     nTheta = 4;//sigma^2, tau^2, phi, nu,
//     sigmaSqIndx = 0; tauSqIndx = 1; phiIndx = 2; nuIndx = 3;
//   }
//   
//   //starting	
//   
//   double *theta = (double *) R_alloc(nTheta, sizeof(double));
//   
//   theta[sigmaSqIndx] = REAL(sigmaSq_r)[0];
//   theta[tauSqIndx] = REAL(tauSq_r)[0];
//   theta[phiIndx] = REAL(phi_r)[0];
//   if(corName == "matern"){
//     theta[nuIndx] = REAL(nu_r)[0];
//   }
//   
//   //allocate for the U index vector that keep track of which locations have the i-th location as a neighbor
//   int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
//   
//   //other stuff
//   int mm = m*m;
//   
//   
//   double *B = (double *) R_alloc(nIndx, sizeof(double));
//   double *F = (double *) R_alloc(n, sizeof(double));
//   
//   double *c =(double *) R_alloc(m*nThreads, sizeof(double));
//   double *C = (double *) R_alloc(mm*nThreads, sizeof(double));
//   int nuUnifb = 0;
//   double *bk = (double *) R_alloc(nThreads*(1.0+static_cast<int>(floor(nuUnifb))), sizeof(double));
//   
//   double nu = 0;
//   if(corName == "matern"){nu = theta[nuIndx];}
//   
//   updateBF1(B, F, c, C, coords, nnIndx, nnIndxLU, n, m, 
//             theta[sigmaSqIndx], theta[phiIndx], nu, 
//             covModel, bk, nuUnifb);
//   
//   
//   double *KerZ = (double *) R_alloc(p*n, sizeof(double));
//   double *Dist = (double *) R_alloc(n, sizeof(double));
//   double *K = (double *) R_alloc(n, sizeof(double));
//   double *KerY = (double *) R_alloc(n, sizeof(double));
//   
//   //double  phi0 = 0.1;//theta[phiIndx];
//   SEXP alpha_r;
//   PROTECT(alpha_r = Rf_allocMatrix(REALSXP, p, nTest)); nProtect++;
//   //double *Alpha = REAL(alpha_r);
//   double *alpha = (double *) R_alloc(p, sizeof(double));
//   // #ifdef _OPENMP
//   // #pragma omp parallel for private(i, j, alpha)
//   // #endif 
//   
//   
//   for(int s = 0; s < nTest; s++){
//     //solve distance from testing location to basis table
//     // #ifdef _OPENMP
//     // #pragma omp parallel 
//     // #endif
//     for(i = 0; i < n; i++){
//       if(GeomVariable == 1){
//         Dist[i] = (Z[i]  - TestZ[s])/h;// pow(Z[i]  - TestZ[s], 2);
//       }else{
//         Dist[i] = dist2(TestCoords[s], TestCoords[nTest + s], coords[i], coords[n + i]);///theta  
//       }
//       K[i] = exp(-abs(Dist[i]))/h;  
//     }
//     
//     for(i = 0; i < n; i++){
//       for(j = 0; j < p; j++){
//         if(j==0){
//           KerZ[n*j + i] = pow(K[i], 0.5);
//         }else{
//           KerZ[n*j + i] = pow(K[i], 0.5)*Dist[i];
//         }
//       }
//       KerY[i] = pow(K[i], 0.5)*y[i];
//     }
//     
//     BetaEst(alpha, B, F, KerZ, KerY, nnIndx, nnIndxLU, n, p);
//     F77_NAME(dcopy)(&p, alpha, &inc, &REAL(alpha_r)[s*p], &inc);
//   }
//   
//   // SEXP Bmat_r, Fmat_r, Qmat_r;
//   // PROTECT(Bmat_r = Rf_allocVector(REALSXP, n*n)); nProtect++;
//   // PROTECT(Fmat_r = Rf_allocVector(REALSXP, n*n)); nProtect++;
//   // PROTECT(Qmat_r = Rf_allocVector(REALSXP, n*n)); nProtect++;
//   // double *Bmat = REAL(Bmat_r);zeros(Bmat, n*n);
//   // double *Fmat = REAL(Fmat_r);zeros(Fmat, n*n);
//   // double *Qmat = REAL(Qmat_r);zeros(Qmat, n*n);
//   //double *Bmat = (double *) R_alloc(n*n, sizeof(double));zeros(Bmat, n*n);
//   //double *Fmat = (double *) R_alloc(n*n, sizeof(double));zeros(Bmat, n*n);
//   // for(i = 0; i < n; i++){
//   // if(nnIndxLU[n + i]>0){		
//   // for(j = 0; j < nnIndxLU[n + i]; j++){
//   // Bmat[i + nnIndx[nnIndxLU[i] + j]*n] = B[nnIndxLU[i] + j]; 		
//   // }
//   // }
//   // for(j = i; j < i + 1; j++){
//   // Fmat[j + n*i] = 1/F[i]; 
//   // Bmat[j + n*i] = 1;		
//   // }
//   // }
//   // char Trans = 'T';
//   // char TransN = 'N';
//   // double *Q_temp = (double *) R_alloc(n*n, sizeof(double));
//   // F77_NAME(dgemm)(&Trans, &TransN, &n, &n, &n, &one, Bmat, &n, Fmat, &n, &zero, Q_temp,&n);
//   // F77_NAME(dgemm)(&TransN, &TransN, &n, &n, &n, &one, Q_temp, &n, Bmat, &n, &zero, Qmat,&n);
//   //F77_NAME(dcopy)(&p, alpha, &inc, &REAL(b_r), &inc);
//   
//   // double *XCondExpZ = (double *) R_alloc(n*p, sizeof(double));zeros(XCondExpZ, n*p);
//   // double DiagKer_Sum, Cxz;
//   
//   // int Px;
//   // int Xp;
//   // double *X = (double *) R_alloc(n*p, sizeof(double));
//   ////double *tildX = (double *) R_alloc(n*p, sizeof(double));
//   // for(int s = 0; s < n; s++){
//   // DiagKer_Sum = 0.0
//   // for(i = 0; i < n; i++){
//   // Dist[i] = dist2(coords[s], coords[n + s], coords[i], coords[n + i]);///theta  
//   // K[i] = exp(-Dist[i]/phi)/phi;  
//   // }
//   // for(Xp = 0; Xp < Px; Xp++){
//   // Cxz = 0.0;
//   // for(i = 0; i < n; i++){
//   // for(j = i; j < i + 1; j++){
//   // DiagKer_Sum += Qmat[j + n*i]*K[i];
//   // Cxz += X[i + n*Xp]*Qmat[j + n*i]*K[i];
//   // }			
//   // }
//   // XCondExpZ[s + n*Xp] = Cxz;
//   // }
//   // for(i = 0; i < n; i++){
//   // for(Xp = 0; Xp < Px; Xp++){
//   // XCondExpZ[i + n*Xp] =  XCondExpZ[i + n*Xp]/DiagKer_Sum;	
//   // }
//   // }
//   
//   
//   // }
//   
//   //make return object
//   SEXP result_r, resultName_r;
//   int nResultListObjs = 1;
//   PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
//   PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
//   
//   SET_VECTOR_ELT(result_r, 0, alpha_r);
//   SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("alpha"));
//   Rf_namesgets(result_r, resultName_r);
//   
//   //unprotect
//   UNPROTECT(nProtect);
//   // SEXP out;
//   // PROTECT(out = Rf_allocVector(REALSXP, nIndx));
//   // B = REAL(out);
//   
//   return(result_r);
// }
