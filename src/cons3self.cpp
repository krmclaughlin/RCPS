/***
 C++ functions for constraints of the 7 types
*/

#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Concat function for two vectors
// [[Rcpp::export]]
NumericVector concat(NumericVector vec1,
                     NumericVector vec2) {

  int n = vec1.size() + vec2.size();
  NumericVector cvec(n);
  for (int ind = 0; ind < vec1.size(); ++ ind) {
    cvec[ind] = vec1[ind];
  }
  int ctr = 0;
  for (int ind = vec1.size(); ind < n; ++ind) {
    cvec[ind] = vec2[ctr];
    ctr += 1;
  }

  return(cvec);

}

/***
 TYPE 1: i>0; j>0; A[i,j]=1
 */

// [[Rcpp::export]]
NumericVector Uab1(int i,
                   int j,
                   IntegerMatrix adjmat,
                   NumericMatrix U,
                   NumericMatrix V) {

  int npeer = V.nrow();

  Function whichC("whichC");

  IntegerMatrix adjmrest = adjmat(Range(0, adjmat.nrow()), Range(npeer, npeer+2) );

  IntegerVector selfm1 = whichC(adjmrest(i-1, _), 0);
  NumericVector Uselfms(selfm1.size());

  for (int nse = 0; nse < selfm1.size(); nse++) {
    Uselfms[nse] = U(i-1, npeer + selfm1[nse] - 1);
  }

  Function optsetU("optsetU");
  NumericVector Uset = optsetU(i, j, adjmat, U, V);

  Function concat("concat");
  NumericVector cvec = as<NumericVector>(concat(Uselfms, Uset));

  NumericVector cvecnoNA = na_omit(cvec);

  double a = max(cvecnoNA);
  double b = arma::math::inf();

  NumericVector bounds = concat(a, b);

  return(bounds);

}

// [[Rcpp::export]]
NumericVector Vab1(int i,
                   int j,
                   IntegerMatrix adjmat,
                   NumericMatrix U,
                   NumericMatrix V) {

  int nrec = U.nrow();

  Function optsetV("optsetV");

  double Vself = V(j-1, nrec);
  NumericVector Vset = optsetV(i, j, adjmat, U, V);

  Function concat("concat");
  NumericVector cvec = concat(Vself, Vset);

  NumericVector cvecnoNA = na_omit(cvec);

  double a = max(cvecnoNA);
  double b = arma::math::inf();

  NumericVector bounds = concat(a, b);

  return(bounds);

}

/***
TYPE 2: i>0; j>0; A[i,j]=0
*/

// [[Rcpp::export]]
NumericVector Uab2(int i,
                   int j,
                   IntegerMatrix adjmat,
                   NumericMatrix U,
                   NumericMatrix V) {

  Function whichC("whichC");

  IntegerVector matchesofi = whichC(adjmat(i-1, _), 1);
  int nmatchi = matchesofi.size();
  NumericVector Umi(nmatchi);
  for (int ind = 0; ind < nmatchi; ++ind) {
    Umi[ind] = U(i-1, matchesofi[ind]-1);
  }

  int matchofj = as<int>(whichC(adjmat(_, j-1), 1));

  double a = -1 * arma::math::inf();
  double b = 0.0;

  if (V(j-1, i-1) > V(j-1, matchofj-1)) {
    b = min(Umi);
  } else {
    b = arma::math::inf();
  }

  Function concat("concat");

  NumericVector bounds = concat(a, b);

  return(bounds);

}

// [[Rcpp::export]]
NumericVector Vab2(int i,
                   int j,
                   IntegerMatrix adjmat,
                   NumericMatrix U,
                   NumericMatrix V) {

  Function whichC("whichC");

  IntegerVector matchesofi = whichC(adjmat(i-1, _), 1);
  int nmatchi = matchesofi.size();
  NumericVector Umi(nmatchi);
  for (int ind = 0; ind < nmatchi; ++ind) {
    Umi[ind] = U(i-1, matchesofi[ind]-1);
  }

  int matchofj = as<int>(whichC(adjmat(_, j-1), 1));

  double a = -1 * arma::math::inf();
  double b = 0.0;

  if (U(i-1, j-1) > min(Umi)) {
    b = V(j-1, matchofj-1);
  } else {
    b = arma::math::inf();
  }

  Function concat("concat");

  NumericVector bounds = concat(a, b);

  return(bounds);

}

/***
TYPE 3: i>0; j=0; A[i,j]=1
*/

// [[Rcpp::export]]
NumericVector Uab3(int i,
                   int j,
                   IntegerMatrix adjmat,
                   NumericMatrix U,
                   NumericMatrix V) {

  int npeer = V.nrow();

  Function whichC("whichC");

  IntegerMatrix adjmrest = adjmat(Range(0, adjmat.nrow()), Range(npeer, npeer+2) );

  IntegerVector selfm1 = whichC(adjmrest(i-1, _), 0);

  Function optsetU("optsetU");
  NumericVector Uset = optsetU(i, j, adjmat, U, V);

  NumericVector cvec(Uset.size() + selfm1.size());
  if (selfm1.size() > 0) {

    NumericVector Uselfms(selfm1.size());
    for (int nse = 0; nse < selfm1.size(); nse++) {
      Uselfms[nse] = U(i-1, npeer + selfm1[nse] - 1);
    }

    Function concat("concat");
    cvec = as<NumericVector>(concat(Uselfms, Uset));

  } else {

    cvec = Uset;

  }

  NumericVector cvecnoNA = na_omit(cvec);

  double a = 0.0;

  if (cvecnoNA.size() == 0) {
    a = -1 * arma::math::inf();
  } else {
    a = max(cvecnoNA);
  }

  double b = arma::math::inf();

  NumericVector bounds(2);
  bounds[0] = a;
  bounds[1] = b;

  return(bounds);

}

/***
TYPE 4: i>0; j=0; A[i,j]=0
*/

// [[Rcpp::export]]
NumericVector Uab4(int i,
                   int j,
                   IntegerMatrix adjmat,
                   NumericMatrix U,
                   NumericMatrix V) {

  Function whichC("whichC");

  IntegerVector matchesofi = whichC(adjmat(i-1, _), 1);
  int nmi = matchesofi.size();
  NumericVector Umi(nmi);
  for (int ind = 0; ind < nmi; ++ind) {
    Umi[ind] = U(i-1, matchesofi[ind]-1);
  }

  double a = -1 * arma::math::inf();
  double b = min(Umi);
  // i does not select self, so Umi is not empty

  Function concat("concat");

  NumericVector bounds = concat(a, b);

  return(bounds);

}

/***
TYPE 5: i=0; j>0; A[i,j]=1
*/

// [[Rcpp::export]]
NumericVector Vab5(int i,
                   int j,
                   IntegerMatrix adjmat,
                   NumericMatrix U,
                   NumericMatrix V) {

  Function optsetV("optsetV");

  NumericVector Vset = optsetV(i, j, adjmat, U, V);

  NumericVector VsetnoNA = na_omit(Vset);
  double a;

  if (VsetnoNA.size() == 0) {
    a = -1 * arma::math::inf();
  } else {
    a = max(VsetnoNA);
  }

  double b = arma::math::inf();

  Function concat("concat");

  NumericVector bounds = concat(a, b);

  return(bounds);

}

/***
TYPE 6: i=0; j>0; A[i,j]=0
*/

// [[Rcpp::export]]
NumericVector Vab6(int i,
                   int j,
                   IntegerMatrix adjmat,
                   NumericMatrix U,
                   NumericMatrix V) {

  Function whichC("whichC");

  int matchofj = as<int>(whichC(adjmat(_, j-1), 1));

  double a = -1 * arma::math::inf();
  double b = V(j-1, matchofj-1);

  Function concat("concat");

  NumericVector bounds = concat(a, b);

  return(bounds);

}
