/***
C++ functions for get_id_type_ur.R
*/

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int getType(int i,
	int j,
	IntegerMatrix adjmat,
	NumericMatrix U,
	NumericMatrix V) {

	int nrec = U.nrow();
	int npeer = V.nrow();

	int idtype = 0;

	if (i < nrec + 1) {

		if (j < npeer + 1) {

			if (adjmat(i-1, j-1) == 1) {
				idtype = 1;
			} else {
				idtype = 2;
			}

		} else if (adjmat(i-1, j-1) == 1) {
			idtype = 3;
		} else {
			idtype = 4;
		}

	} else if (j < npeer + 1) {

		if (adjmat(i-1, j-1) == 1) {
			idtype = 5;
		} else {
			idtype = 6;
		}

	} else {
		idtype = 7;
	}

	return(idtype);

}

// [[Rcpp::export]]
double newUijC(int i,
	int j,
	IntegerMatrix adjmat,
	NumericMatrix U,
	NumericMatrix V,
	NumericVector alpha,
	NumericVector Xij,
	Function rtruncnorm) {

	Function getType("getType");
  int nchar = alpha.size();

	int type = as<int>(getType(i, j, adjmat, U, V));		//cpp function
	double tmean = 0;
	for (int cc = 0; cc < nchar; ++cc) {
	  tmean += alpha[cc] * Xij[cc];
	}

	NumericVector newUbounds(2);

	if (type == 1) {

		Function Uab1("Uab1");
		newUbounds = Uab1(i, j, adjmat, U, V);

	}

	if (type == 2) {

		Function Uab2("Uab2");
		newUbounds = Uab2(i, j, adjmat, U, V);

	}

	if (type == 3) {

		Function Uab3("Uab3");
		newUbounds = Uab3(i, j, adjmat, U, V);

	}

	if (type == 4) {

    Function Uab4("Uab4");
    newUbounds = Uab4(i, j, adjmat, U, V);

	}

	IntegerVector validtypes = seq_len(4);
	double newUij;

	if (is_true(any(validtypes == type))) {

		double lb = newUbounds[0];
		double ub = newUbounds[1];

		// Fix problem with unordered bounds
		if (lb > ub) {

		  Rprintf("Warning: lb > ub\n");
		  newUij = (lb + ub) / 2;

		} else {

  		newUij = as<double>(rtruncnorm(1,
  			lb,
  			ub,
  			tmean,
  			1.0));

		}

	} else {

		newUij = NA_REAL;

	}

	return(newUij);

}

// [[Rcpp::export]]
double newVjiC(int i,
	int j,
	IntegerMatrix adjmat,
	NumericMatrix U,
	NumericMatrix V,
	NumericVector beta,
	NumericVector Yji,
	Function rtruncnorm) {

	Function getType("getType");
  int nchar = beta.size();
  int nrec = U.nrow();

	int type = as<int>(getType(i, j, adjmat, U, V));
	double tmean = 0;
	for (int cc = 0; cc < nchar; ++cc) {
	  tmean += beta[cc] * Yji[cc];
	}

	NumericVector newVbounds(2);

	if (type == 1) {

		Function Vab1("Vab1");
		newVbounds = Vab1(i, j, adjmat, U, V);

	}

	if (type == 2) {

		Function Vab2("Vab2");
		newVbounds = Vab2(i, j, adjmat, U, V);

	}

	if (type == 5) {

		Function Vab5("Vab5");
		newVbounds = Vab5(i, j, adjmat, U, V);

	}

	if (type == 6) {

	  Function Vab6("Vab6");
	  newVbounds = Vab6(i, j, adjmat, U, V);

	}

	IntegerVector validtypes(4);
	validtypes[0] = 1;
	validtypes[1] = 2;
	validtypes[2] = 5;
	validtypes[3] = 6;
	double newVji;

	if (is_true(any(validtypes == type))) {

		double lb = newVbounds[0];
		double ub = newVbounds[1];

		// Fix problem with unordered bounds
		if (lb > ub) {

		  Rprintf("Warning: lb > ub\n");
		  newVji = (lb + ub) / 2;

		} else {

  		newVji = as<double>(rtruncnorm(1,
  			lb,
  			ub,
  			tmean,
  			1.0));

		}

	} else {

		newVji = NA_REAL;

	}

	return(newVji);

}

