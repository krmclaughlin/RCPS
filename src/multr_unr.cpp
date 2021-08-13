/***
C++ functions for multr_unr.R (cpp version of multround_unranked.R)
*/

# include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]
using namespace Rcpp;

// [[Rcpp::export]]
List makeX(IntegerMatrix charrec,
	IntegerMatrix charpeer,
	IntegerMatrix netstage) {

	int nrec = charrec.nrow();
	int npeer = charpeer.nrow();
	List X(nrec);

	int nchar = charrec.ncol();

	for (int i = 0; i < nrec; ++i) {

    // For self-match
		IntegerMatrix temp(nchar+3, npeer+3);

		for (int j = 0; j < npeer; ++j) {

		  for (int cc = 0; cc < nchar; ++cc) {
  			if (charrec(i, cc) == charpeer(j, cc)) {
  				temp(cc, j) = 1;
  			} else {
  				temp(cc, j) = 0;
  			}
		  }

		}

		// Deal with self-matches
		  temp(nchar, npeer) = 1;
		  temp(nchar+1, npeer+1) = 1;
		  temp(nchar+2, npeer+2) = 1;


		X[i] = temp;

	}

	return(X);

}

// [[Rcpp::export]]
List makeY(IntegerMatrix charrec,
	IntegerMatrix charpeer,
	IntegerMatrix netstage) {

	int nrec = charrec.nrow();
	int npeer = charpeer.nrow();
	List Y(npeer);

	int nchar = charrec.ncol();

	for (int j = 0; j < npeer; ++j) {

	  // For self-match
		IntegerMatrix temp(nchar+1, nrec+1);

		for (int i = 0; i < nrec; ++i) {

		  for (int cc = 0; cc < nchar; ++cc) {
  			if (charpeer(j, cc) == charrec(i, cc)) {
  				temp(cc, i) = 1;
  			} else {
  				temp(cc, i) = 0;
  			}
		  }

		}

		temp(nchar, nrec) = 1;

		Y[j] = temp;

	}

	return(Y);

}

// [[Rcpp::export]]
NumericMatrix calcU(NumericVector alpha,
	List X,
	NumericMatrix epsilon) {

	int nrec = X.size();
	int npeer = as<IntegerMatrix>(X[0]).ncol() - 3;
	int nchar = alpha.size();

	NumericMatrix Ustar(nrec, npeer+3);

	for (int i = 0; i < nrec; ++i) {

		NumericMatrix Xi = as<NumericMatrix>(X[i]);

		for (int j = 0; j < npeer+3; ++j) {

		  for (int cc = 0; cc < nchar; ++cc) {
			  Ustar(i, j) += alpha[cc] * Xi(cc, j);
		  }
		  Ustar(i, j) += epsilon(i, j);

		}

	}

	return(Ustar);

}

// [[Rcpp::export]]
NumericMatrix calcV(NumericVector beta,
	List Y,
	NumericMatrix gamma) {

	int nrec = as<IntegerMatrix>(Y[0]).ncol() - 1;
	int npeer = Y.size();
	int nchar = beta.size();

	NumericMatrix Vstar(npeer, nrec+1);

	for (int j = 0; j < npeer; ++j) {

		NumericMatrix Yj = as<NumericMatrix>(Y[j]);

		for (int i = 0; i < nrec+1; ++i) {

		  for (int cc = 0; cc < nchar; ++cc) {
        Vstar(j, i) += beta[cc] * Yj(cc, i);
		  }
		  Vstar(j, i) += gamma(j, i);

		}

	}

	return(Vstar);

}

// Modified to handle missing values
// [[Rcpp::export]]
double ipNoSelf(NumericVector vec1,
	NumericVector vec2) {

	int n = vec1.size() - 1;
  NumericVector temp(n);

	for (int i = 0; i < n; ++i) {
		temp[i] = vec1[i] * vec2[i];
	}

	NumericVector tempNoNA = na_omit(temp);
	double ip = sum(tempNoNA);

	return(ip);

}

// Modified to handle missing values
// [[Rcpp::export]]
double ipSelf(NumericVector vec1,
                NumericVector vec2) {

  int n = vec1.size();
  NumericVector temp(n);

  for (int i = 0; i < n; ++i) {
    temp[i] = vec1[i] * vec2[i];
  }

  NumericVector tempNoNA = na_omit(temp);
  double ip = sum(tempNoNA);

  return(ip);

}

// [[Rcpp::export]]
NumericVector alphaFromU(List X,
	arma::mat U,
	arma::vec ma0,
	arma::mat sa0,
	Function mvrnorm) {

	int nrec = U.n_rows;
  int nchar = sa0.n_rows;

  arma::mat xx(nchar, nchar);
  xx.zeros();
  arma::mat xu(nchar, 1);
  xu.zeros();

	for (int i = 0; i < nrec; ++i) {

		arma::mat Xi = X[i];
		arma::mat Ui = U.row(i);

		arma::mat xt3 = Xi * Xi.t();
		arma::mat xu3 = Xi * Ui.t();

    xx += xt3;
    xu += xu3;

	}

	arma::mat saInvT = (xx + inv(sa0));
	arma::mat maT = inv(saInvT) * (xu + inv(sa0) * ma0);

	NumericVector alphanew = mvrnorm(1, maT, inv(saInvT));

	return(alphanew);

}

// [[Rcpp::export]]
NumericVector betaFromV(List Y,
                        arma::mat V,
                        arma::vec mb0,
                        arma::mat sb0,
                        Function mvrnorm) {

  int npeer = V.n_rows;
  int nchar = sb0.n_rows;

  arma::mat yy(nchar, nchar);
  yy.zeros();
  arma::mat yv(nchar, 1);
  yv.zeros();

  for (int j = 0; j < npeer; ++j) {

    arma::mat Yj = Y[j];
    arma::mat Vj = V.row(j);

    arma::mat yyt = Yj * Yj.t();
    arma::mat yvt = Yj * Vj.t();

      yy += yyt;
      yv += yvt;

  }

  arma::mat sbInvT = (yy + inv(sb0));
  arma::mat mbT = inv(sbInvT) * (yv + inv(sb0) * mb0);

  NumericVector betanew = mvrnorm(1, mbT, inv(sbInvT));

  return(betanew);

}
