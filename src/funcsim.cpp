//
//Writing the sim.rec function in C++
//

#include <Rcpp.h>
using namespace Rcpp;

//' Internal funcion to construct U given alpha, X, epsilon
//'
//' Note that this is the CALCULATION, only used in step 0
//'
//' @param charrec Matrix of covariate values for recruiters.
//' @param charpeer Matrix of covariate values for peers.
//' @param alpha Preference coefficients for recruiters to peers.
//' @param ksi Parameter for number of recruits. Vector length ncoupon.
//' @param epsilon Matrix with random error.  Must have same dimension as U.
//'
//' @return Matrix with calculated values of \eqn{U}, where \eqn{U =\alpha X + \epsilon}
// [[Rcpp::export]]
NumericMatrix makeU(NumericMatrix charrec,
                    NumericMatrix charpeer,
                    NumericVector alpha,
                    NumericVector ksi,
                    NumericMatrix epsilon) {

  int nr = charrec.nrow();
  int np = charpeer.nrow();
  int nc = charrec.ncol();

  NumericMatrix U(nr, np+3);

  for (int i = 0; i < nr; ++i) {
    NumericMatrix Xi(nc, np);

    for (int j = 0; j < np; ++j) {
      for (int cc = 0; cc < nc; ++cc) {
        if (charrec(i, cc) == charpeer(j, cc)) {
          Xi(cc, j) = 1.0;
        } else {
          Xi(cc, j) = 0.0;
        }
      }
    }

    for (int j = 0; j < np; ++j) {
      for (int cc = 0; cc < nc; ++cc) {
        U(i, j) += alpha[cc] * Xi(cc, j);
      }
      U(i, j) += epsilon(i, j);
    }

    U(i, np) = ksi[0] + epsilon(i, np);
    U(i, np+1) = ksi[1] + epsilon(i, np+1);
    U(i, np+2) = ksi[2] + epsilon(i, np+2);

  }

  return(U);

}

//' Internal funcion to construct V given beta, Y, gamma
//'
//' Note that this is the CALCULATION, only used in step 0
//'
//' @param charrec Matrix of covariate values for recruiters.
//' @param charpeer Matrix of covariate values for peers.
//' @param beta Preference coefficients for peers to recruiters. Last value is self.
//' @param zeta Parameter for self match. Vector length 1.
//' @param gamma Matrix with random error.  Must have same dimension as V.
//'
//' @return Matrix with calculated values of \eqn{V}, where \eqn{V =\beta Y + \gamma}
// [[Rcpp::export]]
NumericMatrix makeV(NumericMatrix charrec,
	NumericMatrix charpeer,
	NumericVector beta,
	NumericVector zeta,
	NumericMatrix gamma) {

  int nr = charrec.nrow();
  int np = charpeer.nrow();
  int nc = charrec.ncol();

	NumericMatrix V(np, nr+1);

	for (int j = 0; j < np; ++j) {
		NumericMatrix Yj(nc, nr);

		for (int i = 0; i < nr; ++i) {
		  for (int cc = 0; cc < nc; ++cc) {
  			if (charpeer(j, cc) == charrec(i, cc)) {
  				Yj(cc, i) = 1.0;
  			} else {
  				Yj(cc, i) = 0.0;
  			}
		  }
		}

		for (int i = 0; i < nr; ++i) {
		  for (int cc = 0; cc < nc; ++cc) {
			  V(j, i) += beta[cc] * Yj(cc, i);
		  }
		  V(j, i) += gamma(j, i);
		}

		V(j, nr) = zeta[0] + gamma(j, nr);

	}

	return(V);

}

// [[Rcpp::export]]
DataFrame makeEL(IntegerMatrix adjmat,
	IntegerVector sindex,
	IntegerVector pindex,
	int w) {

	int nrec = adjmat.nrow();
	int npeer = adjmat.ncol();

	int n = sum(adjmat);

	IntegerVector recers(n);
	IntegerVector recd(n);

	int ctr = 0;

	for (int i = 0; i < nrec; ++i) {

		for (int j = 0; j < npeer; ++j) {

			if (adjmat(i, j) == 1) {

				recers[ctr] = sindex[i];
				recd[ctr] = pindex[j];

				ctr += 1;
			}
		}

	}

	IntegerVector wave(n, w);

	return DataFrame::create(_["recruiter"]=recers,
		_["peer"]=recd,
		_["wave"]=wave);

}
