// Rcpp functions for makeNetstage

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix makeNetstageC(IntegerMatrix netstage,
	IntegerVector recvec,
	IntegerVector peervec) {

	int nrec = netstage.nrow() - 1;
	int npeer = netstage.ncol() - 3;

	IntegerMatrix netstageFinal = netstage;

	// Add rec-peer matches
	for (int ind = 0; ind < recvec.size(); ++ind) {
		netstageFinal(recvec[ind]-1, peervec[ind]-1) = 1;
	}

	// Deal with recruiter self-matches
	for (int i = 0; i < nrec; ++i) {
		int rs = sum(netstageFinal(i, _));
		if (rs < 3) {
			netstageFinal(i, npeer) = 1;
		}
		if (rs < 2) {
		  netstageFinal(i, npeer + 1) = 1;
		}
		if (rs < 1) {
		  netstageFinal(i, npeer + 2) = 1;
		}

	}

	// Deal with peer self-matches
	for (int j = 0; j < npeer; ++j) {
		int cs = sum(netstageFinal(_, j));
		if (cs < 1) {
			netstageFinal(nrec, j) = 1;
		}
	}

	return (netstageFinal);

}
