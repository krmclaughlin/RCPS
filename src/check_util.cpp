//
//Function to check if U,V generate observed recruitment tree
//

#include <Rcpp.h>
using namespace Rcpp;

//' Internal function to check if the observed recruitment tree is returned by any
//'   given Ustar and Vstar matrices.
//'
//' @param Ustar Matrix of utilities from recruiters to peers.
//' @param Vstar Matrix of utilities from peers to recruiters.
//' @param netstage Amended adjacency matrix, with observed recruitment information.
//' @param recfromUV3urR Call to another \code{C++} function.
//' @param reslocal Logical: use additional information about the underlying network?
//'    Makes computation faster. If TRUE, assumes that \code{netstage} takes three values:
//'    1 indicates a recruitment tie (and thus a tie in the underlying network); 0 indicates
//'    no recruitment tie, but a tie does exist in the underlying network; NA indicates no
//'    tie in the underlying network, and thus no possibility of recruitment.
//'
//' @return TRUE or FALSE if the given Ustar and Vstar return the observed recruitment tree.
// [[Rcpp::export]]
bool checkutilC(NumericMatrix Ustar,
	NumericMatrix Vstar,
	IntegerMatrix netstage,
	Function recfromUV3urR,
	bool reslocal) {

	int nrec = Ustar.nrow();
	int npeer = Vstar.nrow();

	IntegerMatrix redadj = netstage(Range(0,nrec-1), Range(0,npeer-1));

	IntegerMatrix genadj = recfromUV3urR(Ustar, Vstar, reslocal);

	int matmatch = 0;

	for (int i = 0; i < nrec; ++i) {

		for (int j = 0; j < npeer; ++j) {

			if (redadj(i, j) == genadj(i, j)) {
				matmatch += 1;
			}

		}

	}

	bool allgood = false;
	if (matmatch == nrec * npeer) {
	  allgood = true;
	}

	return(allgood);

}
