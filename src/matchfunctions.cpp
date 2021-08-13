//
//  Functions to carry out parts of the matching algorithm
//


#include <Rcpp.h>
using namespace Rcpp;

//' Internal (cpp) of matching function: peer choices.
//'
//' For each peer, returns their choice from amongst the recruiters
//'   that selected them at this stage and themself. Will always return
//'   a value (self), even if not selected.
//'
//' @param recsels List with the (up to 3) selections of each recruiter.
//' @param Vrank matrix of ranks from each peer to each recruiter and self.
//'
//' @return A vector with the index of each peer's choice.
// [[Rcpp::export]]
IntegerVector peerChoice(List recsels,
                         IntegerMatrix Vrank) {

    int nrec = recsels.size();
    int npeer = Vrank.nrow();

    IntegerVector jchoice(npeer);

    for (int j = 0; j < npeer; ++j) {

        // Recruiters that selected each peer
        IntegerVector whichrecs(nrec + 1);

        for (int i = 0; i < nrec; ++i) {

            IntegerVector tvec = as<IntegerVector>(recsels[i]);
            if (is_true(any(tvec == j + 1))) {
                whichrecs[i] = i + 1;
            }

        }

        // Always have option of choosing self
        whichrecs[nrec] = nrec + 1;

        // Get each peer's choice
        IntegerVector peerfeels(nrec + 1);
        for (int i = 0; i < nrec + 1; ++i) {

            if (whichrecs[i] != 0) {
                peerfeels[i] = Vrank(j, whichrecs[i]-1);
            }

        }

        int maxj = which_max(peerfeels) + 1;
        jchoice[j] = maxj;

    }

    return(jchoice);

}

//' Internal (cpp) of matching function: housekeeping.
//'
//' Function that carries out list reassignment: takes the vector of peer
//'   choices at a given stage of matching and makes a new selection list.
//'
//' @param nrec Number of recruiters.
//' @param npeer Number of peers.
//' @param jchoice Vector of length \code{npeer} with each peer's selection
//'   (can be self).
//' @param RS List.
//'
//' @return A list with the new recruiter selections.
// [[Rcpp::export]]
List listReassign(int nrec,
                  int npeer,
                  IntegerVector jchoice,
                  List RS) {

    Function whichCone("whichCone");

    List selsnew(nrec);
    for (int i = 0; i < nrec; ++i) {

        IntegerVector ips(3);
        int ctr = 0;
        for (int j = 0; j < npeer; ++j) {
            if (jchoice[j] == i + 1) {
                ips[ctr] = j + 1;
                ctr += 1;
            }
        }

        IntegerVector RSrecs = as<IntegerVector>(RS[i]);
        int s1l = as<int>(whichCone(RSrecs, npeer+1));
        int s2l = as<int>(whichCone(RSrecs, npeer+2));
        int s3l = as<int>(whichCone(RSrecs, npeer+3));

        if (s1l > 0) {
          ips[ctr] = npeer + 1;
          ctr += 1;
        }
        if (s2l > 0) {
          ips[ctr] = npeer + 2;
          ctr += 1;
        }
        if (s3l > 0) {
          ips[ctr] = npeer + 3;
          ctr += 1;
        }

        selsnew[i] = ips;
    }

    return(selsnew);

}

// Function that identifies the place of the first 0 in a vector
// NOTE: R-position, so 0 indicates no zeros
//' Internal (cpp) of matching function: finds first zero in a vector.
//'
//' @param x An integer vector.
//'
//' @return Position in vector (R-vector; cannot be 0) where first zero occurs.
//'   Returns 0 if there are no zeros.
// [[Rcpp::export]]
int firstZero(IntegerVector x) {

    int firstzero = 0;

    if (is_true(any(x == 0))) {

        int i = 0;

        while (firstzero == 0) {
            if (x[i] == 0) {
                firstzero = i + 1;
            } else {
                i += 1;
            }
        }

    }

    return(firstzero);

}
