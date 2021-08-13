// Generate adjacency matrix from URANK and VRANK
// C++ version of rec.from.UV.3.ur

#include <Rcpp.h>
using namespace Rcpp;

// My version of the "which" function in R
// Returns only a SINGULAR value (must only have one match)
// [[Rcpp::export]]
int whichCone(NumericVector x,
                     int val) {

    int out = 0;
    for (int i = 0; i < x.size(); ++i) {
        if (x[i] == val) {
            out = i + 1;
        }
    }

    return(out);

}

// My version of the max function that returns max after removing NA values
// This is for integers
// Note that max in Rcpp sugar returns NA
// [[Rcpp::export]]
int maxintNoNA(IntegerVector x) {

  IntegerVector xnoNA = na_omit(x);
  int mx;
  if (xnoNA.size() > 0) {
    mx = max(xnoNA);
  }
  else {
    mx = NA_INTEGER;
  }

  return mx;

}


// [[Rcpp::export]]
List recfromUV3urC(IntegerMatrix Urank,
                            IntegerMatrix Vrank,
                            bool reslocal) {

    int nrec = Urank.nrow();
    int npeer = Vrank.nrow();

    List first3(nrec);
    IntegerVector containsself(nrec);
    IntegerVector minusedid(nrec);

    Function whichCone("whichCone");

    IntegerVector nselfm(nrec);

    // Get first three matches for recruiters (non-self)
    // Set-up step
    for (int i = 0; i < nrec; ++i) {

        // Number of non-NA values in row i
        // (Will be npeer+3 if restrict.local=FALSE)
        Function maxintNoNA("maxintNoNA");
        IntegerVector temprow = Urank(i, _);
        int npotmatch = as<int>(maxintNoNA(temprow));

        // Get first 3 matches (or as many as have peers not NA)
        int ntries;
        if (npotmatch <=3) {
          ntries = npotmatch;
        } else {
          ntries = 3;
        }
        IntegerVector temp(ntries);
        if (npotmatch > 0) {
          for (int npm = 0; npm < ntries; ++npm) {
            temp[npm] = as<int>(whichCone(temprow, npotmatch - npm));
          }
        }

          // Record self match locations and count
          int nselfmi = 0;

          int s1l = as<int>(whichCone(temp, npeer+1));
          if (s1l > 0) {
            nselfmi += 1;
          }
          int s2l = as<int>(whichCone(temp, npeer+2));
          if (s2l > 0) {
            nselfmi += 1;
          }
          int s3l = as<int>(whichCone(temp, npeer+3));
          if (s3l > 0) {
            nselfmi += 1;
          }


          first3(i) = temp;
          nselfm[i] = nselfmi;

        // Record minimum used id
        minusedid[i] = npeer + 1;

    }

    // Peers respond (external function)
    Function peerChoice("peerChoice");

    IntegerVector jchoice = peerChoice(first3, Vrank);

    // List reassignment: list tentative accepts for each recruiter (external function)
    Function listReassign("listReassign");
    List selstemp = listReassign(nrec, npeer, jchoice, first3);


    // Record which recruiters have 3 - n(self so far) peers (temporarily do not recruit)
    IntegerVector hasthree(nrec);

    for (int i = 0; i < nrec; ++i) {

        IntegerVector tv = as<IntegerVector>(selstemp[i]);
        if (is_false(any(tv == 0))) {
            hasthree[i] = 1;
        }

    }

    // Recruiters that do not need to select more this attempt

    Function firstZero("firstZero");

    // Iterate
    while (sum(hasthree) < nrec) {

        List selsnew(nrec);

        // Recruit until select self OR have 3 potential matches
        for (int i = 0; i < nrec; ++i) {

            // Get existing recruits
            IntegerVector rcs = selstemp[i];

            bool addmore = false;
            if (hasthree[i] == 0) {
                addmore = true;
            }

            while (addmore == true) {

                int newid = minusedid[i] - 1;
                int newp = as<int>(whichCone(Urank(i, _), newid));

                // Update minimum used ID
                minusedid[i] -= 1;

                    // If self-match, no more recruits
                if (newp > npeer) {
                    nselfm[i] += 1;
                }

                int fz = as<int>(firstZero(rcs)) - 1;
                rcs[fz] = newp;

                // Check if need to recruit more
                if (is_false(any(rcs == 0))) {
                    addmore = false;
                }

            }

            selsnew[i] = rcs;

        }

        // Peers respond (external function)
        IntegerVector jchoiceup = peerChoice(selsnew, Vrank);

        // List reassignment: list tentative accepts for each recruiter (external function)
        selstemp = listReassign(nrec, npeer, jchoiceup, selsnew);

        // Record which recruiters have 3 peers (temporarily do not recruit)
        for (int i = 0; i < nrec; ++i) {

            IntegerVector tv = selstemp[i];
            if (is_false(any(tv == 0))) {
                hasthree[i] = 1;
            } else {
                hasthree[i] = 0;
            }

        }

    }

    // Create adjacency matrix to return
    IntegerMatrix utiladj(nrec, npeer);

    for (int i = 0; i < nrec; ++i) {
        IntegerVector rcs = as<IntegerVector>(selstemp[i]);

            for (int ind = 0; ind < rcs.size(); ++ind) {

                if (rcs[ind] <= npeer) {
                    utiladj(i, rcs[ind]-1) = 1;
                }

            }

    }

    // If using additional population information, add in NAs
    if (reslocal == true) {

      for (int i = 0; i < nrec; ++ i) {

        for (int j = 0; j < npeer; ++j) {

          IntegerVector tmpval(1);
          tmpval[0] = Urank(i, j);
          LogicalVector tcondv = is_na(tmpval);
          bool tcond = tcondv[0];

          if (tcond == true) {
            utiladj(i, j) = NA_INTEGER;
          }

        }

      }

    }

    return Rcpp::List::create(Rcpp::Named("utiladj") = utiladj,
                              Rcpp::Named("minusedid") = minusedid);

}
