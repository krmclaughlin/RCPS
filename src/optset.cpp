// C++ functions for optset

#include <Rcpp.h>
using namespace Rcpp;

// My version of the "which" function in R
// [[Rcpp::export]]
IntegerVector whichC(NumericVector x,
                     int val) {

  int nmatch = 0;
  for (int i = 0; i < x.size(); ++i) {
    if (x[i] == val) {
      nmatch += 1;
    }
  }

  IntegerVector outvec(nmatch);
  int ctr = 0;
  for (int i = 0; i < x.size(); ++i) {
    if (x[i] == val) {
      outvec[ctr] = i+1;
      ctr += 1;
    }
  }

  return(outvec);

}

// Return set of peers that prefer i over own match
// [[Rcpp::export]]
NumericVector optsetU(int i,
                      int j,
                      IntegerMatrix adjmat,
                      NumericMatrix U,
                      NumericMatrix V) {

  int npeer = V.nrow();

  Function whichC("whichC");

  IntegerVector matchesofi = whichC(adjmat(i-1, _), 1);
  IntegerVector peerset = seq_len(npeer);

  IntegerMatrix adjmrest = adjmat(Range(0,adjmat.nrow()), Range(0, npeer-1) );
  IntegerVector matchesofireal = whichC(adjmrest(i-1, _), 1);

  int nmatchireal = matchesofireal.size();

  int len = npeer - nmatchireal;

  if (len > 0) {

    IntegerVector hlist(len);

    // Vector of peers h not matched to i, h!=0
    int ctr = 0;
    for (int ind = 0; ind < npeer; ++ind) {

      if (is_false(any(matchesofi == peerset[ind]))) {
        hlist[ctr] = ind;
        ctr += 1;
      }

    }

    // Vector of eligible peers that prefer i to own match
    IntegerVector eligh(len);
    for (int ind = 0; ind < len; ++ind) {

      int h = hlist[ind];
      int rh = as<int>(whichC(adjmat(_, h), 1));

      double c1 = V(h, i-1);
      double c2 = V(h, rh-1);

      if (c1 > c2) {
        eligh[ind] = 1;
      } else {
        eligh[ind] = 0;
      }

    }

    IntegerVector chg = whichC(eligh, 1);

    NumericVector Uset(chg.size());
    for (int ind = 0; ind < chg.size(); ++ind) {

      int chgh = hlist[chg[ind] - 1];
      Uset[ind] = U(i-1, chgh);

    }

    // Account for case where no peers prefer i to match
    int Usetsize = Uset.size();
    if (Usetsize == 0) {
      Usetsize = 1;
    }

    NumericVector UsetF(Usetsize);
    if (Uset.size() == 0) {
      UsetF[0] = NA_REAL;
    } else {
      UsetF = Uset;
    }

    return(UsetF);

  } else {
    NumericVector UsetF(1);
    UsetF[0] = NA_REAL;

    return(UsetF);
  }

}

// Return set of recruiters that prefer j over own match
// [[Rcpp::export]]
NumericVector optsetV(int i,
                      int j,
                      IntegerMatrix adjmat,
                      NumericMatrix U,
                      NumericMatrix V) {

  int nrec = U.nrow();

  Function whichC("whichC");

  int matchofj = as<int>(whichC(adjmat(_, j-1), 1));
  IntegerVector recset = seq_len(nrec);

  int len = 0;
  if (matchofj == nrec+1) {
    len = nrec;
  } else {
    len = nrec - 1;
  }

  if (len > 0) {

    IntegerVector llist(len);

    // Vector of recruiters l not matched to j, l!=0
    int ctr = 0;
    for (int ind = 0; ind < nrec; ++ind) {

      if (matchofj != recset[ind]) {
        llist[ctr] = ind;
        ctr += 1;
      }

    }

    // Vector of eligible recruiters that prefer j to own match
    IntegerVector eligl(len);
    for (int ind = 0; ind < len; ++ind) {

      int l = llist[ind];
      IntegerVector pl = whichC(adjmat(l, _), 1);

      NumericVector Ulpl(pl.size());
      for (int ind1 = 0; ind1 < pl.size(); ++ind1) {
        Ulpl[ind1] = U(l, pl[ind1]-1);
      }

      double minUlpl = min(Ulpl);

      if (U(l, j-1) > minUlpl) {
        eligl[ind] = 1;
      } else {
        eligl[ind] = 0;
      }

    }

    IntegerVector chg = whichC(eligl, 1);

    NumericVector Vset(chg.size());
    for (int ind = 0; ind < chg.size(); ++ind) {

      int chgl = llist[chg[ind] - 1];
      Vset[ind] = V(j-1, chgl);

    }

    // Account for case where no recruiters prefer j to match
    int Vsetsize = Vset.size();
    if (Vsetsize == 0) {
      Vsetsize = 1;
    }

    NumericVector VsetF(Vsetsize);
    if (Vset.size() == 0) {
      VsetF = NA_REAL;
    } else {
      VsetF = Vset;
    }

    return(VsetF);

  } else {
    NumericVector VsetF(1);
    VsetF[0] = NA_REAL;

    return(VsetF);
  }

}
