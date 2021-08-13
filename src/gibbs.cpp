//
// C++ functions for Gibbs sampler (combines many other functions)
//

#include <Rcpp.h>
using namespace Rcpp;

//' Internal (cpp) function for the Gibbs sampler iterations.
//'
//' Setup must already have occurred.
//'
//' @param Kcalc The number of iterations to run.  Equal to desired number
//'   of draws times interval, plus the burn-in period.
//' @param alphastart Initialization value for alpha, preference coefficient for recruiters.
//' @param betastart Initialization value for beta, preference coefficient for peers.
//' @param adjmat List of amended adjacency matrices, one for each wave.  Each contains
//'   recruitment information for this wave plus information about self-matches. Assumed to
//'   take two values (0-1) if \code{restrict.local=FALSE}, or three values (0-1-NA) if
//'   \code{restrict.local=TRUE}.
//' @param nrec List of number of recruiters, one value for each wave.
//' @param npeer List of number of peers, one value for each wave.
//' @param X List of list of covariate(s) matches from recruiters to peers.
//' @param Y List of list of covariate(s) matches from peers to recruiters.
//' @param Ustar List of matrices of utilities from recruiters to peers, one for each wave.
//' @param Vstar List of matrices of utilities from peers to recruiters, one for each wave.
//' @param expgrid R function \code{expand.grid} from \code{base}.
//' @param rtruncnorm R function \code{rtruncnorm} from the \code{truncnorm} package.
//' @param sample R function \code{sample} from \code{base}.
//' @param ma0 Prior value \eqn{\mu_\alpha}, mean of normal distribution for alpha.
//' @param sa0 Prior value \eqn{\Sigma_\alpha}, variance of normal distribution for alpha.
//' @param mb0 Prior value \eqn{\mu_\beta}, mean of normal distribution for beta.
//' @param sb0 Prior value \eqn{\Sigma_\beta}, varaince of normal distribution for beta.
//' @param tallowed Number of hours allowed on Hoffman.
//' @param recfromUV3urR R function in this package.
//' @param reslocal Logical: use additional information about the underlying network? Makes
//'   computation faster. If TRUE, assumes that \code{netstage} takes three values: 1 indicates a
//'   recruitment tie (and thus a tie in the underlying network); 0 indicates no recruitment tie,
//'   but a tie does exist in the underlying network; NA indicates no tie in the underlying network,
//'   and thus no possibility of recruitment.
//' @param mvrnorm Function to make multivariate normal draws.
//' @param updlast Logical: flag to stop. Not currently used.
//'
//' @return A list containing
//'   \item{alpha}{Draws of alpha from the posterior distribution.}
//'   \item{beta}{Draws of beta from the posterior distribution.}
//'   \item{rejcounter}{Number of rejections (bc of swapping problem).}
//'   \item{earlystop}{Number of iterations after which Gibbs sampler was forced to abort for time constraints.}
// [[Rcpp::export]]
List gibbsC(int Kcalc,
	NumericVector alphastart,
	NumericVector betastart,
	List adjmat,
	List nrec,
	List npeer,
	List X,
	List Y,
	List Ustar,
	List Vstar,
	Function expgrid,
	Function rtruncnorm,
	Function sample,
	NumericVector ma0,
	NumericMatrix sa0,
	NumericVector mb0,
	NumericMatrix sb0,
	double tallowed,
	Function recfromUV3urR,
	bool reslocal,
	Function mvrnorm,
	bool updlast) {

  // Get initial clock reading for time abort
  double t1 = clock();
  int earlystop = NA_INTEGER;

  // Other cpp functions
  Function newUijC("newUijC");
  Function newVjiC("newVjiC");
  Function checkutilC("checkutilC");
  Function alphaFromU("alphaFromU");
  Function betaFromV("betaFromV");

  //-----------------------------------------//
  //     SET UP TO HANDLE MULTIPLE WAVES     //
  //-----------------------------------------//

  // Determine number of waves
  int nwav = adjmat.size();
  int nchar = alphastart.size();

  List iord(nwav);
  List jord(nwav);
  IntegerVector nind(nwav);
  List nindvec(nwav);

  for (int wav = 0; wav < nwav; ++wav) {

    // Make vector of all i,j index combinations
    // To be used later when pair order chosen randomly at each iteration
  	IntegerVector iind = seq_len(as<int>(nrec[wav])+1);
  	IntegerVector jind = seq_len(as<int>(npeer[wav])+3);

  	DataFrame indcombos = expgrid(iind, jind);
  	   // Remove NA combinations (when additional pop info)
  	   if (reslocal == true) {

  	     int totcombs = iind.size() * jind.size();
  	     IntegerVector iordw = indcombos["Var1"];
  	     IntegerVector jordw = indcombos["Var2"];
  	     IntegerMatrix adjm = as<IntegerMatrix>(adjmat[wav]);

  	     IntegerMatrix NAindcombo(totcombs);

  	     for (int totcombos = 0; totcombos < totcombs; ++totcombos) {
    	     IntegerVector tmpval(1);
    	     tmpval[0] = adjm(iordw[totcombos]-1, jordw[totcombos]-1);
    	     LogicalVector tcondv = is_na(tmpval);
    	     NAindcombo[totcombos] = tcondv[0];
  	     }

  	     int numNotNA = totcombs - sum(NAindcombo);

  	     IntegerVector iordwF(numNotNA);
  	     IntegerVector jordwF(numNotNA);

  	     int ctr = 0;
  	     for (int pairNoNA = 0; pairNoNA < totcombs; ++pairNoNA) {

  	       if (NAindcombo[pairNoNA] == false) {
  	         iordwF[ctr] = iordw[pairNoNA];
  	         jordwF[ctr] = jordw[pairNoNA];
  	         ctr += 1;
  	       }

  	     }

  	     iord[wav] = iordwF;
  	     jord[wav] = jordwF;

  	   } else {

  	     iord[wav] = indcombos["Var1"];
  	     jord[wav] = indcombos["Var2"];

  	   }

    nind[wav] = as<IntegerVector>(iord[wav]).size();
    nindvec[wav] = seq_len(nind[wav]);

  }

  bool sampreplace = false;

  int k = 0;

  NumericVector alphanew = alphastart;
  NumericVector betanew = betastart;

  int rejcounter = 0;

  Rprintf("..finished set up");

  //-----------------------------------------//
  //      ITERATE THROUGH GIBBS SAMPLER      //
  //-----------------------------------------//

  int wav = 0;

  List alphaSIM(Kcalc);
  List betaSIM(Kcalc);

  NumericMatrix alphaSIMtemp(nwav, nchar);
  NumericMatrix betaSIMtemp(nwav, nchar-2);

  int waverejcounter = 0;

  while(k < Kcalc & wav <= nwav-1) {

      if (wav == 0) {

        Rprintf("\nIteration %d: wave ", k+1);

      }

      Rprintf("%d..", wav+1);

      IntegerVector randpairord = sample(nindvec[wav], nind[wav], sampreplace);

      NumericMatrix Ustarw = Ustar[wav];
      NumericMatrix Vstarw = Vstar[wav];

      NumericMatrix Uold = clone(Ustarw);
      NumericMatrix Vold = clone(Vstarw);

      IntegerMatrix adjm = as<IntegerMatrix>(adjmat[wav]);

      if (waverejcounter > 10) {
        Rcpp::print(alphaSIMtemp);
        Rcpp::print(betaSIMtemp);
      }

      List Xw = X[wav];
      List Yw = Y[wav];

      for (int ijpair = 0; ijpair < nind[wav]; ++ijpair) {

        // i and j start at 1
        int ref = randpairord[ijpair] - 1;
        IntegerVector io = as<IntegerVector>(iord[wav]);
        IntegerVector jo = as<IntegerVector>(jord[wav]);
        int i = io[ref];
        int j = jo[ref];

          NumericVector Xij(nchar);
          NumericVector Yji(nchar);

          if (i < as<int>(nrec[wav]) + 1) {
            NumericMatrix tmpx = Xw[i-1];
            Xij = tmpx(_, j-1);
          } else {
            for (int cc = 0; cc < nchar; ++cc) {
              Xij[cc] = NA_REAL;
            }
          }

          if (j < as<int>(npeer[wav]) + 1) {
            NumericMatrix tmpy = Yw[j-1];
            Yji = tmpy(_, i-1);
          } else {
            for (int cc = 0; cc < nchar; ++cc) {
              Yji[cc] = NA_REAL;
            }
          }

          // Draw new Vstar(ji)
            double Vjinew = as<double>(newVjiC(i, j, adjm, Ustarw, Vstarw, betanew, Yji, rtruncnorm));

          NumericVector Vjinewvec(1);
          Vjinewvec[0] = Vjinew;

          if (is_false(any(is_na(Vjinewvec)))) {
            Vstarw(j-1, i-1) = Vjinew;
          }

          // Draw new Ustar(ij)
          double Uijnew = as<double>(newUijC(i, j, adjm, Ustarw, Vstarw, alphanew, Xij, rtruncnorm));

          NumericVector Uijnewvec(1);
          Uijnewvec[0] = Uijnew;

          if (is_false(any(is_na(Uijnewvec)))) {
            Ustarw(i-1, j-1) = Uijnew;
          }

      }

      bool compatiblecheck = as<bool>(checkutilC(Ustarw, Vstarw, adjm, recfromUV3urR, reslocal));

      // If good draw
      if (compatiblecheck == true || k <= 1) {

        waverejcounter = 0;

        // Update alpha and beta
        alphanew = alphaFromU(Xw, Ustarw, ma0, sa0, mvrnorm);

        double betao = betanew[betanew.size()-1];
        betanew = betaFromV(Yw, Vstarw, mb0, sb0, mvrnorm);
        if (wav == nwav - 1 & updlast == false) {
          betanew[betanew.size()-1] = betao;
        }

        alphaSIMtemp(wav, _) = alphanew;
        betaSIMtemp(wav, _) = betanew;

        alphaSIM[k] = clone(alphaSIMtemp);
        betaSIM[k] = clone(betaSIMtemp);

        // If last wave
        if (wav == nwav - 1) {
          k += 1;
          wav = 0;
        } else {
          // Go to next wave
          wav += 1;
        }

      } else {

        // Revert Ustar and Vstar to previous forms
        Ustarw = Uold;
        Vstarw = Vold;

        // Update rejection counter
        waverejcounter += 1;
        rejcounter += 1;
        Rprintf("Reject!\n");

      }

      // Add in time check for Hoffman
      // Allows function to abort when there is less than 1 hour remaining
      //  on the time allotted to the cluster
      double t2 = clock();
      double telapse = (t2 - t1) / CLOCKS_PER_SEC;

      if (telapse > (tallowed-1) * 60 * 60) {

        earlystop = k;
        k = Kcalc;

      }

    }  // END while loop

	return Rcpp::List::create(Rcpp::Named("alpha") = alphaSIM,
                          Rcpp::Named("beta") = betaSIM,
                          Rcpp::Named("rejcounter") = rejcounter,
                          Rcpp::Named("earlystop") = earlystop);

}
