#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace Eigen;
  using namespace density;

  DATA_VECTOR        (N);        /* Observations (counts) */
  DATA_FACTOR        (position); /* Observations */
  DATA_FACTOR        (time);     /* Observations */
  DATA_FACTOR        (sizeGroup);/* Observations */
  DATA_FACTOR        (haulid);   /* Observations */

  DATA_SPARSE_MATRIX (Q0);       /* For random field */
  DATA_SPARSE_MATRIX (I);        /* For random field */
  DATA_SPARSE_MATRIX (A);        /* Design matrix (standardization) */

  /* Prediction */
  DATA_SPARSE_MATRIX (Apredict);
  DATA_FACTOR     (pos_predict); /* Length = fine scale area. nlevels = num grid cells */

  /* Fixed effects */
  PARAMETER          (logdelta); /* For random field (corr) */
  PARAMETER          (logkappa); /* For random field (scale) */
  PARAMETER          (tphi_time);/* One-step time correlation */
  PARAMETER          (tphi_size);/* One-step time correlation */
  PARAMETER          (logsigma); /* Nugget error */
  PARAMETER_VECTOR   (beta);     /* For design matrix */

  /* Random effects */
  PARAMETER_ARRAY(eta);    // 3D: space x time x size
  PARAMETER_ARRAY(etanug); // 3D: size x haul
  PARAMETER_ARRAY(etamean);// 2D: size x time

  /* Parameter transforms */
  Type delta = exp(logdelta);
  Type kappa = exp(logkappa);
  Type sigma = exp(logsigma);
  Type phi_time = tphi_time / sqrt( 1.0 + tphi_time * tphi_time );
  Type phi_size = tphi_size / sqrt( 1.0 + tphi_size * tphi_size );

  /* Random fields */
  SparseMatrix<Type> Q = kappa * (Q0 + delta * I);
  GMRF_t<Type>      nldens_spatial = GMRF(Q);
  AR1_t<N01<Type> > nldens_time    = AR1(phi_time);
  AR1_t<N01<Type> > nldens_size    = AR1(phi_size);

  Type nll = 0; // Negative log likelhood

  // Flat beta prior
  DATA_SCALAR(huge_sd);
  nll -= dnorm(beta, Type(0), huge_sd, true).sum();
  nll -= dnorm(vector<Type>(etamean), Type(0), huge_sd, true).sum();

  // Process likelihood
  nll += SEPARABLE(nldens_size,SEPARABLE(nldens_time, nldens_spatial))(eta);

  // Nugget
  SCALE_t< AR1_t<N01<Type> > >  scaled_nldens_size = SCALE(nldens_size, sigma);
  for(int i=0; i < NLEVELS(haulid); i++){
    nll += scaled_nldens_size(vector<Type>(etanug.col(i)));
  }

  // Measurement likelihood
  vector<Type> predictor = A*beta;
  for(int i=0; i<N.size(); i++){
    nll -= dpois(N(i), 
		 exp(predictor(i) +
		     etamean(sizeGroup(i), time(i)) +
		     eta(position(i), time(i), sizeGroup(i)) +
		     etanug(sizeGroup(i), haulid(i))
		     ),
		 true);
  }

  // Spatial abundance
  vector<Type> lin_predictor = Apredict * beta;
  matrix<Type> dens(Apredict.rows(), NLEVELS(time));
  dens.setZero();
  for(int j=0; j<dens.cols(); j++){ // Time
    for(int i=0; i<dens.rows(); i++){ // Space
      for(int k=0; k<NLEVELS(sizeGroup); k++){ // Size
        dens(i, j) += exp(lin_predictor(i) + etamean(k, j) + eta(pos_predict(i), j, k));
      }
      dens(i, j) /= NLEVELS(sizeGroup);
    }
  }

  vector<Type> logindex = dens.colwise().mean().array().log();

  ADREPORT(logindex);
  REPORT(logindex);
  REPORT(dens);

  return nll;

}
