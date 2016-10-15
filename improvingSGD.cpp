#include <RcppEigen.h>
#include <algorithm>    // std::max

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::LLT;
using Eigen::Lower;
using Eigen::MatrixXi;
using Eigen::Upper;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::SparseVector;
typedef Eigen::MappedSparseMatrix<double>  MapMatd;
typedef Map<MatrixXi>  MapMati;
typedef Map<VectorXd>  MapVecd;
typedef Map<VectorXi>  MapVeci;

typedef MapMatd::InnerIterator InIterMat;
typedef Eigen::SparseVector<double> SpVec;
typedef SpVec::InnerIterator InIterVec;

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

// Use this function to calculate inverse square root.
inline double invSqrt( const double& x ) {
double y = x;
double xhalf = ( double )0.5 * y;
long long i = *( long long* )( &y );
i = 0x5fe6ec85e7de30daLL - ( i >> 1 ); //LL suffix for (long long) type for GCC
y = *( double* )( &i );
y = y * ( ( double )1.5 - xhalf * y * y );

return y;
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

List improvingSGD(MapMatd X, VectorXd Y, VectorXd m, double step, VectorXd beta0, double lambda=0.0, int npass=1){
  // X is the design matrix stored in column-major format
  // i.e. with features for case i stores in column i
  // Y is the vector of counts
  // M is the vector of sample sizes
  // Thus Y[i] ~ Binomial( M[i], w[i])  )
  // w[i] = 1/{1+exp(- x[i] dot beta)}
  // where Beta is the regression vector we want to estimate
  // lambda is the l1 regularization parameter

  //INITIALIZE VALUES:
  int numobs = X.cols(); 
  int numfearures = X.rows();

  int iter = 0;
  int j = 0;
  double epsi = 1e-4;
  double Yi = Y(0);
  double mi = m(0);
  double grad = 0;
  double wi = .5;
  double wi_exponent = 0;

  VectorXd beta(numfearures);			//Initialize Beta_hat vector
  VectorXd hist_grad(numfearures);		//Initialize vector to store running hist_grad
  VectorXd adj_grad(numfearures);		//Initialize vector to store adj_grad
  for (int i=0;i<numfearures;i++){
    hist_grad(i) = 1e-3;
    beta(i)=beta0(i);
    adj_grad(i) = 0;
  }

  //Initialize elements to hold X, Y and m for a single observation (column).
  SparseVector<double> Xi(numobs);
  Xi=X.innerVector(0);

  // Bookkeeping: how long has it been since the last update of each feature?
  NumericVector last_updated(numfearures,0.0);
  
  //Initialize vectors to hold log-likelihood and running avg neg log-likelihood.
  double neglik = 0;
  NumericVector neglik_ra(numobs*npass,0.0);
  
  //Initialize variable to calculate penalty.
  double skip = 1;
  double accum_l2_penalty = 0;
  double gam = 0;
  
  // Outer loop: number of passes over data set
  double k = 0; // global interation counter
  for (int pass=0; pass < npass; pass++){

    //Loop over each observation (columns of x)
    for (int i=0; i<numobs; i++){

      //For linear predictor and E(Y[i]) from features
      Xi = X.innerVector(i);
      Yi = Y(i);
      mi = m(i);			

      //Update wi probability value
      wi_exponent = Xi.dot(beta);
      wi = 1 / (1 + exp(-wi_exponent));

      //Update neglik
      neglik = -(Yi * log(wi + epsi) + (mi - Yi) * log(1 - wi + epsi));
      if(iter > 0) {
        neglik_ra(iter) = (neglik_ra(iter-1) * (iter-1) + neglik) / iter;
      }			

      //Iterate over the active features for this instance
      for (SparseVector<double>::InnerIterator it(Xi); it; ++it){
        j = it.index();
        skip = k - last_updated(j);
        if (skip > 5){  skip = 5;} //Cap skip at 5
        last_updated(j) = k;

        //Update penalty and beta
        gam = step * adj_grad(j);
        accum_l2_penalty = beta(j) * ((1 - pow(1 + lambda * gam,skip)) / (1 - lambda * gam) );
        beta(j) -= accum_l2_penalty; 

        //Calculate l2 norm penalty
        double l2penalty = 2*lambda*beta(j);

        //Update hist_grad and gradient
        grad = (mi * wi - Yi) * it.value() + l2penalty;  
        hist_grad(j) += grad * grad;

        //Calculate adj_grad and beta
        adj_grad(j) = grad * invSqrt(hist_grad(j) + epsi);
        beta(j) -= step*adj_grad(j);	
      }
      k++; //increment global counter
    }	
  }
  
  // At the very end, apply the accumulated penalty for the variables we haven't touched recently
  for (int j=0; j< numfearures; ++j){
    skip = (iter - 1) - last_updated(j); 
    if (skip > 5){  skip = 5;} //Cap skip at 5
    gam = step * adj_grad(j);	
    accum_l2_penalty = beta(j) * ((1 - pow(1 + lambda * gam,skip)) / (1 - lambda * gam) );
    beta(j) -= accum_l2_penalty; 
  }	

  return List::create(Named("beta") = beta,
                      Named("numobs") = numobs,
                      Named("numfearures") = numfearures,
                      Named("iter") = iter,
                      Named("loglik") = neglik_ra);
  
}