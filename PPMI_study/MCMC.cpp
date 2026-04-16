#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::depends("RcppArmadillo")]]

///////////////////////////////////////////////////////////////////
///////////// Tools: Generate variable from distribution  ////
///////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::mat rmvnorm(unsigned int n, const arma::vec& mu, const arma::mat& S) {
  unsigned int ncols = S.n_cols;
  if(arma::symmatu(S).is_sympd()){
    arma::mat Y(n, ncols);
    Y.imbue( norm_rand ) ;
    return arma::repmat(mu, 1, n).t() + Y * arma::chol(S);
  }else{
    arma::mat out(n,ncols);
    for(unsigned int j=0; j<n; ++j){
      for(unsigned int i=0; i<ncols; ++i){
        out(j,i) = rnorm(1, mu(i),sqrt(S(i,i)))[0];
      }
    }
    return out;
  }
}

// norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to be in the interval
// (a,b) via rejection sampling.
// ======================================================================
// [[Rcpp::export]]

double norm_rs(double a, double b)
{
   double  x;
   x = Rf_rnorm(0.0, 1.0);
   while( (x < a) || (x > b) ) x = norm_rand();
   return x;
}

// half_norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) (with a > 0) using half normal rejection sampling.
// ======================================================================

// [[Rcpp::export]]
double half_norm_rs(double a, double b)
{
   double   x;
   x = fabs(norm_rand());
   while( (x<a) || (x>b) ) x = fabs(norm_rand());
   return x;
}

// unif_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using uniform rejection sampling. 
// ======================================================================

// [[Rcpp::export]]
double unif_rs(double a, double b)
{
   double xstar, logphixstar, x, logu;

   // Find the argmax (b is always >= 0)
   // This works because we want to sample from N(0,1)
   if(a <= 0.0) xstar = 0.0;
   else xstar = a;
   logphixstar = R::dnorm(xstar, 0.0, 1.0, 1.0);

   x = R::runif(a, b);
   logu = log(R::runif(0.0, 1.0));
   while( logu > (R::dnorm(x, 0.0, 1.0,1.0) - logphixstar))
   {
      x = R::runif(a, b);
      logu = log(R::runif(0.0, 1.0));
   }
   return x;
}

// exp_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using exponential rejection sampling.
// ======================================================================

// [[Rcpp::export]]
double exp_rs(double a, double b)
{
  double  z, u, rate;

//  Rprintf("in exp_rs");
  rate = 1/a;
//1/a

   // Generate a proposal on (0, b-a)
   z = R::rexp(rate);
   while(z > (b-a)) z = R::rexp(rate);
   u = R::runif(0.0, 1.0);

   while( log(u) > (-0.5*z*z))
   {
      z = R::rexp(rate);
      while(z > (b-a)) z = R::rexp(rate);
      u = R::runif(0.0,1.0);
   }
   return(z+a);
}




// rnorm_trunc( mu, sigma, lower, upper)
//
// generates one random normal RVs with mean 'mu' and standard
// deviation 'sigma', truncated to the interval (lower,upper), where
// lower can be -Inf and upper can be Inf.
//======================================================================

// [[Rcpp::export]]
double sim_rtnorm(double mu, double sigma, double lower, double upper)
{
int change;
 double a, b;
 double logt1 = log(0.150), logt2 = log(2.18), t3 = 0.725;
 double z, tmp, lograt;

 change = 0;
 a = (lower - mu)/sigma;
 b = (upper - mu)/sigma;

 // First scenario
 if( (a == R_NegInf) || (b == R_PosInf))
   {
     if(a == R_NegInf)
       {
     change = 1;
     a = -b;
     b = R_PosInf;
       }

     // The two possibilities for this scenario
     if(a <= 0.45) z = norm_rs(a, b);
     else z = exp_rs(a, b);
     if(change) z = -z;
   }
 // Second scenario
 else if((a * b) <= 0.0)
   {
     // The two possibilities for this scenario
     if((R::dnorm(a, 0.0, 1.0,1.0) <= logt1) || (R::dnorm(b, 0.0, 1.0, 1.0) <= logt1))
       {
     z = norm_rs(a, b);
       }
     else z = unif_rs(a,b);
   }
 // Third scenario
 else
   {
     if(b < 0)
       {
     tmp = b; b = -a; a = -tmp; change = 1;
       }

     lograt = R::dnorm(a, 0.0, 1.0, 1.0) - R::dnorm(b, 0.0, 1.0, 1.0);
     if(lograt <= logt2) z = unif_rs(a,b);
     else if((lograt > logt1) && (a < t3)) z = half_norm_rs(a,b);
     else z = exp_rs(a,b);
     if(change) z = -z;
   }
   double output;
   output = sigma*z + mu;
 return (output);
}

// [[Rcpp::export]]
arma::vec sim_dirichlet(arma::vec a){
  unsigned int K = a.n_elem;
  
  arma::vec x_temp(K);
  for(unsigned k = 0; k<K; ++k){
    x_temp(k) = rgamma(1, a(k), 1)[0];
  }
  double x_sum = arma::sum(x_temp);
  arma::vec y_out = x_temp/x_sum;
  
  return y_out;
}

// [[Rcpp::export]]
double sim_invgauss(double a, double b){
// here "a" denotes the mean parameter, "b" denotes the shape parameter
    double result = 0;
    double v = rnorm(1,0,1)[0];
    double y = v*v;
    double x = a + a * a * y / (2 * b) - a * sqrt(4 * a * b * y + a * a * y * y)/ (2 * b);

    double test = runif(1,0,1)[0];
    if (test < a / (a + x)) {
        result = x;
    }
    else
    {
        result = a * a / x;
    }
    return result;
}

// [[Rcpp::export]]
bool containsNaN(NumericVector A) {
  for (int i = 0; i < A.size(); i++) {
    if(NumericVector::is_na(A[i])){
      return true;  // A contains NaN
    }else if(Rcpp::traits::is_infinite<REALSXP>(A[i])){
      return true;  // A contains NaN
    }
  }
  return false;  // A does not contain NaN
}

// [[Rcpp::export]]
List f_judge(NumericVector mu, int M) {
  
  if (is_true(any(duplicated(mu)))) {
    Rf_warning("There are duplicates in 'mu'; order not guaranteed to match that of R's base::order");
  }
  NumericVector sorted = clone(mu).sort();
  IntegerVector norder = match(sorted, mu);
  
  bool judge = FALSE;
  for(int m=0; m<M; ++m){
    if(norder[m] != m+1){
      judge = TRUE;
    }
  }
  return List::create(Rcpp::Named("judge") = judge,
                      Rcpp::Named("norder") = norder);
}

// [[Rcpp::export]]
arma::uvec my_setdiff(arma::uvec& x, const arma::uvec& y){
  for (size_t j = 0; j < y.n_elem; j++) {
    arma::uword q1 = arma::conv_to<arma::uword>::from(arma::find(x == y(j)));
    x.shed_row(q1);
  }
  return x;
}

///////////////////////////////////////////////////////////////////
////////////////// Parameter Update  //////////////////////////////
///////////////////////////////////////////////////////////////////


// [[Rcpp::export]]
void update_omega_matrix(const arma::mat& Y,
                         const arma::mat& X,
                         const arma::mat& I_K,
                         arma::mat& omega,
                         arma::mat& b_temp,
                         const arma::vec& psi_eps_inv,
                         const arma::mat& zeta_mat,
                         const arma::mat& Z_mat,
                         const arma::vec& psi_e_inv,
                         const arma::vec& Queue,
                         const arma::vec& NT) {

  unsigned int SNT = arma::sum(NT);
  unsigned int Q = omega.n_cols;
  unsigned int R = X.n_cols;

  // Update omega
  // EFA part
  arma::mat Psi_inv = arma::diagmat(psi_eps_inv);
  arma::mat Psi_inv_Lam = Psi_inv * b_temp;
  arma::mat V_EFA_omega = I_K + b_temp.t() * Psi_inv_Lam;
  arma::mat m_EFA_omega = Y * Psi_inv_Lam;

  // MMM part
  arma::mat fzi(SNT*R, Q);
  arma::mat azi(SNT,R);
  for (unsigned int i = 0; i < SNT; ++i) {
    unsigned int index = Z_mat(i,0) - 1;
    for(unsigned int r = 0; r < R; ++r){
      fzi.row(r*SNT+i) = zeta_mat.submat((index*R+r), 1, (index*R+r), Q);
      azi(i,r) = zeta_mat((index*R+r), 0);
    }
  }
  //Rprintf("Passed MMM\n");
  arma::mat V_a = X - azi;

  for (unsigned int i = 0; i < SNT; ++i) {
    arma::mat V_MMM_omega(Q,Q);
    arma::mat m_MMM_omega(Q,1);
    for(unsigned int r = 0; r < R; ++r){
      // Sigma omega
      arma::mat FOF = fzi.row(r*SNT+i).t() * fzi.row(r*SNT+i);
      V_MMM_omega += psi_e_inv(Z_mat(i,0) - 1) * FOF;
      m_MMM_omega += V_a(i,r) * psi_e_inv(Z_mat(i,0) - 1) * fzi.row(r*SNT+i).t();
    }
    // Sigma and mean
    arma::mat Sigma = arma::inv(V_EFA_omega + V_MMM_omega);
    arma::mat m_omega = Sigma * (m_EFA_omega.row(i).t() + m_MMM_omega);
    omega.row(i) = m_omega.t() + arma::randn<arma::rowvec>(1, Q) * arma::chol(Sigma);
  }
}

// [[Rcpp::export]]
void update_psi_eps(const arma::mat& Y,
                         arma::mat& omega,
                         arma::mat& b_temp,
                         arma::vec& psi_eps_inv) {
  // One distinct uniqueness variance per variable

  int P = b_temp.n_rows;
  unsigned int N = omega.n_rows;
  
  double b1 = 1;
  double b2 = 1;
  
  arma::mat Y_hat = omega * b_temp.t();
  
  arma::mat Y_resids = Y - Y_hat;
  arma::colvec d2 = trans( sum( pow(Y_resids, 2), 0) );

  double u;
  int n_obs;
  for (int j = 0; j < P; j++) {
    
    u = R::runif(0.0,1.0);
    n_obs = N;
    
    // qgamma is parameterized as shape, scale rather than shape and rate
    psi_eps_inv(j) = (R::qgamma(u,
                   double(n_obs + 2.0*b1)/2.0,
                   2.0/(2.0*b2 + d2[j]),1,0));
    //Rprintf("j: %d; alpha: %f; beta: %f; psi_eps_inv: %f\n", (j+1), (n_obs + 2.0*b1),  (2.0*b2 + d2[j]), psi_eps_inv(j));
  }
}


// [[Rcpp::export]]
void update_b_temp(const arma::mat& Y,
                arma::uvec& r_idx,
                arma::mat& omega,
                arma::mat& b_temp,
                arma::vec& psi_eps_inv,
                arma::vec& gammaj, 
                arma::mat& ukj,
                arma::mat& tau2, 
                double my_gamma) {

  unsigned int P = b_temp.n_rows;
  unsigned int Q = omega.n_cols;
  
  //update b_temp
  arma::mat OpO = omega.t()*omega;
  arma::mat OpY = omega.t()*Y;
  arma::rowvec b_tempj(Q);

  // Constrained case
  unsigned int j;
  for(unsigned int q=0;q<Q;++q){
    j = r_idx(q);
    arma::mat OpO_1toj = OpO(arma::span(0,q),arma::span(0,q));
    arma::vec OpYj = OpY(arma::span(0,q),j);
    arma::mat invClamj((q+1), (q+1), arma::fill::zeros);
    for (unsigned s = 0; s < (q+1); ++s) {
      double c_jq = ukj(j,s) * gammaj(s) + (1 - ukj(j,s) * gammaj(s)) * 0.001;
      /*if(s == q){
        Rprintf("j: %d; c_jq: %f; tau2: %f\n", (j+1), c_jq, tau2(j,s));
      }*/
      invClamj(s,s) = 1/(c_jq*tau2(j,s));
    }
    arma::mat C_1toj = inv(psi_eps_inv(j)*OpO_1toj+invClamj);
    arma::vec m_1toj = psi_eps_inv(j)*C_1toj*OpYj;
    
    // update b_temp_jj
    // my_gamma = 0.5
    /*double b_tempjj=sim_gamma_type(b_temp(j,q),my_gamma+1.0,m_1toj(q),2.*C_1toj(q,q));
    if (!arma::is_finite(b_tempjj)) {
      b_tempjj = 0;
      Rprintf("j: %d; b_tempjj is not finite\n", (j+1));
      Rprintf("m_1toj: %f; C_1toj: %f\n", m_1toj(q), C_1toj(q,q));
      Rprintf("invClamj: %f; tau2: %f\n", invClamj(q,q), tau2(j,q));
    }*/
    double b_tempjj = sim_rtnorm(m_1toj(q),sqrt(C_1toj(q,q)),0,R_PosInf);
    
    // update b_tempj conditioned upon b_tempjj
    if (q > 0) {
      arma::vec m1 = m_1toj(arma::span(0,q-1));
      arma::mat S1 = C_1toj(arma::span(0,q-1),arma::span(0,q-1));
      arma::vec S12= C_1toj(arma::span(0,q-1),q);
      arma::vec mb = m1+S12*(b_tempjj-m_1toj(q))/C_1toj(q,q);
      arma::mat Sb = S1-S12*S12.t()/C_1toj(q,q);
      arma::rowvec b_temp_1tjm1 = rmvnorm(1,mb,Sb);
      b_temp(j,arma::span(0,q-1))=b_temp_1tjm1;
    }
    b_temp(j,q) = b_tempjj;
    
    if ((q+1) < Q) {
      arma::rowvec fill_zeros = arma::zeros<arma::rowvec>(Q-(q+1));
      b_temp(j,arma::span(q+1,Q-1)) = fill_zeros;
    }
    
  }
  
  // Unconstrained case
  arma::uvec one_to_J = arma::linspace<arma::uvec>(0, P-1, P);
  arma::uvec not_r_idx = my_setdiff(one_to_J, r_idx);
  unsigned int I = not_r_idx.n_elem;
  for(unsigned int i=0; i<I; ++i) {
    j = not_r_idx(i);
    arma::mat invClamj(Q, Q, arma::fill::zeros);
    for (unsigned s = 0; s < Q; ++s) {
      double c_jq = ukj(j,s) * gammaj(s) + (1 - ukj(j,s) * gammaj(s)) * 0.001;
      invClamj(s,s) = 1/(c_jq*tau2(j,s));
    }
    arma::mat Cj = inv(psi_eps_inv(j)*OpO+invClamj);
    arma::vec mj = psi_eps_inv(j)*Cj*OpY.col(j);
    b_tempj = rmvnorm(1,mj,Cj);
    b_temp.row(j) = b_tempj;
  }
}

// [[Rcpp::export]]
void update_tau(arma::mat& b_temp,
                  arma::vec& gammaj, 
                  arma::mat& ukj,
                  arma::mat& tau2){

  unsigned int P = b_temp.n_rows;
  unsigned int Q = b_temp.n_cols;

  double a_tau_temp = 5 + 0.5;

  for(unsigned int j=0; j<Q; ++j){
    for(unsigned int k=0; k<P; ++k){
      double ckj = ukj(k,j)*gammaj(j)+(1-ukj(k,j)*gammaj(j))*0.001;
      double b_tau_temp = 25 + 0.5 * pow(b_temp(k,j),2)/(2.0*ckj);;
      tau2(k,j) = 1/rgamma(1,a_tau_temp,(1/b_tau_temp))[0];
    }
  }
}

// [[Rcpp::export]]
void update_ukj(arma::mat& b_temp,
                arma::vec& gammaj, 
                arma::mat& ukj,
                arma::mat& p_ukj,
                arma::mat& tau2){

  unsigned int P = b_temp.n_rows;
  unsigned int Q = b_temp.n_cols;

  NumericVector samp(2);
  for(int m=0; m<2; ++m){
    samp[m] = m;
  }

  for(unsigned int j=0; j<Q; ++j){
    if (gammaj(j)==0){
      ukj.row(j).zeros();
    }else{
      for(unsigned int k=0; k<P; ++k){
        NumericVector prob(2);
        double nominator = exp(-pow(b_temp(k,j),2)/(2.0*tau2(k,j)*0.001))/sqrt(2.0*M_PI*tau2(k,j)*0.001);
        double denominator = exp(-pow(b_temp(k,j),2)/(2.0*tau2(k,j)))/sqrt(2.0*M_PI*tau2(k,j));
        double A = 0.0;
        if(denominator != 0.0){
          A = (1-p_ukj(k,j))*nominator/(p_ukj(k,j)*denominator);
        }
        prob[0] = 1-1/(1+A);
        prob[1] = 1/(1+A);
        //Rprintf("j:%d; k:%d; prob1: %f; ukj:%f\n", (j+1), (k+1), prob[1], ukj_out(k,j));
        ukj(k,j) = Rcpp::sample(samp, 1, false, prob)[0];
      }
    }
  }
}

// [[Rcpp::export]]
void update_p_ukj(arma::mat& ukj,
                  arma::mat& p_ukj){
              
  unsigned int P = ukj.n_rows;
  unsigned int Q = ukj.n_cols;

  arma::mat prior_p_ukj(P, Q, arma::fill::ones);
  arma::mat a_p_ukj = prior_p_ukj + ukj;
  arma::mat b_p_ukj = prior_p_ukj + prior_p_ukj - ukj;

  for(unsigned int j=0; j<Q; ++j){
    for(unsigned int k=0; k<P; ++k){
      p_ukj(k,j) = rbeta(1,a_p_ukj(k,j),b_p_ukj(k,j))[0];
    }
  }
}

// [[Rcpp::export]]
void update_gammaj(arma::mat& b_temp,
              arma::vec& gammaj,
              arma::vec& p_gammaj,
              arma::mat& ukj,
              arma::mat& p_ukj,
              arma::mat& tau2) {
  unsigned int P = b_temp.n_rows;
  unsigned int Q = b_temp.n_cols;

  NumericVector samp(2);
  for (unsigned int m = 0; m < 2; ++m) {
    samp[m] = m;
  }

  double sqrt2Pi = sqrt(2.0 * M_PI);
  double sqrt2PiTauFactor = sqrt(2.0 * M_PI * 0.001);

  for (unsigned int j = 0; j < Q; ++j) {
    NumericVector prob(2);
    double nominator = (1 - p_gammaj(j));
    double denominator = p_gammaj(j);

    for (unsigned int k = 0; k < P; ++k) {
      double tau2Factor = 2.0 * tau2(k, j) * 0.001;
      double nominatorExp = exp(-(b_temp(k, j) * b_temp(k, j)) / tau2Factor) / (sqrt2PiTauFactor * sqrt(tau2(k, j)));
      double denominatorExp = exp(-(b_temp(k, j) * b_temp(k, j)) / (2.0 * tau2(k, j))) / (sqrt2Pi * sqrt(tau2(k, j)));

      nominator *= nominatorExp;
      denominator *= ((1 - p_ukj(k, j)) * nominatorExp + p_ukj(k, j) * denominatorExp);
    }

    double A = nominator / denominator;
    prob[0] = 1 - 1 / (1 + A);
    prob[1] = 1 / (1 + A);
    gammaj(j) = Rcpp::sample(samp, 1, false, prob)[0];
  }
}

// [[Rcpp::export]]
void update_p_gammaj(arma::vec& gammaj,
                    arma::vec& p_gammaj){

  unsigned int Q = gammaj.n_elem;

  for(unsigned int j=0; j<Q; ++j){
    double a = 1.0 + gammaj(j);
    double b = 2.0 - gammaj(j);
    p_gammaj(j) = rbeta(1,a,b)[0];
  }
}


// [[Rcpp::export]]
void update_alpha(const arma::mat& W,
                  const arma::vec& beta,
                  arma::mat& omega,
                  arma::mat& b_temp,
                  arma::vec& psi_eps_inv,
                  arma::vec& alpha,
                  const arma::vec& sigma_alpha,
                  const arma::vec& beta_cumsum,
                  const unsigned int PO,
                  arma::vec& at_alpha,
                  arma::mat& Y) {
  
  unsigned int N = omega.n_rows;
  unsigned int Q = omega.n_cols;
  
  
  arma::mat OpB = omega * b_temp.t();
  arma::mat mu = omega*b_temp.t();
  for(unsigned int p=0; p<PO; ++p){
    unsigned int beta_p = beta(p);
    if(beta_p > 2){
      // sample new alpha_star
      arma::vec alpha_star_p(beta_p);
      arma::vec alpha_p = alpha.subvec(beta_cumsum[p], beta_cumsum[p+1]-1);
      alpha_star_p(0) = alpha_p(0);
      for(unsigned int w = 1; w<beta_p-1; ++w){
        alpha_star_p(w) = sim_rtnorm(alpha_p(w), sigma_alpha(p), alpha_star_p(w-1), alpha_p(w+1));
        //Rprintf("p: %d; w: %d; alpha_star_p: %f; alpha_p: %f\n", (p+1), w, alpha_star_p(w-1), alpha_p(w+1));
      }
      //alpha_star_p(beta_p-1) = sim_rtnorm(alpha_p(beta_p-1), sigma_alpha(p), alpha_star_p(beta_p-2), R_PosInf);
      alpha_star_p(beta_p-1) = alpha_p(beta_p-1);

      // Compute ratio 1
      arma::vec deno_diff1 = alpha_star_p.subvec(2, beta_p-1) - alpha_star_p.subvec(1, beta_p-2);
      arma::vec deno_diff2 = alpha_p.subvec(0, beta_p-3) - alpha_star_p.subvec(1, beta_p-2);
      arma::vec nomi_diff1 = alpha_p.subvec(2, beta_p-1) - alpha_p.subvec(1, beta_p-2);
      arma::vec nomi_diff2 = alpha_star_p.subvec(0, beta_p-3) - alpha_p.subvec(1, beta_p-2);

      arma::vec nomi_diff_cdf = arma::normcdf(nomi_diff1 / sigma_alpha(p)) - arma::normcdf(nomi_diff2 / sigma_alpha(p));
      arma::vec deno_diff_cdf1 = arma::normcdf(deno_diff1 / sigma_alpha(p)) - arma::normcdf(deno_diff2 / sigma_alpha(p));
      arma::vec deno_diff_cdf = 1.0 / deno_diff_cdf1;
      arma::vec ratio_beta = nomi_diff_cdf % deno_diff_cdf;

      // Compute ratio 2
      arma::vec alpha_tilde(beta_p + 2);
      alpha_tilde(0) = R_NegInf;
      alpha_tilde(beta_p + 1) = R_PosInf;
      alpha_tilde.subvec(1, beta_p) = alpha_p;

      arma::vec alpha_star_tilde(beta_p + 2);
      alpha_star_tilde(0) = R_NegInf;
      alpha_star_tilde(beta_p + 1) = R_PosInf;
      alpha_star_tilde.subvec(1, beta_p) = alpha_star_p;

      arma::vec alpha_w_1(N);
      arma::vec alpha_w(N);
      arma::vec alpha_w_star_1(N);
      arma::vec alpha_w_star(N);
      //Rprintf("Passed3\n");
      for(unsigned int i=0;i<N;++i){
        alpha_w_1(i) = alpha_tilde(W(i,p)+1);
        alpha_w(i) = alpha_tilde(W(i,p));
        alpha_w_star_1(i) = alpha_star_tilde(W(i,p)+1);
        alpha_w_star(i) = alpha_star_tilde(W(i,p));  
      }  

      double psi_eps_inv_p = std::sqrt(psi_eps_inv(p));
      
      arma::vec deno_diff3 = (alpha_w_1 - OpB.col(p)) * psi_eps_inv_p;
      arma::vec deno_diff4 = (alpha_w - OpB.col(p)) * psi_eps_inv_p;
      arma::vec nomi_diff3 = (alpha_w_star_1 - OpB.col(p)) * psi_eps_inv_p;
      arma::vec nomi_diff4 = (alpha_w_star - OpB.col(p)) * psi_eps_inv_p;

      arma::vec nomi_N = arma::normcdf(nomi_diff3) - arma::normcdf(nomi_diff4);
      arma::vec deno_N1 = arma::normcdf(deno_diff3) - arma::normcdf(deno_diff4);
      arma::vec deno_N = 1.0 / deno_N1;
      arma::vec ratio_N = nomi_N % deno_N;
      //Rprintf("i: %d; nomi_diff3: %f; nomi_diff4: %f; deno_diff3: %f; deno_diff4: %f; nomi_N: %f; deno_N1: %f; deno_N: %f; ratio_N: %f\n", 1, nomi_diff3(0), nomi_diff4(0), deno_diff3(0), deno_diff4(0), nomi_N(0), deno_N1(0), deno_N(0), ratio_N(0));

      // Compute ratio
      double ratio1 = arma::prod(ratio_beta);
      double ratio2 = arma::prod(ratio_N);
      double ratio = std::min(ratio1*ratio2,1.0);

      // Judge
      double randnum = arma::randu();

      /*Rprintf("p: %d; ratio1: %f; ratio2: %f; ratio: %f; randnum: %f\n", (p+1), ratio1, ratio2, ratio, randnum);
      arma::vec ratio_prod(10);
      for(unsigned int i=0; i<10; ++i){
        ratio_prod(i) = arma::prod(ratio_N.subvec(i*300+0, i*300+299));
      }
      Rcpp::Rcout << ratio_prod.t() << " ";
      Rcpp::Rcout << std::endl;
      Rcpp::Rcout << alpha_tilde.t() << " ";
      Rcpp::Rcout << std::endl;
      Rcpp::Rcout << alpha_star_tilde.t() << " ";
      Rcpp::Rcout << std::endl;*/

      /*if(beta_p == 4){
        Rprintf("p: %d; ratio1: %f; ratio2: %f; ratio: %f; randnum: %f\n", (p+1), ratio1, ratio2, ratio, randnum);
        arma::vec ratio_prod(10);
        for(unsigned int i=0; i<10; ++i){
          ratio_prod(i) = arma::prod(ratio_N.subvec(i*10+0, i*10+9));
        }
        Rcpp::Rcout << ratio_prod.t() << " ";
        Rcpp::Rcout << std::endl;
        Rcpp::Rcout << alpha_tilde.t() << " ";
        Rcpp::Rcout << std::endl;
        Rcpp::Rcout << alpha_star_tilde.t() << " ";
        Rcpp::Rcout << std::endl;
      }*/
    
      if (randnum < ratio) {
        alpha.subvec(beta_cumsum(p), beta_cumsum(p+1)-1) = alpha_star_p;
        at_alpha(p) += 1;
        // --------- update_Y -------------------
        double sigma = sqrt(1/psi_eps_inv(p));
        arma::vec alpha_tilde_p(beta(p) + 2);
        alpha_tilde_p(0) = R_NegInf;
        alpha_tilde_p(beta(p) + 1) = R_PosInf;
        alpha_tilde_p.subvec(1, beta(p)) = alpha.subvec(beta_cumsum(p), beta_cumsum(p+1)-1);
        for (unsigned int i = 0; i < N; ++i) {
          Y(i, p) = sim_rtnorm(mu(i, p), sigma, alpha_tilde_p(W(i,p)), alpha_tilde_p(W(i,p) + 1));
        }
      }
    }else{
      // --------- update_Y -------------------
      double sigma = sqrt(1/psi_eps_inv(p));
      arma::vec alpha_tilde_p(beta(p) + 2);
      alpha_tilde_p(0) = R_NegInf;
      alpha_tilde_p(beta(p) + 1) = R_PosInf;
      alpha_tilde_p.subvec(1, beta(p)) = alpha.subvec(beta_cumsum(p), beta_cumsum(p+1)-1);
      for (unsigned int i = 0; i < N; ++i) {
        Y(i, p) = sim_rtnorm(mu(i, p), sigma, alpha_tilde_p(W(i,p)), alpha_tilde_p(W(i,p) + 1));
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////
////////////////// MMM part ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
void update_delta(double& delta_0, 
            arma::vec& pi, 
            arma::mat& g_mat,
            double& at_delta,
            double& c_delta){

  unsigned int N = g_mat.n_rows;
  unsigned int M = g_mat.n_cols;
  
  arma::vec delta = pi*delta_0;
  arma::vec delta_star(M);
  
  // sample new delta_star
  for(unsigned int m=0; m<M; ++m){
    delta_star(m) = rlnorm(1, log(delta(m)), c_delta)[0];
  }
  double delta_star_0 = arma::sum(delta_star);
  arma::vec pi_star = delta_star/delta_star_0;
  arma::vec delta_diff = delta_star - delta;
  arma::mat log_g_mat = arma::log(g_mat);
  arma::mat g1_col_sum = arma::sum(log_g_mat, 0);

  // compute ratio
  double ratio_1 = 5.0 * (delta_0 - delta_star_0);
  double ratio_2 = sum(log(delta_star) - log(delta));
  double ratio_3 = (double(N)) * (lgamma(delta_star_0) - lgamma(delta_0));
  double ratio_4 = 0.0;
  for(unsigned int m=0; m<M; ++m){
    ratio_4 += (lgamma(delta(m)) - lgamma(delta_star(m)));
  }
  ratio_4 = (double(N)) *ratio_4;
  double ratio_5 = arma::dot(delta_diff, g1_col_sum);
  double ratio_6 = (1.0 - 1.0) * (log(delta_star_0) - log(delta_0));
  double ratio_delta = ratio_1 + ratio_2 + ratio_3 + ratio_4 + ratio_5 + ratio_6;

  // judge
  ratio_delta = exp(ratio_delta);
  
  double randnum = R::runif(0.0,1.0);
  if(randnum<ratio_delta){
    delta_0 = delta_star_0;
    pi = delta_star / delta_star_0;
    at_delta += 1;
  }
}

// [[Rcpp::export]]
void update_gi(arma::vec& pi,
              double& delta_0,
              arma::mat& Z_mat,
              arma::mat& g_mat,
              const arma::vec& Queue,
              const arma::vec& NT){

  unsigned int N = g_mat.n_rows;
  unsigned int M = g_mat.n_cols;      

  arma::mat delta(N, M);    
  delta.each_row() = pi.t() * delta_0;
  
  for(unsigned int i=0; i<N; ++i){
    unsigned int TT = NT(i);
    for(unsigned int t = 0; t<TT; ++t){
      unsigned int index = Z_mat(Queue(i)+t,0)-1;
      delta(i,index) += 1.0;
    }
    g_mat.row(i) = sim_dirichlet(delta.row(i).t()).t();
    for(unsigned int m=0; m<M; ++m){
      if(g_mat(i,m) == 0){
        g_mat(i,m) = 1e-323;
      }
    }
  }
}

// [[Rcpp::export]]
void update_Z(const arma::mat& V,
              arma::mat& omega,
              arma::mat& zeta_mat,
              arma::mat& Z_mat,
              arma::mat& g_mat,
              arma::vec& psi_e_inv,
              const arma::vec& Queue,
              const arma::vec& NT,
              const arma::vec& MMM_threshold){

  unsigned int N = g_mat.n_rows;
  unsigned int Q = omega.n_cols;
  unsigned int M = g_mat.n_cols;
  unsigned int SNT = arma::sum(NT);
  unsigned int R = V.n_cols;

  NumericVector samp(M);
  for(unsigned int m=0; m<M; ++m){
    samp[m] = m+1;
  }

  arma::mat xi(SNT, (Q+1));
  xi.col(0).fill(1.0); 
  xi.submat(0,1,SNT-1,Q) = omega;
  
  arma::mat E_v(SNT, M*R);
  for (unsigned int i = 0; i < SNT; ++i) {
    for(unsigned r=0; r<R; ++r){
      for (unsigned int m = 0; m < M; ++m) {
        E_v(i, m*R+r) = arma::dot(zeta_mat.row(m*R+r),xi.row(i).t()) - MMM_threshold(0);
      }
    }
  }
  
  arma::mat V_rep = arma::repmat(V, 1, M);    
  arma::mat prob1 = arma::normcdf(E_v);
  arma::mat prob2 = prob1 % V_rep + (1 - prob1) % (1 - V_rep);
  arma::mat prob3(SNT,M, arma::fill::ones);
  for (unsigned int i = 0; i < SNT; ++i) {
    for(unsigned int m=0; m<M; ++m){
      for (unsigned int r = 0; r < R; ++r){
        prob3(i,m) *= prob2(i,m*R+r);
      }
    }
  }
  arma::mat gV_exp(SNT, M);
  for (unsigned int i = 0; i < N; ++i) {
    // ----------------- probability -------------------
    unsigned int TT = NT(i);
    for(unsigned int t = 0; t<TT; ++t){
      gV_exp.row((Queue(i)+t)) = (g_mat.row(i) % prob3.row((Queue(i)+t)));
    }
  }

  double threshold = 1.0 / double(M);
  arma::vec prob_sum = arma::sum(gV_exp, 1);
  arma::mat prob(SNT,M);
  for (unsigned int i = 0; i < SNT; ++i) {
    prob.row(i) = gV_exp.row(i) / prob_sum(i);
    double prob_max = arma::max(prob.row(i));
    bool condition1 = containsNaN(Rcpp::wrap(arma::conv_to<arma::mat>::from(prob.row(i))));
    bool condition = prob_max> 1.0||prob_max < threshold|| condition1;
    if (condition) {
        prob.row(i).fill(threshold);
    }
    NumericVector prob_num = Rcpp::wrap(arma::conv_to<arma::vec>::from(prob.row(i)));
    //Rprintf("i:%d; prob_num1:%f; prob_num:%f\n", (i+1), prob_num[0], prob_num[1]);
    Z_mat(i,0) = Rcpp::sample(samp, 1, false, prob_num)[0];
  }
}


// [[Rcpp::export]]
void update_psi_e_inv(const arma::mat& X,
                  arma::mat& omega,
                  arma::mat& zeta_mat,
                  arma::vec& Z_mat,
                  arma::vec& psi_e_inv){
  // One distinct uniqueness variance per variable

  unsigned int SNT = X.n_rows;
  unsigned int Q = omega.n_cols;
  unsigned int M = zeta_mat.n_rows;

  arma::mat xi(SNT, (Q+1));
  xi.col(0).fill(1.0); 
  xi.submat(0,1,SNT-1,Q) = omega;

  arma::mat E_v(SNT,1);
  for (unsigned int i = 0; i < SNT; ++i) {
    unsigned int index = Z_mat(i) - 1;
    E_v(i,0) = arma::dot(zeta_mat.row(index),xi.row(i).t());
  }
  
  arma::mat eta_resids = X - E_v;
  arma::mat d2_temp = eta_resids % eta_resids;

  double u;
  for (unsigned int m=0; m<M; m++) {
    u = R::runif(0.0,1.0);
    unsigned int n_m = arma::accu(Z_mat == (m+1));
    double d2 = arma::accu(d2_temp % (Z_mat == (m+1)));
    // qgamma is parameterized as shape, scale rather than shape and rate
    psi_e_inv(m) = (R::qgamma(u,
                   double(n_m + 9.0)/2.0,
                   2.0/(4.0 + d2),1,0));
  }
}


// [[Rcpp::export]]
void f_rho_temp(arma::mat& zeta_mat, 
                arma::mat& rho_temp,
                const unsigned int M,
                const unsigned int R) {

  unsigned int Q = zeta_mat.n_cols;
  arma::mat new_zeta(M*R,Q);
  if(M > 1){
    for(unsigned int r = 0; r<R; ++r){
      // ---------------------- give new order --------------------------
      NumericVector mua(M);
      for (unsigned int m = 0; m < M; ++m){
        mua[m] = zeta_mat(m*R+r,0);
      }
      List all = f_judge(mua, M);
      IntegerVector norder = all[1];
      // ---------------------- assign new order to rho --------------------------
      for(unsigned int m=0; m<M; ++m){
        unsigned int new_m = norder[m]-1;
        for(unsigned int j=0; j<Q; ++j){
          new_zeta(m*R+r,j) = zeta_mat(new_m*R+r,j);
        }
      }
    }
      // ----------------------- compute rho -------------------------
    for(unsigned int r = 0; r<R; ++r){
      for(unsigned int j=0; j<Q; ++j){
        rho_temp(r,j) = new_zeta(r,j);
        for(unsigned int m=1; m<M; ++m){
          rho_temp(m*R+r,j) = new_zeta(m*R+r,j) - new_zeta(((m-1)*R+r),j);
        }
      }
    }
  }else{
    for(unsigned int r = 0; r<R; ++r){
      for(unsigned int j=0; j<Q; ++j){
        rho_temp(r,j) = new_zeta(r,j);
      }
    }
  }  
}

// [[Rcpp::export]]
void update_rho(const arma::mat& X,
                arma::mat& omega,
                arma::mat& zeta_mat,
                arma::mat& rho_temp,
                arma::mat& Z_mat,
                arma::vec& psi_e_inv,
                arma::mat& tau2_eta,
                const arma::vec& Queue,
                const arma::vec& NT){

  unsigned int SNT = X.n_rows;
  unsigned int Q = omega.n_cols;
  unsigned int R = X.n_cols;
  unsigned int M = tau2_eta.n_rows;
  
  NumericVector stat_vec(M);

  arma::mat xi(SNT, (Q+1));
  xi.col(0).fill(1.0); 
  xi.submat(0,1,SNT-1,Q) = omega;
  for(unsigned int m=0; m<M; ++m){
    for(unsigned int r=0; r<R; ++r){
      // ------------------------------------- covariace -------------------------------------
      arma::mat Sigma_rho((Q+1),(Q+1));
      arma::mat Mean_rho((Q+1),1);
      for(unsigned int j=0; j<(Q+1); ++j){
        Sigma_rho(j,j) = 1/tau2_eta(m,r);
        //Sigma_rho(j,j) = 1/tau2_eta(m);
      }
      // ----------------------------------------- difference --------------------------------------------------------
      if(m>0){
        arma::mat zeta_temp(1,(Q+1));
        for(unsigned int j=0; j<(Q+1); ++j){
          for(unsigned int s=0; s<m; ++s){
            zeta_temp(0,j) += rho_temp(s*R+r,j);
          }
        }
        //arma::mat zeta_temp = zeta_mat.row((m-1)*R+r);
        /*Rprintf("m: %d; r: %d; zeta_temp\n",(m+1),(r+1));
        Rcpp::Rcout << zeta_temp.row(0) << " ";
        Rcpp::Rcout << std::endl;*/
        for(unsigned int i=0; i<SNT; ++i){
          unsigned int index = Z_mat(i,0);
          if(index == (m+1)){
            arma::mat xx = xi.row(i).t() * xi.row(i);
            Sigma_rho = Sigma_rho + xx*1;
            // ------------------------------------- mean -------------------------------------
            arma::mat x3 = xi.row(i).t()*(X(i,r) - arma::dot(zeta_temp,xi.row(i)))*1;
            Mean_rho = Mean_rho + x3;
          }
        }
      // ----------------------------------------- 1st class --------------------------------------------------------
      }else{
        for(unsigned int i=0; i<SNT; ++i){
          unsigned int index = Z_mat(i,0);
          if(index == (m+1)){
            arma::mat xx = xi.row(i).t() * xi.row(i);
            Sigma_rho = Sigma_rho + xx*1;
            // ------------------------------------- mean -------------------------------------
            arma::mat x3 = xi.row(i).t() *X(i,r)* 1;          
            Mean_rho = Mean_rho + x3;
          }
        }
      }
      
      // ------------------------------------- sample -------------------------------------
      Sigma_rho = inv(Sigma_rho);
      Mean_rho = Sigma_rho*Mean_rho;
      rho_temp.row(m*R+r) = Mean_rho.t() + arma::randn<arma::rowvec>(1, (Q+1)) * arma::chol(Sigma_rho);
      /*for (unsigned int q = 0; q < (Q+1); ++q){
        rho_temp(m,q) = R::rnorm(Mean_rho(q,0), sqrt(Sigma_rho(q,q)));
      }*/

      /*Rprintf("m: %d; r: %d; Sigma_rho\n",(m+1),(r+1));
      for(unsigned int j=0; j<(Q+1); ++j){
        Rcpp::Rcout << Sigma_rho.row(j) << " ";
        Rcpp::Rcout << std::endl;
      }
      Rprintf("Mean_rho\n");
      Rcpp::Rcout << Mean_rho.t() << " ";
      Rcpp::Rcout << std::endl;
      Rprintf("rho\n");
      Rcpp::Rcout << rho_temp.row(m*R+r) << " ";
      Rcpp::Rcout << std::endl;*/
      /*for(unsigned int j=0; j<(Q+1); ++j){
        rho_length[m] += rho_temp((m*q1+k),j)*rho_temp((m*q1+k),j);
      }*/
    }
  }
}


// [[Rcpp::export]]
void f_rho_to_zeta(arma::mat& zeta_mat, 
                arma::mat& rho_temp,
                arma::vec& rho_length,
                const unsigned int M,
                const unsigned int R) {

  unsigned int Q = rho_temp.n_cols;
  
  if(M > 1){
    rho_length.fill(0.0);
    arma::mat new_zeta(M*R,Q);
    arma::mat new_rho(M*R,Q);
    for(unsigned int r = 0; r<R; ++r){
      for(unsigned int m=0; m<M; ++m){
        for (unsigned int j=0; j<Q; ++j){
          new_zeta(m*R+r,j) = rho_temp(r,j);
        }
        if(m>0){
          for(unsigned int s=1; s<=m; ++s){
            for (unsigned int j=0; j<Q; ++j){
              new_zeta(m*R+r,j) += rho_temp(s*R+r,j);
            }
          }
        }
      }
    }

    for(unsigned int r = 0; r<R; ++r){
      // ---------------------- give new order --------------------------
      NumericVector mua(M);
      for (unsigned int m = 0; m < M; ++m){
        mua[m] = new_zeta(m*R+r,0);
      }
      List all = f_judge(mua, M);
      IntegerVector norder = all[1];
      // ---------------------- assign new order to rho --------------------------
      for(unsigned int m=0; m<M; ++m){
        unsigned int new_m = norder[m]-1;
        for(unsigned int j=0; j<Q; ++j){
          zeta_mat(m*R+r,j) = new_zeta(new_m*R+r,j);
          new_rho(m*R+r,j) = rho_temp(new_m*R+r,j);
        }
        rho_length(m) += arma::dot(new_rho.row(m*R+r),new_rho.row(m*R+r).t());
      }
    }
  }else{
    for(unsigned int r = 0; r<R; ++r){
      for(unsigned int j=0; j<Q; ++j){
        zeta_mat(r,j) = rho_temp(r,j);
      }
      rho_length(0) += arma::dot(rho_temp.row(r),rho_temp.row(r).t());
    }
  }

  
}

// [[Rcpp::export]]
void update_tau2_eta(arma::mat& tau2_eta,
                     arma::mat& rho_temp,
                     arma::vec& psi_e_inv,
                     arma::vec& gamma2_eta) {
  unsigned int M = tau2_eta.n_rows;
  unsigned int R = tau2_eta.n_cols;
  
  arma::mat rho_norm_temp = rho_temp % rho_temp;

  for (unsigned int m = 0; m < M; ++m) {
    for(unsigned int r=0; r<R; ++r){
      double rho_norm = arma::sum(rho_norm_temp.row(m*R+r));
      double mean_tau = sqrt(gamma2_eta(m)/(rho_norm*psi_e_inv(m)));
      tau2_eta(m,r) = 1 / sim_invgauss(mean_tau, 1/gamma2_eta(m));
    }
  }
}

// [[Rcpp::export]]
void update_gamma2_eta(arma::vec& gamma2_eta,
                        arma::mat& tau2_eta,
                        const unsigned int& Q) {
  unsigned int M = tau2_eta.n_rows;
  unsigned int R = tau2_eta.n_cols;

  double a = 1.0 + double(Q+1.0)*double(R)/2.0;
  for(unsigned int m=0; m<M; ++m){
    double b = arma::sum(tau2_eta.row(m))/2.0 + 0.1;
    gamma2_eta(m) = rgamma(1,a, 1/b)[0];
  }
}


// [[Rcpp::export]]
void update_X(const arma::mat& V,
              arma::mat& X,
              arma::mat& omega,
              arma::mat& zeta_mat,
              arma::mat& Z_mat,
              const arma::vec& MMM_threshold) {

  unsigned int N = omega.n_rows;
  unsigned int Q = omega.n_cols;
  unsigned int R = X.n_cols;
  unsigned int M = (zeta_mat.n_rows)/R;

  
  //compute mean
  arma::mat xi(N, (Q+1));
  xi.col(0).fill(1.0); 
  xi.submat(0,1,N-1,Q) = omega;

  arma::mat mean_X(N,R);
  for (unsigned int i = 0; i < N; ++i) {
    unsigned int index = Z_mat(i,0) - 1;
    for(unsigned int r=0; r<R; ++r){
      mean_X(i,r) = arma::dot(zeta_mat.row(index*R+r),xi.row(i).t());
      if(V(i,r) == 1.0){
        X(i, r) = sim_rtnorm(mean_X(i, r), 1.0, MMM_threshold(0), R_PosInf);
      }else{
        X(i, r) = sim_rtnorm(mean_X(i, r), 1.0, R_NegInf, MMM_threshold(0));
      }
      //X(i, 0) = mean_X(i, 0) + rnorm(1,0,1)[0];
      //X.row(i) = mean_X.row(i) + arma::randn<arma::rowvec>(1, 1);
    }
  }
}


// [[Rcpp::export]]
double AIC(const arma::mat& V,
                const arma::mat& omega,
                const arma::mat& zeta_mat,
                const arma::mat& Z_mat,
                const arma::mat& g_mat,
                const arma::vec& pi,
                double delta_0,
                const arma::vec& psi_e_inv,
                const arma::vec& rho_length,
                double AIC_diff){
  unsigned int SNT = omega.n_rows;
  unsigned int N = g_mat.n_rows;
  unsigned int Q = omega.n_cols;
  unsigned int R = V.n_cols;
  unsigned int M = (zeta_mat.n_rows)/R;
  unsigned int index_smallest = arma::index_min(rho_length);
  
  if(index_smallest == 0){
    AIC_diff = 100.0;
  }else{
    //compute mean
    arma::mat xi(SNT, (Q+1));
    xi.col(0).fill(1.0); 
    xi.submat(0,1,SNT-1,Q) = omega;
    
    arma::mat mean_X1(SNT,R);
    arma::mat mean_X2(SNT,R);
    for (unsigned int i = 0; i < SNT; ++i) {
      unsigned int index = Z_mat(i,0) - 1;
      if(index == index_smallest){
        for(unsigned int r=0; r<R; ++r){
          mean_X1(i,r) = arma::dot(zeta_mat.row((index-1)*R+r),xi.row(i).t());
          mean_X2(i,r) = arma::dot(zeta_mat.row(index*R+r),xi.row(i).t());
        }
      }else{
        for(unsigned int r=0; r<R; ++r){
          mean_X1(i,r) = arma::dot(zeta_mat.row(index*R+r),xi.row(i).t());
          mean_X2(i,r) = arma::dot(zeta_mat.row(index*R+r),xi.row(i).t());
        }
      }
    }
    //compute likelihood
    arma::mat prob1 = arma::normcdf(mean_X1);
    arma::mat prob2 = prob1 % V + (1 - prob1) % (1 - V);
    arma::mat logProb_0 = arma::log(prob2);
    double likelihood_01 = arma::accu(logProb_0);
    
    arma::mat prob3 = arma::normcdf(mean_X2);
    arma::mat prob4 = prob3 % V + (1 - prob3) % (1 - V);
    arma::mat logProb_1 = arma::log(prob4);
    double likelihood_11 = arma::accu(logProb_1);
    double likelihood0 = (double((M-1)*R*(Q+1)))*2.0 - 2.0*likelihood_01;
    double likelihood1 = (double(M*R*(Q+1)))*2.0 - 2.0*likelihood_11;
    AIC_diff =  likelihood0 - likelihood1;
  }
  return AIC_diff;
  
}

///////////////////////////////////////////////////////////////////
////////////////// Main function: MCMC  ///////////////////////////
///////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
Rcpp::List mcmc(List observed_data,
                unsigned int PO,
                unsigned int Q,
                unsigned int M,
                List mh,
                List para,
                List latent,
                List acc_list){

  // MH
  double c_delta = mh[0];
  double my_gamma = mh[1];
  const arma::vec sigma_alpha = mh[2];

  // Data
  const arma::mat W = Rcpp::as<arma::mat>(observed_data[0]);
  const arma::mat V = Rcpp::as<arma::mat>(observed_data[1]);
  const arma::vec beta = Rcpp::as<arma::vec>(observed_data[2]);
  const arma::vec Queue = Rcpp::as<arma::vec>(observed_data[3]);
  const arma::vec NT = Rcpp::as<arma::vec>(observed_data[4]);
  const arma::vec beta_cumsum = Rcpp::as<arma::vec>(observed_data[5]);
  const arma::vec MMM_threshold = Rcpp::as<arma::vec>(observed_data[6]);
  unsigned int P = W.n_cols;
  
  arma::mat I_K = arma::eye(Q,Q);

  // Parameters
  // -------------- EFA ----------------
  arma::mat b_temp = Rcpp::as<arma::mat>(para[0]);
  arma::vec psi_eps_inv = Rcpp::as<arma::vec>(para[1]);
  arma::mat tau2 = Rcpp::as<arma::mat>(para[2]);
  arma::mat ukj = Rcpp::as<arma::mat>(para[3]);
  arma::mat p_ukj = Rcpp::as<arma::mat>(para[4]);
  arma::vec gammaj = Rcpp::as<arma::vec>(para[5]);
  arma::vec p_gammaj = Rcpp::as<arma::vec>(para[6]);
  arma::vec alpha = Rcpp::as<arma::vec>(para[7]);
  // ------------ MMM ----------------
  arma::mat zeta_mat = Rcpp::as<arma::mat>(para[8]);
  arma::vec psi_e_inv = Rcpp::as<arma::vec>(para[9]);
  double delta_0 = para[10];
  arma::vec pi = Rcpp::as<arma::vec>(para[11]);
  arma::mat g_mat = Rcpp::as<arma::mat>(para[12]);
  arma::mat tau2_eta = Rcpp::as<arma::mat>(para[13]);
  arma::vec gamma2_eta = Rcpp::as<arma::vec>(para[14]);
  

  // latent variables
  arma::mat omega = Rcpp::as<arma::mat>(latent[0]);
  arma::mat Y = Rcpp::as<arma::mat>(latent[1]);
  arma::mat Z_mat = Rcpp::as<arma::mat>(latent[2]);
  arma::mat X = Rcpp::as<arma::mat>(latent[3]);

  // Acceptance
  double at_delta = acc_list[0];
  arma::vec at_alpha = Rcpp::as<arma::vec>(acc_list[1]);

  arma::uvec r_idx;
  arma::uvec one_to_J = arma::linspace<arma::uvec>(0, P-1, P);
  r_idx = one_to_J.subvec(0, Q-1);

  const unsigned int R = X.n_cols;
  // Main Algorithm Loop
  arma::mat rho_temp(M*R,(Q+1));
  arma::vec rho_length(M, arma::fill::zeros);
  double AIC_diff = 0.0;

  update_omega_matrix(Y, X, I_K, omega, b_temp, psi_eps_inv, zeta_mat, Z_mat, psi_e_inv, Queue, NT);

  // -------------  MMM part ----------------
  f_rho_temp(zeta_mat, rho_temp, M, R);

  update_rho(X, omega, zeta_mat, rho_temp, Z_mat, psi_e_inv, tau2_eta, Queue, NT);

  f_rho_to_zeta(zeta_mat, rho_temp, rho_length, M, R);

  update_delta(delta_0, pi, g_mat, at_delta, c_delta); 
  
  update_gi(pi, delta_0, Z_mat, g_mat, Queue, NT);

  update_Z(V, omega, zeta_mat, Z_mat, g_mat, psi_e_inv, Queue, NT, MMM_threshold);

  update_tau2_eta(tau2_eta, rho_temp, psi_e_inv, gamma2_eta);

  update_gamma2_eta(gamma2_eta, tau2_eta, Q);

  update_X(V, X, omega, zeta_mat, Z_mat, MMM_threshold);

  // -------------  EFA part ----------------
  update_alpha(W, beta, omega, b_temp, psi_eps_inv, alpha, sigma_alpha, beta_cumsum, PO, at_alpha, Y);

  update_b_temp(Y, r_idx, omega, b_temp, psi_eps_inv, gammaj, ukj, tau2, my_gamma);

  update_psi_eps(Y, omega, b_temp, psi_eps_inv);

  update_tau(b_temp, gammaj, ukj, tau2);

  update_ukj(b_temp, gammaj, ukj, p_ukj, tau2);

  update_p_ukj(ukj, p_ukj);

  update_gammaj(b_temp, gammaj, p_gammaj, ukj, p_ukj, tau2);

  update_p_gammaj(gammaj, p_gammaj);
  
  AIC_diff = AIC(V, omega, zeta_mat, Z_mat, g_mat, pi, delta_0, psi_e_inv, rho_length, AIC_diff);

  Rcpp::List para_list = Rcpp::List::create(Rcpp::Named("b_temp",b_temp),
                            Rcpp::Named("psi_eps_inv",psi_eps_inv),
                            Rcpp::Named("tau2", tau2),
                            Rcpp::Named("ukj", ukj),
                            Rcpp::Named("p_ukj", p_ukj),
                            Rcpp::Named("gammaj", gammaj),
                            Rcpp::Named("p_gammaj", p_gammaj),
                            Rcpp::Named("alpha",alpha),
                            Rcpp::Named("zeta_mat", zeta_mat),
                            Rcpp::Named("psi_e_inv", psi_e_inv),
                            Rcpp::Named("delta_0", delta_0),
                            Rcpp::Named("pi", pi),
                            Rcpp::Named("g_mat", g_mat),
                            Rcpp::Named("tau2_eta", tau2_eta),
                            Rcpp::Named("gamma2_eta", gamma2_eta));
  
  Rcpp::List M_list = Rcpp::List::create(
    Rcpp::Named("M") = M,
    Rcpp::Named("AIC_diff") = AIC_diff
  );

  Rcpp::List latent_list = Rcpp::List::create(
    Rcpp::Named("omega") = omega,
    Rcpp::Named("Y") = Y,
    Rcpp::Named("Z_mat") = Z_mat,
    Rcpp::Named("X") = X
  );

  Rcpp::List acc_list_new = Rcpp::List::create(
    Rcpp::Named("at_delta") = at_delta,
    Rcpp::Named("at_alpha") = at_alpha
  );

  return Rcpp::List::create(
    Rcpp::Named("para_list") = para_list,
    Rcpp::Named("latent_list") = latent_list,
    Rcpp::Named("M_list") = M_list,
    Rcpp::Named("acc_list") = acc_list_new
  );
}
