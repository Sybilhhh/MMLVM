############ Generate data ############
library(MASS)
library(Rcpp)
library(progress)
library(bit)
library(psych)
library(lavaan)
library(LaplacesDemon)

set.seed(7293)

############## MCMC ###############
REP = 1    ## number of replication
alg_iter = 40000 # Algorithm MCMC iterations
alg_burn_in = 20000

observed_data = readRDS(file = "observed_data.RData")
sourceCpp("MCMC.cpp")
sink("output.txt")
start_time <- Sys.time()
# Create a progress bar object
pb <- progress_bar$new(total = REP * alg_iter, format = "[:bar] :percent Remaining time: :eta (Rep: :replication)")
for(iter in 1:REP){
  num_list = 7
  
  ##################### Get rough initial value #########################
  Q_max = 15
  M_max = 6
  M_min = 2
  Q_max0 = 15
  M_max0 = 6
  
  save_result = list()
  acc_list = list()
  save_latent = list()
  
  W = as.data.frame(observed_data[[(iter-1)*num_list+1]])
  V = observed_data[[(iter-1)*num_list+2]]
  NT = observed_data[[(iter-1)*num_list+5]]
  SNT = sum(NT)
  N = length(NT)
  P = ncol(W)
  R = ncol(V)
  PO = P
  beta = observed_data[[(iter-1)*num_list+3]]
  beta_sum = observed_data[[(iter-1)*num_list+6]]
  Queue = observed_data[[(iter-1)*num_list+4]]
  
  c_delta = 0.05
  my_gamma = 0.5
  sigma_alpha = rep(0.05,P)
  
  mh = list()
  mh[[1]] = c_delta
  mh[[2]] = my_gamma
  mh[[3]] = sigma_alpha
  
  # Used in real data
  alpha_cpp = rep(NA, beta_sum[P+1])
  for(p in 1:P){
    beta_temp = beta[p]
    for(s in 1:beta[p]){
      alpha_cpp[beta_sum[p]+s] = qnorm(sum(W[,p]<=(s-1))/SNT)
    }
  }
  
  for(i in 1:Q_max0){
    if(i == 1){
      b_temp_cpp = matrix(1,nrow = P,ncol = i)
      omega_cpp = as.matrix(rnorm(SNT,0,1))
      psi_eps_inv_cpp = 1.0/rep(0.1,P)
      Y_cpp = omega_cpp %*% t(b_temp_cpp)
      
      tau2_cpp = matrix(1,nrow = P,ncol = i)
      gammaj_cpp = rep(1,i)
      p_gammaj_cpp = rep(0.5,i)
      ukj_cpp = matrix(0,nrow = P,ncol = i)
      p_ukj_cpp = matrix(0.1,nrow = P,ncol = i)
      for(k in 1:P){
        for(h in 1:i){
          if(b_temp_cpp[k,h]>0.2){
            ukj_cpp[k,h] = 1
            p_ukj_cpp[k,h] = 0.9
          }
        }
      }
    }else{
      'b_temp1 = diag(1,nrow = i)
      b_temp2 = matrix(1,nrow = (P-i),ncol = i)
      b_temp_cpp = rbind(b_temp1, b_temp2)
      omega_cpp = as.matrix(factor.scores(W,b_temp_cpp)$scores)
      psi_eps_inv_cpp = 1.0/rep(0.1,P)
      Y_cpp = omega_cpp %*% t(b_temp_cpp) + rmvnorm(SNT,rep(0,P),diag(1,P))'
      poly_model = efa(data = W, ordered = TRUE, nfactors = i,output = "efa")
      QR = qr(t(poly_model$loadings))
      L <- qr.R(QR)
      Q_mat <- qr.Q(QR)
      # Adjust the diagonal of L to be positive
      diagL <- diag(L)
      sign_diagL <- sign(diagL)
      Q_modified <- Q_mat %*% diag(sign_diagL)
      # Compute Q using modified L and R
      b_temp_cpp = poly_model$loadings %*%Q_modified
      b_temp_cpp = round(b_temp_cpp,3)
      for(s in 1:i){
        for(h in 1:P){
          if(b_temp_cpp[h,s] > 1.0){
            b_temp_cpp[h,s] = 1.0
          }
          if(b_temp_cpp[h,s] < -1.0){
            b_temp_cpp[h,s] = -1.0
          }
        }
        if(abs(b_temp_cpp[s,s])<0.1){
          b_temp_cpp[s,s] = 0.1
        }
      }
      omega_cpp = as.matrix(factor.scores(W,b_temp_cpp)$scores)
      Y_cpp = omega_cpp %*% t(b_temp_cpp) + rmvnorm(SNT,rep(0,P),diag(1,P))
      
      tau2_cpp = matrix(1,nrow = P,ncol = i)
      gammaj_cpp = rep(1,i)
      p_gammaj_cpp = rep(0.5,i)
      ukj_cpp = matrix(0,nrow = P,ncol = i)
      p_ukj_cpp = matrix(0.1,nrow = P,ncol = i)
      for(k in 1:P){
        for(h in 1:i){
          if(b_temp_cpp[k,h]>0.2){
            ukj_cpp[k,h] = 1
            p_ukj_cpp[k,h] = 0.9
          }
        }
      }
    }
    
    for(j in 1:M_max0){
      temp_a = sort(sample(seq(-0.5,0.5,0.1),j))
      temp_b = sort(sample(seq(-0.5,0.5,0.1),j))
      zeta_cpp = matrix(NA,nrow = j*R,ncol = (i+1), byrow = T)
      for(r in 1:R){
        for (h in 1:j) {
          zeta_cpp[((h-1)*R+r),] = c(temp_a[h],rep(temp_b[h],i))
        }
      }
      psi_ecpp = rep(1.0,j)
      delta_0_cpp = 1
      Pi_cpp = rep(1,j)/j
      if(j>1){
        g_cpp = rdirichlet(N,delta_0_cpp*Pi_cpp)
        Z_cpp = matrix(1,nrow = SNT, ncol = 1)
        for(s in 1:N){
          Z_cpp[(Queue[s]+1):(Queue[s]+NT[s]),] = sample(1:j, NT[s], replace = TRUE, prob = g_cpp[s,])
        }
      }else{
        g_cpp = matrix(1,nrow = N, ncol = 1)
        Z_cpp = matrix(1,nrow = SNT, ncol = 1)
      }
      tau2_eta_cpp = matrix(30,nrow = j, ncol = R)
      #tau2_eta_cpp = rep(30,j)
      gamma2_eta_cpp = rep(0.1,j)
      #X_cpp = true_X[[iter]]
      X_cpp = matrix(NA, nrow = SNT, ncol = R)
      for (s in 1:SNT) {
        index = Z_cpp[s,1] - 1
        for(r in 1:R){
          X_cpp[s,r] = zeta_cpp[(index*R+r),]%*%c(1,omega_cpp[s,])+rnorm(1,0,1)
        }
      }
      at_delta = 0
      at_alpha_cpp = rep(0,P)
      
      acc_list_temp = list()
      acc_list_temp[[1]] = at_delta
      acc_list_temp[[2]] = at_alpha_cpp
      
      para = list()
      para[[1]] = b_temp_cpp
      para[[2]] = psi_eps_inv_cpp
      para[[3]] = tau2_cpp
      para[[4]] = ukj_cpp
      para[[5]] = p_ukj_cpp
      para[[6]] = gammaj_cpp
      para[[7]] = p_gammaj_cpp
      para[[8]] = alpha_cpp
      
      para[[9]] = zeta_cpp
      para[[10]] = psi_ecpp
      para[[11]] = delta_0_cpp
      para[[12]] = Pi_cpp
      para[[13]] = g_cpp
      para[[14]] = tau2_eta_cpp
      para[[15]] = gamma2_eta_cpp
      
      latent = list()
      latent[[1]] = omega_cpp
      latent[[2]] = Y_cpp
      latent[[3]] = Z_cpp
      latent[[4]] = X_cpp
      
      save_latent[[(j-1)*Q_max0+i]] = latent
      save_result[[(j-1)*Q_max0+i]] = list(para)
      acc_list[[(j-1)*Q_max0+i]] = acc_list_temp
    }
  }
  
  # ---------------------- MCMC running -------------------------------
  iter_times = rep(1,Q_max0*M_max0)
  M_str_mcmc = rep(M_max,alg_iter)
  q_str_mcmc = rep(Q_max,alg_iter)
  for (i in 1:alg_iter) {
    #sourceCpp("MCMC.cpp")
    q_mcmc = q_str_mcmc[i]
    M_mcmc = M_str_mcmc[i]
    s = iter_times[(M_mcmc-1)*Q_max0+q_mcmc]
    update_para = clone(save_result[[(M_mcmc-1)*Q_max0+q_mcmc]][[s]])
    update_latent = clone(save_latent[[(M_mcmc-1)*Q_max0+q_mcmc]])
    acc_list_temp = clone(acc_list[[(M_mcmc-1)*Q_max0+q_mcmc]])
    result = mcmc(observed_data[((iter-1)*num_list+1):((iter-1)*num_list+7)], PO, q_mcmc, M_mcmc, mh, update_para, update_latent, acc_list_temp)
    # ------------ save parameters -------------------
    save_result[[(M_mcmc-1)*Q_max0+q_mcmc]][[s+1]] = result$para_list
    iter_times[(M_mcmc-1)*Q_max0+q_mcmc] = s+1
    acc_list[[(M_mcmc-1)*Q_max0+q_mcmc]] = result$acc_list
    save_latent[[(M_mcmc-1)*Q_max0+q_mcmc]] = result$latent_list
    # ------------ judge the iteration of q -------------
    q_temp = sum(result$para_list$gammaj == 1)
    if(q_temp < Q_max){
      Q_max = min(q_mcmc,Q_max)
      q_str_mcmc[i+1] = max(q_temp,1)
      # ------------ judge the iteration of M -------------
      M_list = result$M_list
      M_str_mcmc[i+1] = M_list$M
    }else{
      q_str_mcmc[i+1] = q_mcmc
      # ------------ judge the iteration of M -------------
      M_list = result$M_list
      M_str_mcmc[i+1] = M_list$M
      AIC_diff = M_list$AIC_diff
      if(AIC_diff>=0){
        M_str_mcmc[i+1] = min(M_mcmc+1,M_max-1)
      }else if(AIC_diff<0){
        M_max = min(M_mcmc,M_max)
        M_str_mcmc[i+1] = max(M_mcmc-1,M_min)
      }
    }
    
    # Update the progress bar
    pb$tick(tokens = list(replication = iter))
  }
  
  q_final = which.max(tabulate(q_str_mcmc[alg_burn_in:alg_iter]))
  M_final = which.max(tabulate(M_str_mcmc[alg_burn_in:alg_iter]))
}
end_time <- Sys.time()
run_time = end_time - start_time

run_time

q_final
M_final

sink()

save.image(file = "real_LME.RData")
