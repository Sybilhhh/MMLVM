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
REP = 100    ## number of replication
alg_iter = 6000 # Algorithm MCMC iterations
alg_burn_in = 3000
Q_max0 = 7
M_max0 = 7

source("generate_data.r")
sourceCpp("MCMC.cpp")
#sink("output.txt")
start_time <- Sys.time()
# Create a progress bar object
pb <- progress_bar$new(total = REP * alg_iter, format = "[:bar] :percent Remaining time: :eta (Rep: :replication)")
for(iter in 1:REP){
  num_list = 7
  
  ##################### Get rough initial value #########################
  Q_max = 7
  M_max = 7
  M_min = 2
  save_result = list()
  acc_list = list()
  save_latent = list()
  
  W0 = as.data.frame(observed_data[[(iter-1)*num_list+1]])
  Queue = observed_data[[(iter-1)*num_list+4]]
  NT = observed_data[[(iter-1)*num_list+5]]
  SNT = sum(NT)
  # Used in simulation
  alpha_cpp = rep(NA, beta_sum[P+1])
  for(p in 1:P){
    beta_temp = beta[p]
    alpha_cpp[beta_sum[p]+1] = -1.5
    alpha_cpp[beta_sum[p+1]] = 1.5
    if(beta_temp > 2){
      for(s in 2:(beta_temp-1)){
        alpha_cpp[beta_sum[p]+s] = qnorm(sum(W[,p]<=(s-1))/SNT)
      }
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
      poly_model = efa(data = W0, ordered = TRUE, nfactors = i,output = "efa")
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
      omega_cpp = as.matrix(factor.scores(W0,b_temp_cpp)$scores)
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
    'poly_model = efa(data = W0, ordered = TRUE, nfactors = i,output = "efa")
    #Y_cpp = lavPredict(poly_model,type = "yhat")
    QR = qr(t(poly_model$loadings))
    L <- qr.R(QR)
    Q_mat <- qr.Q(QR)
    # Adjust the diagonal of L to be positive
    diagL <- diag(L)
    sign_diagL <- sign(diagL)
    Q_modified <- Q_mat %*% diag(sign_diagL)
    # Compute Q using modified L and R
    if(i != 3){
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
    }else{
      b_temp_cpp = matrix(c(
        1, 0, 0, 
        0, 1, 0, 
        0, 0, 1, 
        0.8, 0, 0, 
        0.8, 0, 0, 
        0.8, 0, 0, 
        0, 0.8, 0,
        0, 0.8, 0, 
        0, 0.8, 0, 
        0, 0, 0.8, 
        0, 0, 0.8, 
        0, 0, 0.8), nrow = P, ncol = Q, byrow = T)
    }
    
    if(i < 5){
      omega_cpp_temp = lavPredict(poly_model)
      omega_cpp = as.matrix(omega_cpp_temp%*%Q_modified)
    }else{
      omega_cpp = matrix(rnorm(SNT*i,0,1),nrow = SNT,ncol = i)
      #omega_cpp = as.matrix(factor.scores(W,b_temp_cpp)$scores)
    }
   
    Y_cpp = omega_cpp %*% t(b_temp_cpp) + rmvnorm(SNT,rep(0,P),diag(0.1,P))
    psi_eps_inv_cpp = rep(1/0.1,P)
    tau2_cpp = matrix(1,nrow = P,ncol = i)
    ukj_cpp = matrix(c(
      1, 0, 0, 
      0, 1, 0, 
      0, 0, 1, 
      1, 0, 0, 
      1, 0, 0, 
      1, 0, 0, 
      0, 1, 0, 
      0, 1, 0, 
      0, 1, 0,
      0, 0, 1,
      0, 0, 1,
      0, 0, 1), nrow = P, ncol = Q, byrow = T)
    p_ukj_cpp = matrix(0.5,nrow = P,ncol = i)
    gammaj_cpp = rep(1,i)
    p_gammaj_cpp = rep(1,i)
    }'
    for(j in 1:M_max0){
      'g_cpp = replicate(1, true_g[[iter]], simplify=FALSE)[[1]]
      Z_cpp = replicate(1, true_Z[[iter]], simplify=FALSE)[[1]]
      X_cpp = true_X[[iter]]
      zeta_cpp = matrix(c(-1,-0.5,-0.5,-0.5,
                          -1,-0.5,-0.5,-0.5,
                          -1,-0.5,-0.5,-0.5,
                          0,0.5,-0.5,-0.5,
                          0,0.5,-0.5,-0.5,
                          0,0.5,-0.5,-0.5,
                          1,0.5,0.5,0.5,
                          1,0.5,0.5,0.5,
                          1,0.5,0.5,0.5),
                        nrow = M*R,ncol = (Q+1), byrow = T)
      psi_e_cpp = rep(1,M)
      delta_0_cpp = 1
      Pi_cpp = c(0.3,0.4,0.3)
      tau2_eta_cpp = matrix(30,nrow = j, ncol = R)
      gamma2_eta_cpp = rep(0.1,j)'
      if(j < 4){
        temp_a = sort(sample(seq(-1,1,1),j))
      }else{
        temp_a = sort(sample(seq(-0.5,0.5,0.1),j))
      }
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
  
  q_str_rep[iter] = q_temp
  M_final = which.max(tabulate(M_str_mcmc[alg_iter:alg_burn_in]))
  M_str_rep[iter] = M_final
  
  if(q_temp == Q){
    ns=iter_times[(M_final-1)*Q_max0+Q]
    #----------------- Save data ----------------------#
    #---------- Mixed Membership model---------------#
    if(ns<alg_burn_in){
      if(ns%%2==0){
        cut_off = 0.5*ns
      }else{
        cut_off = 0.5*(ns+1)
      }
    }else{
      cut_off = alg_burn_in
    }
    b_mcmc = array(0, dim = c(ns,P,Q))
    psi_epsmcmc = array(0,dim=c(ns,P))
    tau2_mcmc = array(0,dim = c(ns,P,Q))
    gammaj_mcmc = array(0,dim = c(ns,Q))
    p_gammaj_mcmc = array(0,dim = c(ns,Q))
    ukj_mcmc = array(0,dim = c(ns,P,Q))
    p_ukj_mcmc = array(0,dim = c(ns,P,Q))
    alpha_mcmc = array(0,dim = c(ns,beta_sum[PO+1]))
    
    for(i in 1:ns){
      b_mcmc[i,,] = save_result[[(M_final-1)*Q_max0+Q]][[i]][[1]]
      psi_epsmcmc[i,] = 1.0/save_result[[(M_final-1)*Q_max0+Q]][[i]][[2]]
      tau2_mcmc[i,,]= save_result[[(M_final-1)*Q_max0+Q]][[i]][[3]]
      ukj_mcmc[i,,]= save_result[[(M_final-1)*Q_max0+Q]][[i]][[4]]
      p_ukj_mcmc[i,,]= save_result[[(M_final-1)*Q_max0+Q]][[i]][[5]]
      gammaj_mcmc[i,]= save_result[[(M_final-1)*Q_max0+Q]][[i]][[6]]
      p_gammaj_mcmc[i,]= save_result[[(M_final-1)*Q_max0+Q]][[i]][[7]]
      alpha_mcmc[i,]= save_result[[(M_final-1)*Q_max0+Q]][[i]][[8]]
    }
    
    b_rep[iter,,] = apply(b_mcmc[cut_off:ns,,],2:3,mean)
    psi_epsrep[iter,] = apply(psi_epsmcmc[cut_off:ns,],2,mean)
    tau2_rep[iter,,] = apply(tau2_mcmc[cut_off:ns,,],2:3,mean)
    gammaj_rep[iter,] = apply(gammaj_mcmc[cut_off:ns,],2,mean)
    p_gammaj_rep[iter,] = apply(p_gammaj_mcmc[cut_off:ns,],2,mean)
    ukj_rep[iter,,] = apply(ukj_mcmc[cut_off:ns,,],2:3,mean)
    p_ukj_rep[iter,,] = apply(p_ukj_mcmc[cut_off:ns,,],2:3,mean)
    alpha_rep[iter,] = apply(alpha_mcmc[cut_off:ns,],2,mean)
    
    if(M_final == M){
      ns=iter_times[(M-1)*Q_max0+Q]
      #----------------- Save data ----------------------#
      #---------- Mixed Membership model---------------#
      if(ns<alg_burn_in){
        if(ns%%2==0){
          cut_off = 0.5*ns
        }else{
          cut_off = 0.5*(ns+1)
        }
      }else{
        cut_off = alg_burn_in
      }
      zeta_mcmc = array(0, dim = c(ns,M*R,(Q+1)))
      g_mcmc = array(0,dim=c(ns,N,M))
      delta_0_mcmc = rep(0,ns)
      Pi_mcmc = array(0, dim = c(ns,M))
      psi_e_mcmc = array(0,dim=c(ns,M))
      tau2_eta_mcmc = array(0,dim = c(ns,M,R))
      gamma2_eta_mcmc = array(0,dim = c(ns,M))
      
      for(i in 1:ns){
        zeta_mcmc[i,,] = save_result[[(M-1)*Q_max0+Q]][[i]][[9]]
        delta_0_mcmc[i] = save_result[[(M-1)*Q_max0+Q]][[i]][[11]]
        Pi_mcmc[i,] = save_result[[(M-1)*Q_max0+Q]][[i]][[12]]
        g_mcmc[i,,] = save_result[[(M-1)*Q_max0+Q]][[i]][[13]]
        tau2_eta_mcmc[i,,]= save_result[[(M-1)*Q_max0+Q]][[i]][[14]]
        gamma2_eta_mcmc[i,]= save_result[[(M-1)*Q_max0+Q]][[i]][[15]]
      }
      
      at_delta = acc_list[[(M-1)*Q_max0+Q]][[1]]/ns
      at_alpha = acc_list[[(M-1)*Q_max0+Q]][[2]]/ns
      
      zeta_rep[iter,,] = apply(zeta_mcmc[cut_off:ns,,],2:3,mean)
      delta_0_rep[iter] = mean(delta_0_mcmc[cut_off:ns])
      Pi_rep[iter,] = apply(Pi_mcmc[cut_off:ns,],2,mean)
      g_rep[iter,,] = apply(g_mcmc[cut_off:ns,,],2:3,mean)
      tau2_eta_rep[iter,,] = apply(tau2_eta_mcmc[cut_off:ns,,],2:3,mean)
      gamma2_eta_rep[iter,] = apply(gamma2_eta_mcmc[cut_off:ns,],2,mean)
    }
  }
}
end_time <- Sys.time()
run_time = end_time - start_time

################################# analysis #######################
b_mean = apply(b_rep,2:3,mean)
psi_eps_mean = apply(psi_epsrep,2,mean)
tau2_mean = apply(tau2_rep, 2:3, mean)
ukj_mean = apply(ukj_rep, 2:3, mean)
p_ukj_mean = apply(p_ukj_rep, 2:3, mean)
gammaj_mean = apply(p_gammaj_rep, 2, mean)
p_gammaj_mean = apply(gammaj_rep, 2, mean)
alpha_mean = apply(alpha_rep, 2, mean)

zeta_mean = apply(zeta_rep,2:3,mean)
delta_mean = mean(delta_0_rep) 
Pi_mean = apply(Pi_rep,2,mean)
tau2_eta_mean = apply(tau2_eta_rep, 2:3, mean)
gamma2_eta_mean = apply(gamma2_eta_rep, 2, mean)

round(b_mean-b_matrix,3)
round(psi_eps_mean - rep(0.3,P),3)
round(alpha_mean - alpha,3)

tau2_mean
ukj_mean
p_ukj_mean
gammaj_mean
p_gammaj_mean
  
round(zeta_mean - zeta,3)
round(delta_mean - delta_0,3)
round(Pi_mean - Pi,3)
apply(abs(g_rep[iter,,] - true_g[[iter]]), 2, mean)
round(tau2_eta_mean,3)
round(gamma2_eta_mean,5)

#sink()

run_time

save.image(file = "m3_100rep.RData")

write.table(q_str_mcmc, file = "q_mcmc1.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(M_str_mcmc, file = "M_mcmc1.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

q_str_rep
M_str_rep

