N = 500 # sample size
P = 12 # dimension of y/observed variables
PO = 12 # dimension of y/observed variables of ordinal data
PC = 0# dimension of y/observed variables of continuous data
Q = 3 # number of latent variables
M = 3 # number of extreme profile
MT = 12 # maximum number of repeated measurement of each patient
#beta = sample(2:4, PO, replace = TRUE) # Number of threshold for each observed variables
beta = c(rep(c(3,4),6))
beta_sum = c(0,cumsum(beta))
R = 3 # dimension of y/observed variables for MMM

#################### 1.True value and Parameters #########################################################
# ------ measurement models --------
psi_eps = diag(rep(0.3,P+1))[1:P,1:P] # 1/psi
b_label =  as.numeric(matrix(c(
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
  0, 0, 1), nrow = P, ncol = Q, byrow = T))
b_matrix = matrix(c(
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
# Threshold for each observed variables
alpha = rep(0,beta_sum[PO+1])
for(p in 1:PO){
  if(beta[p] == 4){
    alpha[(beta_sum[p] + 1):(beta_sum[p] + beta[p])] = c(-1.5,-0.5,0.5,1.5)
  }
  if(beta[p] == 3){
    alpha[(beta_sum[p] + 1):(beta_sum[p] + beta[p])] = c(-1.5,0,1.5)
  }
  if(beta[p] == 2){
    alpha[(beta_sum[p] + 1):(beta_sum[p] + beta[p])] = c(-1.5,1.5)
  }
}

# --------- structural equation models ----------
zeta = matrix(c(-1,-0.5,-0.5,-0.5,
                -1,-0.5,-0.5,-0.5,
                -1,-0.5,-0.5,-0.5,
                0,0.5,-0.5,-0.5,
                0,0.5,-0.5,-0.5,
                0,0.5,-0.5,-0.5,
                1,0.5,0.5,0.5,
                1,0.5,0.5,0.5,
                1,0.5,0.5,0.5),
              nrow = M*R,ncol = (Q+1), byrow = T)
psi_e = diag(rep(0.3,M+1))[1:M,1:M] # 1/psi
phi = diag(Q)
delta_0 = 0.5 # spread parameter
Pi = c(0.3,0.4,0.3) # average proportion
delta = delta_0 * Pi # concentration parameter


#################### 2.Saving parameters #########################################################
b_rep = array(0,dim = c(REP,P,Q))
psi_epsrep = array(0,dim=c(REP,P))
tau2_rep = array(0,dim = c(REP,P,Q))
gammaj_rep = array(0,dim = c(REP,Q))
p_gammaj_rep = array(0,dim = c(REP,Q))
ukj_rep = array(0,dim = c(REP,P,Q))
p_ukj_rep = array(0,dim = c(REP,P,Q))
alpha_rep = array(0,dim = c(REP,beta_sum[PO+1]))
q_str_rep = rep(Q,REP)

zeta_rep = array(0, dim = c(REP,M*R,(Q+1)))
delta_0_rep = rep(0,REP)
Pi_rep = array(0,dim = c(REP,M))
g_rep = array(0,dim = c(REP,N,M))
tau2_eta_rep = array(0,dim = c(REP,M,R))
gamma2_eta_rep = array(0,dim = c(REP,M))
M_str_rep = rep(2,REP)


c_delta = 0.1
my_gamma = 0.5
sigma_alpha = rep(c(0.05,0.03),6)
c_omega = 0.3

mh = list()
mh[[1]] = c_delta
mh[[2]] = my_gamma
mh[[3]] = sigma_alpha
mh[[4]] = c_omega

######################## generate data ##################################################
observed_data = list()
true_omega = list()
true_g = list()
true_Z = list()
true_Y = list()
true_X = list()

for(iter in 1:REP){
  # Generate different time for each individuals
  Queue = c(0)
  NT = numeric(N)
  for(i in 1:N){
    NT[i] = sample(8:MT,1)
    #NT[i] = MT
    Queue=c(Queue,sum(NT))
  }
  SNT = sum(NT)
  
  # Initialize data
  g = rdirichlet(N,delta)
  Z = matrix(NA,nrow = SNT, ncol = 1)
  e_test = matrix(NA,nrow = SNT, ncol = R)
  V = matrix(NA,nrow = SNT, ncol = R)
  X = matrix(NA,nrow = SNT, ncol = R)
  omega = matrix(rnorm(SNT*Q, 0,1), nrow = SNT, ncol = Q)
  
  mu = omega %*% t(b_matrix)
  eps_test = mvrnorm(SNT,rep(0,P),psi_eps)
  Y = mu + eps_test
  W = matrix(NA, nrow = SNT, ncol = P)
  ############  3. Sample data  #############
  for(i in 1:N){
    #---------- MMM ---------------
    # generate latent variable
    Z[(Queue[i]+1):(Queue[i+1]),1] = sample(1:M, NT[i], replace = TRUE, prob = g[i,])
    
    for (t in 1:NT[i]) {
      zi = Z[Queue[i]+t,1]
      for(r in 1:R){
        #Z[i,] = rep(1,Q1)
        e_test[Queue[i]+t,r] = rnorm(1,0,1)
        X[Queue[i]+t,r] = t(zeta[((zi-1)*R+r),])%*%c(1,omega[Queue[i]+t,]) + e_test[Queue[i]+t,r]
      }
    }
  }
  MMM_thre = c(median(X))
  
  for(i in 1:N){
    for (t in 1:NT[i]) {
      for(r in 1:R){
        if (X[Queue[i]+t,r]>MMM_thre[1]) V[Queue[i]+t,r]<-1 else V[Queue[i]+t,r]<-0 
      }
      # --------- EFA ---------------
      # generate response variable
      for (k in 1:PO) {
        alpha_temp = alpha[(beta_sum[k] + 1):(beta_sum[k] + beta[k])]
        W[Queue[i]+t,k] <- sum(Y[Queue[i]+t,k]>= alpha_temp)
      }
      'for (k in ((PO+1):P)) {
        W[Queue[i]+t,k] <- Y[Queue[i]+t,k]
      }'
    }
  }
  'for (k in 1:PO) {
    W[,k] = factor(W[,k],ordered = TRUE)
  }'
  colnames(W) <- paste0("V", 1:P)
  ############  4. Save data  #############
  num_list = 7
  
  observed_data[[(iter-1)*num_list+1]] = W
  observed_data[[(iter-1)*num_list+2]] = V
  observed_data[[(iter-1)*num_list+3]] = beta
  observed_data[[(iter-1)*num_list+4]] = Queue
  observed_data[[(iter-1)*num_list+5]] = NT
  observed_data[[(iter-1)*num_list+6]] = beta_sum
  observed_data[[(iter-1)*num_list+7]] = MMM_thre
  
  true_omega[[iter]] = omega
  true_g[[iter]] = g
  true_Z[[iter]] = Z
  true_Y[[iter]] = Y
  true_X[[iter]] = X
}

saveRDS(observed_data, file="observed_data.RData")
saveRDS(true_omega, file="true_omega.RData")
saveRDS(true_g, file="true_g.RData")
saveRDS(true_Z, file="true_Z.RData")
saveRDS(true_Y, file="true_Y.RData")
saveRDS(true_X, file="true_X.RData")