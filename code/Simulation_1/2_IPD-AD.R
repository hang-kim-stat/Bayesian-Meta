rm(list = ls())

rep_no = 1 
# Note: This number was changed from 1 to 300, for getting 300 repeated simulation results. The authors used a batch script to run the 300 repeated simulations in multiple cores parallel. 

set.seed(rep_no+1000)

FolderName = "2_IPD-AD"
if (!file.exists(FolderName)){ dir.create(FolderName, showWarnings = TRUE, recursive = FALSE, mode = "0777") }

library(ModelMetrics) ; library(mvtnorm) ; library(invgamma) ; library(MCMCpack)

load("0_SimulData.RData")

#####################################
# Functions                         # 
#####################################

logit = function(p) return( log(p)-log(1-p) ) ;
logit_inv = function(alpha) return( exp(alpha)/(exp(alpha)+1) ) ;

# Equation (10), integral of q_l for Type 1 AD
MB.est.ADtype1 = function(theta, x, w, alpha)
{  
 Wn = diag(w); Xn = as.matrix(x); 
 subgroup.Xn = as.matrix(x[,-4]); theta=as.vector(theta); 
 alpha=as.vector(alpha);  w.DRM=diag(as.vector(exp(as.matrix(x[,1:2])%*%alpha))) # added for the density ratio model
 Wn = Wn * w.DRM # added for the density ratio model 
 Mn1 = t(subgroup.Xn)%*%Wn%*%subgroup.Xn
 Mn2 = t(subgroup.Xn)%*%Wn%*%Xn%*%theta
 param = solve(Mn1)%*%Mn2
 return(param)
}

# Equation (10), integral of q_l for Type 2 AD
MB.est.ADtype2 = function(theta, x, w, subgroup.x, alpha)
{ 
 Wn = diag(w); Xn = as.matrix(x); 
 subgroup.Xn = as.matrix(subgroup.x); theta=as.vector(theta); 
 alpha=as.vector(alpha);  w.DRM=diag(as.vector(exp(as.matrix(x[,1:2])%*%alpha))) # added for the density ratio model 
 Wn = Wn * w.DRM # added for the density ratio model 
 Mn1 = t(subgroup.Xn)%*%Wn%*%subgroup.Xn
 Mn2 = t(subgroup.Xn)%*%Wn%*%Xn%*%theta
 param = solve(Mn1)%*%Mn2
 return(param)
}

# Equation (10), integral of q_l for Type 3 AD
MB.est.ADtype3 = function(theta, x, w, alpha)
{  
 Wn = diag(w); Xn = as.matrix(x); 
 theta=as.vector(theta); 
 alpha=as.vector(alpha);  w.DRM=diag(as.vector(exp(as.matrix(x[,1:2])%*%alpha))) # added for the density ratio model 
 Wn = Wn * w.DRM # added for the density ratio model 
 Mn1 = t(Xn)%*%Wn%*%Xn
 Mn2 = t(Xn)%*%Wn%*%Xn%*%theta
 param = solve(Mn1)%*%Mn2
 return(param)
}

L = dim(SimulData[[1]]$X_cube)[[1]] # total number of studies (IPD + AD)
p_theta = dim(SimulData[[1]]$X_cube)[[3]] ; p_beta = p_theta ; p_alpha = 2

###############################################
# Setting Hyperparameters and starting values # 
###############################################

# Hyperparameters of prior distributions 
invLambda_theta = diag(1/10^4,p_theta) # mu ~ N(0, Lambda_theta)
nu0 = 0.1 ; Phi0 = diag(0.1,p_theta) # Sigma ~ InvWishart(0.1, 0.1 I)

# Divide simulation data into IPD dataset and AD dataset 
SEQ_IPD = 31:40 ; J = length(SEQ_IPD)
X_IPD = SimulData[[rep_no]]$X_cube[SEQ_IPD,,]
Y_mat = SimulData[[rep_no]]$Y_mat[SEQ_IPD,]
SEQ_AD = 1:30 ; K = length(SEQ_AD)
beta_tilde_mat = SimulData[[rep_no]]$beta_mat[SEQ_AD,] 
V_tilde_cube = SimulData[[rep_no]]$V_beta_cube[SEQ_AD,,]
for (kk in 1:30) { 
 V_tilde_cube[kk,,] = diag(diag(V_tilde_cube[kk,,])) 
} # for (kk)

# Reference data for the density ratio model
tilde_D_x = tilde_D_x_subgroup = NULL ; hat_tau_vec = hat_Gamma_tau_vec = rep(0,length(SEQ_AD))
mean_X1_IPD = rep(0,J)
for (j in 1:J){
 mean_X1_IPD[j] = mean(X_IPD[j,,"X1"])
} # 
for (k in SEQ_AD) {	
 X_AD_l = SimulData[[rep_no]]$X_cube[k,,]	  
 hat_tau_vec[k] = mean( X_AD_l[,"X1"] )
 hat_Gamma_tau_vec[k] = var( X_AD_l[,"X1"] ) / length(X_AD_l[,"X1"])
 WHICH = which.min( abs( hat_tau_vec[k] - mean_X1_IPD ) )
 tilde_D_x_l = tilde_D_x[[k]] = X_IPD[WHICH,,]
 phi = rep(0, nrow(tilde_D_x_l))
 for (i in 1:nrow(tilde_D_x_l)) {
  if (tilde_D_x_l[i,"X1"] > 0 && tilde_D_x_l[i,"X2"] == 0) { phi[i] = 1 }
  else if (tilde_D_x_l[i,"X1"] > 0 && tilde_D_x_l[i,"X2"] == 1) { phi[i] = 2 }
  else if (tilde_D_x_l[i,"X1"] <= 0 && tilde_D_x_l[i,"X2"] == 0) { phi[i] = 3 } 
  else { phi[i] = 4 }
 } # for (i)
 ind.1 = 1*(phi == 1); ind.2 = 1*(phi == 2); ind.3 = 1*(phi == 3); ind.4 = 1*(phi == 4)
 tilde_D_x_subgroup[[k]] = as.data.frame(cbind(ind.1, ind.2, ind.3, ind.4))
} # for (k)

# Starting values 
theta_mat_IPD = SimulData[[rep_no]]$theta_l_mat[SEQ_IPD,] 
mu_vec = true_mu
Sigma_theta_mat = diag(1,p_theta) 
sig2 = true_kappa	
tau_vec = hat_Gamma_tau_vec # for the density ratio model
alpha_mat = array(0.1,c(30,p_alpha)) # for the density ratio model 
theta_mat_AD = SimulData[[rep_no]]$theta_l_mat[SEQ_AD,] 		
beta_mat_AD = beta_tilde_mat
for (k in 1:10) {	# type 1 AD
 beta_mat_AD[k,1:3] = MB.est.ADtype1(theta=theta_mat_AD[k,], x=tilde_D_x[[k]], w=rep(1,nrow(tilde_D_x[[k]])), alpha=alpha_mat[k,])
} 
for (k in 11:20) {	# type 2 AD
 beta_mat_AD[k,] = MB.est.ADtype2(theta=theta_mat_AD[k,], x=tilde_D_x[[k]], w=rep(1,nrow(tilde_D_x[[k]])), subgroup.x=tilde_D_x_subgroup[[k]], alpha=alpha_mat[k,])
}
for (k in 21:30) {	# type 3 AD
 beta_mat_AD[k,] = MB.est.ADtype3(theta=theta_mat_AD[k,], x=tilde_D_x[[k]], w=rep(1,nrow(tilde_D_x[[k]])), alpha=alpha_mat[k,])
} 

###############################################
# Running MCMC
###############################################

# MCMC setting 
burnin = 10000 ; mainrun = 10000
n_iter = burnin + mainrun 
stepsize_theta = 0.2 ; stepsize_alpha = 0.01 ; stepsize_tau = 0.02

draw_mu = array(0,c(n_iter,p_theta))
draw_theta = array(0,c(n_iter,(J+K),p_theta))
draw_sig2 = rep(0,n_iter)
draw_Sigma_theta = array(0,c(n_iter,p_theta)) 
draw_alpha = array(0,c(n_iter,K,p_alpha))
draw_tau = array(0,c(n_iter,K))

is_acc_theta_AD = array(0,c(n_iter,K))
is_acc_alpha = array(0,c(n_iter,K))
is_acc_theta_IPD = array(0,c(n_iter,J))
is_acc_tau = array(0,c(n_iter,K))

no_solve_failed = array(0,n_iter)

######
program_start_time = Sys.time()
format(Sys.time(), "%b %d %Y, %a, %H:%M:%S ")
Prevtime = proc.time()[3]

# iterations 
for (i_iter in 1:n_iter) {

 ##################################
 # Update theta, beta, and alpha 
 ##################################
 
 # Type 1 AD
 
 for (k in 1:10){
  
  # theta and beta
  
  theta_vec_q = rnorm(n=p_theta, mean=theta_mat_AD[k,], sd=stepsize_theta)	
  ww = rnorm(nrow(tilde_D_x[[k]]),mean=1,sd=1) # multiplier bootstrap 	
  
  beta_vec_q = tryCatch({
   MB.est.ADtype1(theta=theta_vec_q, x=tilde_D_x[[k]], w=ww, alpha=alpha_mat[k,])
  }, error = function(e) {
   NULL  # Return NULL if an error occurs
  })
  
  if (is.null(beta_vec_q)){
   
   no_solve_failed[i_iter] = no_solve_failed[i_iter] + 1
   
  } else {
   
   logNum = dmvnorm(beta_tilde_mat[k,2:3], mean=beta_vec_q[2:3], sigma=V_tilde_cube[k,2:3,2:3], log=TRUE)
   logNum = logNum + dmvnorm(theta_vec_q, mean=mu_vec, sigma=Sigma_theta_mat, log=TRUE) 
   
   logDen = dmvnorm(beta_tilde_mat[k,2:3], mean=beta_mat_AD[k,2:3], sigma=V_tilde_cube[k,2:3,2:3], log=TRUE)
   logDen = logDen + dmvnorm(theta_mat_AD[k,], mean=mu_vec, sigma=Sigma_theta_mat, log=TRUE)
   logAcc = logNum - logDen 
   
   if ( runif(n=1) < exp(logAcc) ) {
    is_acc_theta_AD[i_iter,k] = 1
    theta_mat_AD[k,] = theta_vec_q
    beta_mat_AD[k,1:3] = beta_vec_q
   }
   
  } # if (is.null(beta_vec_q))
  
  # alpha and beta
  
  alpha_vec_q = rnorm(n=p_alpha, mean=alpha_mat[k,], sd=stepsize_alpha)	# density ratio
  ww = rnorm(nrow(tilde_D_x[[k]]),mean=1,sd=1) # multiplier bootstrap 			
  
  beta_vec_q = tryCatch({
   MB.est.ADtype1(theta=theta_mat_AD[k,], x=tilde_D_x[[k]], w=ww, alpha=alpha_vec_q)
  }, error = function(e) {
   NULL  # Return NULL if an error occurs
  })
  
  if (is.null(beta_vec_q)){
   
   no_solve_failed[i_iter] = no_solve_failed[i_iter] + 1
   
  } else {
   
   Exp_alpha_psi_q = exp( tilde_D_x[[k]][,1:2]%*%alpha_vec_q )
   q_l_mat_q = cbind(Exp_alpha_psi_q - 1, tilde_D_x[[k]][,"X1"] * Exp_alpha_psi_q - tau_vec[k])
   bar_q_l_q = apply(q_l_mat_q,2,mean) ; Sigma_q_l_q = var(q_l_mat_q) / nrow(q_l_mat_q)
   
   logNum = dmvnorm(beta_tilde_mat[k,2:3], mean=beta_vec_q[2:3], sigma=V_tilde_cube[k,2:3,2:3], log=TRUE)
   logNum = logNum + dmvnorm(bar_q_l_q, mean=rep(0,2), sigma=Sigma_q_l_q, log=TRUE) 
   
   Exp_alpha_psi_k = exp( tilde_D_x[[k]][,1:2]%*%alpha_mat[k,] )
   q_l_mat = cbind(Exp_alpha_psi_k - 1, tilde_D_x[[k]][,"X1"] * Exp_alpha_psi_k - tau_vec[k])
   bar_q_l = apply(q_l_mat,2,mean) ; Sigma_q_l = var(q_l_mat) / nrow(q_l_mat)
   
   logDen = dmvnorm(beta_tilde_mat[k,2:3], mean=beta_mat_AD[k,2:3], sigma=V_tilde_cube[k,2:3,2:3], log=TRUE)
   logDen = logDen + dmvnorm(bar_q_l, mean=rep(0,2), sigma=Sigma_q_l, log=TRUE) 
   
   logAcc = logNum - logDen 
   
   if ( runif(n=1) < exp(logAcc) ) {
    is_acc_alpha[i_iter,k] = 1
    alpha_mat[k,] = alpha_vec_q # density ratio
    beta_mat_AD[k,1:3] = beta_vec_q
   }
   
  } # if (is.null(beta_vec_q))
  
 } # for (k)
 
 # Type 2 AD
 
 for (k in 11:20){
  
  # theta and beta
  
  theta_vec_q = rnorm(n=p_theta, mean=theta_mat_AD[k,], sd=stepsize_theta)	
  ww = rnorm(nrow(tilde_D_x[[k]]),mean=1,sd=1) # multiplier bootstrap 
  
  beta_vec_q = tryCatch({
   MB.est.ADtype2(theta=theta_vec_q, x=tilde_D_x[[k]], w=ww, subgroup.x=tilde_D_x_subgroup[[k]], alpha=alpha_mat[k,])
  }, error = function(e) {
   NULL  # Return NULL if an error occurs
  })
  
  if (is.null(beta_vec_q)){
   
   no_solve_failed[i_iter] = no_solve_failed[i_iter] + 1
   
  } else {
   
   logNum = dmvnorm(beta_tilde_mat[k,], mean=beta_vec_q, sigma=V_tilde_cube[k,,], log=TRUE)
   logNum = logNum + dmvnorm(theta_vec_q, mean=mu_vec, sigma=Sigma_theta_mat, log=TRUE) 
   
   logDen = dmvnorm(beta_tilde_mat[k,], mean=beta_mat_AD[k,], sigma=V_tilde_cube[k,,], log=TRUE)
   logDen = logDen + dmvnorm(theta_mat_AD[k,], mean=mu_vec, sigma=Sigma_theta_mat, log=TRUE)
   
   logAcc = logNum - logDen 
   
   if ( runif(n=1) < exp(logAcc) ) {
    is_acc_theta_AD[i_iter,k] = 1
    theta_mat_AD[k,] = theta_vec_q
    beta_mat_AD[k,] = beta_vec_q
   }
   
  } # if (is.null(beta_vec_q))
  
  # alpha and beta
  
  alpha_vec_q = rnorm(n=p_alpha, mean=alpha_mat[k,], sd=stepsize_alpha)	# density ratio
  ww = rnorm(nrow(tilde_D_x[[k]]),mean=1,sd=1) # multiplier bootstrap 	
  
  beta_vec_q = tryCatch({
   MB.est.ADtype2(theta=theta_mat_AD[k,], x=tilde_D_x[[k]], w=ww, subgroup.x=tilde_D_x_subgroup[[k]], alpha=alpha_vec_q)
  }, error = function(e) {
   NULL  # Return NULL if an error occurs
  })
  
  if (is.null(beta_vec_q)){
   
   no_solve_failed[i_iter] = no_solve_failed[i_iter] + 1
   
  } else {
   
   Exp_alpha_psi_q = exp( tilde_D_x[[k]][,1:2]%*%alpha_vec_q )
   q_l_mat_q = cbind(Exp_alpha_psi_q - 1, tilde_D_x[[k]][,"X1"] * Exp_alpha_psi_q - tau_vec[k])
   bar_q_l_q = apply(q_l_mat_q,2,mean) ; Sigma_q_l_q = var(q_l_mat_q) / nrow(q_l_mat_q)
   
   logNum = dmvnorm(beta_tilde_mat[k,], mean=beta_vec_q, sigma=V_tilde_cube[k,,], log=TRUE)
   logNum = logNum + dmvnorm(bar_q_l_q, mean=rep(0,2), sigma=Sigma_q_l_q, log=TRUE) 
   
   Exp_alpha_psi_k = exp( tilde_D_x[[k]][,1:2]%*%alpha_mat[k,] )
   q_l_mat = cbind(Exp_alpha_psi_k - 1, tilde_D_x[[k]][,"X1"] * Exp_alpha_psi_k - tau_vec[k])
   bar_q_l = apply(q_l_mat,2,mean) ; Sigma_q_l = var(q_l_mat) / nrow(q_l_mat)
   
   logDen = dmvnorm(beta_tilde_mat[k,], mean=beta_mat_AD[k,], sigma=V_tilde_cube[k,,], log=TRUE)
   logDen = logDen + dmvnorm(bar_q_l, mean=rep(0,2), sigma=Sigma_q_l, log=TRUE) 
   
   logAcc = logNum - logDen 
   
   if ( runif(n=1) < exp(logAcc) ) {
    is_acc_alpha[i_iter,k] = 1
    alpha_mat[k,] = alpha_vec_q # density ratio
    beta_mat_AD[k,] = beta_vec_q
   }
   
  } # if (is.null(beta_vec_q))
  
 } # for (k)
 
 # Type 3 AD
 
 for (k in 21:30){
  
  # theta and beta
  
  theta_vec_q = rnorm(n=p_theta, mean=theta_mat_AD[k,], sd=stepsize_theta)	
  ww = rnorm(nrow(tilde_D_x[[k]]),mean=1,sd=1) # multiplier bootstrap 
  
  beta_vec_q = tryCatch({
   beta_vec_q=MB.est.ADtype3(theta=theta_vec_q, x=tilde_D_x[[k]], w=ww, alpha=alpha_mat[k,]) 
  }, error = function(e) {
   NULL  # Return NULL if an error occurs
  })
  
  if (is.null(beta_vec_q)){
   
   no_solve_failed[i_iter] = no_solve_failed[i_iter] + 1
   
  } else {
   
   logNum = dmvnorm(beta_tilde_mat[k,3:4], mean=beta_vec_q[3:4], sigma=V_tilde_cube[k,3:4,3:4], log=TRUE)
   logNum = logNum + dmvnorm(theta_vec_q, mean=mu_vec, sigma=Sigma_theta_mat, log=TRUE)
   
   logDen = dmvnorm(beta_tilde_mat[k,3:4], mean=beta_mat_AD[k,3:4], sigma=V_tilde_cube[k,3:4,3:4], log=TRUE)
   logDen = logDen + dmvnorm(theta_mat_AD[k,], mean=mu_vec, sigma=Sigma_theta_mat, log=TRUE)
   
   logAcc = logNum - logDen 
   
   if ( runif(n=1) < exp(logAcc) ) {
    is_acc_theta_AD[i_iter,k] = 1
    theta_mat_AD[k,] = theta_vec_q
    beta_mat_AD[k,] = beta_vec_q
   }
   
  } # if (is.null(beta_vec_q))
  
  # alpha and beta
  
  alpha_vec_q = rnorm(n=p_alpha, mean=alpha_mat[k,], sd=stepsize_alpha)	# density ratio
  ww = rnorm(nrow(tilde_D_x[[k]]),mean=1,sd=1) # multiplier bootstrap 
  
  beta_vec_q = tryCatch({
   beta_vec_q = MB.est.ADtype3(theta=theta_mat_AD[k,], x=tilde_D_x[[k]], w=ww, alpha=alpha_vec_q) 
  }, error = function(e) {
   NULL  # Return NULL if an error occurs
  })
  
  if (is.null(beta_vec_q)){
   
   no_solve_failed[i_iter] = no_solve_failed[i_iter] + 1
   
  } else {
   
   Exp_alpha_psi_q = exp( tilde_D_x[[k]][,1:2]%*%alpha_vec_q )
   q_l_mat_q = cbind(Exp_alpha_psi_q - 1, tilde_D_x[[k]][,"X1"] * Exp_alpha_psi_q - tau_vec[k])
   bar_q_l_q = apply(q_l_mat_q,2,mean) ; Sigma_q_l_q = var(q_l_mat_q) / nrow(q_l_mat_q)
   
   logNum = dmvnorm(beta_tilde_mat[k,3:4], mean=beta_vec_q[3:4], sigma=V_tilde_cube[k,3:4,3:4], log=TRUE)
   logNum = logNum + dmvnorm(theta_mat_AD[k,], mean=mu_vec, sigma=Sigma_theta_mat, log=TRUE)
   logNum = logNum + dmvnorm(bar_q_l_q, mean=rep(0,2), sigma=Sigma_q_l_q, log=TRUE) 
   
   Exp_alpha_psi_k = exp( tilde_D_x[[k]][,1:2]%*%alpha_mat[k,] )
   q_l_mat = cbind(Exp_alpha_psi_k - 1, tilde_D_x[[k]][,"X1"] * Exp_alpha_psi_k - tau_vec[k])
   bar_q_l = apply(q_l_mat,2,mean) ; Sigma_q_l = var(q_l_mat) / nrow(q_l_mat)
   
   logDen = dmvnorm(beta_tilde_mat[k,3:4], mean=beta_mat_AD[k,3:4], sigma=V_tilde_cube[k,3:4,3:4], log=TRUE)
   logDen = logDen + dmvnorm(theta_mat_AD[k,], mean=mu_vec, sigma=Sigma_theta_mat, log=TRUE)
   logDen = logDen + dmvnorm(bar_q_l, mean=rep(0,2), sigma=Sigma_q_l, log=TRUE) 
   
   logAcc = logNum - logDen 
   
   if ( runif(n=1) < exp(logAcc) ) {
    is_acc_alpha[i_iter,k] = 1
    alpha_mat[k,] = alpha_vec_q # density ratio
    beta_mat_AD[k,] = beta_vec_q
   }
   
  } # if (is.null(beta_vec_q))
  
 } # for (k)
 
 # Type 4 - ID
 
 for (j in 1:J){
  
  theta_vec_q = rnorm(n=p_theta, mean=theta_mat_IPD[j,], sd=stepsize_theta)	
  x.the_q = X_IPD[j,,1:4] %*% theta_vec_q
  x.the = X_IPD[j,,1:4] %*% theta_mat_IPD[j,]
  logAcc = sum( dnorm(Y_mat[j,], x.the_q, sqrt(sig2), log=TRUE) - dnorm(Y_mat[j,], x.the, sqrt(sig2), log=TRUE) )
  logAcc = logAcc + dmvnorm(theta_vec_q, mean=mu_vec, sigma=Sigma_theta_mat, log=TRUE)
  logAcc = logAcc - dmvnorm(theta_mat_IPD[j,], mean=mu_vec, sigma=Sigma_theta_mat, log=TRUE)
  
  if (runif(n=1)<exp(logAcc)){
   is_acc_theta_IPD[i_iter,j] = 1
   theta_mat_IPD[j,] = theta_vec_q
  }
  
 } # for (j)	
 
 theta_merge = rbind(theta_mat_IPD, theta_mat_AD)
 
 ##################################
 # Update mu_vec 
 ##################################
 
 inv_Sigma = solve(Sigma_theta_mat)
 inv_Var = invLambda_theta + (J+K) * inv_Sigma
 Var = solve(inv_Var)
 Mean_2ndpart = inv_Sigma %*% apply(theta_merge,2,sum)
 Mean = Var %*% Mean_2ndpart
 
 mu_vec = rmvnorm(n=1, mean=Mean, sigma=Var) 
 
 ##################################
 # Update Sigma_theta_mat 
 ##################################
 
 SS = array(0,c(p_theta,p_theta))
 for (l in 1:(J+K)){
  SS = SS + t(theta_merge[l,]-mu_vec) %*% t(t(theta_merge[l,]-mu_vec))
 } # 
 
 Sigma_theta_mat = riwish((nu0+J+K),(Phi0+SS))	
 
 ##################################
 # Update sig2 (IPD only for j=1,...,J)
 ##################################
 
 llik = rep(0,J)
 for (j in 1:J){
  x.the = X_IPD[j,,1:4] %*% theta_mat_IPD[j,]
  llik[j] = sum( (Y_mat[j,] - x.the)^2 )
 }
 sum_llik = sum( llik )
 
 sig2 = 1.0 / rgamma(1, shape = (1+prod(dim(Y_mat))/2), rate = (1+sum_llik/2) ) # rate !!!
 
 ##################################
 # Update tau_l (AD only)
 ##################################
 
 for (k in 1:30){
  
  Exp_alpha_psi_k = exp( tilde_D_x[[k]][,1:2]%*%alpha_mat[k,] ) # common for q and current
  
  tau_q = rnorm(n=1, mean=tau_vec[k], sd=stepsize_tau)	
  
  q_l_mat_q = cbind(Exp_alpha_psi_k - 1, tilde_D_x[[k]][,"X1"] * Exp_alpha_psi_k - tau_q)
  bar_q_l_q = apply(q_l_mat_q,2,mean) ; Sigma_q_l_q = var(q_l_mat_q) / nrow(q_l_mat_q)
  logNum = dnorm(hat_tau_vec[k], mean=tau_q, sd=sqrt(hat_Gamma_tau_vec[k]), log=TRUE)
  logNum = logNum + dmvnorm(bar_q_l_q, mean=rep(0,2), sigma=Sigma_q_l_q, log=TRUE) 
  q_l_mat = cbind(Exp_alpha_psi_k - 1, tilde_D_x[[k]][,"X1"] * Exp_alpha_psi_k - tau_vec[k])
  bar_q_l = apply(q_l_mat,2,mean) ; Sigma_q_l = var(q_l_mat) / nrow(q_l_mat)
  logDen = dnorm(hat_tau_vec[k], mean=tau_vec[k], sd=sqrt(hat_Gamma_tau_vec[k]), log=TRUE)
  logDen = logDen + dmvnorm(bar_q_l, mean=rep(0,2), sigma=Sigma_q_l, log=TRUE) 
  
  logAcc = logNum - logDen 
  
  if ( runif(n=1) < exp(logAcc) ) {
   is_acc_tau[i_iter,k] = 1
   tau_vec[k] = tau_q
  }
  
 } # 
 
 ################
 # Store 
 
 draw_mu[i_iter,] = mu_vec
 draw_Sigma_theta[i_iter,] = diag(Sigma_theta_mat)
 draw_sig2[i_iter] = sig2
 
 # Print iteration numbers and plot mu's  
 
 if (i_iter%%1000==0) {
  
  print( paste0( "rep_no = ",rep_no, ", iteration: ",i_iter," / ",n_iter) )
  Currenttime = proc.time()[3]
  LastBatch = Currenttime-Prevtime ; Time_to_Go = (n_iter-i_iter)*(LastBatch/1000)
  Prevtime = Currenttime
  print( paste("The last 1000 iter=",round(LastBatch/60,1),"min, Est. Time to go=",round(Time_to_Go/60,1),"min" ))
  
  png(file=paste0(FolderName,"/W_rep_",rep_no,"_mu.png"),width=1000,height=1800,pointsize=40)
  par(mfrow=c(4,1),mai=c(1.4,1.1,0.6,0.4),family="serif",mgp = c(1.5, 0.5, 0)) # b l t r 
  for (jj in 1:p_theta){
   plot(draw_mu[1:i_iter,jj], type="l", xlab="Iteration", ylab=paste0("theta",jj), main="mu")  
   abline(h=true_mu[jj], col="red", lwd=3)
   abline(v=burnin, col="blue", lty="dotted", lwd=3)
  }
  dev.off()
  
  save(draw_mu,draw_sig2,draw_Sigma_theta,file=paste0(FolderName,"/Rep_",rep_no,".RData"))
  
 } # if (i_iter%%1000)
 
} # for (i_iter)

############## Iteration ends ##############
########################################################
