rm(list = ls())

rep_no = 1 
# Note: This number was changed from 1 to 300, for getting 300 repeated simulation results. The authors used a batch script to run the 300 repeated simulations in multiple cores parallel. 
set.seed(rep_no+1000)

# create output folder if it does not exist
OutputFolder = "../../output/Simulation_2/3_IPD-AD-pooled"
if (!file.exists(OutputFolder)){ dir.create(OutputFolder, showWarnings = TRUE, recursive = FALSE, mode = "0777") }
RDataFolder = paste0(OutputFolder,"/RData")
if (!file.exists(RDataFolder)){ dir.create(RDataFolder, showWarnings = TRUE, recursive = FALSE, mode = "0777") }

# Ask if diagnostic plots will be produced
DrawDiagnostics = TRUE # TRUE / FALSE
if (DrawDiagnostics==T){
  PlotFolder = paste0(OutputFolder,"/DiagnosticPlot")
  if (!file.exists(PlotFolder)){ dir.create(PlotFolder, showWarnings = TRUE, recursive = FALSE, mode = "0777") }
} # 

# libraries 
library(ModelMetrics) ; library(mvtnorm) ; library(invgamma) ; library(MCMCpack)

#####################################
# Load 300 repeated datasets of Simulatino Study 2
#####################################

load("../../data/SimulationData_2.RData")
ls() 
# For the meaning of each object, refer to README in "data" folder of the github repository

#####################################
# Functions                         # 
#####################################

logit_inv = function(a) return( exp(a)/(exp(a)+1) ) ;

MB.logit.est=function(param, theta, x, w, alpha)
{  
  x = as.matrix(x); theta = as.vector(theta); param = as.vector(param)  
  # added for the density ratio model 
  alpha = as.vector(alpha) ; w.DRM= exp(x[,1:2]%*%alpha)
  w = w.DRM * w
  # added for the density ratio model 
  resid = w * (logit_inv(x%*%theta)-logit_inv(x[,-4]%*%param))
  g.i = cbind( x[,1]*resid, x[,2]*resid, x[,3]*resid ) 
  g.sum = apply( g.i,2,sum )
  return( sum(g.sum^2) )
}

L = dim(SimulData[[1]]$X_cube)[[1]] # total number of studies (IPD + AD) 
p_theta = dim(SimulData[[1]]$X_cube)[[3]] ; p_beta = p_theta - 1 ; p_alpha = 2

###############################################
# Setting Hyperparameters and starting values # 
###############################################

# Hyperparameters of prior distributions 
invLambda_theta = diag(1/10^4,p_theta)  # mu ~ N(0, Lambda_theta)
nu0 = 0.1 ; Phi0 = diag(0.1,p_theta)    # Sigma ~ InvWishart(0.1, 0.1 I)

# Divide simulation data into IPD dataset and AD dataset 
SEQ_IPD = which(SimulData[[rep_no]]$delta_biased_access==1)	
J = length(SEQ_IPD)
X_IPD = SimulData[[rep_no]]$X_cube[SEQ_IPD,,]
Y_mat = SimulData[[rep_no]]$Y_mat[SEQ_IPD,]
SEQ_AD = which(SimulData[[rep_no]]$delta_biased_access==0)
K = length(SEQ_AD)
beta_tilde_mat = SimulData[[rep_no]]$beta_mat[SEQ_AD,] 
V_tilde_cube = SimulData[[rep_no]]$V_beta_cube[SEQ_AD,,] 
for (kk in 1:K) { 
  V_tilde_cube[kk,,]=diag(diag(V_tilde_cube[kk,,])) 
} # for (kk)

# Made X_pooled
X_pooled = NULL
for (j in 1:J) {
  X_pooled = rbind(X_pooled, X_IPD[j,,])
}

# Reference data for the density ratio model
tilde_D_x = NULL ; hat_tau_vec = hat_Gamma_tau_vec = rep(0,K)
for (k in 1:K){	
  i_study = SEQ_AD[k]
  X_AD_l = SimulData[[ rep_no ]]$X_cube[i_study,,]
  hat_tau_vec[k] = mean( X_AD_l[,"X1"] )
  hat_Gamma_tau_vec[k] = var( X_AD_l[,"X1"] ) / length( X_AD_l[,"X1"] )
  SEQ = sample(c(1:nrow(X_pooled)), size=nrow(X_AD_l), replace=T) # from pooled X
  tilde_D_x[[k]] = X_pooled[SEQ,]
} # for (k)

# Starting values 
theta_mat_IPD = SimulData[[rep_no]]$theta_l_mat[SEQ_IPD,] 
mu_vec = true_mu
Sigma_theta_mat = diag(1,p_theta) 	
tau_vec = hat_Gamma_tau_vec # for the density ratio model
alpha_mat = array(0.1,c(K,p_alpha)) # for the density ratio model 
theta_mat_AD = SimulData[[rep_no]]$theta_l_mat[SEQ_AD,] 		
beta_mat_AD = beta_tilde_mat

###############################################
# Running MCMC
###############################################

# MCMC setting 
burnin = 10000 ; mainrun = 10000
n_iter = burnin + mainrun 
stepsize_theta_AD = c(0.2, 0.1, 0.1, 0.2) * 1.5
stepsize_theta_IPD = 0.15 ; stepsize_alpha = 0.005 ; stepsize_tau = 0.1

# Prepare repository for posterior draws 
draw_mu = array(0,c(n_iter,p_theta))
draw_Sigma_theta = array(0,c(n_iter,p_theta)) 

######
program_start_time = Sys.time()
format(Sys.time(), "%b %d %Y, %a, %H:%M:%S ")
Prevtime = proc.time()[3]

###### Start of MCMC 

for (i_iter in 1:n_iter) {
  
  ##################################
  # Update theta_AD, beta, and alpha 
  ##################################
  
  for (k in 1:K){
    
    # update theta and beta (Supplementary Material, Section 2, Step 2)
    
    theta_vec_q = rnorm(n=p_theta, mean=theta_mat_AD[k,], sd=stepsize_theta_AD)	
    ww = rnorm(nrow(tilde_D_x[[k]]),mean=1,sd=1) # multiplier bootstrap 			
    temp.optim = optim(beta_mat_AD[k,1:3], fn=MB.logit.est, theta=theta_vec_q, x=tilde_D_x[[k]], w=ww, control=list(reltol=1e-7, maxit=1000), alpha=alpha_mat[k,])		

    if ( temp.optim$convergence==0 ){
      
      beta_vec_q = temp.optim$par
      
      logNum = dmvnorm(beta_tilde_mat[k,2:3], mean=beta_vec_q[2:3], sigma=V_tilde_cube[k,2:3,2:3], log=TRUE)
      logNum = logNum + dmvnorm(theta_vec_q, mean=mu_vec, sigma=Sigma_theta_mat, log=TRUE) 
      
      logDen = dmvnorm(beta_tilde_mat[k,2:3], mean=beta_mat_AD[k,2:3], sigma=V_tilde_cube[k,2:3,2:3], log=TRUE)
      logDen = logDen + dmvnorm(theta_mat_AD[k,], mean=mu_vec, sigma=Sigma_theta_mat, log=TRUE)
      
      logAcc = logNum - logDen 
      
      if ( runif(n=1) < exp(logAcc) ) {
        theta_mat_AD[k,] = theta_vec_q
        beta_mat_AD[k,1:3] = beta_vec_q
      } # if ( runif(n=1) < exp(logAcc) )
      
    } # if ( temp.optim$convergence==0 )
    
    # alpha and beta
    
    alpha_vec_q = rnorm(n=p_alpha, mean=alpha_mat[k,], sd=stepsize_alpha)	# density ratio
    ww = rnorm(nrow(tilde_D_x[[k]]),mean=1,sd=1) # multiplier bootstrap 			
    temp.optim = optim(beta_mat_AD[k,1:3], fn=MB.logit.est, theta=theta_mat_AD[k,], x=tilde_D_x[[k]], w=ww, control=list(reltol=1e-7, maxit=1000), alpha=alpha_vec_q)		

    if ( temp.optim$convergence==0 ){
      
      beta_vec_q = temp.optim$par
      
      g_l1_vec_q = exp( tilde_D_x[[k]][,1:2]%*%alpha_vec_q )
      g_l_mat_q = cbind(g_l1_vec_q - 1, tilde_D_x[[k]][,"X1"] * g_l1_vec_q - tau_vec[k])
      bar_g_l_q = apply(g_l_mat_q,2,mean) ; Sigma_g_l_q = var(g_l_mat_q) / nrow(g_l_mat_q)

      logNum = dmvnorm(beta_tilde_mat[k,2:3], mean=beta_vec_q[2:3], sigma=V_tilde_cube[k,2:3,2:3], log=TRUE)
      logNum = logNum + dmvnorm(bar_g_l_q, mean=rep(0,2), sigma=Sigma_g_l_q, log=TRUE)
      
      g_l1_vec = exp( tilde_D_x[[k]][,1:2]%*%alpha_mat[k,] )
      g_l_mat = cbind(g_l1_vec - 1, tilde_D_x[[k]][,"X1"] * g_l1_vec - tau_vec[k])
      bar_g_l = apply(g_l_mat,2,mean) ; Sigma_g_l = var(g_l_mat) / nrow(g_l_mat)
      
      logDen = dmvnorm(beta_tilde_mat[k,2:3], mean=beta_mat_AD[k,2:3], sigma=V_tilde_cube[k,2:3,2:3], log=TRUE)
      logDen = logDen + dmvnorm(bar_g_l, mean=rep(0,2), sigma=Sigma_g_l, log=TRUE)
      
      logAcc = logNum - logDen 
      
      if ( runif(n=1) < exp(logAcc) ) {
        alpha_mat[k,] = alpha_vec_q # density ratio
        beta_mat_AD[k,1:3] = beta_vec_q
      }
      
    } # if ( temp.optim$convergence==0 )
    
    g_l1_vec = exp( tilde_D_x[[k]][,1:2]%*%alpha_mat[k,] )
    g_l_mat = cbind(g_l1_vec - 1, tilde_D_x[[k]][,"X1"] * g_l1_vec - tau_vec[k])
    bar_g_l = apply(g_l_mat,2,mean) ; Sigma_g_l = var(g_l_mat) / nrow(g_l_mat)
    
  } # for (k)
  
  ##################################
  # Update theta_IPD (Supplementary Material, Section 2, Step 1)
  ##################################
  
  for (j in 1:J){
    
    theta_vec_q = rnorm(n=p_theta, mean=theta_mat_IPD[j,], sd=stepsize_theta_IPD)	
    x.the_q = X_IPD[j,,] %*% theta_vec_q
    x.the = X_IPD[j,,] %*% theta_mat_IPD[j,]
    logAcc = sum( x.the_q * Y_mat[j,] - log( 1 + exp(x.the_q) ) - x.the * Y_mat[j,] + log( 1 + exp(x.the) ) )
    logAcc = logAcc + dmvnorm(theta_vec_q, mean=mu_vec, sigma=Sigma_theta_mat, log=TRUE)
    logAcc = logAcc - dmvnorm(theta_mat_IPD[j,], mean=mu_vec, sigma=Sigma_theta_mat, log=TRUE)
    
    if (runif(n=1)<exp(logAcc)){
      theta_mat_IPD[j,] = theta_vec_q
    }
    
  } # for (j)	
  
  theta_merge = rbind(theta_mat_IPD, theta_mat_AD)
  
  ##################################
  # Update mu_vec (Supplementary Material, Section 2, Step 5)
  ##################################
  
  inv_Sigma = solve(Sigma_theta_mat)
  inv_Var = invLambda_theta + L * inv_Sigma
  Var = solve(inv_Var)
  Mean_2ndpart = inv_Sigma %*% apply(theta_merge,2,sum)
  Mean = Var %*% Mean_2ndpart
  
  mu_vec = rmvnorm(n=1, mean=Mean, sigma=Var) 
  
  ##################################
  # Update Sigma_theta_mat (Supplementary Material, Section 2, Step 5)
  ##################################
  
  SS = array(0,c(p_theta,p_theta))
  for (l in 1:L){
    SS = SS + t(theta_merge[l,]-mu_vec) %*% t(t(theta_merge[l,]-mu_vec))
  } # 
  
  Sigma_theta_mat = riwish((nu0+L),(Phi0+SS))	
  
  ##################################
  # Update tau_l (Supplementary Material, Section 2, Step 4)
  ##################################
  
  for (k in 1:K){
    
    g_l1_vec = exp( tilde_D_x[[k]][,1:2]%*%alpha_mat[k,] )

    tau_q = rnorm(n=1, mean=tau_vec[k], sd=stepsize_tau)	
    
    g_l_mat_q = cbind(g_l1_vec - 1, tilde_D_x[[k]][,"X1"] * g_l1_vec - tau_q)
    bar_g_l_q = apply(g_l_mat_q,2,mean) ; Sigma_g_l_q = var(g_l_mat_q) / nrow(g_l_mat_q)
    logNum = dnorm(hat_tau_vec[k], mean=tau_q, sd=sqrt(hat_Gamma_tau_vec[k]), log=TRUE)
    logNum = logNum + dmvnorm(bar_g_l_q, mean=rep(0,2), sigma=Sigma_g_l_q, log=TRUE) 
    g_l_mat = cbind(g_l1_vec - 1, tilde_D_x[[k]][,"X1"] * g_l1_vec - tau_vec[k])
    bar_g_l = apply(g_l_mat,2,mean) ; Sigma_g_l = var(g_l_mat) / nrow(g_l_mat)
    logDen = dnorm(hat_tau_vec[k], mean=tau_vec[k], sd=sqrt(hat_Gamma_tau_vec[k]), log=TRUE)
    logDen = logDen + dmvnorm(bar_g_l, mean=rep(0,2), sigma=Sigma_g_l, log=TRUE) 
    
    logAcc = logNum - logDen 
    
    if ( runif(n=1) < exp(logAcc) ) {
      tau_vec[k] = tau_q
    }
    
  } # 
  
  ##################################
  # Store posterior draws (every iteration)
  ##################################
  
  draw_mu[i_iter,] = mu_vec
  draw_Sigma_theta[i_iter,] = diag(Sigma_theta_mat)
  
  ################################
  # Print iteration numbers (every 1K iterations) 
  #  and plot mu's (if DrawDiagnostics==TRUE)  
  ################################
  
  if (i_iter%%1000==0) {
    
    print( paste0( "i_rep = ",i_rep, ", iteration: ",i_iter," / ",n_iter) )
    Currenttime = proc.time()[3]
    LastBatch = Currenttime-Prevtime ; Time_to_Go = (n_iter-i_iter)*(LastBatch/1000)
    Prevtime = Currenttime
    print( paste("The last 1000 iter=",round(LastBatch/60,1),"min, Est. Time to go=",round(Time_to_Go/60,1),"min" ))
    
    if (DrawDiagnostics==T){
      
      png(file=paste0(PlotFolder,"/rep_",rep_no,"_mu.png"),width=1000,height=1800,pointsize=40)
      par(mfrow=c(4,1),mai=c(1.4,1.1,0.6,0.4),family="serif",mgp = c(1.5, 0.5, 0)) # b l t r 
      for (jj in 1:p_theta){
        plot(draw_mu[1:i_iter,jj], type="l", xlab="Iteration", ylab=paste0("mu",jj), main="mu")  
        abline(h=true_mu[jj], col="red", lwd=3)
        abline(v=burnin, col="blue", lty="dotted", lwd=3)
      }
      dev.off()
      
    } # if (DrawDiagnostics==T)
    
  } # if (i_iter%%1000)
  
} # for (i_iter)

###### End of MCMC 

# Save the posterior draws after burn-in
SEQ = (burnin+1):(burnin+mainrun)
posterior_mu = draw_mu[SEQ,]
posterior_Sigma_theta = draw_Sigma_theta[SEQ,]
save(posterior_mu,posterior_Sigma_theta,file=paste0(RDataFolder,"/rep_",rep_no,".RData"))