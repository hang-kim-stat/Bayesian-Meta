rm(list = ls())

rep_no = 1 
# Note: This number was changed from 1 to 300, for getting 300 repeated simulation results. The authors used a batch script to run the 300 repeated simulations in multiple cores parallel. 
set.seed(rep_no+1000)

# create output folder if it does not exist
OutputFolder = "../../Output/Simulation_1/4_IPD_only"
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
# Load 300 repeated datasets of Simulation Study 1
#####################################

load("../../Data/SimulationData_1.RData")
ls() 
# For the meaning of each object, refer to README_data in "Data" folder of the github repository

L = dim(SimulData[[1]]$X_cube)[[1]] # total number of studies (IPD + AD)
p_theta = dim(SimulData[[1]]$X_cube)[[3]] ; p_beta = p_theta 

###############################################
# Setting Hyperparameters and starting values # 
###############################################

invLambda_theta = diag(1/10^4,p_theta)  # mu ~ N(0, Lambda_theta)
nu0 = 0.1 ; Phi0 = diag(0.1,p_theta)  # Sigma ~ InvWishart(0.1, 0.1 I)

# Extract IPD dataset from the (stacked) simulation dataset
SEQ_IPD = 31:40 ; J = length(SEQ_IPD) # 31:40
X_IPD = SimulData[[rep_no]]$X_cube[SEQ_IPD,,]
Y_mat = SimulData[[rep_no]]$Y_mat[SEQ_IPD,]

# Starting values 
theta_mat_IPD = SimulData[[rep_no]]$theta_l_mat[SEQ_IPD,] 
mu_vec = true_mu
Sigma_theta_mat = diag(1,p_theta) 
sig2 = true_kappa	

###############################################
# Running MCMC
###############################################

# MCMC setting 
burnin = 10000 ; mainrun = 10000
n_iter = burnin + mainrun 
stepsize_theta = 0.15

# Prepare repository for posterior draws 
draw_mu = array(0,c(n_iter,p_theta))
draw_Sigma_theta = array(0,c(n_iter,p_theta)) 
draw_sig2 = rep(0,n_iter)

######
program_start_time = Sys.time()
format(Sys.time(), "%b %d %Y, %a, %H:%M:%S ")
Prevtime = proc.time()[3]

###### Start of MCMC 

for (i_iter in 1:n_iter) {
  
  ##################################
  # Update theta_mat_IPD
  ##################################
  
  for (j in 1:J){
    
    theta_vec_q = rnorm(n=p_theta, mean=theta_mat_IPD[j,], sd=stepsize_theta)	
    x.the_q = X_IPD[j,,1:4] %*% theta_vec_q
    x.the = X_IPD[j,,1:4] %*% theta_mat_IPD[j,]
    logAcc = sum( dnorm(Y_mat[j,], x.the_q, sqrt(sig2), log=TRUE) - dnorm(Y_mat[j,], x.the, sqrt(sig2), log=TRUE) )
    logAcc = logAcc + dmvnorm(theta_vec_q, mean=mu_vec, sigma=Sigma_theta_mat, log=TRUE)
    logAcc = logAcc - dmvnorm(theta_mat_IPD[j,], mean=mu_vec, sigma=Sigma_theta_mat, log=TRUE)
    
    if (runif(n=1)<exp(logAcc)){
      theta_mat_IPD[j,] = theta_vec_q
    }
    
  } # for (j)	
  
  ##################################
  # Update mu_vec 
  ##################################
  
  inv_Sigma = solve(Sigma_theta_mat)
  inv_Var = invLambda_theta + J * inv_Sigma
  Var = solve(inv_Var)
  Mean_2ndpart = inv_Sigma %*% apply(theta_mat_IPD,2,sum)
  Mean = Var %*% Mean_2ndpart
  
  mu_vec = rmvnorm(n=1, mean=Mean, sigma=Var) 
  
  ##################################
  # Update Sigma_theta_mat 
  ##################################
  
  SS = array(0,c(p_theta,p_theta))
  for (l in 1:J){
    SS = SS + t(theta_mat_IPD[l,]-mu_vec) %*% t(t(theta_mat_IPD[l,]-mu_vec))
  } # 
  
  Sigma_theta_mat = riwish((nu0+J),(Phi0+SS))	
  
  ##################################
  # Update sig2
  ##################################
  
  llik = rep(0,J)
  for (j in 1:J){
    x.the = X_IPD[j,,1:4] %*% theta_mat_IPD[j,]
    llik[j] = sum( (Y_mat[j,] - x.the)^2 )
  }
  sum_llik = sum( llik )
  
  sig2 = 1.0 / rgamma(1, shape = (1+prod(dim(Y_mat))/2), rate = (1+sum_llik/2) ) 
  
  ##################################
  # Store posterior draws (every iteration)
  ##################################
  
  draw_mu[i_iter,] = mu_vec
  draw_Sigma_theta[i_iter,] = diag(Sigma_theta_mat)
  draw_sig2[i_iter] = sig2
  
  ################################
  # Print iteration numbers (every 1K iterations) 
  #  and plot mu's (if DrawDiagnostics==TRUE)  
  ################################
  
  if (i_iter%%1000==0) {
    
    print( paste0( "rep_no = ",rep_no, ", iteration: ",i_iter," / ",n_iter) )
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
posterior_sig2 = draw_sig2[SEQ]
save(posterior_mu,posterior_Sigma_theta,posterior_sig2,file=paste0(RDataFolder,"/rep_",rep_no,".RData"))