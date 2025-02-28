rm(list=ls())

# create output folder if it does not exist
Folder = paste0("IPD_AD_meta")
if (!file.exists(Folder)){ dir.create(Folder, showWarnings = TRUE, recursive = FALSE, mode = "0777") }

# libraries 
library(mvtnorm) ; library(MCMCpack) ; library(truncnorm) ; library(gplots) ; library(MASS)

#####################################
# Import the individual participant data (IPD) 
#####################################

IPD_import = read.csv(file="S_RealData/11_Import.csv") 
Unique_study_name = sort(unique(IPD_import[,"study_name"]))
J_IPD = length(Unique_study_name) 

IPD_listobj = NULL 
for (l in 1:J_IPD){
  SEQ = which(IPD_import[,"study_name"]==Unique_study_name[[l]])
  groupData = IPD_import[SEQ,]
  y = groupData[,"w_gain"]
  x11 = groupData[,"bmi_cat1"] ; x12 = groupData[,"bmi_cat2"] ; x13 = groupData[,"bmi_cat3"]
  x2 = groupData[,"b_wt"]
  x31 = groupData[,"bmi_cat1.trt"] ; x32 = groupData[,"bmi_cat2.trt"] ; x33 = groupData[,"bmi_cat3.trt"]
  X = cbind(x11,x12,x13,x2,x31,x32,x33)
  n_l = length(y)
  IPD_listobj[[l]] = list(y=y,X=X,n_l=n_l)
} # for (l)

#####################################
# Import the aggregate data (AD)
#####################################

AD_type1 = read.csv( file="S_RealData/AD_type1.csv", header=T ) ; J_type1 = nrow(AD_type1)
AD_type2 = suppressWarnings( read.csv( file="S_RealData/AD_type2.csv", header=T ) ) ; J_type2 = nrow(AD_type2)
AD_type3_betahat = read.csv( file="S_RealData/AD_type3_a_mean.csv", header=T ) ; J_type3 = nrow(AD_type3_betahat)
AD_type3_SE_betahat = read.csv( file="S_RealData/AD_type3_b_SE.csv", header=T) 
AD_type3_DRM = read.csv( file="S_RealData/AD_type3_c_drm.csv", header=T)
AD_type3_betahat = AD_type3_betahat[,1:6]
AD_type3_SE_betahat = AD_type3_SE_betahat[,1:6]

L = J_IPD + J_type1 + J_type2 + J_type3

#####################################
# Make reference X 
#####################################

REF_X_type3 = NULL

i_REF_X_type3 = 1 # Spain
ref_studyname = c("Barakat 2008","Barakat 2012a","Perales 2014","Ruiz 2013")
WHICH = which(IPD_import[,"study_name"] %in% ref_studyname)
REF_X_type3[[i_REF_X_type3]] = IPD_import[WHICH,]

i_REF_X_type3 = i_REF_X_type3 + 1 # Australia
ref_studyname = c("Dodd 2014","Harrison 2013","Ong 2009")
WHICH = which(IPD_import[,"study_name"] %in% ref_studyname)
REF_X_type3[[i_REF_X_type3]] = IPD_import[WHICH,]

i_REF_X_type3 = i_REF_X_type3 + 1 # UK / Ireland
ref_studyname = c("Walsh 2012") # for Daley 2015 from UK
WHICH = which(IPD_import[,"study_name"] %in% ref_studyname)
REF_X_type3[[i_REF_X_type3]] = IPD_import[WHICH,]

i_REF_X_type3 = i_REF_X_type3 + 1 # USA
ref_studyname = c("Phelan 2011")
WHICH = which(IPD_import[,"study_name"] %in% ref_studyname)
REF_X_type3[[i_REF_X_type3]] = IPD_import[WHICH,]

i_REF_X_type3 = i_REF_X_type3 + 1 # Canada
ref_studyname = c("Hui 2011")
WHICH = which(IPD_import[,"study_name"] %in% ref_studyname)
REF_X_type3[[i_REF_X_type3]] = IPD_import[WHICH,]

######

tilde_D_x_type3 = NULL 

for (k in 1:J_type3){
  ref_no = AD_type3_DRM[k,"reference_no"]
  tilde_D_x_type3[[k]] = REF_X_type3[[ref_no]]
} # for (k) 

# Remove of IPD whose subgroup does not exist in AD 

WHICH = which(AD_type3_DRM[,"study_name"]=="Barakat 2015") # use all 
WHICH = which(AD_type3_DRM[,"study_name"]=="Brownfoot 2016") # use all 

WHICH = which(AD_type3_DRM[,"study_name"]=="Daly 2017")
SEQ = which(tilde_D_x_type3[[WHICH]][,"bmi_cat3"]==1)
tilde_D_x_type3[[WHICH]] = tilde_D_x_type3[[WHICH]][SEQ,]

WHICH = which(AD_type3_DRM[,"study_name"]=="Dekker 2015")
SEQ = which(tilde_D_x_type3[[WHICH]][,"bmi_cat3"]==1)
tilde_D_x_type3[[WHICH]] = tilde_D_x_type3[[WHICH]][SEQ,]

WHICH = which(AD_type3_DRM[,"study_name"]=="Hui 2014")
SEQ1 = which(tilde_D_x_type3[[WHICH]][,"bmi_cat1"]==1)
SEQ2 = which(tilde_D_x_type3[[WHICH]][,"bmi_cat2"]==1)
tilde_D_x_type3[[WHICH]] = tilde_D_x_type3[[WHICH]][c(SEQ1,SEQ2),]

WHICH = which(AD_type3_DRM[,"study_name"]=="Kong 2014")
SEQ1 = which(tilde_D_x_type3[[WHICH]][,"bmi_cat2"]==1)
SEQ2 = which(tilde_D_x_type3[[WHICH]][,"bmi_cat3"]==1)
tilde_D_x_type3[[WHICH]] = tilde_D_x_type3[[WHICH]][c(SEQ1,SEQ2),]

WHICH = which(AD_type3_DRM[,"study_name"]=="Polley 2002")
SEQ1 = which(tilde_D_x_type3[[WHICH]][,"bmi_cat1"]==1)
SEQ2 = which(tilde_D_x_type3[[WHICH]][,"bmi_cat2"]==1)
tilde_D_x_type3[[WHICH]] = tilde_D_x_type3[[WHICH]][c(SEQ1,SEQ2),]

######
Design_X_type3 = nu_X_type3 = threshold_type3 = NULL 

for (k in 1:J_type3){
  
  tempX = tilde_D_x_type3[[k]]
  Design_X_type3[[k]] = tempX[,c("bmi_cat1","bmi_cat2","bmi_cat3","b_wt","bmi_cat1.trt","bmi_cat2.trt","bmi_cat3.trt")]
  
  s1 = (tempX[,"new_trt"]==0) * (tempX[,"bmi_cat1"]==1) * 1
  s2 = (tempX[,"new_trt"]==0) * (tempX[,"bmi_cat2"]==1) * 1
  s3 = (tempX[,"new_trt"]==0) * (tempX[,"bmi_cat3"]==1) * 1
  s4 = (tempX[,"new_trt"]==1) * (tempX[,"bmi_cat1"]==1) * 1
  s5 = (tempX[,"new_trt"]==1) * (tempX[,"bmi_cat2"]==1) * 1
  s6 = (tempX[,"new_trt"]==1) * (tempX[,"bmi_cat3"]==1) * 1
  nu_X_type3[[k]] = cbind(s1,s2,s3,s4,s5,s6)
  
  WHICH = which( apply(nu_X_type3[[k]],1,sum)!=1 )
  if (length(WHICH)>0){
    Design_X_type3[[k]] = Design_X_type3[[k]][-WHICH,]
    nu_X_type3[[k]] = nu_X_type3[[k]][-WHICH,]
    tilde_D_x_type3[[k]] = tilde_D_x_type3[[k]][-WHICH,]
  } # if 
  
  WHICH = which( apply(is.na(Design_X_type3[[k]]),1,sum)>0 )
  if (length(WHICH)>0){
    Design_X_type3[[k]] = Design_X_type3[[k]][-WHICH,]
    nu_X_type3[[k]] = nu_X_type3[[k]][-WHICH,]
    tilde_D_x_type3[[k]] = tilde_D_x_type3[[k]][-WHICH,]
    
  } # if 
  
  WHICH = which( apply(is.na(nu_X_type3[[k]]),1,sum)>0 )
  if (length(WHICH)>0){
    Design_X_type3[[k]] = Design_X_type3[[k]][-WHICH,]
    nu_X_type3[[k]] = nu_X_type3[[k]][-WHICH,]
    tilde_D_x_type3[[k]] = tilde_D_x_type3[[k]][-WHICH,]
  } # if 
  
  threshold_type3[[k]] = rep(0,dim(nu_X_type3[[k]])[[1]])
  WHICH = which(Design_X_type3[[k]][,"bmi_cat1"]==1) 	
  threshold_type3[[k]][WHICH] = 15.9 
  WHICH = which(Design_X_type3[[k]][,"bmi_cat2"]==1) 	
  threshold_type3[[k]][WHICH] = 11.3 
  WHICH = which(Design_X_type3[[k]][,"bmi_cat3"]==1) 	
  threshold_type3[[k]][WHICH] = 9.1 
  
} # k 

#####################################
# Density ratio model
#####################################

hat_tau_type3_vec = hat_Gamma_tau_type3_vec = rep(0,J_type3)

n_AD_mat = AD_type3_DRM[,c("n_AD","n_AD1","n_AD2","n_AD3")]

DRM_X_type3 = NULL

for (k in 1:J_type3){
  if (AD_type3_DRM[k,"USE"]=="bmi"){
    hat_tau_type3_vec[k] = AD_type3_DRM[k,"b_bmi_mean"]
    hat_Gamma_tau_type3_vec[k] = AD_type3_DRM[k,"b_bmi_sd"]^2 / ( AD_type3_DRM[k,"n_AD"]-1 )
    temp_X = tilde_D_x_type3[[k]]
    temp_X = temp_X[,c("bmi_cat1","bmi_cat2","bmi_cat3","b_bmi")]
    colnames(temp_X)[[4]] = "DRM_X_type3"
    int = rep(1,nrow(temp_X)) ; temp_X = cbind(temp_X,int)
    
    a_n = rep(0,nrow(temp_X))
    temp_X_IndOnly = temp_X[,c("bmi_cat1","bmi_cat2","bmi_cat3")]
    proportion_IPD = apply( temp_X_IndOnly, 2, mean )
    proportion_AD = n_AD_mat[k,c("n_AD1","n_AD2","n_AD3")] / n_AD_mat[k,"n_AD"]
    for (i_temp in 1:3){
      SEQ = which(temp_X_IndOnly[i_temp]==1)
      if ( length(SEQ)>0 ){
        a_n[SEQ] = as.numeric( proportion_AD[i_temp] / proportion_IPD[i_temp] )
      } #
    } # for (i_temp)
    temp_X = cbind(temp_X,a_n)
    
    DRM_X_type3[[k]] = temp_X
  } # if (AD_type3_DRM[k,"USE"]=="bmi")
  if (AD_type3_DRM[k,"USE"]=="wt"){
    hat_tau_type3_vec[k] = AD_type3_DRM[k,"b_wt_mean"]
    hat_Gamma_tau_type3_vec[k] = AD_type3_DRM[k,"b_wt_sd"]^2 / ( AD_type3_DRM[k,"n_AD"]-1 )
    temp_X = tilde_D_x_type3[[k]]
    temp_X = tilde_D_x_type3[[k]]
    temp_X = temp_X[,c("bmi_cat1","bmi_cat2","bmi_cat3","b_wt")]
    colnames(temp_X)[[4]] = "DRM_X_type3"
    int = rep(1,nrow(temp_X)) ; temp_X = cbind(temp_X,int)
    
    a_n = rep(0,nrow(temp_X))
    temp_X_IndOnly = temp_X[,c("bmi_cat1","bmi_cat2","bmi_cat3")]
    proportion_IPD = apply( temp_X_IndOnly, 2, mean )
    proportion_AD = n_AD_mat[k,c("n_AD1","n_AD2","n_AD3")] / n_AD_mat[k,"n_AD"]
    for (i_temp in 1:3){
      SEQ = which(temp_X_IndOnly[i_temp]==1)
      if ( length(SEQ)>0 ){
        a_n[SEQ] = as.numeric( proportion_AD[i_temp] / proportion_IPD[i_temp] )
      } #
    } # for (i_temp)
    temp_X = cbind(temp_X,a_n)
    
    DRM_X_type3[[k]] = temp_X
  } # if (AD_type3_DRM[k,"USE"]=="wt")

} # for (k)

##################################
# Functions 
##################################

MB.logit.est = function(param, theta_vec, Design_mat, bootstrap_w, nu_x_mat, threshold_type3, sig2_IPD_type3, UseDRM, alpha, REF_X_type3_input)
{
	Design_mat = as.matrix(Design_mat); theta_vec=as.vector(theta_vec); param=as.vector(param); nu_x_mat=as.matrix(nu_x_mat);
	std_dev = (threshold_type3-Design_mat%*%theta_vec)/sqrt(sig2_IPD_type3)
	exp_y_above_tau = 1-pnorm(std_dev, mean=0, sd=1)
	if (UseDRM == TRUE){
	  REF_X_type3_mat = REF_X_type3_input[,c("int","DRM_X_type3")]
	  w.DRM = exp(as.matrix(REF_X_type3_mat)%*%alpha) # alpha_l * psi(x)
	  bootstrap_w = w.DRM * bootstrap_w
	}
	exp_minus_fitted = bootstrap_w * as.vector( exp_y_above_tau - nu_x_mat%*%param )
	g.i = cbind(nu_x_mat[,1]*exp_minus_fitted,nu_x_mat[,2]*exp_minus_fitted,nu_x_mat[,3]*exp_minus_fitted, nu_x_mat[,4]*exp_minus_fitted, nu_x_mat[,5]*exp_minus_fitted, nu_x_mat[,6]*exp_minus_fitted) 
	g.sum = apply(g.i,2,sum)
	return( sum(g.sum^2) )
}

###############################################
# Setting Hyperparameters and starting values # 
###############################################

p_theta = 7 ; p_alpha = 2

# Prior distributions 
mu_sig2 = 20 ; tau_sig2 = 20 
mu_0_vec = rep(0,p_theta) ; Sigma_0_mat = diag(100,p_theta) ; inv_Sigma_0_mat = solve(Sigma_0_mat)
nu_0 = 0.1 ; Psi_mat = diag(0.1,p_theta)

# Starting values 
theta_l_mat = array(0.1,c(L,p_theta)) 
sig2_IPD_type3 = rep(10.0,(J_IPD+J_type3)) ; mu_vec = rep(0,p_theta) ; Sigma_mat = diag(1,p_theta) ; inv_Sigma_mat = solve(Sigma_mat)
tau_vec = hat_tau_type3_vec                       # for the density ratio model
alpha_mat_type3 = array(0.01,c(J_type3,p_alpha))  # for the density ratio model 

p_beta_type3 = 6
beta_l_type3_mat = array(0.1,c(J_type3,p_beta_type3)) 
for (k in 1:J_type3){
    l = J_IPD + J_type1 + J_type2 + k
		n_l = dim(Design_X_type3[[k]])[[1]]
		ww = rnorm(n_l,mean=1,sd=1) # multiplier bootstrap 
		temp.optim = optim(beta_l_type3_mat[k,], fn=MB.logit.est, theta_vec=theta_l_mat[l,], Design_mat=Design_X_type3[[k]], bootstrap_w=ww, nu_x_mat=nu_X_type3[[k]], threshold_type3=threshold_type3[[k]], sig2_IPD_type3=sig2_IPD_type3[J_IPD+k], UseDRM=AD_type3_DRM[k,"DRM"], alpha=alpha_mat_type3[k,], REF_X_type3_input=DRM_X_type3[[k]], control=list(pgtol=1e-7, maxit=1000), lower=0, upper=1, method="L-BFGS-B")
		beta_l_type3_mat[k,] = temp.optim$par
} # for (k)

stepsize_type1 = stepsize_type2 = 0.2 
stepsize_type3 = 0.1 ; stepsize_sigma2 = 2.0 ; stepsize_alpha = 0.001 ; stepsize_tau = 0.2

###############################################
# Running MCMC
###############################################

# MCMC setting 
n_burnin = n_main = 10000
n_iter = n_burnin+n_main 

# Prepare repository for posterior draws 
Draw_mu_vec = array(0,c(n_iter,p_theta))
Draw_Sigma_mat = array(0,c(n_iter,p_theta,p_theta))
set.seed(117)

######
program_start_time = Sys.time()
format(Sys.time(), "%b %d %Y, %a, %H:%M:%S ")
Prevtime = proc.time()[3]

###### Start of MCMC 

for (i_iter in 1:n_iter){
	
  ##################################
  # Update theta, beta, and alpha
  ##################################
  
  # update theta for IPD studies (Supplementary Material, Section 6, Step 1)
  
	for (l in 1:J_IPD){
		y = IPD_listobj[[l]]$y ; X = IPD_listobj[[l]]$X 
		inv_Var = t(X)%*%X / sig2_IPD_type3[l] + inv_Sigma_mat 
		Var = solve(inv_Var)
		Mean_temp_part = t(X)%*%y / sig2_IPD_type3[l] + inv_Sigma_mat %*% t(t(mu_vec))
		Mean = Var %*% Mean_temp_part
		theta_l_mat[l,] = rmvnorm(n=1, mean=Mean, sigma=Var)
	} # for (l)
	
  # update theta and beta for Type 1 AD studies (Supplementary Material, Section 6, Step 2)
  
	for (k in 1:J_type1){
		l = J_IPD + k
		theta_cur = theta_l_mat[l,] ; beta_cur = theta_cur[7]
		theta_q = rnorm(n=p_theta,theta_cur,stepsize_type1) ; beta_q = theta_q[7]
		logNUM = dmvnorm(theta_q, mu_vec, Sigma_mat, log=T) + dnorm(AD_type1[k,"beta_hat"], beta_q, AD_type1[k,"sqrt_hat_V"],log=T)
		logDEN = dmvnorm(theta_cur, mu_vec, Sigma_mat, log=T) + dnorm(AD_type1[k,"beta_hat"], beta_cur, AD_type1[k,"sqrt_hat_V"],log=T)	
		if (runif(1)<=exp(logNUM-logDEN)){
			theta_l_mat[l,] = theta_q 
		}	# 
	} # for (l)
	
  # update theta and beta for Type 2 AD studies (Supplementary Material, Section 6, Step 2)
  
	for (k in 1:J_type2){
		l = J_IPD + J_type1 + k
		theta_cur = theta_l_mat[l,]
		p_vec = AD_type2[k,1:6]
		beta_cur = (p_vec["p_11"]-p_vec["p_01"])*theta_cur[1] + p_vec["p_11"]*theta_cur[5]
		beta_cur = beta_cur + (p_vec["p_12"]-p_vec["p_02"])*theta_cur[2] + p_vec["p_12"]*theta_cur[6]
		beta_cur = beta_cur + (p_vec["p_13"]-p_vec["p_03"])*theta_cur[3] + p_vec["p_13"]*theta_cur[7]
		beta_cur = as.numeric(beta_cur)
		theta_q = rnorm(n=p_theta,theta_cur,stepsize_type2)
		beta_q = (p_vec["p_11"]-p_vec["p_01"])*theta_q[1] + p_vec["p_11"]*theta_q[5]
		beta_q = beta_q + (p_vec["p_12"]-p_vec["p_02"])*theta_q[2] + p_vec["p_12"]*theta_q[6]
		beta_q = beta_q + (p_vec["p_13"]-p_vec["p_03"])*theta_q[3] + p_vec["p_13"]*theta_q[7]
		beta_q = as.numeric(beta_q)
		logNUM = dmvnorm(theta_q, mu_vec, Sigma_mat, log=T) + dnorm(AD_type2[k,"beta_hat"], beta_q, AD_type2[k,"sqrt_hat_V"], log=T)
		logDEN = dmvnorm(theta_cur, mu_vec, Sigma_mat, log=T) + dnorm(AD_type2[k,"beta_hat"], beta_cur, AD_type2[k,"sqrt_hat_V"], log=T)
		if (runif(1)<=exp(logNUM-logDEN)){
			theta_l_mat[l,] = theta_q
		}	#
	} # for (l)
	
  # update theta, beta, and alpha for Type 3 AD studies (Supplementary Material, Section 6, Step 3)
  
	for (k in 1:J_type3){
		
	  n_l = dim(Design_X_type3[[k]])[[1]]
	  
	  # theta 
		l = J_IPD + J_type1 + J_type2 + k
		
		theta_type3_vec_q = rnorm(n=p_theta, mean=theta_l_mat[l,], sd=stepsize_type3)	
		logNUMDEN1 = dmvnorm(theta_type3_vec_q, mu_vec, Sigma_mat, log=T) - dmvnorm(theta_l_mat[l,], mu_vec, Sigma_mat, log=T)
				
		ww = rnorm(n_l,mean=1,sd=1) # multiplier bootstrap 
		temp.optim = optim(beta_l_type3_mat[k,], fn=MB.logit.est, theta_vec=theta_type3_vec_q, Design_mat=Design_X_type3[[k]], bootstrap_w=ww, nu_x_mat=nu_X_type3[[k]], threshold_type3=threshold_type3[[k]], sig2_IPD_type3=sig2_IPD_type3[J_IPD+k],  UseDRM=AD_type3_DRM[k,"DRM"], alpha=alpha_mat_type3[k,], REF_X_type3_input=DRM_X_type3[[k]], control=list(pgtol=1e-7, maxit=1000), lower=0, upper=1, method="L-BFGS-B")
		
		if ( temp.optim$convergence==0 ){
		  
		  beta_type3_vec_q = temp.optim$par
		  
		  hat_beta_temp = AD_type3_betahat[k,]
		  SEQ = which(is.na(hat_beta_temp)==F)
		  hat_beta_temp = as.numeric(hat_beta_temp[SEQ])
		  se_hat_beta_temp = as.numeric(AD_type3_SE_betahat[k,SEQ])
		  logNUMDEN2 = sum(dnorm(hat_beta_temp, beta_type3_vec_q[SEQ], se_hat_beta_temp, log=T) - dnorm(hat_beta_temp, beta_l_type3_mat[k,SEQ], se_hat_beta_temp, log=T))
		  
		  logAcc = logNUMDEN1 + logNUMDEN2
		  if ( runif(n=1) < exp(logAcc) ) {
		    theta_l_mat[l,] = theta_type3_vec_q
		    beta_l_type3_mat[k,] = beta_type3_vec_q
		  } # if 	
		  
		} # if ( temp.optim$convergence==0 )
		
		# alpha 
		if (AD_type3_DRM[k,"DRM"]==TRUE){
		  
		  alpha_vec_q = rnorm(n=p_alpha, mean=alpha_mat_type3[k,], sd=stepsize_alpha)	# density ratio
		  
		  a_n = DRM_X_type3[[k]][,"a_n"] #### design effect
		  
		  Exp_alpha_psi_q = exp( as.matrix(DRM_X_type3[[k]][,c("int","DRM_X_type3")])%*%alpha_vec_q )
		  g_l_mat_q = cbind( a_n * (Exp_alpha_psi_q - 1), a_n * (DRM_X_type3[[k]][,"DRM_X_type3"] * Exp_alpha_psi_q - tau_vec[k]) )  #### design effect 
		  bar_g_l_q = apply(g_l_mat_q,2,mean) ; Sigma_g_l_q = var(g_l_mat_q) / nrow(g_l_mat_q)
		  
		  Exp_alpha_psi = exp( as.matrix(DRM_X_type3[[k]][,c("int","DRM_X_type3")])%*%alpha_mat_type3[k,] )
		  g_l_mat = cbind( a_n * (Exp_alpha_psi - 1), a_n * (DRM_X_type3[[k]][,"DRM_X_type3"] * Exp_alpha_psi - tau_vec[k]) )  #### design effect  
		  bar_g_l = apply(g_l_mat,2,mean) ; Sigma_g_l = var(g_l_mat)  / nrow(g_l_mat)
		  
		  logNUMDEN1 = dmvnorm(bar_g_l_q, mean=rep(0,2), sigma=Sigma_g_l_q, log=TRUE) - dmvnorm(bar_g_l, mean=rep(0,2), sigma=Sigma_g_l, log=TRUE)
		  
		  ww = rnorm(n_l,mean=1,sd=1) # multiplier bootstrap 
		  temp.optim = optim(beta_l_type3_mat[k,], fn=MB.logit.est, theta_vec=theta_l_mat[l,], Design_mat=Design_X_type3[[k]], bootstrap_w=ww, nu_x_mat=nu_X_type3[[k]], threshold_type3=threshold_type3[[k]], sig2_IPD_type3=sig2_IPD_type3[J_IPD+k],  UseDRM=AD_type3_DRM[k,"DRM"], alpha=alpha_vec_q, REF_X_type3_input=DRM_X_type3[[k]], control=list(pgtol=1e-7, maxit=1000), lower=0, upper=1, method="L-BFGS-B")
		  
		  if ( temp.optim$convergence==0 ){
		    
		    beta_type3_vec_q = temp.optim$par
		    
		    hat_beta_temp = AD_type3_betahat[k,]
		    SEQ = which(is.na(hat_beta_temp)==F)
		    hat_beta_temp = as.numeric(hat_beta_temp[SEQ])
		    se_hat_beta_temp = as.numeric(AD_type3_SE_betahat[k,SEQ])
		    logNUMDEN2 = sum(dnorm(hat_beta_temp, beta_type3_vec_q[SEQ], se_hat_beta_temp, log=T) - dnorm(hat_beta_temp, beta_l_type3_mat[k,SEQ], se_hat_beta_temp, log=T))
		    
		    logAcc = logNUMDEN1 + logNUMDEN2
		    if ( runif(n=1) < exp(logAcc) ) {
		      alpha_mat_type3[k,] = alpha_vec_q # density ratio 
		      beta_l_type3_mat[k,] = beta_type3_vec_q
		    } # if 
		    
		  } # if ( temp.optim$convergence==0 )
		  
		} # if (AD_type3_DRM[k,"DRM"]==TRUE)
	
	} # for (k)
	
  ##################################
  # Update sig2 
  ##################################
  
  # update sig2 for IPD studies (Supplementary Material, Section 6, Step 4)

	for (l in 1:J_IPD){
		
		sig2_q = rnorm(n=1, mean=sig2_IPD_type3[l], sd=stepsize_sigma2)
		
		if (sig2_q > 0){
		  
		  y = IPD_listobj[[l]]$y ; X = IPD_listobj[[l]]$X ; n_l = IPD_listobj[[l]]$n_l
		  theta_l_vec = t(t(theta_l_mat[l,]))
		  logNUMDEN1 = sum(dnorm(y,X%*%theta_l_vec,sqrt(sig2_q),log=T)) - sum(dnorm(y,X%*%theta_l_vec,sqrt(sig2_IPD_type3[l]),log=T))
		  logSig2 = log( dtruncnorm(sig2_q, a=0, b=Inf, mean=mu_sig2, sd=tau_sig2) ) - log( dtruncnorm(sig2_IPD_type3[l], a=0, b=Inf, mean=mu_sig2, sd=tau_sig2) )
		  
		  logAcc = logNUMDEN1 + logSig2
		  if ( runif(n=1) < exp(logAcc) ) {
		    sig2_IPD_type3[l] = sig2_q 
		  } # if 		
		  
		} # if (sig2_q > 0)
		
	} # for (l in 1:J_IPD)
 	
  # update sig2 for Type 3 AD studies (Supplementary Material, Section 6, Step 5)
  
	for (k in 1:J_type3){
		
	  l = J_IPD + J_type1 + J_type2 + k
	  
	  n_l = dim(Design_X_type3[[k]])[[1]]
	  
		sig2_q = rnorm(n=1, mean=sig2_IPD_type3[J_IPD+k], sd=stepsize_sigma2) # USE	n_l = dim(Design_X_type3[[k]])[[1]]
		
		if ( sig2_q>0 ){
		  
		  ww = rnorm(n_l,mean=1,sd=1) # multiplier bootstrap 
		  temp.optim = optim(beta_l_type3_mat[k,], fn=MB.logit.est, theta_vec=theta_l_mat[l,], Design_mat=Design_X_type3[[k]], bootstrap_w=ww, nu_x_mat=nu_X_type3[[k]], threshold_type3=threshold_type3[[k]], sig2_IPD_type3=sig2_q, UseDRM=AD_type3_DRM[k,"DRM"], alpha=alpha_mat_type3[k,], REF_X_type3_input=DRM_X_type3[[k]], control=list(pgtol=1e-7, maxit=1000), lower=0, upper=1, method="L-BFGS-B")
		  
		  if ( temp.optim$convergence==0 ){
		    
		    beta_type3_vec_q = temp.optim$par
		    
		    hat_beta_temp = AD_type3_betahat[k,]
		    SEQ = which(is.na(hat_beta_temp)==F)
		    hat_beta_temp = as.numeric(hat_beta_temp[SEQ])
		    se_hat_beta_temp = as.numeric(AD_type3_SE_betahat[k,SEQ])
		    
		    logNUMDEN2 = sum(dnorm(hat_beta_temp, beta_type3_vec_q[SEQ], se_hat_beta_temp, log=T) - dnorm(hat_beta_temp, beta_l_type3_mat[k,SEQ], se_hat_beta_temp, log=T))
		    logSig2 = log( dtruncnorm(sig2_q, a=0, b=Inf, mean=mu_sig2, sd=tau_sig2) ) - log( dtruncnorm(sig2_IPD_type3[J_IPD+k], a=0, b=Inf, mean=mu_sig2, sd=tau_sig2) )
		    
		    logAcc = logNUMDEN2 + logSig2 
		    
		    if ( runif(n=1) < exp(logAcc) ) {
		      beta_l_type3_mat[k,] = beta_type3_vec_q
		      sig2_IPD_type3[J_IPD+k] = sig2_q 
		    } # if 	
		    
		  } # if ( temp.optim$convergence==0 )
		  
		} # 
		
	} # for (k in 1:J_type3)
	
  ##################################
  # Update mu_vec (Supplementary Material, Section 6, Step 7)
  ##################################
  
	inv_Var = inv_Sigma_0_mat + L * inv_Sigma_mat
	Var = solve(inv_Var)
	theta_sum_vec = apply(theta_l_mat,2,sum)
	Mean_temp_part = inv_Sigma_mat %*% theta_sum_vec
	Mean = Var %*% Mean_temp_part
	mu_vec = t(rmvnorm(n=1, mean=Mean, sigma=Var))
	
	##################################
	# Update Sigma_theta_mat (Supplementary Material, Section 6, Step 8)
	##################################
	
	df1 = nu_0 + L
	SqSum = array(0,c(p_theta,p_theta))
	for (l in 1:L){
		theta_l_vec = t(t(theta_l_mat[l,]))
		SqSum = SqSum + (theta_l_vec-mu_vec)%*%t(theta_l_vec-mu_vec)
	} # (l)
	df2 = Psi_mat + SqSum 
	Sigma_mat = riwish(v=df1, S=df2) ; inv_Sigma_mat = solve(Sigma_mat)
	
	##################################
	# Update tau_l (Supplementary Material, Section 6, Step 6)
	##################################
	
	if ( i_iter > 1000 ){
	  
	  for (k in 1:J_type3){
	    
	    if (AD_type3_DRM[k,"DRM"]==TRUE){
	      
	      Exp_alpha_psi = exp( as.matrix(DRM_X_type3[[k]][,c("int","DRM_X_type3")])%*%alpha_mat_type3[k,] )                 # common for q and current
	      a_n = DRM_X_type3[[k]][,"a_n"]  #### design effect
	      
	      tau_q = rnorm(n=1, mean=tau_vec[k], sd=stepsize_tau)	
	      
	      g_l_mat_q = cbind( a_n * (Exp_alpha_psi - 1), a_n * (DRM_X_type3[[k]][,"DRM_X_type3"] * Exp_alpha_psi - tau_q) ) #### design effect 
	      bar_g_l_q = apply(g_l_mat_q,2,mean) ; Sigma_g_l_q = var(g_l_mat_q) / nrow(g_l_mat_q)
	      logNum = dnorm(hat_tau_type3_vec[k], mean=tau_q, sd=sqrt(hat_Gamma_tau_type3_vec[k]), log=TRUE)
	      logNum = logNum + dmvnorm(bar_g_l_q, mean=rep(0,2), sigma=Sigma_g_l_q, log=TRUE) 
	      
	      g_l_mat = cbind( a_n * (Exp_alpha_psi - 1), a_n * (DRM_X_type3[[k]][,"DRM_X_type3"] * Exp_alpha_psi - tau_vec[k]) ) #### design effect 
	      bar_g_l = apply(g_l_mat,2,mean) ; Sigma_g_l = var(g_l_mat) / nrow(g_l_mat)
	      logDen = dnorm(hat_tau_type3_vec[k], mean=tau_vec[k], sd=sqrt(hat_Gamma_tau_type3_vec[k]), log=TRUE)
	      logDen = logDen + dmvnorm(bar_g_l, mean=rep(0,2), sigma=Sigma_g_l, log=TRUE) 
	      
	      logAcc = logNum - logDen 
	      
	      if ( runif(n=1) < exp(logAcc) ) {
	        tau_vec[k] = tau_q
	      }
	      
	    } # if (AD_type3_DRM[k,"DRM"]==TRUE)
	    
	  } # for (k)
	  
	} # 
	
	##################################
	# Store posterior draws (every iteration)
	##################################
	
	Draw_mu_vec[i_iter,] = mu_vec
	Draw_Sigma_mat[i_iter,,] = Sigma_mat
	
	################################
	# Print iteration numbers (every 100 iterations) 
	#  and plot mu's 
	################################
	
	if (i_iter%%100==0){
		
	  print( paste0( "rep_no = ",rep_no, ", iteration: ",i_iter," / ",n_iter) )
	  Currenttime = proc.time()[3]
	  LastBatch = Currenttime-Prevtime ; Time_to_Go = (n_iter-i_iter)*(LastBatch/100)
	  Prevtime = Currenttime
	  print( paste("The last 1000 iter=",round(LastBatch/60,1),"min, Est. Time to go=",round(Time_to_Go/60,1),"min" ))
	  
		png(file=paste0(Folder,"/W_a_mu_progress.png"),width=2000,height=2000,pointsize=30)
		par(mfrow=c(3,3),mai=c(0.8,0.8,0.6,0.6),family="AppleGothic",mgp = c(1.5, 0.5, 0)) # b l t r
		for (k in 1:7){
		  plot(Draw_mu_vec[1:i_iter,k], type="l", main=paste0("mu_",k))
		} # for 
		dev.off()

	} # if (i_iter%%100==0)

	if (i_iter%%1000==0){
		 save.image(paste0(Folder,"_temp.RData"))
	} # if (i_iter%%1000)

} # for (i_iter)

###### End of MCMC 

# Save the posterior draws after burn-in
SEQ = (burnin+1):(burnin+mainrun)
posterior_mu = Draw_mu_vec[SEQ,]
posterior_Sigma_mat = Draw_Sigma_mat[SEQ,]
save(posterior_mu,posterior_Sigma_mat,file=paste0(RDataFolder,"/rep_",rep_no,".RData"))

save.image(paste0(Folder,".RData"))