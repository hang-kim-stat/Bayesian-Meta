###################################################
# Generating simulation data for Section 4.2
###################################################

rm(list=ls())

library(MASS)

######################################## 
n_rep = 300     # number of replicates 
L = 40          # total number of studies 
n_sample = 200  # sample size of each study 

# Model setting
true_mu = c(log(0.4/(1-0.4)), 0.5, -0.3, -0.6) 
true_Sigma_theta = rep(0.1,4)

# biased access
alpha_biased_access = c(-0.75,0.95) 
fr_biased = rep(0, n_rep)
RR = NULL

p_theta = length(true_mu) ; p_beta = p_theta 

# function 
logit_inv = function(alpha) return( exp(alpha)/(exp(alpha)+1) ) ;

########################################
SimulData = NULL

# for each simulation 
for (i_rep in 1:n_rep){
  
	print(" ============================================== ")	
	print("                                                ")
  
	set.seed(i_rep+33) ; print(i_rep)
  
	theta_l_mat = array(0,c(L,p_theta)) ; dimnames(theta_l_mat)[[2]] = c("theta0","theta1","theta2","theta3")
	theta_hat_l_mat = array(0,c(L,p_theta)) ; dimnames(theta_hat_l_mat)[[2]] = c("theta0","theta1","theta2","theta3")
	X_cube = array(0,c(L,n_sample,p_theta)) ; dimnames(X_cube)[[3]] = c("Int","X1","X2","X1X2") 
	A_cube = array(0,c(L,p_theta,p_theta))
	Y_mat = array(0,c(L,n_sample)) 
	beta_mat = array(0,c(L,p_beta)) ; dimnames(beta_mat)[[2]] = c("beta0","beta1","beta2","beta3")
	V_beta_cube = array(0,c(L,p_beta,p_beta)) ; dimnames(V_beta_cube)[[2]] = dimnames(V_beta_cube)[[3]] = c("beta0","beta1","beta2","beta3")     
  
	for (i_study in 1:L) {
		
	  # theta_l ~ N( mu, Sigma )
		theta_l_mat[i_study,] = rnorm(n=p_theta, mean=true_mu, sd=sqrt(true_Sigma_theta)) 
    
		X_cube[i_study,,1] = rep(1, n_sample)      # dim(X_cube) is (L, n_sample, p_theta)
		upsilon_l = rnorm(1, mean=0, sd=sqrt(0.1)) # upsilon_l ~ N(0,0.1) added for heterogeneous X
		X_cube[i_study,,2] = rnorm(n_sample, mean=upsilon_l, sd=sqrt(0.5))
		X_cube[i_study,,3] = (runif(n_sample) <= 0.5) * 1.0  # treatment allocation (binary)
		X_cube[i_study,,4] = X_cube[i_study,,2] * X_cube[i_study,,3]
		
		logit_p_Y = X_cube[i_study,,] %*% theta_l_mat[i_study,] # logit( p(Y_li=1|X_li1,X_li2) )
		Y_mat[i_study,] = ( runif(n=n_sample) <= logit_inv(logit_p_Y) ) * 1.0  
		
		tempData = data.frame(Y=Y_mat[i_study,], X1=X_cube[i_study,,2], X2=X_cube[i_study,,3], X1X2=X_cube[i_study,,4])		
 	  theta_glm = glm(Y ~ X1 * X2, data = tempData, family = "binomial")
		theta_hat_l_mat[i_study,] = summary(theta_glm)$coefficients[,1]
    
		# beta for AD 
		beta_glm = glm(Y ~ X1 + X2, data = tempData, family = "binomial")		
		beta_mat[i_study,2:3] = summary(beta_glm)$coefficients[2:3,1]
		V_beta_cube[i_study,2:3,2:3] = vcov( summary(beta_glm) )[2:3,2:3]
			
	} # for (i_study)
	
	temp_standardized_beta2 = beta_mat[,3] / sqrt(V_beta_cube[,3,3])
	logit_p_delta_biased = alpha_biased_access[1] + alpha_biased_access[2]*abs(temp_standardized_beta2)
	delta_biased_access = ( runif(n=L) <= logit_inv(logit_p_delta_biased) ) * 1.0
	
	fr_biased[i_rep] = sum(delta_biased_access)/length(delta_biased_access)
	  
	SimulData[[i_rep]] = list(theta_l_mat = theta_l_mat, theta_hat_l_mat = theta_hat_l_mat, X_cube = X_cube, Y_mat = Y_mat, beta_mat = beta_mat, V_beta_cube = V_beta_cube, true_mu = true_mu, delta_biased_access = delta_biased_access)
  
} # for (i_rep in 1:n_rep)

mean(fr_biased) # [1] 0.62375

save.image(file="SimulationData_2.RData")
