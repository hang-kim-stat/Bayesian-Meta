rm(list = ls())

#####################################
# Load simulation truth of Simulation Study 1 for comparison 
#####################################

load("../../data/SimulationData_1.RData")
p_mu = length(true_mu)

#####################################
# Compare methods 
#####################################

Method_list = c("1_Benchmark","2_IPD-AD","3_IPD-AD-pooled","4_IPD_only")
Table4 = NULL

for (i_method in 1:length(Method_list)){
  
  PostMean = PostVar = isCover1 = isCover2 = array(0,c(n_rep,2*p_mu))
  
  for (rep_no in 1:n_rep){
    
    load(paste0("../../output/Simulation_1/",Method_list[i_method],"/RData/rep_",rep_no,".RData"))
    
    PostMean[rep_no,1:p_mu] = apply(posterior_mu,2,mean)
    PostVar[rep_no,1:p_mu] = apply(posterior_mu,2,var)
    for (i_var in 1:p_mu){
      Post2_5 = quantile(posterior_mu[,i_var],probs=0.025)
      Post97_5 = quantile(posterior_mu[,i_var],probs=0.975)
      isCover1[rep_no,i_var] = (Post2_5<=true_mu[i_var]) && (true_mu[i_var]<=Post97_5)
    } # for (i_var in 1:p_mu)
    
    PostMean[rep_no,p_mu+(1:p_mu)] = apply(posterior_Sigma_theta,2,mean)
    PostVar[rep_no,p_mu+(1:p_mu)] = apply(posterior_Sigma_theta,2,var)
    for (i_var in 1:p_mu){
      Post2_5 = quantile(posterior_Sigma_theta[,i_var],probs=0.025)
      Post97_5 = quantile(posterior_Sigma_theta[,i_var],probs=0.975)
      isCover1[rep_no,p_mu+i_var] = (Post2_5<=true_Sigma_theta[i_var]) && (true_Sigma_theta[i_var]<=Post97_5)
    } # for (i_var in 1:p_mu)
    
  } # for (rep_no in 1:n_rep)
  
  RESULT = array(0,c(8,4))
  dimnames(RESULT)[[1]] = c("mu1","mu2","mu3","mu4","Sig11","Sig22","Sig33","Sig44")
  dimnames(RESULT)[[2]] = c("Est","Bias","MSE","CredCov95")
  RESULT[,"Est"] = apply(PostMean,2,mean)
  RESULT[,"Bias"] = apply(PostMean,2,mean)-c(true_mu,true_Sigma_theta)
  for (i_var in 1:p_theta){
    RESULT[i_var,"MSE"] = mean((PostMean[,i_var]-true_mu[i_var])^2)
  }
  for (i_var in 1:p_theta){
    RESULT[p_theta+i_var,"MSE"] = mean((PostMean[,i_var]-true_Sigma_theta[i_var])^2)
  }
  RESULT[,"CredCov95"] = apply(isCover1,2,mean)
  
  Table4 = cbind(Table4,RESULT)
  
} # for (i_method in 1:4)

##############################
write.csv(Table4,file="../../Output/Figure_Table/Supplement_Table4.csv")
##############################
