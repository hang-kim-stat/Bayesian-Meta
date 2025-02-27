rm(list = ls())

load("0_SimulData.RData")

simulSEQ = 1:300
n_simul = length(simulSEQ)
p_mu = length(true_mu)
AfterBurnin = 10001:20000

Method_list = c("1_Benchmark","2_IPD-AD","3_IPD-AD-pooled","4_IPD-only")

Table1 = NULL

for (i_method in 1:length(Method_list)){
  
  PostMean = PostVar = isCover1 = isCover2 = array(0,c(n_simul,2*p_mu))
  
  for (i_simul in 1:n_simul){
    
    load(paste0(Method_list[i_method],"/Rep_",simulSEQ[i_simul],".RData"))
    
    PostMean[i_simul,1:p_mu] = apply(draw_mu[AfterBurnin,],2,mean)
    PostVar[i_simul,1:p_mu] = apply(draw_mu[AfterBurnin,],2,var)
    for (i_var in 1:p_mu){
      Post2_5 = quantile(draw_mu[AfterBurnin,i_var],probs=0.025)
      Post97_5 = quantile(draw_mu[AfterBurnin,i_var],probs=0.975)
      isCover1[i_simul,i_var] = (Post2_5<=true_mu[i_var]) && (true_mu[i_var]<=Post97_5)
    } # for (i_var in 1:p_mu)
    
    PostMean[i_simul,p_mu+(1:p_mu)] = apply(draw_Sigma_theta[AfterBurnin,],2,mean)
    PostVar[i_simul,p_mu+(1:p_mu)] = apply(draw_Sigma_theta[AfterBurnin,],2,var)
    for (i_var in 1:p_mu){
      Post2_5 = quantile(draw_Sigma_theta[AfterBurnin,i_var],probs=0.025)
      Post97_5 = quantile(draw_Sigma_theta[AfterBurnin,i_var],probs=0.975)
      isCover1[i_simul,p_mu+i_var] = (Post2_5<=true_Sigma_theta[i_var]) && (true_Sigma_theta[i_var]<=Post97_5)
    } # for (i_var in 1:p_mu)
    
  } # for (i_simul)
  
  #
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
  
  #
  
  Table1 = cbind(Table1,RESULT)
  
} # for (i_method in 1:4)


##############################
write.csv(Table1,file="W_Table1.csv")
##############################

Table2 = array(NA,c(16,5))
colnames(Table2) = c("Method","Parameter","Bias","MSE","Coverage")
Table2[,1] = rep(c("Benchmark","IPD-AD","IPD-AD(pooled)","IPD only"),4)
Table2[,2] = rep(c("mu1","mu2","mu3","mu4"),each=4)
for (i_method in 1:4){
  SEQ1 = c(0,4,8,12)+i_method
  SEQ2 = (2:4)+(i_method-1)*4
  Table2[SEQ1,3:5] = Table1[1:4,SEQ2]
} # 

write.csv(Table2,row.names=F,file="W_Table2.csv")
