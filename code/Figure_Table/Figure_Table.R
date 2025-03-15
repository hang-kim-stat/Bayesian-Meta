rm(list=ls())

# library(gplots) ; library(MASS) ; library(xlsx) ; library(Hmisc) ; library(rlist) ; library(meta)

# if (!file.exists("Plot_Reference")){ dir.create("Plot_Reference", showWarnings = TRUE, recursive = FALSE, mode = "0777") }
# if (!file.exists("Plot_Manuscript")){ dir.create("Plot_Manuscript", showWarnings = TRUE, recursive = FALSE, mode = "0777") }













#############################

load("2_IPD_only.RData")

SEQ = (n_burnin+1):(n_burnin+n_main)

PostMean_IPD = apply(Draw_mu_vec[SEQ,],2,mean)
L95_IPD = U95_IPD = rep(0,7)
for (j in 1:7){
  L95_IPD[j] = quantile(Draw_mu_vec[SEQ,j],prob=0.025)
  U95_IPD[j] = quantile(Draw_mu_vec[SEQ,j],prob=0.975)
}

ipd.x <- Draw_mu_vec[SEQ,6]-Draw_mu_vec[SEQ,5] 
ipd.y <- Draw_mu_vec[SEQ,6]-Draw_mu_vec[SEQ,7]

dens.ipd.x <- density(ipd.x)
dens.ipd.y <- density(ipd.y)

ProbComparison = rep(0,3)
ProbComparison[1] = mean( Draw_mu_vec[SEQ,6] > Draw_mu_vec[SEQ,5] )  
ProbComparison[2] = mean( Draw_mu_vec[SEQ,6] > Draw_mu_vec[SEQ,7] )  
ProbComparison[3] = mean( ( Draw_mu_vec[SEQ,6] > Draw_mu_vec[SEQ,5] ) & ( Draw_mu_vec[SEQ,6] > Draw_mu_vec[SEQ,7] ) )
ProbComparison
# [1] 0.6107 0.5619 0.4468

#############################

ADD_NAME = "3_IPD_AD_meta_seed117_temp"
DATA_NAME = paste0(ADD_NAME,".RData")
RESULT_NAME = "final" 

##
load(DATA_NAME)

AfterBurn = (i_iter-10000+1):i_iter
SEQ = (i_iter-10000+1):i_iter

PostMean_IPD_AD = apply(Draw_mu_vec[AfterBurn,],2,mean)
L95_IPD_AD = U95_IPD_AD = rep(0,7)
for (j in 1:7){
	L95_IPD_AD[j] = quantile(Draw_mu_vec[AfterBurn,j],prob=0.025)
	U95_IPD_AD[j] = quantile(Draw_mu_vec[AfterBurn,j],prob=0.975)
}

ad.x <- Draw_mu_vec[AfterBurn,6]-Draw_mu_vec[AfterBurn,5] 
ad.y <- Draw_mu_vec[AfterBurn,6]-Draw_mu_vec[AfterBurn,7]

dens.ad.x <- density(ad.x)
dens.ad.y <- density(ad.y)

ProbComparison = rep(0,3)
ProbComparison[1] = mean( Draw_mu_vec[AfterBurn,6] > Draw_mu_vec[AfterBurn,5] )  # normal to overweight 
ProbComparison[2] = mean( Draw_mu_vec[AfterBurn,6] > Draw_mu_vec[AfterBurn,7] )  # obese to overweight
ProbComparison[3] = mean( ( Draw_mu_vec[AfterBurn,6] > Draw_mu_vec[AfterBurn,5] ) & ( Draw_mu_vec[AfterBurn,6] > Draw_mu_vec[AfterBurn,7] ) )
ProbComparison
# [1] 0.9470 0.9374 0.9090

##

TABLE = rbind(L95_IPD,PostMean_IPD,U95_IPD,L95_IPD_AD,PostMean_IPD_AD,U95_IPD_AD)
rownames(TABLE) = c("L95_IPD","PostMean_IPD","U95_IPD","L95_IPD_AD","PostMean_IPD_AD","U95_IPD_AD")
write.xlsx(TABLE, file=paste0("Plot_Reference/Result_",RESULT_NAME,".xlsx"))

##

FixedResult = read.xlsx(file="FixedEffectAnalysis/HangCleaning.xlsx", sheetName="Sheet1")
row.names(FixedResult) = FixedResult[,1]
FixedResult = FixedResult[,2:4]
dimnames(FixedResult)

####################################################################################################
# ForestPlot_pop.png
####################################################################################################

PlusMinus = 0.45

png(file=paste0("Plot_Manuscript/ForestPlot_pop_",RESULT_NAME,".png"),width=2000,height=800,pointsize=9,res=400)

par(mfrow=c(1,1),family="AppleGothic",mgp = c(2, 0.5, 0),mai=c(0.45,1,0.1,0.2)) # b l t r

plot(c(-1.4, 0.2), c(0, 6), type = "n", xlab = "Weight Gain Reduction (kg)", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
box()
# Add horizontal lines for effect sizes
segments(L95_IPD[5:7], c(5,3,1)+PlusMinus, U95_IPD[5:7], c(5,3,1)+PlusMinus, lty = "dashed", col="blue",lwd=1.2) # IPD-only # lty = 2 for dashed lines
segments(L95_IPD_AD[5:7], c(5,3,1), U95_IPD_AD[5:7], c(5,3,1), col="red",lwd=1.2)  # IPD-AD
segments(FixedResult[,"X95.cr.l"], c(5,3,1)-PlusMinus, FixedResult[,"X95.cr.r"], c(5,3,1)-PlusMinus, lty = "dotted", col="brown",lwd=1.2)  # Fixed effect model
# Add points for effect sizes
points(PostMean_IPD[5:7], c(5,3,1)+PlusMinus, pch = 17, cex=0.9, col="blue") # IPD-only
points(PostMean_IPD_AD[5:7], c(5,3,1), pch = 20, cex=1.2, col="red") # IPD-AD
points(FixedResult[,"Est"], c(5,3,1)-PlusMinus, pch = 15, cex=1, col="brown") # Fixed effect model
# Add a vertical line at effect size = 0
abline(v = 0, lty = "dotted", col="gray")  # lty = 2 for dashed line
# Add axes
axis(side = 1, at = seq(-1.6,0.4, by = 0.4), las = 1, cex.axis=0.8) # cex.axis=0.8
axis(side = 2, at = c(5,3,1), cex.axis = 1, labels=c(expression( BMI < 25),expression(paste(25 <=  BM, "",I < 30)),expression(BMI >= 30)), las=1)

dev.off()

####################################################################################################
# ContourPlot.png
####################################################################################################

NBINS = 31 

png(file=paste0("Plot_Manuscript/ContourPlot_",RESULT_NAME,".png"),width=2000,height=900,pointsize=10,res=400)
par(mfrow=c(1,2), mgp = c(2, 0.5, 0))

#
par(mai=c(0.5,0.6,0.1,0.1), cex.lab=1.1) # b l t r

ci2d( ipd.x, ipd.y, nbins=rep(NBINS,2), method="bkde2D", factor=1.0, ci.levels=c(0.50,0.75,0.85,0.90,0.95), show.points=FALSE, show = "contour", range.x=list(c(-2,2), c(-2,2)), col = "blue", xlab = expression(E(theta[l32])-E(theta[l31])), ylab = expression(E(theta[l32])-E(theta[l33])), pch=par("pch") )  
abline(h=0, v=0, lty="dotted", col="gray",lwd=1)

dens.ipd.x <- density(ipd.x)
dens.ipd.y <- density(ipd.y)
par(new = TRUE)
plot(dens.ipd.x$x, dens.ipd.x$y - 2, type = "l",  xlim=c(-2,2),ylim=c(-2,2) , xlab = "", ylab = "",xaxt='n',yaxt='n')
par(new = TRUE)
plot(-dens.ipd.y$y + 2, dens.ipd.y$x , type = "l",  xlim=c(-2,2),ylim=c(-2,2) , xlab = "", ylab = "",xaxt='n',yaxt='n')

#
par(mai=c(0.5,0.5,0.1,0.1), cex.lab=1.1) # b l t r

ci2d(ad.x, ad.y, nbins=rep(NBINS,2), method="bkde2D", factor=1.0, ci.levels=c(0.50,0.75,0.85,0.90,0.95), show.points=FALSE, show = "contour",  range.x=list(c(-2,2), c(-2,2)), col = "red", xlab =expression(E(theta[l32])-E(theta[l31])), ylab = "", pch=par("pch"))

abline(h=0, v=0, lty="dotted", col="gray",lwd=1)
dens.ad.x <- density(ad.x)
dens.ad.y <- density(ad.y)
par(new = TRUE)
plot(dens.ad.x$x, dens.ad.x$y - 2, type = "l",  xlim=c(-2,2),ylim=c(-2,2) , xlab = "", ylab = "", xaxt='n',yaxt='n')
par(new = TRUE)
plot(-dens.ad.y$y + 2, dens.ad.y$x , type = "l",  xlim=c(-2,2),ylim=c(-2,2) , xlab = "", ylab = "",xaxt='n',yaxt='n')
#
dev.off()

####################################################################################################
# FunnelPlot_Studylevel_Difference.png
####################################################################################################

### Nadine's log file
### 4 diet, 14 Excercise, 14 mixed methods,
study_id<-c("Khoury 2005", "Vitolo 2011", "Walsh 2012", "Wolff 2008",
            "Baciuk 2008","Barakat 2008", "Barakat 2012a", "Haakstad 2011","Khaledan 2010",
            "Nascimento 2011", "Ong 2009", "Oostdam 2012", "Perales 2014",
            "Prevedel 2003", "Ruiz 2013", "Stafne 2012",
            "Althuizen 2012", "Bogaerts 2012", "Dodd 2014", "Guelinckx 2010", "Harrison 2013",
            "Hui 2011", "Jeffries 2009", "Luoto 2011", "Petrella 2013", "Phelan 2011",
            "Poston unpub", "Rauh 2013", "Sagedal unpub", "Vinter 2011", "Renault 2013")

# data <- read.table("ipdLME2.csv", sep=",", header=T)
data <- read.table("S_RealData/IPD.csv", sep=",", header=T)

ipd.data <- data
study_id2 <- study_id
ipd.data$totgwg2 <- ipd.data$f_wt - ipd.data$b_wt

thres1 <- c(25,30)
thres2 <- c(15.9, 11.3, 9.1)
for(i in c(1: nrow(ipd.data))){
  if(ipd.data$b_bmi[i] < thres1[1]){
    if(ipd.data$totgwg2[i] > thres2[1]){
      ipd.data$wg_indi[i] <- 1
    } else{
      ipd.data$wg_indi[i] <- 0
    }
  } else if (ipd.data$b_bmi[i] > thres1[1]){
    if(ipd.data$totgwg2[i] > thres2[2]){
      ipd.data$wg_indi[i] <- 1
    } else{
      ipd.data$wg_indi[i] <- 0
    }
  } else if (ipd.data$b_bmi[i] >= thres1[2]){
    if(ipd.data$totgwg2[i] > thres2[3]){
      ipd.data$wg_indi[i] <- 1
    } else{
      ipd.data$wg_indi[i] <- 0
    }
  }
}

############
# nml
############

sum.data <- matrix(NA, ncol = 2, nrow = length(study_id2))
wg.data <- matrix(NA, ncol = 2, nrow = length(study_id2))
nsample.data <- sum.data
rownames(sum.data) <- study_id2
colnames(sum.data) <- c("diff_GWG_rate", "SE_diff_GWG_rate")
for (i in c(1: length(study_id2))){
  dat.tmp <- ipd.data[ipd.data$study_name == study_id2[i] & ipd.data$b_bmi < thres1[1], ]  # different by group
  n0 <- sum(dat.tmp$new_trt == 0)
  n1 <- sum(!dat.tmp$new_trt == 0)
  ctrl.sum.indi <- dat.tmp[dat.tmp$new_trt == 0,]$wg_indi
  inv.sum.indi <- dat.tmp[!dat.tmp$new_trt == 0,]$wg_indi
  ctrl.p <- sum(ctrl.sum.indi) / n0
  inv.p <-  sum(inv.sum.indi)  / n1

  ctrl.wg.tmp <- dat.tmp[dat.tmp$new_trt == 0,]$totgwg2
  inv.wg.tmp <- dat.tmp[!dat.tmp$new_trt == 0,]$totgwg2
  diff.wg <-  mean( inv.wg.tmp ) - mean( ctrl.wg.tmp)
  sd.wg <- sqrt(var(ctrl.wg.tmp) / n0 + var(inv.wg.tmp ) / n1)
  wg.data[i,] <- c( diff.wg, sd.wg)

  diff <-  inv.p - ctrl.p
  var <- ctrl.p * (1 - ctrl.p) / n0 + inv.p * (1 - inv.p) / n1
  sd <- sqrt(var)
  sum.data[i,] <- c(diff , sd)
  nsample.data[i,] <- c(n1, n0)
}
sum.data2_nml <- sum.data[!is.na(wg.data[,2]),]
nsample.data2_nml <- nsample.data[!is.na(wg.data[,2]),]
wg.data2_nml <- wg.data[!is.na(wg.data[,2]),]
ad.data_nml <- read.csv("S_RealData/ForGraph/ad.rate.bmi1.csv", header = TRUE)  # different by group
sizes_nml <- sqrt(rowSums(nsample.data2_nml)) / 1000

############
# over
############

sum.data <- matrix(NA, ncol = 2, nrow = length(study_id2))
wg.data <- matrix(NA, ncol = 2, nrow = length(study_id2))
nsample.data <- sum.data
rownames(sum.data) <- study_id2
colnames(sum.data) <- c("diff_GWG_rate", "SE_diff_GWG_rate")
for (i in c(1: length(study_id2))){
  dat.tmp <- ipd.data[ipd.data$study_name == study_id2[i] & ipd.data$b_bmi >= thres1[1] & ipd.data$b_bmi < thres1[2],] # different by group
  n0 <- sum(dat.tmp$new_trt == 0)
  n1 <- sum(!dat.tmp$new_trt == 0)
  ctrl.sum.indi <- dat.tmp[dat.tmp$new_trt == 0,]$wg_indi
  inv.sum.indi <- dat.tmp[!dat.tmp$new_trt == 0,]$wg_indi
  ctrl.p <- sum(ctrl.sum.indi) / n0
  inv.p <-  sum(inv.sum.indi)  / n1

  ctrl.wg.tmp <- dat.tmp[dat.tmp$new_trt == 0,]$totgwg2
  inv.wg.tmp <- dat.tmp[!dat.tmp$new_trt == 0,]$totgwg2
  diff.wg <-  mean( inv.wg.tmp ) - mean( ctrl.wg.tmp)
  sd.wg <- sqrt(var(ctrl.wg.tmp) / n0 + var(inv.wg.tmp ) / n1)
  wg.data[i,] <- c( diff.wg, sd.wg)

  diff <-  inv.p - ctrl.p
  var <- ctrl.p * (1 - ctrl.p) / n0 + inv.p * (1 - inv.p) / n1
  sd <- sqrt(var)
  sum.data[i,] <- c(diff , sd)
  nsample.data[i,] <- c(n1, n0)
}
sum.data2_over <- sum.data[!is.na(wg.data[,2]),]
nsample.data2_over <- nsample.data[!is.na(wg.data[,2]),]
wg.data2_over <- wg.data[!is.na(wg.data[,2]),]
ad.data_over <- read.csv("S_RealData/ForGraph/ad.rate.bmi2.csv", header = TRUE) # different by group
sizes_over <- sqrt(rowSums(nsample.data2_over)) / 1000

############
# obs
############

sum.data <- matrix(NA, ncol = 2, nrow = length(study_id2))
wg.data <- matrix(NA, ncol = 2, nrow = length(study_id2))
nsample.data <- sum.data
rownames(sum.data) <- study_id2
colnames(sum.data) <- c("diff_GWG_rate", "SE_diff_GWG_rate")
for (i in c(1: length(study_id2))){
  dat.tmp <- ipd.data[ipd.data$study_name == study_id2[i] & ipd.data$b_bmi >= thres1[2] ,] # different by group
  n0 <- sum(dat.tmp$new_trt == 0)
  n1 <- sum(!dat.tmp$new_trt == 0)
  ctrl.sum.indi <- dat.tmp[dat.tmp$new_trt == 0,]$wg_indi
  inv.sum.indi <- dat.tmp[!dat.tmp$new_trt == 0,]$wg_indi
  ctrl.p <- sum(ctrl.sum.indi) / n0
  inv.p <-  sum(inv.sum.indi)  / n1

  ctrl.wg.tmp <- dat.tmp[dat.tmp$new_trt == 0,]$totgwg2
  inv.wg.tmp <- dat.tmp[!dat.tmp$new_trt == 0,]$totgwg2
  diff.wg <-  mean( inv.wg.tmp ) - mean( ctrl.wg.tmp)
  sd.wg <- sqrt(var(ctrl.wg.tmp) / n0 + var(inv.wg.tmp ) / n1)
  wg.data[i,] <- c( diff.wg, sd.wg)

  diff <-  inv.p - ctrl.p
  var <- ctrl.p * (1 - ctrl.p) / n0 + inv.p * (1 - inv.p) / n1
  sd <- sqrt(var)
  sum.data[i,] <- c(diff , sd)
  nsample.data[i,] <- c(n1, n0)
}
sum.data2_obs <- sum.data[!is.na(wg.data[,2]),]
nsample.data2_obs <- nsample.data[!is.na(wg.data[,2]),]
wg.data2_obs <- wg.data[!is.na(wg.data[,2]),]
ad.data_obs <- read.csv("S_RealData/ForGraph/ad.rate.bmi3.csv", header = TRUE) # different by group
sizes_obs <- sqrt(rowSums(nsample.data2_obs)) / 1000

############
# graph 
############

png(file="Plot_Manuscript/FunnelPlot_Studylevel_Difference.png",width=2000,height=650,pointsize=10,res=400)

par(mfrow=c(1,3),family="AppleGothic",mgp = c(2, 0.5, 0)) # b l t r

c_lab = 1.2
xlimplot <- c(-0.6, 0.6) ; ylimplot <-  c(0, 0.4)

#
par(mai=c(0.35,0.4,0.15,0.0)) # b l t r

plot(sum.data2_nml, xlim = xlimplot, ylim = ylimplot, col = "white", xlab = "Difference in proportions", ylab = "Standard error", xaxt = "n", yaxt = "n", cex.lab=c_lab, main=expression( BMI < 25))
axis(side = 1, at = seq(-0.6, 0.6, 0.3), label = FALSE)
axis(side = 1, at = seq(-0.6, 0.6, 0.3), cex.axis = 1, line = 0, lty = "blank")
axis(side = 2, at = seq(0, 0.4, 0.2), cex.axis = 1)
polygon(c(xlimplot[1], mean(c(sum.data2_nml[,1],ad.data_nml[,1])), xlimplot[1]),
        c(0,0,1/qnorm(0.975) * (mean(c(sum.data2_nml[,1],ad.data_nml[,1])) -xlimplot[1])),col = "gray80")
polygon(c(mean(c(sum.data2_nml[,1],ad.data_nml[,1])), xlimplot[2],xlimplot[2]),
        c(0,0,1/qnorm(0.975) * (-mean(c(sum.data2_nml[,1],ad.data_nml[,1])) + xlimplot[2])),col = "gray80")
points(sum.data2_nml,cex = 1,col = "blue")
points(ad.data_nml, pch = 4,cex = 1,col = "red")

#
par(mai=c(0.35,0.35,0.15,0.05)) # b l t r

plot(sum.data2_over, xlim = xlimplot, ylim = ylimplot, col = "white", xlab = "Difference in proportions", ylab = "",xaxt = "n", yaxt = "n", cex.lab=c_lab, main=expression(paste(25 <=  BM, "",I < 30)))
axis(side = 1, at = seq(-0.6, 0.6, 0.3), label = FALSE)
axis(side = 1, at = seq(-0.6, 0.6, 0.3), cex.axis = 1, line = 0, lty = "blank")
axis(side = 2, at = seq(0, 0.4, 0.2), cex.axis = 1)
polygon(c(xlimplot[1], mean(c(sum.data2_over[,1],ad.data_over[,1])), xlimplot[1]),
        c(0,0,1/qnorm(0.975) * (mean(c(sum.data2_over[,1],ad.data_over[,1])) -xlimplot[1])),col = "gray80")
polygon(c(mean(c(sum.data2_over[,1],ad.data_over[,1])), xlimplot[2],xlimplot[2]),
        c(0,0,1/qnorm(0.975) * (-mean(c(sum.data2_over[,1],ad.data_over[,1])) + xlimplot[2])),col = "gray80")
points(sum.data2_over,cex = 1,col = "blue")
points(ad.data_over, pch = 4,cex = 1, col = "red")

#
par(mai=c(0.35,0.3,0.15,0.1)) # b l t r

plot(sum.data2_obs, xlim = xlimplot, ylim = ylimplot, col = "white", xlab = "Difference in proportions", ylab = "",xaxt = "n", yaxt = "n", cex.lab=c_lab, main=expression(BMI >= 30))
axis(side = 1, at = seq(-0.6, 0.6, 0.3), label = FALSE)
axis(side = 1, at = seq(-0.6, 0.6, 0.3), cex.axis = 1, line = 0, lty = "blank")
axis(side = 2, at = seq(0, 0.4, 0.2), cex.axis = 1)
polygon(c(xlimplot[1], mean(c(sum.data2_obs[,1],ad.data_obs[,1])), xlimplot[1]),
        c(0,0,1/qnorm(0.975) * (mean(c(sum.data2_obs[,1],ad.data_obs[,1])) -xlimplot[1])),col = "gray80")
polygon(c(mean(c(sum.data2_obs[,1],ad.data_obs[,1])), xlimplot[2],xlimplot[2]),
        c(0,0,1/qnorm(0.975) * (-mean(c(sum.data2_obs[,1],ad.data_obs[,1])) + xlimplot[2])),col = "gray80")
points(sum.data2_obs,cex = 1,col = "blue")
points(ad.data_obs, pch = 4,cex = 1,col = "red")

#
dev.off()
