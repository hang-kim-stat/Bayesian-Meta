rm(list=ls())

library(xlsx) ; library(gplots)

#############################

# Load IPD-only analysis results 
load("../../Output/Application_Study/IPDonly_Posterior.RData")

PostMean_IPD = apply(posterior_mu,2,mean)
L95_IPD = U95_IPD = rep(0,7)
for (j in 1:7){
  L95_IPD[j] = quantile(posterior_mu[,j],prob=0.025)
  U95_IPD[j] = quantile(posterior_mu[,j],prob=0.975)
}
ipd.x <- posterior_mu[,6]-posterior_mu[,5] 
ipd.y <- posterior_mu[,6]-posterior_mu[,7]
dens.ipd.x <- density(ipd.x)
dens.ipd.y <- density(ipd.y)

ProbComparison = rep(0,3)
ProbComparison[1] = mean( posterior_mu[,6] > posterior_mu[,5] )  
ProbComparison[2] = mean( posterior_mu[,6] > posterior_mu[,7] )  
ProbComparison[3] = mean( ( posterior_mu[,6] > posterior_mu[,5] ) & ( posterior_mu[,6] > posterior_mu[,7] ) )
ProbComparison # [1] 0.6107 0.5619 0.4468 

#############################

# Load IPD-AD integrated analysis results 
load("../../Output/Application_Study/IPD_AD_Posterior.RData")
     
PostMean_IPD_AD = apply(posterior_mu,2,mean)
L95_IPD_AD = U95_IPD_AD = rep(0,7)
for (j in 1:7){
	L95_IPD_AD[j] = quantile(posterior_mu[,j],prob=0.025)
	U95_IPD_AD[j] = quantile(posterior_mu[,j],prob=0.975)
}
ipd_ad.x <- posterior_mu[,6]-posterior_mu[,5] 
ipd_ad.y <- posterior_mu[,6]-posterior_mu[,7]
dens.ipd_ad.x <- density(ipd_ad.x)
dens.ipd_ad.y <- density(ipd_ad.y)

ProbComparison = rep(0,3)
ProbComparison[1] = mean( posterior_mu[,6] > posterior_mu[,5] )  # normal to overweight 
ProbComparison[2] = mean( posterior_mu[,6] > posterior_mu[,7] )  # obese to overweight
ProbComparison[3] = mean( ( posterior_mu[,6] > posterior_mu[,5] ) & ( posterior_mu[,6] > posterior_mu[,7] ) )
ProbComparison # [1] 0.9470 0.9374 0.9090

#############################

# Load fixed-effect analysis results 

FixedResult = read.xlsx(file="../../Output/Application_Study/FixedEffect.xlsx", sheetName="Sheet1")
row.names(FixedResult) = FixedResult[,1] ; FixedResult = FixedResult[,2:4]

#############################

# Figure 4 of the main text

png(file=paste0("../../Output/Figure_Table/Maintext_Figure4.png"),width=2000,height=800,pointsize=9,res=400)

PlusMinus = 0.45

par(mfrow=c(1,1),family="AppleGothic",mgp = c(2, 0.5, 0),mai=c(0.45,1,0.1,0.2)) # b l t r

plot(c(-1.4, 0.2), c(0, 6), type = "n", xlab = "Weight Gain Reduction (kg)", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
box()
segments(L95_IPD[5:7], c(5,3,1)+PlusMinus, U95_IPD[5:7], c(5,3,1)+PlusMinus, lty = "dashed", col="blue",lwd=1.2) # IPD-only 
segments(L95_IPD_AD[5:7], c(5,3,1), U95_IPD_AD[5:7], c(5,3,1), col="red",lwd=1.2)  # IPD-AD
segments(FixedResult[,"X95.cr.l"], c(5,3,1)-PlusMinus, FixedResult[,"X95.cr.r"], c(5,3,1)-PlusMinus, lty = "dotted", col="brown",lwd=1.2)  # Fixed effect model
# Add points for effect sizes
points(PostMean_IPD[5:7], c(5,3,1)+PlusMinus, pch = 17, cex=0.9, col="blue") # IPD-only
points(PostMean_IPD_AD[5:7], c(5,3,1), pch = 20, cex=1.2, col="red") # IPD-AD
points(FixedResult[,"Est"], c(5,3,1)-PlusMinus, pch = 15, cex=1, col="brown") # Fixed effect model
abline(v = 0, lty = "dotted", col="gray")  # Add a vertical line at effect size = 0
axis(side = 1, at = seq(-1.6,0.4, by = 0.4), las = 1, cex.axis=0.8) # Add axes
axis(side = 2, at = c(5,3,1), cex.axis = 1, labels=c(expression( BMI < 25),expression(paste(25 <=  BM, "",I < 30)),expression(BMI >= 30)), las=1) # Add axes

dev.off()

# Figure 5 of the main text

png(file=paste0("../../Output/Figure_Table/Maintext_Figure5.png"),width=2000,height=900,pointsize=10,res=400)

par(mfrow=c(1,2), mgp = c(2, 0.5, 0), mai=c(0.5,0.6,0.1,0.1), cex.lab=1.1) # b l t r
NBINS = 31 

# IPD-only result 
ci2d( ipd.x, ipd.y, nbins=rep(NBINS,2), method="bkde2D", factor=1.0, ci.levels=c(0.50,0.75,0.85,0.90,0.95), show.points=FALSE, show = "contour", range.x=list(c(-2,2), c(-2,2)), col = "blue", xlab = expression(E(theta[l32])-E(theta[l31])), ylab = expression(E(theta[l32])-E(theta[l33])), pch=par("pch") )  
abline(h=0, v=0, lty="dotted", col="gray",lwd=1)

par(new = TRUE)
plot(dens.ipd.x$x, dens.ipd.x$y - 2, type = "l",  xlim=c(-2,2),ylim=c(-2,2) , xlab = "", ylab = "",xaxt='n',yaxt='n')

par(new = TRUE)
plot(-dens.ipd.y$y + 2, dens.ipd.y$x , type = "l",  xlim=c(-2,2),ylim=c(-2,2) , xlab = "", ylab = "",xaxt='n',yaxt='n')

# IPD-AD result 
par(mai=c(0.5,0.5,0.1,0.1), cex.lab=1.1) # b l t r

ci2d(ipd_ad.x, ipd_ad.y, nbins=rep(NBINS,2), method="bkde2D", factor=1.0, ci.levels=c(0.50,0.75,0.85,0.90,0.95), show.points=FALSE, show = "contour",  range.x=list(c(-2,2), c(-2,2)), col = "red", xlab =expression(E(theta[l32])-E(theta[l31])), ylab = "", pch=par("pch"))
abline(h=0, v=0, lty="dotted", col="gray",lwd=1)

par(new = TRUE)
plot(dens.ipd_ad.x$x, dens.ipd_ad.x$y - 2, type = "l",  xlim=c(-2,2),ylim=c(-2,2) , xlab = "", ylab = "", xaxt='n',yaxt='n')

par(new = TRUE)
plot(-dens.ipd_ad.y$y + 2, dens.ipd_ad.y$x , type = "l",  xlim=c(-2,2),ylim=c(-2,2) , xlab = "", ylab = "",xaxt='n',yaxt='n')

dev.off()