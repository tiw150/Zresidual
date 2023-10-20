library(Zresidual)
library(survival)
data("BreastCancer")
BreastCancer<-data
head(BreastCancer)
BreastCancer$grade <- as.factor(BreastCancer$grade)

#############fit the coxph model with frailty#################
fit_AFTlognormal_BreastCancer<- survreg(Surv(BreastCancer$time,BreastCancer$status)~.,
                                        data=BreastCancer[,-c(1:3)],dist="lognormal")

#############Calculate residuals###############################
allresid.BreastCancer<-Zresidual(fit.object = fit_AFTlognormal_BreastCancer,
                                 data=BreastCancer)

#Z-residual
Zresid.BreastCancer<-allresid.BreastCancer$Zresid
#censored Z-residual
censored.Zresid.BreastCancer<-allresid.BreastCancer$censored.Zresid
#Survival Probabilities
sp.BreastCancer<-allresid.BreastCancer$SP
#unmodified CS residuals
ucs.BreastCancer<-allresid.BreastCancer$ucs
#modified CS residuals
mcs.BreastCancer<-allresid.BreastCancer$mcs
#Martingale residuals
martg.BreastCancer<-allresid.BreastCancer$martg
#Deviance residuals
dev.BreastCancer<-allresid.BreastCancer$dev

########################### Z-residual Graphical Diagnosis ##################
plot.zresid(Zresid.BreastCancer)

qqnorm.zresid (Zresid.BreastCancer)
boxplot.zresid(Zresid.BreastCancer,
               fitted.values=fit_AFTlognormal_BreastCancer$linear.predictors)

plot.zresid.fitted.value(Zresid.BreastCancer,
                         fitted.values=fit_AFTlognormal_BreastCancer$linear.predictors,
                         xlab="Linear Predictor",
                         cenered=BreastCancer$status)

boxplot.zresid(Zresid.BreastCancer,fitted.values=BreastCancer$age)
boxplot.zresid(Zresid.BreastCancer,fitted.values=BreastCancer$grade)

########################### other residuals Graphical Diagnosis ##################
##unmodified CS residuals
km.ln.BreastCancer <- survfit(Surv(ucs.BreastCancer, BreastCancer$status)~1,type='fleming')
id.ln.BreastCancer<-order(ucs.BreastCancer)

plot(km.ln.BreastCancer, fun="cumhaz", xlab=("Cox-Snell Residuals"),
     ylab=("Cumulative Hazard Function"),
     main="Log-normal, Cum. Hazard of CS Residuals",
     ylim= c(0,4),xlim=c(0,4))
abline(0, 1, col="red", lty=2)
points(km.ln.BreastCancer$time, -log(km.ln.BreastCancer$surv),
       col=c("blue","darkolivegreen4")[BreastCancer$status[id.ln.BreastCancer]+1],
       pch=c(3,2)[BreastCancer$status[id.ln.BreastCancer]+1] )
legend(x = "topleft",
       legend = c("Uncensored", "Censored"), col=c("darkolivegreen4","blue"),
       pch=c(2,3),cex=1,xpd = TRUE,bty="L")


#Martingale residuals
plot(BreastCancer$size,martg.BreastCancer,ylab="Martingale",
     xlab="size",
     main="Martingale Residual Plot",
     col=c("blue","darkolivegreen4")[BreastCancer$status+1],
     pch=c(3,2)[BreastCancer$status+1])
#abline(h=c(3,-3),col="grey")
lines(lowess(martg.BreastCancer~ BreastCancer$size),
      col = "red",lwd = 3)
legend(x = "bottomright",
       legend = c("Uncensored", "Censored"),
       col=c("darkolivegreen4","blue"),
       pch=c(2,3),cex=1,xpd = TRUE,bty="l")

#Deviance residuals
plot(BreastCancer$size ,dev.BreastCancer,ylab="Deviance",
     xlab="size",
     main="Deviance Residual Plot",
     col=c("blue","darkolivegreen4")[BreastCancer$status+1],
     pch=c(3,2)[BreastCancer$status+1])
#abline(h=c(3,-3),col="grey")
lines(lowess(dev.BreastCancer~ BreastCancer$size),
      col = "red",lwd = 3)
legend(x = "bottomright",
       legend = c("Uncensored", "Censored"),
       col=c("darkolivegreen4","blue"),
       pch=c(2,3),cex=1,xpd = TRUE,bty="l")


########################### Statistical test Diagnosis##################
sw.test.zresid(Zresid.BreastCancer)
sf.test.zresid(Zresid.BreastCancer)
gof.censore.zresid(censored.Zresidual=censored.Zresid.BreastCancer,censored=BreastCancer$status)

anov.test.zresid(Zresid.BreastCancer,fitted.values=fit_AFTlognormal_BreastCancer$linear.predictors, k.anova=10)
bartlett.test.zresid(Zresid.BreastCancer,fitted.values=fit_AFTlognormal_BreastCancer$linear.predictors, k.bl=10)

anov.test.zresid(Zresid.BreastCancer,fitted.values=BreastCancer$age, k.anova=10)
bartlett.test.zresid(Zresid.BreastCancer,fitted.values=BreastCancer$age, k.bl=10)

anov.test.zresid(Zresid.BreastCancer,fitted.values=BreastCancer$grade, k.anova=10)
bartlett.test.zresid(Zresid.BreastCancer,fitted.values=BreastCancer$grade, k.bl=10)

################pmin value############################

n_sims<-1000
cur_time = proc.time()
sw.BreastCancer<- rep(0,n_sims)
sf.BreastCancer<- rep(0,n_sims)
anov.BreastCancer.lp<- rep(0,n_sims)

for(j in 1:n_sims ){
  cat(paste('Simulation ',j,' out of ',n_sims,'\n'))
  if(j ==2){
    elapsed=as.numeric(proc.time()-cur_time)[3]
    cat(paste("Time for 1 simulation: ",elapsed/3600," hours \n"))
    cat(paste("Estimated time remaining: ",elapsed/3600*(n_sims-1)," hours \n"))
  }
  allresid.BreastCancer<-Zresidual(fit.object = fit_AFTlognormal_BreastCancer,data= BreastCancer)
  Zresid.BreastCancer<-allresid.BreastCancer$Zresid

  sw.BreastCancer[j]<-sw.test.zresid(Zresid.BreastCancer)
  sf.BreastCancer[j]<-sf.test.zresid(Zresid.BreastCancer)

  anov.BreastCancer.lp[j]<-anov.test.zresid(Zresid.BreastCancer,
                                      fitted.values=fit_AFTlognormal_BreastCancer$linear.predictors,
                                      k.anova=10)
}
pmin.sw.BreastCancer<-bounds_pvalues(pv=sw.BreastCancer);pmin.sw.BreastCancer
pmin.sf.BreastCancer<-bounds_pvalues(pv=sf.BreastCancer);pmin.sf.BreastCancer
pmin.aov.lp.BreastCancer<-bounds_pvalues(pv=anov.BreastCancer.lp);pmin.aov.lp.BreastCancer

