library(Zresidual)
library(survival)
data("LeukSurv")
head(LeukSurv)
LeukSurv<-LeukSurv[LeukSurv$age<60,]
is.factor(LeukSurv$district)
is.factor(LeukSurv$sex)
LeukSurv$district<-as.factor(LeukSurv$district)
LeukSurv$sex<-as.factor(LeukSurv$sex)

#############fit the coxph model with frailty#################
fit_LeukSurv <- tryCatch(
  coxph(Surv(time, cens) ~ age  +sex+ wbc +tpi  +
          frailty(district, distribution="gamma"), data= LeukSurv),
  error = function(e) NA,
  warning = function(w) NA
)
#############Calculate residuals###############################
allresid.LeukSurv<-Zresidual(fit.object = fit_LeukSurv,data= LeukSurv)

#Z-residual
Zresid.LeukSurv<-allresid.LeukSurv$Zresid
#censored Z-residual
censored.Zresid.LeukSurv<-allresid.LeukSurv$censored.Zresid
#Survival Probabilities
sp.LeukSurv<-allresid.LeukSurv$SP
#unmodified CS residuals
ucs.LeukSurv<-allresid.LeukSurv$ucs
#modified CS residuals
mcs.LeukSurv<-allresid.LeukSurv$mcs
#Martingale residuals
martg.LeukSurv<-allresid.LeukSurv$martg
#Deviance residuals
dev.LeukSurv<-allresid.LeukSurv$dev


########################### Graphical Diagnosis##################
plot.zresid(Zresid.LeukSurv)
plot.zresid.fitted.value(Zresid.LeukSurv,
                         fitted.values=fit_LeukSurv$linear.predictors,
                         xlab="Linear Predictor",
                         cenered=LeukSurv$cens)
qqnorm.zresid (Zresid.LeukSurv)
boxplot.zresid(Zresid.LeukSurv,fitted.values=LeukSurv$wbc)
boxplot.zresid(Zresid.LeukSurv,fitted.values=LeukSurv$sex)

########################### other residuals Graphical Diagnosis ##################
##unmodified CS residuals
km.ln.LeukSurv <- survfit(Surv(ucs.LeukSurv, LeukSurv$cens)~1,type='fleming')
id.ln.LeukSurv<-order(ucs.LeukSurv)

plot(km.ln.LeukSurv, fun="cumhaz", xlab=("Cox-Snell Residuals"),
     ylab=("Cumulative Hazard Function"),
     main="Log-normal, Cum. Hazard of CS Residuals",
     ylim= c(0,4),xlim=c(0,4))
abline(0, 1, col="red", lty=2)
points(km.ln.LeukSurv$time, -log(km.ln.LeukSurv$surv),
       col=c("blue","darkolivegreen4")[LeukSurv$cens[id.ln.LeukSurv]+1],
       pch=c(3,2)[LeukSurv$cens[id.ln.LeukSurv]+1] )
legend(x = "topleft",
       legend = c("Uncensored", "Censored"), col=c("darkolivegreen4","blue"),
       pch=c(2,3),cex=1,xpd = TRUE,bty="L")


#Martingale residuals
plot(LeukSurv$age,martg.LeukSurv,ylab="Martingale",
     xlab="age",
     main="Martingale Residual Plot",
     col=c("blue","darkolivegreen4")[LeukSurv$cens+1],
     pch=c(3,2)[LeukSurv$cens+1])
#abline(h=c(3,-3),col="grey")
lines(lowess(martg.LeukSurv~ LeukSurv$age),
      col = "red",lwd = 3)
legend(x = "bottomright",
       legend = c("Uncensored", "Censored"),
       col=c("darkolivegreen4","blue"),
       pch=c(2,3),cex=1,xpd = TRUE,bty="l")

#Deviance residuals
plot(LeukSurv$age ,dev.LeukSurv,ylab="Deviance",
     xlab="size",
     main="Deviance Residual Plot",
     col=c("blue","darkolivegreen4")[LeukSurv$cens+1],
     pch=c(3,2)[LeukSurv$cens+1])
#abline(h=c(3,-3),col="grey")
lines(lowess(dev.LeukSurv~ LeukSurv$age),
      col = "red",lwd = 3)
legend(x = "bottomright",
       legend = c("Uncensored", "Censored"),
       col=c("darkolivegreen4","blue"),
       pch=c(2,3),cex=1,xpd = TRUE,bty="l")


########################### Statistical test Diagnosis##################
sw.test.zresid(Zresid.LeukSurv)
sf.test.zresid(Zresid.LeukSurv)
gof.censore.zresid(censored.Zresidual=censored.Zresid.LeukSurv,censored=LeukSurv$cens)

anov.test.zresid(Zresid.LeukSurv,fitted.values=fit_LeukSurv$linear.predictors, k.anova=10)
bartlett.test.zresid(Zresid.LeukSurv,fitted.values=fit_LeukSurv$linear.predictors, k.bl=10)

anov.test.zresid(Zresid.LeukSurv,fitted.values=LeukSurv$wbc, k.anova=10)
bartlett.test.zresid(Zresid.LeukSurv,fitted.values=LeukSurv$wbc, k.bl=10)

anov.test.zresid(Zresid.LeukSurv,fitted.values=LeukSurv$sex, k.anova=10)
bartlett.test.zresid(Zresid.LeukSurv,fitted.values=LeukSurv$sex, k.bl=10)

################pmin value############################

n_sims<-1000
cur_time = proc.time()
sw.LeukSurv<- rep(0,n_sims)
sf.LeukSurv<- rep(0,n_sims)
anov.LeukSurv.lp<- rep(0,n_sims)

for(j in 1:n_sims ){
  cat(paste('Simulation ',j,' out of ',n_sims,'\n'))
  if(j ==2){
    elapsed=as.numeric(proc.time()-cur_time)[3]
    cat(paste("Time for 1 simulation: ",elapsed/3600," hours \n"))
    cat(paste("Estimated time remaining: ",elapsed/3600*(n_sims-1)," hours \n"))
  }
  allresid.LeukSurv<-Zresidual(fit.object = fit_LeukSurv,data= LeukSurv)
  Zresid.LeukSurv<-allresid.LeukSurv$Zresid

  sw.LeukSurv[j]<-sw.test.zresid(Zresid.LeukSurv)
  sf.LeukSurv[j]<-sf.test.zresid(Zresid.LeukSurv)

  anov.LeukSurv.lp[j]<-anov.test.zresid(Zresid.LeukSurv,
                                   fitted.values=fit_LeukSurv$linear.predictors,
                                   k.anova=10)
}
pmin.sw.LeukSurv<-bounds_pvalues(pv=sw.LeukSurv);pmin.sw.LeukSurv
pmin.sf.LeukSurv<-bounds_pvalues(pv=sf.LeukSurv);pmin.sf.LeukSurv
pmin.aov.lp.LeukSurv<-bounds_pvalues(pv=anov.LeukSurv.lp);pmin.aov.lp.LeukSurv


