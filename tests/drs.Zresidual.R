library(Zresidual)
library(survival)
data("drs")
head(drs)
drs$subject_id<-as.factor(drs$subject_id)
drs$treated<-as.factor(drs$treated)
drs$laser_type<-as.factor(drs$laser_type)
drs$diabetes_type<-as.factor(drs$diabetes_type)


#############fit the coxph model with frailty#################
fit_drs  <- tryCatch(
  coxph(Surv(time, status) ~ treated + age_at_onset +
          laser_type+ diabetes_type+
          frailty(subject_id, distribution="gamma"), data= drs),
  error = function(e) NA,
  warning = function(w) NA
)
#############Calculate residuals###############################
allresid.drs<-Zresidual(fit.object = fit_drs,data= drs)
#Z-residual
Zresid.drs<-allresid.drs$Zresid
#censored Z-residual
censored.Zresid.drs<-allresid.drs$censored.Zresid
#Survival Probabilities
sp.drs<-allresid.drs$SP
#unmodified CS residuals
ucs.drs<-allresid.drs$ucs
#modified CS residuals
mcs.drs<-allresid.drs$mcs
#Martingale residuals
martg.drs<-allresid.drs$martg
#Deviance residuals
dev.drs<-allresid.drs$dev

########################### Graphical Diagnosis##################
plot.zresid(Zresid.drs)
plot.zresid.fitted.value(Zresid.drs,
                         fitted.values=fit_drs$linear.predictors,
                         xlab="Linear Predictor",
                         cenered=drs$status)
qqnorm.zresid (Zresid.drs)
boxplot.zresid(Zresid.drs,fitted.values=drs$age_at_onset)
boxplot.zresid(Zresid.drs,fitted.values=drs$treated)

########################### other residuals Graphical Diagnosis ##################
##unmodified CS residuals
km.ln.drs <- survfit(Surv(ucs.drs, drs$status)~1,type='fleming')
id.ln.drs<-order(ucs.drs)

plot(km.ln.drs, fun="cumhaz", xlab=("Cox-Snell Residuals"),
     ylab=("Cumulative Hazard Function"),
     main="Log-normal, Cum. Hazard of CS Residuals",
     ylim= c(0,4),xlim=c(0,4))
abline(0, 1, col="red", lty=2)
points(km.ln.drs$time, -log(km.ln.drs$surv),
       col=c("blue","darkolivegreen4")[drs$status[id.ln.drs]+1],
       pch=c(3,2)[drs$status[id.ln.drs]+1] )
legend(x = "topleft",
       legend = c("Uncensored", "Censored"), col=c("darkolivegreen4","blue"),
       pch=c(2,3),cex=1,xpd = TRUE,bty="L")


#Martingale residuals
plot(drs$age_at_onset,martg.drs,ylab="Martingale",
     xlab="age",
     main="Martingale Residual Plot",
     col=c("blue","darkolivegreen4")[drs$status+1],
     pch=c(3,2)[drs$status+1])
#abline(h=c(3,-3),col="grey")
lines(lowess(martg.drs~ drs$age_at_onset),
      col = "red",lwd = 3)
legend(x = "bottomright",
       legend = c("Uncensored", "Censored"),
       col=c("darkolivegreen4","blue"),
       pch=c(2,3),cex=1,xpd = TRUE,bty="l")

#Deviance residuals
plot(drs$age_at_onset ,dev.drs,ylab="Deviance",
     xlab="size",
     main="Deviance Residual Plot",
     col=c("blue","darkolivegreen4")[drs$status+1],
     pch=c(3,2)[drs$status+1])
#abline(h=c(3,-3),col="grey")
lines(lowess(dev.drs~ drs$age_at_onset),
      col = "red",lwd = 3)
legend(x = "bottomright",
       legend = c("Uncensored", "Censored"),
       col=c("darkolivegreen4","blue"),
       pch=c(2,3),cex=1,xpd = TRUE,bty="l")

########################### Statistical test Diagnosis##################
sw.test.zresid(Zresid.drs)
sf.test.zresid(Zresid.drs)
gof.censore.zresid(censored.Zresidual=censored.Zresid.drs,censored=drs$status)

anov.test.zresid(Zresid.drs,fitted.values=fit_drs$linear.predictors, k.anova=10)
bartlett.test.zresid(Zresid.drs,fitted.values=fit_drs$linear.predictors, k.bl=10)

anov.test.zresid(Zresid.drs,fitted.values=drs$age_at_onset, k.anova=10)
bartlett.test.zresid(Zresid.drs,fitted.values=drs$age_at_onset, k.bl=10)

anov.test.zresid(Zresid.drs,fitted.values=drs$treated, k.anova=10)
bartlett.test.zresid(Zresid.drs,fitted.values=drs$treated, k.bl=10)

################pmin value############################

n_sims<-1000
cur_time = proc.time()
sw.drs<- rep(0,n_sims)
sf.drs<- rep(0,n_sims)
anov.drs.lp<- rep(0,n_sims)

for(j in 1:n_sims ){
  cat(paste('Simulation ',j,' out of ',n_sims,'\n'))
  if(j ==2){
    elapsed=as.numeric(proc.time()-cur_time)[3]
    cat(paste("Time for 1 simulation: ",elapsed/3600," hours \n"))
    cat(paste("Estimated time remaining: ",elapsed/3600*(n_sims-1)," hours \n"))
  }
  allresid.drs<-Zresidual(fit.object = fit_drs,data= drs)
  Zresid.drs<-allresid.drs$Zresid

  sw.drs[j]<-sw.test.zresid(Zresid.drs)
  sf.drs[j]<-sf.test.zresid(Zresid.drs)

  anov.drs.lp[j]<-anov.test.zresid(Zresid.drs,
                                      fitted.values=fit_drs$linear.predictors,
                                      k.anova=10)
}
pmin.sw.drs<-bounds_pvalues(pv=sw.drs);pmin.sw.drs
pmin.sf.drs<-bounds_pvalues(pv=sf.drs);pmin.sf.drs
pmin.aov.lp.drs<-bounds_pvalues(pv=anov.drs.lp);pmin.aov.lp.drs

