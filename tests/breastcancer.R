library(Zresidual)
library(survival)
data("BreastCancer")
BreastCancer<-data
head(BreastCancer)
BreastCancer$grade <- as.factor(BreastCancer$grade)

#############fit the coxph model with frailty#################
fit_AFTlognormal_BreastCancer<- survreg(Surv(time,status)~treat+age+men+
                                          size+grade+nodes+prog+oest,
                                        data=BreastCancer,dist="lognormal")

#############Calculate residuals###############################
Zresid.BreastCancer<-Zresidual(fit.object = fit_AFTlognormal_BreastCancer)

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
gof.censored.zresidual(censored.Zresidual=censored.Zresid.BreastCancer,censored=BreastCancer$status)

aov.test.zresid(Zresid.BreastCancer,fitted.values=fit_AFTlognormal_BreastCancer$linear.predictors, k.anova=10)
bartlett.test.zresid(Zresid.BreastCancer,fitted.values=fit_AFTlognormal_BreastCancer$linear.predictors, k.bl=10)

aov.test.zresid(Zresid.BreastCancer,fitted.values=BreastCancer$age, k.anova=10)
bartlett.test.zresid(Zresid.BreastCancer,fitted.values=BreastCancer$age, k.bl=10)

aov.test.zresid(Zresid.BreastCancer,fitted.values=BreastCancer$grade, k.anova=10)
bartlett.test.zresid(Zresid.BreastCancer,fitted.values=BreastCancer$grade, k.bl=10)

\
