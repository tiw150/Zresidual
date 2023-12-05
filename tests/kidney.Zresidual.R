library(Zresidual)
library(survival)
data("kidney")
head(kidney)
kidney$sex <- ifelse(kidney$sex == 1, "male", "female")
kidney$sex<-as.factor(kidney$sex)
kidney$id<-as.factor(kidney$id)

#############fit the coxph model with frailty#################
fit_kidney <- tryCatch(
  coxph(Surv(time, status) ~ age + sex + disease+
          frailty(id, distribution="gamma"), data= kidney),
  error = function(e) NA,
  warning = function(w) NA
)
#############Calculate residuals###############################

allresid.kidney<-Zresidual(fit.object = fit_kidney,data= kidney)

#Z-residual
Zresid.kidney<-allresid.kidney$Zresid
#censored Z-residual
censored.Zresid.kidney<-allresid.kidney$censored.Zresid
#Survival Probabilities
sp.kidney<-allresid.kidney$SP
#unmodified CS residuals
ucs.kidney<-allresid.kidney$ucs
#modified CS residuals
mcs.kidney<-allresid.kidney$mcs
#Martingale residuals
martg.kidney<-allresid.kidney$martg
#Deviance residuals
dev.kidney<-allresid.kidney$dev

########################### Graphical Diagnosis##################
plot.zresid(Zresid.kidney)
plot.zresid.fitted.value(Zresid.kidney,
                         fitted.values=fit_kidney$linear.predictors,
                         xlab="Linear Predictor",
                         cenered=kidney$status)
qqnorm.zresid (Zresid.kidney)
boxplot.zresid(Zresid.kidney,fitted.values=kidney$age)
boxplot.zresid(Zresid.kidney,fitted.values=kidney$sex)

########################### other residuals Graphical Diagnosis ##################
##unmodified CS residuals
km.ln.kidney <- survfit(Surv(ucs.kidney, kidney$status)~1,type='fleming')
id.ln.kidney<-order(ucs.kidney)

plot(km.ln.kidney, fun="cumhaz", xlab=("Cox-Snell Residuals"),
     ylab=("Cumulative Hazard Function"),
     main="Log-normal, Cum. Hazard of CS Residuals",
     ylim= c(0,4),xlim=c(0,4))
abline(0, 1, col="red", lty=2)
points(km.ln.kidney$time, -log(km.ln.kidney$surv),
       col=c("blue","darkolivegreen4")[kidney$status[id.ln.kidney]+1],
       pch=c(3,2)[kidney$status[id.ln.kidney]+1] )
legend(x = "topleft",
       legend = c("Uncensored", "Censored"), col=c("darkolivegreen4","blue"),
       pch=c(2,3),cex=1,xpd = TRUE,bty="L")


#Martingale residuals
plot(kidney$age,martg.kidney,ylab="Martingale",
     xlab="age",
     main="Martingale Residual Plot",
     col=c("blue","darkolivegreen4")[kidney$status+1],
     pch=c(3,2)[kidney$status+1])
#abline(h=c(3,-3),col="grey")
lines(lowess(martg.kidney~ kidney$age),
      col = "red",lwd = 3)
legend(x = "bottomright",
       legend = c("Uncensored", "Censored"),
       col=c("darkolivegreen4","blue"),
       pch=c(2,3),cex=1,xpd = TRUE,bty="l")

#Deviance residuals
plot(kidney$age ,dev.kidney,ylab="Deviance",
     xlab="size",
     main="Deviance Residual Plot",
     col=c("blue","darkolivegreen4")[kidney$status+1],
     pch=c(3,2)[kidney$status+1])
#abline(h=c(3,-3),col="grey")
lines(lowess(dev.kidney~ kidney$age),
      col = "red",lwd = 3)
legend(x = "bottomright",
       legend = c("Uncensored", "Censored"),
       col=c("darkolivegreen4","blue"),
       pch=c(2,3),cex=1,xpd = TRUE,bty="l")

########################### Statistical test Diagnosis##################
sw.test.zresid(Zresid.kidney)
sf.test.zresid(Zresid.kidney)
gof.censored.zresidual(censored.Zresidual=censored.Zresid.kidney,censored=kidney$status)

aov.test.zresid(Zresid.kidney,fitted.values=fit_kidney$linear.predictors, k.anova=10)
bartlett.test.zresid(Zresid.kidney,fitted.values=fit_kidney$linear.predictors, k.bl=10)

aov.test.zresid(Zresid.kidney,fitted.values=kidney$age, k.anova=10)
bartlett.test.zresid(Zresid.kidney,fitted.values=kidney$age, k.bl=10)

aov.test.zresid(Zresid.kidney,fitted.values=kidney$sex, k.anova=10)
bartlett.test.zresid(Zresid.kidney,fitted.values=kidney$sex, k.bl=10)


