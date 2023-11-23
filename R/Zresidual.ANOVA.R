test.nl.aov <- function(Zresidual, fitted.value, k.anova=10)
{
  if(is.factor(fitted.value)){
    fitted.value<-as.numeric(fitted.value)-1
    lpred.bin <- fitted.value
    anova(lm(Zresidual ~ lpred.bin))$`Pr(>F)`[1]
  }

 if(!is.factor(fitted.value)){
   lpred.bin <- cut(fitted.value, k.anova)
   less2_factor<-which(tapply(lpred.bin,lpred.bin,length)<= 2)
   if(rlang::is_empty(names(less2_factor))){
     anova(lm(Zresidual ~ lpred.bin))$`Pr(>F)`[1]
   }else{
     list_less2_factor<-list()
     for(j in 1:length(less2_factor)){
       list_less2_factor[[j]]<-which(lpred.bin==names(less2_factor[j]))
     }
     vector_less2_factor<-unlist(list_less2_factor, use.names = FALSE)
     new.lpred.bin<- lpred.bin[-vector_less2_factor]
     new.Zresidual<-Zresidual[-vector_less2_factor]
     anova(lm(new.Zresidual ~ new.lpred.bin))$`Pr(>F)`[1]
   }

 }

}

