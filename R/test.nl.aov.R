#' A function to calculate ANOVA
#'
#' @param Zresidual A z-residual.
#' @param fitted.value Fitted values
#' @param k.anova Number of bins if applicable

test.nl.aov <- function(Zresidual, fitted.value, k.anova=10)
{
  unique.vals <- unique(fitted.value)
  n.unique <- length(unique.vals)

  if (is.factor(fitted.value) || n.unique <= k.anova) {
    # treat as categorical / discrete variable
    if (!is.factor(fitted.value)) {
      fitted.value <- factor(fitted.value)
    }
    lpred.bin <- as.numeric(fitted.value) - 1
    anova(lm(Zresidual ~ lpred.bin))$`Pr(>F)`[1]
  } else if(!is.factor(fitted.value)){
    lpred.bin <- droplevels(cut(fitted.value, k.anova))
    less2_factor<-which(tapply(lpred.bin,lpred.bin,length)<= 2)
    is.bins2 <- (nlevels(lpred.bin) - length(less2_factor))>2

    if(!is.bins2) {
    fitted.value <- log(fitted.value)
    message("Contrasts can be applied only to factors with 2 or more levels. Fitted values converted to log.")
    lpred.bin <- droplevels(cut(fitted.value, k.anova))
    less2_factor<-which(tapply(lpred.bin,lpred.bin,length)<= 2)
    }

    if(rlang::is_empty(names(less2_factor))){
      anova(lm(Zresidual ~ lpred.bin))$`Pr(>F)`[1]
    }else {
      list_less2_factor<-list()
      for(j in 1:length(less2_factor)){
        list_less2_factor[[j]]<-which(lpred.bin==names(less2_factor[j]))
      }
      vector_less2_factor<-unlist(list_less2_factor, use.names = FALSE)
      new.lpred.bin<- lpred.bin[-vector_less2_factor]
      new.Zresidual<-Zresidual[-vector_less2_factor]
      anova(lm(new.Zresidual ~ new.lpred.bin))$`Pr(>F)`[1]
    }
  } else {
    warning("Problem in ANOVA.")
  }
}
