#' Truncation of probability
#' @param p a numeric value between 0 and 1 to be truncated
truncProb <- function(p) {
  return(pmin(pmax(p, 10^{-9}), 1-10^{-9}))
}

#' Generate the test statistic
#' @param st a ranger forest object.
#' @param lab an integer value containing the class labels
genU <- function(st, lab) {

  # OOB preds
  preds <- st$predictions

  # number of actual classes
  ui <- rep(NA, ncol(preds))

  for (i in 1:length(ui)) {

    ind <- colnames(preds)[i]
    p <- preds[lab==ind,i]
    p <- truncProb(p)
    lab2 <- lab[is.finite(p)]
    p <- p[is.finite(p)]

    p.a <- preds[lab!=ind,i]
    p.a <- truncProb(p.a)
    lab2.a <- lab[is.finite(p.a)]
    p.a <- p.a[is.finite(p.a)]

    #class.freq <- mean(lab==ind)
    #ui[i] <- mean(log(p/(1-p)) - log(class.freq/(1-class.freq )))

    ui[i] <- mean(log(p/(1-p))) - mean(log(p.a/(1-p.a)))
  }

  return(mean(ui, na.rm=TRUE))
  #return(sum(ui, na.rm=TRUE))
}
