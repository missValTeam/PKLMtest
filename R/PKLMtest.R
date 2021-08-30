#' PKLMtest: compute a p-value for testing MCAR
#'
#' @param X a matrix containing missing values, the data.
#' @param num.proj an integer specifying the number of projections to consider for the score.
#' @param num.trees.per.proj an integer, the number of trees per projection.
#' @param nrep an integer, the number of permutations.
#' @param min.node.size the minimum number of nodes in a tree.
#' @param size.resp.set an integer (>= 2), maximum number of classes allowed to be compared in each projection.
#' @param compute.partial.pvals boolean, indicate if partial p-values shopuld be computed as well.
#' @param ... additional parameters.
#' @import ranger
#' @import parallel
#'
#' @examples
#' n <- 500
#' X <- cbind(rnorm(n),rnorm(n))
#' X.NA <- X
#' X.NA[,1] <- ifelse(stats::runif(n)<=0.2, NA, X[,1])
#'
#' pval <- PKLMtest(X)
#'
#' @return a numeric value, the p-value(s) for the MCAR test.
#'
#' @export
PKLMtest <- function(X,
                    num.proj = 300,
                    num.trees.per.proj = 10,
                    nrep = 500,
                    min.node.size=10,
                    size.resp.set = 2,
                    compute.partial.pvals = FALSE,
                    ...) {


  ## Combining p-values in Nicolais way without refitting the forest each time

  size.resp.set <- size.resp.set - 1

  # checks
  if (!is.matrix(X)) {
    stop("X should be a matrix.")
  }

  # get the mask for the data provided in X
  M <- is.na(X)

  Mperm<-lapply(1:nrep, function(j) M[sample(nrow(M)),])

  # retrieve the global unique patterns (not really needed)
  #unique.pat <- unique(M)
  #ids.patterns <- sapply(1:nrow(unique.pat),
  #                       function(i) which(apply(M,
  #                                               1,
  #                                               function(m) identical(m, unique.pat[i,]))))

  # list of results
  list.obs <- list()
  list.null <- list()
  list.var.cov <- list()
  list.var.resp <- list()

  # main loop
  while (length(list.obs) != num.proj) {


    # print(length(list.obs))

    # sample a projection for the response
    new.var <- c()
    num.classes <- 1
    while (length(new.var)==0 ||  num.classes <= 1 ) {

      # sample a projection for the covariates
      size.cov <- sample(1:ncol(X), size = 1, replace = FALSE)
      var.cov <- sample(1:ncol(X), size = size.cov, replace = FALSE)

      new.var <- setdiff(1:ncol(X), var.cov)

      if (length(new.var)>1) {

        if (size.resp.set=="random"){
          size.resp <- sample(1:length(new.var), size = 1, replace = FALSE)
        }else{
          size.resp<- sample(1:min(length(new.var), size.resp.set), size=1, replace=FALSE)
        }

        var.resp <- sample(new.var, size = size.resp, replace = FALSE)
      } else if (length(new.var)==1) {
        var.resp <- new.var
      } else if (length(new.var)==0){
        var.resp <- var.cov
      }



      # which one are complete there
      ids.keep <- which(apply(X[,var.cov,drop=FALSE],
                              1, function(x) !any(is.na(x))))

      # define the new response
      M.resp <- M[ids.keep,var.resp,drop=F]
      unique.pat.resp <- unique(M.resp)
      num.classes<- nrow(unique(M.resp))

    }

    list.var.cov[[length(list.var.cov)+1]] <- var.cov
    list.var.resp[[length(list.var.resp)+1]] <- var.resp


    ids.patterns.resp <- sapply(1:nrow(unique.pat.resp),
                                function(i) which(apply(M.resp,
                                                        1,
                                                        function(m) identical(m,
                                                                              unique.pat.resp[i,]))), simplify = F)

    # get the data for the projection
    X.proj <- X[ids.keep,var.cov,drop=F]

    # create the class for the projection
    class.resp <- rep(NA, length(ids.keep))

    for (i in 1:length(ids.patterns.resp)) {
      class.resp[ids.patterns.resp[[i]]] <- i
    }


    class.resp.perm<-lapply(1:nrep, FUN=function(j){


      M.respj <- Mperm[[j]][ids.keep,var.resp,drop=F]
      unique.pat.respj <- unique(M.respj)
      num.classesj<- nrow(unique(M.respj))


      ids.patterns.resp <- sapply(1:nrow(unique.pat.respj),
                                  function(i) which(apply(M.respj,
                                                          1,
                                                          function(m) identical(m,
                                                                                unique.pat.respj[i,]))), simplify = F)

      # create the class for the projection
      class.respj <- rep(NA, length(ids.keep))

      for (i in 1:length(ids.patterns.resp)) {
        class.respj[ids.patterns.resp[[i]]] <- i
      }

      return(class.respj)


    } )

    # permute under the null
    #class.resp.perm <- lapply(1:nrep, FUN = function(x) sample(class.resp))

    # if more than one class
    #if (length(unique(class.resp))!=1) {

    d <- data.frame(y = factor(class.resp),
                    X = X.proj)

    st <- tryCatch({ranger::ranger(y~., data = d,
                                   num.trees = num.trees.per.proj,
                                   classification = TRUE,
                                   min.node.size = min.node.size,
                                   probability = TRUE)
    }, error = function(e) NA)

    # st<-ranger::ranger(y~., data = d,
    #                num.trees = num.trees.per.proj,
    #                classification = TRUE,
    #                min.node.size = min.node.size,
    #                probability = TRUE)


    if(!any(is.na(st))) {

      # generate the observed stat
      obs <- genU(st, lab = class.resp)
      # generate the stat under H0
      null <- sapply(class.resp.perm, function(p) genU(st, lab = p))
      list.obs[[length(list.obs)+1]] <- obs
      list.null[[length(list.null)+1]] <- null

    }

  }
  print(paste0("number of NAs in given projection/permutation: ", mean(is.na(unlist(list.null)))))


  #lapply(1:length(list.null), function(i) (sum(list.obs[[i]] > list.null[[i]])+1)/(length(list.null[[i]])+1) )


  #print(dim(isIn))
  #print(dim(isIn))
  stat.perm.raw <- Reduce(list.null,f = rbind)



  ##### partial p-values #########

  isIn <- sapply(1:ncol(X), function(i) sapply(list.var.resp, function(l) !(i %in% l))) # list.var.resp
  isIn <- cbind(rep(TRUE, num.proj ), isIn)
  if (!compute.partial.pvals) {
    isIn <- isIn[,1,drop=FALSE]
  }

  partial.perm <- apply(isIn ,2, function(x) {
    ids <- which(x)
    #apply(stat.perm.raw[ids,,drop=FALSE], 2, function(x) mean(x, na.rm=T ))  # mean(x[is.finite(x)])
    ##


    colMeans(stat.perm.raw[ids,,drop=FALSE], na.rm=T)


  })
  if (ncol(isIn)==1) {
    partial.perm <- matrix(partial.perm, nrow=1)
  }

  partial.obs <- apply(isIn ,2, function(x) {
    ids <- which(x)

    if (kind=="mean"){
      ## mean
      mean(unlist(list.obs[ids]),
           na.rm=T)
    } else if (kind=="sup"){
      ## max
      max(unlist(list.obs[ids])[is.finite(unlist(list.obs[ids]))],
          na.rm=T)
    }


  })

  if (compute.partial.pvals){
    pvals <- sapply(1:(ncol(X)+1), function(i) (sum(partial.perm[,i] >= partial.obs[i])+1)/(length(partial.perm[,i])+1))
    names(pvals) <- c("all", paste0("-",1:ncol(X)))
  }else{
    pvals<-(sum(partial.perm >= partial.obs)+1)/(length(partial.perm)+1)
    names(pvals) <- "all"
  }

  ##############################



  # old
  ####pval <- (sum(stat.null>obs)+1)/(length(stat.null)+1)

  return(pvals)
  # isIn = isIn))

}
