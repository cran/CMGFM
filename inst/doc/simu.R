## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  library(CMGFM)

## ----eval = FALSE-------------------------------------------------------------
#  pveclist <- list('gaussian'=c(50, 150),'poisson'=c(50, 150),
#                   'binomial'=c(100,60))
#  q <- 6
#  sigmavec <- rep(1,3)
#  pvec <- unlist(pveclist)
#  methodNames <- c("CMGFM", "GFM", "MRRR", "LFR" , "GPCA", "COAP")
#  
#  datlist <- gendata_cmgfm(pveclist = pveclist, seed = 1, n = 300,d = 3,
#                             q = q, rho = rep(1,length(pveclist)), rho_z=0.2,
#                             sigmavec=sigmavec, sigma_eps=1)
#  str(datlist)
#  XList <- datlist$XList
#  Z <- datlist$Z
#  numvarmat <- datlist$numvarmat
#  Uplist <- list()
#  k <- 1
#  for(id in 1:length(datlist$Uplist)){
#      for(im in 1:length(datlist$Uplist[[id]])){
#        Uplist[[k]] <- datlist$Uplist[[id]][[im]]
#        k <- k + 1
#      }
#  }
#  types <- datlist$types

## ----eval = FALSE-------------------------------------------------------------
#  system.time({
#    tic <- proc.time()
#    rlist <- CMGFM(XList, Z, types=types, q=q, numvarmat=numvarmat)
#    toc <- proc.time()
#    time_smgfm <- toc[3] - tic[3]
#  })

## ----eval = FALSE-------------------------------------------------------------
#  library(ggplot2)
#  dat_iter <- data.frame(iter=1:length(rlist$ELBO_seq), ELBO=rlist$ELBO_seq)
#  ggplot(data=dat_iter, aes(x=iter, y=ELBO)) + geom_line() + geom_point() + theme_bw(base_size = 20)
#  

## ----eval = FALSE-------------------------------------------------------------
#  library(GFM)
#  hUplist <- lapply(seq_along(rlist$Bf), function(m) cbind(rlist$muf[[m]], rlist$Bf[[m]]))
#  metricList <- list()
#  metricList$CMGFM <- list()
#  meanTr <- function(hBlist, Blist,  type='trace_statistic'){
#    ###It is noted that the trace statistics is not symmetric, the true value must be in last
#    trvec <- sapply(1:length(Blist), function(j) measurefun(hBlist[[j]], Blist[[j]], type = type))
#    return(mean(trvec))
#  }
#  normvec <- function(x) sqrt(sum(x^2/ length(x)))
#  metricList$CMGFM$Tr_F <- measurefun(rlist$M, datlist$F0)
#  metricList$CMGFM$Tr_Upsilon <- meanTr(hUplist, Uplist)
#  metricList$CMGFM$Tr_U <- measurefun(Reduce(cbind,rlist$Xif), Reduce(cbind, datlist$U0))
#  metricList$CMGFM$beta_norm <- normvec(as.vector(Reduce(cbind,rlist$betaf)- Reduce(cbind,datlist$beta0List)))
#  metricList$CMGFM$Time <- rlist$time_use

## ----eval = FALSE-------------------------------------------------------------
#  metricList$GFM <- list()
#  tic <- proc.time()
#  res_gfm <- gfm(XList, types=types, q=q)
#  toc <- proc.time()
#  time_gfm <- toc[3] - tic[3]
#  mat2list <- function(B, pvec, by_row=TRUE){
#    Blist <- list()
#    pcum = c(0, cumsum(pvec))
#    for(i in 1:length(pvec)){
#      if(by_row){
#        Blist[[i]] <- B[(pcum[i]+1):pcum[i+1],]
#      }else{
#        Blist[[i]] <- B[, (pcum[i]+1):pcum[i+1]]
#      }
#    }
#    return(Blist)
#  }
#  metricList$GFM$Tr_F <- measurefun(res_gfm$hH, datlist$F0)
#  metricList$GFM$Tr_Upsilon <- meanTr(mat2list(cbind(res_gfm$hmu,res_gfm$hB), pvec), Uplist)
#  metricList$GFM$Tr_U <- NA
#  metricList$GFM$beta_norm <- NA
#  metricList$GFM$Time <- time_gfm

## ----eval = FALSE-------------------------------------------------------------
#  mrrr_run <- function(Y, Z, numvarmat, rank0,family=list(poisson()),
#                       familygroup,  epsilon = 1e-4, sv.tol = 1e-2,
#                       lambdaSVD=0.1, maxIter = 2000, trace=TRUE, trunc=500){
#    # epsilon = 1e-4; sv.tol = 1e-2; maxIter = 30; trace=TRUE,lambdaSVD=0.1
#    Diag <- function(vec){
#    q <- length(vec)
#    if(q > 1){
#      y <- diag(vec)
#    }else{
#      y <- matrix(vec, 1,1)
#    }
#    return(y)
#  }
#    require(rrpack)
#    q <- rank0
#    n <- nrow(Y); p <- ncol(Y)
#    X <- cbind(cbind(1, Z),  diag(n))
#    d <- ncol(Z)
#    ## Trunction
#    Y[Y>trunc] <- trunc
#    Y[Y< -trunc] <- -trunc
#  
#    tic <- proc.time()
#    pvec <- as.vector(t(numvarmat))
#    pvec <- pvec[pvec>0]
#    pcums <- cumsum(pvec)
#    idxlist <- list()
#    idxlist[[1]] <- 1:pcums[1]
#    if(length(pvec)>1){
#      for(i in 2:length(pvec)){
#        idxlist[[i]] <- (pcums[i-1]+1):pcums[i]
#      }
#    }
#  
#    svdX0d1 <- svd(X)$d[1]
#    init1 = list(kappaC0 = svdX0d1 * 5)
#    offset = NULL
#    control = list(epsilon = epsilon, sv.tol = sv.tol, maxit = maxIter,
#                   trace = trace, gammaC0 = 1.1, plot.cv = TRUE,
#                   conv.obj = TRUE)
#    res_mrrr <- mrrr(Y=Y, X=X[,-1], family = family, familygroup = familygroup,
#                     penstr = list(penaltySVD = "rankCon", lambdaSVD = lambdaSVD),
#                     control = control, init = init1, maxrank = rank0+d) #
#  
#    hmu <- res_mrrr$coef[1,]
#    hbeta <- t(res_mrrr$coef[2:(d+1),])
#    hTheta <- res_mrrr$coef[-c(1:(d+1)),]
#  
#    hbeta_rf <- NULL
#    for(i in seq_along(pvec)){
#      hbeta_rf <- cbind(hbeta_rf, colMeans(hbeta[idxlist[[i]],]))
#    }
#  
#    # Matrix::rankMatrix(hTheta)
#    svd_Theta <- svd(hTheta, nu=q, nv=q)
#    hH <- svd_Theta$u
#    hB <- svd_Theta$v %*% Diag(svd_Theta$d[1:q])
#    toc <- proc.time()
#    time_mrrr <- toc[3] - tic[3]
#  
#    return(list(hH=hH, hB=hB, hmu= hmu, beta=hbeta_rf, time_use=time_mrrr))
#  }
#  
#  family_use <- list(gaussian(), poisson(), binomial())
#  familygroup <- sapply(seq_along(datlist$XList), function(j) rep(j, ncol(datlist$XList[[j]])))
#  res_mrrr <- mrrr_run(Reduce(cbind, datlist$XList),  Z = Z, numvarmat, rank0=q, family=family_use,
#                         familygroup = unlist(familygroup), maxIter=100) #
#  
#  
#  
#  

## ----eval = FALSE-------------------------------------------------------------
#  metricList$MRRR <- list()
#  metricList$MRRR$Tr_F <- measurefun(res_mrrr$hH, datlist$F0)
#  metricList$MRRR$Tr_Upsilon <-meanTr(mat2list(cbind(res_mrrr$hmu,res_mrrr$hB), pvec), Uplist)
#  metricList$MRRR$Tr_U <-NA
#  metricList$MRRR$beta_norm <- normvec(as.vector(res_mrrr$beta- Reduce(cbind,datlist$beta0List)))
#  metricList$MRRR$Time <- res_mrrr$time_use

## ----eval = FALSE-------------------------------------------------------------
#  list2vec <- function(xlist){
#    nn <- length(xlist)
#    me <- rep(NA, nn)
#    idx_noNA <- which(sapply(xlist, function(x) !is.null(x)))
#    for(r in idx_noNA) me[r] <- xlist[[r]]
#    return(me)
#  }
#  
#  dat_metric <- data.frame(Tr_F = sapply(metricList, function(x) x$Tr_F),
#                           Tr_Upsilon = sapply(metricList, function(x) x$Tr_Upsilon),
#                           Tr_U = sapply(metricList, function(x) x$Tr_U),
#                           beta_norm =sapply(metricList, function(x) x$beta_norm),
#                           Time = sapply(metricList, function(x) x$Time),
#                           Method = names(metricList))
#  dat_metric

## ----eval = FALSE, fig.width=9, fig.height=6----------------------------------
#  library(cowplot)
#  p1 <- ggplot(data=subset(dat_metric, !is.na(Tr_F)), aes(x= Method, y=Tr_F, fill=Method)) + geom_bar(stat="identity") + xlab(NULL) + scale_x_discrete(breaks=NULL) + theme_bw(base_size = 16)
#  p2 <- ggplot(data=subset(dat_metric, !is.na(Tr_Upsilon)), aes(x= Method, y=Tr_Upsilon, fill=Method)) + geom_bar(stat="identity") + xlab(NULL) + scale_x_discrete(breaks=NULL)+ theme_bw(base_size = 16)
#  p3 <- ggplot(data=subset(dat_metric, !is.na(beta_norm)), aes(x= Method, y=beta_norm, fill=Method)) + geom_bar(stat="identity") + xlab(NULL) + scale_x_discrete(breaks=NULL)+ theme_bw(base_size = 16)
#  p4 <- ggplot(data=subset(dat_metric, !is.na(Time)), aes(x= Method, y=Time, fill=Method)) + geom_bar(stat="identity") + xlab(NULL) + scale_x_discrete(breaks=NULL)+ theme_bw(base_size = 16)
#  plot_grid(p1,p2,p3, p4, nrow=2, ncol=2)

## ----eval = FALSE-------------------------------------------------------------
#  
#  hq <- MSVR(XList, Z, types=types, numvarmat, q_max=20)
#  
#  print(c(q_true=q, q_est=hq))
#  

## -----------------------------------------------------------------------------
sessionInfo()

