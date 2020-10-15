#### MFM for matrix normal
#################################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("matrixNormal","mniw","MCMCpack","mvtnorm","sna","MixMatrix")
ipak(packages)
#################################
getDahl <- function(MFMfit, burn)
{
  ################################################################
  
  ## Input: MFMfit = the result from CDMFM_new ##
  ##        burn = the number of burn-in iterations ##
  
  ## Output: 
  ##         zout = estimated clustering configuration, a n by 1 vector##
  ##         Qout = estimated probability matrix, a k by k matrix ##
  
  #################################################################
  iters <- MFMfit$Iterates[-(1:burn)]
  n <- length(iters[[1]][[1]])
  niters <- length(iters)
  membershipMatrices <- lapply(iters, function(x){
    clusterAssign <- x[[1]]
    outer(clusterAssign, clusterAssign, FUN = "==")
  })
  membershipAverage <- Reduce("+", membershipMatrices)/niters
  SqError <- sapply(membershipMatrices, function(x, av) sum((x - av)^2),
                    av = membershipAverage)
  DahlIndex <- which.min(SqError)
  DahlAns <- iters[[DahlIndex]]
  attr(DahlAns, "iterIndex") <- burn + DahlIndex
  attr(DahlAns, "burnin") <- burn
  DahlAns
}

## function for Collapsed sampler
CDMFM_new1 <- function(data, niterations, alpha, beta, psi, rho, GAMMA, 
                       M0, Sigma0, Omega0, initNClusters, VN, MLE.initial = FALSE)
{
  ## Model: Y_{i}|z,theta \sim MN(theta_{z_i}) ##
  ##        M_{r} \sim MN(M0, Sigma0, Omega0), r = 1,...,k ##
  ##        U ~ IW(2alpha, (2beta)^{-1})
  ##        V ~ IW(2psi, (2rho)^{-1})
  ##        P(z_i = j) = \pi_j, j = 1,...,k ##
  ##        \pi \sim Dirichlet_k(GAMMA,...,GAMMA) ##
  ##        k-1 \sim possion(1) ##
  
  # define a helper function to calculate the marginal probability Y|U, V
  dMNorm_marginal <- function(Y, U, V, M0, Sigma0, Omega0, log=FALSE){
    p <- nrow(Y)
    q <- ncol(Y)
    solve_V = solve(V)
    solve_U = solve(U)
    solve_Sigma0 = solve(Sigma0)
    solve_Omega0 = solve(Omega0)
    Sigma_tilde <- solve( solve_V %x% solve_U + solve_Omega0 %x% solve_Sigma0 )
    Mu_tilde <- Sigma_tilde %*% ( (solve_V %x% solve_U ) %*% c(Y) + (solve_Omega0 %x% solve_Sigma0) %*% c(M0) )
    log_numerator <- (p*q/2) * log(det(Sigma_tilde))  + 1/2 * t(Mu_tilde) %*% solve(Sigma_tilde) %*% Mu_tilde - 1/2 * ( t(c(Y)) %*% (solve_V %x% solve_U) %*% c(Y) + t(c(M0)) %*% (solve_Omega0 %x% solve_Sigma0) %*% c(M0) ) 
    log_denominator <- (p*q/2) * log(2*pi) + p/2 * log(det(V)) + q/2 * log(det(U)) + p/2 * log(det(Omega0)) + q/2 * log(det(Sigma0))
    if(log == TRUE){
      return(log_numerator-log_denominator)
    }
    else{
      return( exp(log_numerator - log_denominator) )
    }
  }
  
  ################################################################
  
  ## Input: data = the vector of responses ##
  ##        niterations = the total number of iterations in MFM-SBM ##
  ##        alpha, beta = hyperparameters (shape, rate) for the prior on elements in lambda in Gamma distribution ##
  ##        GAMMA = the parameter in Dirichlet distribution that controls the relative size of clusters ##
  ##        LAMBDA = the parameter for Poisson distrition ##
  ##        initNClusters = the initial number of clusters ##
  
  
  ## Output: 
  ##         zout = clustering configuration, a n^2 by 1 vector##
  ##         phiout = possion parameters, a k by 1 vector ##
  
  #################################################################
  # n = length(data)
  n = dim(data)[3]
  #precomputation for prespecified coefficient VN
  # lambda <- LAMBDA
  #gamma <- GAMMA
  N=n ## n is the number of oberservations
  # dimension of the each data observation
  p <- nrow(data[,,1])
  q <- ncol(data[,,1])
  solve_Omega0 <- solve(Omega0)
  solve_Sigma0 <- solve(Sigma0)
  # VN<-0
  # tmax = n+10
  # for (t in 1:tmax)
  # {
  #   r = log(0)
  #   for (k in t:500)
  #   {
  #     b = sum(log((k-t+1):k))-sum(log((k*gamma):(k*gamma+N-1))) + dpois(k-1, lambda, log = TRUE)
  #     m = max(b,r)
  #     r = log(exp(r-m) + exp(b-m)) + m
  #   }
  #   VN[t] = r
  # }
  #initialization of clustering configuration
  clusterAssign <- c(sample(1:initNClusters, size = initNClusters, replace = FALSE),
                     sample(1:initNClusters, size = n-initNClusters, replace = TRUE))
  # print(clusterAssign)
  #
  # phi<-rgamma(initNClusters, shape = alpha, rate = beta)
  # need to intialize M, U, V
  # M <- array(0, dim=c(p,q,initNClusters))
  # M <- array(mniw::rMNorm(initNClusters, 
  #                         Lambda = M0, 
  #            SigmaR = Sigma0, SigmaC = Omega0), 
  #            dim=c(p,q,initNClusters)) # three dimensions
  if(MLE.initial){
    MLE_data <- MixMatrix::MLmatrixnorm(data)
    M <- array(0, dim=c(p,q,initNClusters))
    for(r in 1:initNClusters){
      M[,,r] <- MLE_data$mean
    }
    # the last dimension correspons to the sample size, initNClusters
    U <- (sum(diag(MLE_data$V))/nrow(MLE_data$V))*MLE_data$U
    V <- nrow(MLE_data$V)*((MLE_data$V)/sum(diag(MLE_data$V)))
  }
  else{
    M <- array(mniw::rMNorm(initNClusters,
                            Lambda = M0,
               SigmaR = Sigma0, SigmaC = Omega0),
               dim=c(p,q,initNClusters)) # three dimensions
    U <- MCMCpack::riwish(2*alpha, 2*beta)
    V <- MCMCpack::riwish(2*psi, 2*rho) 
  }
  # U <- MCMCpack::riwish(2*alpha, 2*beta)
  # V <- MCMCpack::riwish(2*psi, 2*rho) 
  solve_V = solve(V)
  solve_U = solve(U)
  # 
  # print(M)
  # MCMCpack::InvWishart
  History <- vector("list", niterations)
  
  ##start Gibb's sampling
  for (iter in 1:niterations)
  {
    print(iter)
    ## update z ##
    clusterSizes = table(as.factor(clusterAssign))
    nClusters = length(clusterSizes)
    
    # preparation to avoid extra calculation
    solve_Sigma_tilde <- solve_V %x% solve_U + solve_Omega0 %x% solve_Sigma0
    Sigma_tilde <- solve( solve_Sigma_tilde )
    for (i in 1:n)
    { #determine whether ith component is a singleton 
      cur.cluster.i = clusterAssign[i]
      if (clusterSizes[clusterAssign[i]] > 1){
        # not a singleton, have |C|+1 choices
        c.counts.noi = clusterSizes  #c.counts.noi corresponds to |C|
        c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1
        #finding the probs for sampling process
        # browser()
        if(length(dim(M))==2){
          M = array(M, dim=c(p,q,1))
        }
        clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          (GAMMA+c.counts.noi[x])*mniw::dMNorm(data[,,i], 
                                               Lambda = M[,,x], 
                                               SigmaR = U, 
                                               SigmaC = V)  # dpois(data[i],phi[x])
        })
        clusterAssign_1 = clusterAssign
        clusterAssign_1[i] = nClusters+1
        
        Mu_tilde <- Sigma_tilde %*% ( (solve_V %x% solve_U ) %*% c(data[,,i]) + (solve_Omega0 %x% solve_Sigma0) %*% c(M0) )
        log_numerator <- (p*q/2) * log(det(Sigma_tilde))  + 1/2 * t(Mu_tilde) %*% solve_Sigma_tilde %*% Mu_tilde - 1/2 * ( t(c(data[,,i])) %*% (solve_V %x% solve_U) %*% c(data[,,i]) + t(c(M0)) %*% (solve_Omega0 %x% solve_Sigma0) %*% c(M0) ) 
        log_denominator <- (p*q/2) * log(2*pi) + p/2 * log(det(V)) + q/2 * log(det(U)) + p/2 * log(det(Omega0)) + q/2 * log(det(Sigma0))
        
        clusterProbs[nClusters+1]<-GAMMA*exp(VN[nClusters+1]-VN[nClusters])*exp(log_numerator - log_denominator)
        
        # clusterProbs[nClusters+1]<- GAMMA*exp(VN[nClusters+1]-VN[nClusters])*dMNorm_marginal(Y=data[,,i], 
        #                                                                                      U=U, 
        #                                                                                      V=V, 
        #                                                                                      M0=M0, 
        #                                                                                      Sigma0=Sigma0, 
        #                                                                                      Omega0=Omega0)
        # print(clusterProbs)
        #browser()
        # GAMMA*beta^alpha/gamma(alpha)*gamma(alpha+data[i])/(beta+1)^(data[i]+alpha)/factorial(data[i])*exp(VN[nClusters+1]-VN[nClusters])
        
        #choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1,
                            prob = clusterProbs)
        clusterAssign[i] <- cluster.i
        
        if (cluster.i > nClusters)
        {
          # phinew = rep(0,nClusters+1)
          # phinew[1:nClusters] = phi
          # phinew[nClusters+1] = rgamma(1, shape = alpha, rate = beta)
          # phi = phinew
          Mnew = array(0, dim=c(p,q,nClusters+1))
          Mnew[,,1:nClusters] = M
          Mnew[,,nClusters+1] = mniw::rMNorm(1, Lambda = M0, SigmaR = Sigma0, SigmaC = Omega0) # simulate from prior
          M = Mnew
          # dimension of M
          if(length(dim(M))==2){
            M = array(M, dim=c(p,q,1))
          }
          clusterSizes <- table(as.factor(clusterAssign)) # sorts according to labels
          nClusters <- length(clusterSizes)} else
          {# phi = phi
           M = M
          clusterSizes <- table(as.factor(clusterAssign))
          nClusters <- length(clusterSizes)}
      } else {
        # force the dimension to be 3D
        if(length(dim(M))==2){
          M = array(M, dim=c(p,q,1))
        }
        # a singleton, have |C| choices
        c.counts.noi = clusterSizes
        c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1 - GAMMA# can offset the gamma adding later
        #finding the probs for sampling process
        clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          (GAMMA+c.counts.noi[x])*mniw::dMNorm(data[,,i], Lambda = M[,,x], SigmaR = U, SigmaC = V)
        })
        clusterAssign_1 = clusterAssign
        clusterAssign_1[i] = nClusters+1
        
        Mu_tilde <- Sigma_tilde %*% ( (solve_V %x% solve_U ) %*% c(data[,,i]) + (solve_Omega0 %x% solve_Sigma0) %*% c(M0) )
        log_numerator <- (p*q/2) * log(det(Sigma_tilde))  + 1/2 * t(Mu_tilde) %*% solve_Sigma_tilde %*% Mu_tilde - 1/2 * ( t(c(data[,,i])) %*% (solve_V %x% solve_U) %*% c(data[,,i]) + t(c(M0)) %*% (solve_Omega0 %x% solve_Sigma0) %*% c(M0) ) 
        log_denominator <- (p*q/2) * log(2*pi) + p/2 * log(det(V)) + q/2 * log(det(U)) + p/2 * log(det(Omega0)) + q/2 * log(det(Sigma0))
        
        clusterProbs[nClusters+1]<-GAMMA*exp(VN[nClusters+1]-VN[nClusters])*exp(log_numerator - log_denominator)
        # clusterProbs[nClusters+1]<-GAMMA*exp(VN[nClusters+1]-VN[nClusters])*dMNorm_marginal(Y=data[,,i], U=U, 
        #                                                                                     V=V, M0=M0, 
        #                                                                                     Sigma0=Sigma0, 
        #                                                                                     Omega0=Omega0)
        #choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1,
                            prob = clusterProbs)
        # remove the empty cluster
        if (cluster.i > nClusters)
        {      clusterAssign[i] <- cur.cluster.i #put the new cluster in the place of the only singleten one
        clusterSizes <- table(as.factor(clusterAssign)) # sorts according to labels
        } else {      
          clusterAssign[i] <- cluster.i
          clusterAssign <- ifelse(clusterAssign > cur.cluster.i, clusterAssign-1, clusterAssign) # to delete the previous group index
          clusterSizes <- table(as.factor(clusterAssign))
          nClusters <- length(clusterSizes) 
          # phi = phi[-cur.cluster.i]}
          M = M[,,-cur.cluster.i] }
      }
     }
    # end for loop over subjects i
    ## update phi ##
    ## update theta ##
    # M = array(0, dim=c(p,q,nClusters))
    # U = matrix(0, nrow = p, ncol = p)
    # V = matrix(0, nrow = q, ncol = q)
    # Z, matrix n by nClusters
    # Z = matrix(0, nrow = n, ncol = nClusters)
    # for(ii in 1:n){
    #   # print(clusterAssign[ii])
    #   Z[ii,clusterAssign[ii]] <- 1
    # }
    # n_Cluster_size <- apply(Z,2,sum)
    # print(nClusters)
    #print(n_Cluster_size)
    # V <- cov2cor(V) # for identifiability purpose
    # rep(0, nClusters)
    # AA = rep(0,nClusters)
    # NN = rep(0,nClusters)
    n_Cluster_size <- table( factor(clusterAssign, levels=1:max(clusterAssign)) )
    M <- array(M, dim = c(p,q,nClusters))
    for (r in 1:nClusters){
      temp_zY_sum <- matrix(0, nrow=p, ncol=q)
      for(ii in 1:n){
        if(clusterAssign[ii] == r){
          temp_zY_sum <- temp_zY_sum + data[,,ii]
        }
        # temp_zY_sum <- temp_zY_sum + Z[ii,r] * data[,,i]
      }
      # print(temp_zY_sum)
      temp_zeta <- n_Cluster_size[r] * (solve_V %x% solve_U) + (solve_Omega0 %x% solve_Sigma0)
      temp_xi <- c( solve_Sigma0 %*% M0 %*% solve_Omega0 + solve_U %*% temp_zY_sum %*% solve_V )
      # first generate a vector
      solve_temp_zeta <- solve(temp_zeta)
      # print(solve_temp_zeta %*% temp_xi)
      temp_M_vec <- mvtnorm::rmvnorm(1, mean = solve_temp_zeta %*% temp_xi, sigma = solve_temp_zeta)
      # browser()
      M[,,r] <- matrix(temp_M_vec, nrow = p, ncol = q) # convert a vector to matrix by column
      # clusterAssign is a vector, convert it to Z_ij
      # rgamma(1,alpha + sum(data[clusterAssign == r]), beta + sum(clusterAssign == r))
    }
    # update U
    U_scale <- matrix(0, nrow=p, ncol=p) + 2*beta 
    #
    for(ii in 1:n){
      U_scale <- U_scale + 1 * (data[,,ii] - M[,,clusterAssign[ii]]) %*% solve_V %*% t(data[,,ii] - M[,,clusterAssign[ii]])
    }
    U <- MCMCpack::riwish(2*alpha+n*q, U_scale)
    solve_U = solve(U)
    # Update V
    V_scale <- matrix(0, nrow=q, ncol=q) + 2*rho
    for(ii in 1:n){
      V_scale <- V_scale + 1 * t(data[,,ii] - M[,,clusterAssign[ii]]) %*% solve_U %*% (data[,,ii] - M[,,clusterAssign[ii]])
    }
    V <- MCMCpack::riwish(2*psi+n*p, V_scale)
    # fix the trace of V to be q
    V <- q * V / sum(diag(V))
    solve_V = solve(V)
    # History[[iter]] <- list(zout = clusterAssign, phiout = phi)
    History[[iter]] <- list(zout = clusterAssign, Mout = M, Uout = U, Vout = V)
    cat(" iteration:", iter,"\n",clusterAssign,"\n")
  }# for loop over iterations
  
  list(Iterates = History)
}

################
## function for Collapsed sampler
CDMFM_new1_log <- function(data, niterations, alpha, beta, psi, rho, GAMMA, 
                       M0, Sigma0, Omega0, initNClusters, VN, MLE.initial=FALSE)
{
  ## Model: Y_{i}|z,theta \sim MN(theta_{z_i}) ##
  ##        M_{r} \sim MN(M0, Sigma0, Omega0), r = 1,...,k ##
  ##        U ~ IW(2alpha, (2beta)^{-1})
  ##        V ~ IW(2psi, (2rho)^{-1})
  ##        P(z_i = j) = \pi_j, j = 1,...,k ##
  ##        \pi \sim Dirichlet_k(GAMMA,...,GAMMA) ##
  ##        k-1 \sim possion(1) ##
  
  # define a helper function to calculate the marginal probability Y|U, V
  dMNorm_marginal <- function(Y, U, V, M0, Sigma0, Omega0, log=FALSE){
    p <- nrow(Y)
    q <- ncol(Y)
    solve_V = solve(V)
    solve_U = solve(U)
    solve_Sigma0 = solve(Sigma0)
    solve_Omega0 = solve(Omega0)
    Sigma_tilde <- solve( solve_V %x% solve_U + solve_Omega0 %x% solve_Sigma0 )
    Mu_tilde <- Sigma_tilde %*% ( (solve_V %x% solve_U ) %*% c(Y) + (solve_Omega0 %x% solve_Sigma0) %*% c(M0) )
    log_numerator <- (p*q/2) * log(det(Sigma_tilde))  + 1/2 * t(Mu_tilde) %*% solve(Sigma_tilde) %*% Mu_tilde - 1/2 * ( t(c(Y)) %*% (solve_V %x% solve_U) %*% c(Y) + t(c(M0)) %*% (solve_Omega0 %x% solve_Sigma0) %*% c(M0) ) 
    log_denominator <- (p*q/2) * log(2*pi) + p/2 * log(det(V)) + q/2 * log(det(U)) + p/2 * log(det(Omega0)) + q/2 * log(det(Sigma0))
    if(log == TRUE){
      return(log_numerator-log_denominator)
    }
    else{
      return( exp(log_numerator - log_denominator) )
    }
  }
  
  ################################################################
  
  ## Input: data = the vector of responses ##
  ##        niterations = the total number of iterations in MFM-SBM ##
  ##        alpha, beta = hyperparameters (shape, rate) for the prior on elements in lambda in Gamma distribution ##
  ##        GAMMA = the parameter in Dirichlet distribution that controls the relative size of clusters ##
  ##        LAMBDA = the parameter for Poisson distrition ##
  ##        initNClusters = the initial number of clusters ##
  
  
  ## Output: 
  ##         zout = clustering configuration, a n^2 by 1 vector##
  ##         phiout = possion parameters, a k by 1 vector ##
  
  #################################################################
  # n = length(data)
  n = dim(data)[3]
  #precomputation for prespecified coefficient VN
  # lambda <- LAMBDA
  #gamma <- GAMMA
  N=n ## n is the number of oberservations
  # dimension of the each data observation
  p <- nrow(data[,,1])
  q <- ncol(data[,,1])
  solve_Omega0 <- solve(Omega0)
  solve_Sigma0 <- solve(Sigma0)
  # VN<-0
  # tmax = n+10
  # for (t in 1:tmax)
  # {
  #   r = log(0)
  #   for (k in t:500)
  #   {
  #     b = sum(log((k-t+1):k))-sum(log((k*gamma):(k*gamma+N-1))) + dpois(k-1, lambda, log = TRUE)
  #     m = max(b,r)
  #     r = log(exp(r-m) + exp(b-m)) + m
  #   }
  #   VN[t] = r
  # }
  #initialization of clustering configuration
  clusterAssign <- c(sample(1:initNClusters, size = initNClusters, replace = FALSE),
                     sample(1:initNClusters, size = n-initNClusters, replace = TRUE))
  # print(clusterAssign)
  #
  # phi<-rgamma(initNClusters, shape = alpha, rate = beta)
  # need to intialize M, U, V
  if(MLE.initial){
    MLE_data <- MixMatrix::MLmatrixnorm(data)
    M <- array(0, dim=c(p,q,initNClusters))
    for(r in 1:initNClusters){
      M[,,r] <- MLE_data$mean
    }
    # the last dimension correspons to the sample size, initNClusters
    U <- (sum(diag(MLE_data$V))/nrow(MLE_data$V))*MLE_data$U
    V <- nrow(MLE_data$V)*((MLE_data$V)/sum(diag(MLE_data$V)))
  }
  else{
    M <- array(mniw::rMNorm(initNClusters,
                            Lambda = M0,
                            SigmaR = Sigma0, SigmaC = Omega0),
               dim=c(p,q,initNClusters)) # three dimensions
    U <- MCMCpack::riwish(2*alpha, 2*beta)
    V <- MCMCpack::riwish(2*psi, 2*rho) 
  }
  # U <- MCMCpack::riwish(2*alpha, 2*beta)
  # V <- MCMCpack::riwish(2*psi, 2*rho) 
  solve_V = solve(V)
  solve_U = solve(U)
  # print(M)
  # MCMCpack::InvWishart
  History <- vector("list", niterations)
  
  ##start Gibb's sampling
  for (iter in 1:niterations)
  {
    print(iter)
    ## update z ##
    clusterSizes = table(as.factor(clusterAssign))
    nClusters = length(clusterSizes)
    # avoid repeated calculations in marginal likelihood part, pre-compute these common numbers in advance
    solve_Sigma_tilde <- solve_V %x% solve_U + solve_Omega0 %x% solve_Sigma0
    Sigma_tilde <- solve( solve_Sigma_tilde )
    for (i in 1:n)
    { #determine whether ith component is a singleton 
      cur.cluster.i = clusterAssign[i]
      if (clusterSizes[clusterAssign[i]] > 1){
        # not a singleton, have |C|+1 choices
        c.counts.noi = clusterSizes  #c.counts.noi corresponds to |C|
        c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1
        #finding the probs for sampling process
        # browser()
        # clusterProbs = sapply(1:nClusters, function(x) {
        #   clusterAssign_temp = clusterAssign
        #   clusterAssign_temp[i] = x
        #   (GAMMA+c.counts.noi[x])*mniw::dMNorm(data[,,i], 
        #                                        Lambda = M[,,x], 
        #                                        SigmaR = U, 
        #                                        SigmaC = V)  # dpois(data[i],phi[x])
        # })
        ###
        if(length(dim(M))==2){
          M = array(M, dim=c(p,q,1))
        }
        log_clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          log(GAMMA+c.counts.noi[x]) + mniw::dMNorm(data[,,i], 
                                               Lambda = M[,,x], 
                                               SigmaR = U, 
                                               SigmaC = V, log=TRUE)  # dpois(data[i],phi[x])
        })
        clusterAssign_1 = clusterAssign
        clusterAssign_1[i] = nClusters+1
        
        Mu_tilde <- Sigma_tilde %*% ( (solve_V %x% solve_U ) %*% c(data[,,i]) + (solve_Omega0 %x% solve_Sigma0) %*% c(M0) )
        log_numerator <- (p*q/2) * log(det(Sigma_tilde))  + 1/2 * t(Mu_tilde) %*% solve_Sigma_tilde %*% Mu_tilde - 1/2 * ( t(c(data[,,i])) %*% (solve_V %x% solve_U) %*% c(data[,,i]) + t(c(M0)) %*% (solve_Omega0 %x% solve_Sigma0) %*% c(M0) ) 
        log_denominator <- (p*q/2) * log(2*pi) + p/2 * log(det(V)) + q/2 * log(det(U)) + p/2 * log(det(Omega0)) + q/2 * log(det(Sigma0))
        
        log_clusterProbs[nClusters+1] <- log(GAMMA) + (VN[nClusters+1]-VN[nClusters]) + (log_numerator - log_denominator)
        # log_clusterProbs[nClusters+1]<- log(GAMMA) + (VN[nClusters+1]-VN[nClusters]) + dMNorm_marginal(Y=data[,,i], 
        #                                                                                      U=U, 
        #                                                                                      V=V, 
        #                                                                                      M0=M0, 
        #                                                                                      Sigma0=Sigma0, 
        #                                                                                      Omega0=Omega0,
        #                                                                                      log=TRUE)
        #print(log_clusterProbs)
        #browser()
        # if(iter==2){
        #   browser()
        # }
        # GAMMA*beta^alpha/gamma(alpha)*gamma(alpha+data[i])/(beta+1)^(data[i]+alpha)/factorial(data[i])*exp(VN[nClusters+1]-VN[nClusters])
        
        #choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1,
                            prob = exp(log_clusterProbs - sna::logSum(log_clusterProbs) ) )
        clusterAssign[i] <- cluster.i
        
        if (cluster.i > nClusters)
        {
          # phinew = rep(0,nClusters+1)
          # phinew[1:nClusters] = phi
          # phinew[nClusters+1] = rgamma(1, shape = alpha, rate = beta)
          # phi = phinew
          Mnew = array(0, dim=c(p,q,nClusters+1))
          Mnew[,,1:nClusters] = M
          Mnew[,,nClusters+1] = mniw::rMNorm(1, Lambda = M0, SigmaR = Sigma0, SigmaC = Omega0) # simulate from prior
          M = Mnew
          clusterSizes <- table(as.factor(clusterAssign)) # sorts according to labels
          nClusters <- length(clusterSizes)} else
          {# phi = phi
            M = M
            clusterSizes <- table(as.factor(clusterAssign))
            nClusters <- length(clusterSizes)}
      } else {
        if(length(dim(M))==2){
          M = array(M, dim=c(p,q,1))
        }
        # a singleton, have |C| choices
        c.counts.noi = clusterSizes
        c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1 - GAMMA# can offset the gamma adding later
        #finding the probs for sampling process
        log_clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          log(GAMMA+c.counts.noi[x]) + mniw::dMNorm(data[,,i], Lambda = M[,,x], SigmaR = U, SigmaC = V, log=TRUE)
        })
        clusterAssign_1 = clusterAssign
        clusterAssign_1[i] = nClusters+1
        
        Mu_tilde <- Sigma_tilde %*% ( (solve_V %x% solve_U ) %*% c(data[,,i]) + (solve_Omega0 %x% solve_Sigma0) %*% c(M0) )
        log_numerator <- (p*q/2) * log(det(Sigma_tilde))  + 1/2 * t(Mu_tilde) %*% solve_Sigma_tilde %*% Mu_tilde - 1/2 * ( t(c(data[,,i])) %*% (solve_V %x% solve_U) %*% c(data[,,i]) + t(c(M0)) %*% (solve_Omega0 %x% solve_Sigma0) %*% c(M0) ) 
        log_denominator <- (p*q/2) * log(2*pi) + p/2 * log(det(V)) + q/2 * log(det(U)) + p/2 * log(det(Omega0)) + q/2 * log(det(Sigma0))
        
        log_clusterProbs[nClusters+1] <- log(GAMMA) + (VN[nClusters+1]-VN[nClusters]) + (log_numerator - log_denominator)
        # log_clusterProbs[nClusters+1]<-log(GAMMA) + VN[nClusters+1]-VN[nClusters] + dMNorm_marginal(Y=data[,,i], U=U, 
        #                                                                                     V=V, M0=M0, 
        #                                                                                     Sigma0=Sigma0, 
        #                                                                                     Omega0=Omega0,
        #                                                                                     log=TRUE)
        #choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1,
                            prob = exp(log_clusterProbs - sna::logSum(log_clusterProbs) ))
        # remove the empty cluster
        if (cluster.i > nClusters)
        {      clusterAssign[i] <- cur.cluster.i #put the new cluster in the place of the only singleten one
        clusterSizes <- table(as.factor(clusterAssign)) # sorts according to labels
        } else {      
          clusterAssign[i] <- cluster.i
          clusterAssign <- ifelse(clusterAssign > cur.cluster.i, clusterAssign-1, clusterAssign) # to delete the previous group index
          clusterSizes <- table(as.factor(clusterAssign))
          nClusters <- length(clusterSizes) 
          # phi = phi[-cur.cluster.i]}
          M = M[,,-cur.cluster.i] }
      }
    }
    # end for loop over subjects i
    ## update phi ##
    ## update theta ##
    # M = array(0, dim=c(p,q,nClusters))
    # U = matrix(0, nrow = p, ncol = p)
    # V = matrix(0, nrow = q, ncol = q)
    # Z, matrix n by nClusters
    # Z = matrix(0, nrow = n, ncol = nClusters)
    # for(ii in 1:n){
    #   # print(clusterAssign[ii])
    #   Z[ii,clusterAssign[ii]] <- 1
    # }
    # n_Cluster_size <- apply(Z,2,sum)
    # print(nClusters)
    #print(n_Cluster_size)
    # V <- cov2cor(V) # for identifiability purpose
    # rep(0, nClusters)
    # AA = rep(0,nClusters)
    # NN = rep(0,nClusters)
    n_Cluster_size <- table( factor(clusterAssign, levels=1:max(clusterAssign)) )
    # ensure M is a 3-D array
    M <- array(M, dim = c(p,q,nClusters))
    for (r in 1:nClusters){
      temp_zY_sum <- matrix(0, nrow=p, ncol=q)
      for(ii in 1:n){
        if(clusterAssign[ii] == r){
          temp_zY_sum <- temp_zY_sum + data[,,ii]
        }
        # temp_zY_sum <- temp_zY_sum + Z[ii,r] * data[,,i]
      }
      # print(temp_zY_sum)
      temp_zeta <- n_Cluster_size[r] * (solve_V %x% solve_U) + (solve_Omega0 %x% solve_Sigma0)
      temp_xi <- c( solve_Sigma0 %*% M0 %*% solve_Omega0 + solve_U %*% temp_zY_sum %*% solve_V )
      # first generate a vector
      solve_temp_zeta <- solve(temp_zeta)
      # print(solve_temp_zeta %*% temp_xi)
      temp_M_vec <- mvtnorm::rmvnorm(1, mean = solve_temp_zeta %*% temp_xi, sigma = solve_temp_zeta)
      # browser()
      #print(temp_M_vec)
      #print(r)
      #print(nClusters)
      #print(dim(M))
      #print(p)
      #print(q)
      M[,,r] <- matrix(temp_M_vec, nrow = p, ncol = q) # convert a vector to matrix by column
      # clusterAssign is a vector, convert it to Z_ij
      # rgamma(1,alpha + sum(data[clusterAssign == r]), beta + sum(clusterAssign == r))
    }
    # update U
    U_scale <- matrix(0, nrow=p, ncol=p) + 2*beta 
    #
    for(ii in 1:n){
      U_scale <- U_scale + 1 * (data[,,ii] - M[,,clusterAssign[ii]]) %*% solve_V %*% t(data[,,ii] - M[,,clusterAssign[ii]])
    }
    U <- MCMCpack::riwish(2*alpha+n*q, U_scale)
    solve_U = solve(U)
    # Update V
    V_scale <- matrix(0, nrow=q, ncol=q) + 2*rho
    for(ii in 1:n){
      V_scale <- V_scale + 1 * t(data[,,ii] - M[,,clusterAssign[ii]]) %*% solve_U %*% (data[,,ii] - M[,,clusterAssign[ii]])
    }
    V <- MCMCpack::riwish(2*psi+n*p, V_scale)
    # fix the trace of V to be q
    V <- q * V / sum(diag(V))
    solve_V = solve(V)
    print(U)
    print(V)
    # History[[iter]] <- list(zout = clusterAssign, phiout = phi)
    History[[iter]] <- list(zout = clusterAssign, Mout = M, Uout = U, Vout = V)
    cat(" iteration:", iter,"\n",clusterAssign,"\n")
  }# for loop over iterations
  
  list(Iterates = History)
}

### VN only depends on the prior
gamma = 1; lambda = 1; N = 400
VN<-0
tmax = 400+10
for (t in 1:tmax)
{
  r = log(0)
  for (k in t:500)
  {
    b = sum(log((k-t+1):k))-sum(log((k*gamma):(k*gamma+N-1))) + dpois(k-1, lambda, log = TRUE)
    m = max(b,r)
    r = log(exp(r-m) + exp(b-m)) + m
  }
  VN[t] = r
}
