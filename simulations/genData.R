# Generate a sample of n brownian motions
# Input:
#   N: Defines 0:N/N, the points at which the brownian motions are calculated
#      in [0,1]
#   n: The sample size
# Output:
#   process: Functional data object containing a sample of n standard brownian
#            motions.
BM <- function(N = 300, n = 100) {
  
  data = matrix(rnorm(N*n), ncol = n)
  values = matrix(nrow = N+1, ncol = n)
  values[1,] = 0
  for(i in 1:n) {
    values[2:(N+1),i] = 1/sqrt(N) * cumsum(data[,i])
  }
  
  basis = create.bspline.basis(norder = 2, breaks = 0:N/N)
  process = Data2fd(argvals = 0:N/N, y = values, basisobj = basis, lambda = 0)
  
  return(process)
}

# Generate a sample of brownian bridges
# Input:
#   dataInfo: A list which contains all the necessary Info:
#     $n: Sample size
#     $nArgvals: Defines 0:N/N, the points at which the brownian motions are
#                calculated in [0,1]
#     $basisType: One of c("bspline", "fourier"), determines the basis functions
#     $muInfo: A list containing the Info about the expectation:
#       $type: One of c("sin", "horv") which determines the shape of the
#              expectation function
#       $a: The expectation function is multiplied with this factor
BB <- function(dataInfo) {
  
  n = dataInfo$n
  N = dataInfo$nArgvals
  muInfo = dataInfo$muInfo
  
  data = matrix(rnorm(N*n), ncol = n)
  values = matrix(nrow = N+1, ncol = n)
  values[1,] = 0
  t = 1:N / N
  for(i in 1:n) {
    cumulative = 1/sqrt(N) * cumsum(data[,i])
    values[2:(N+1),i] =  cumulative - t * cumulative[N]
  }
  
  t = 0:N / N
  # define the mu function which is added to the sample
  switch(muInfo$type,
         sin = {
           mu = muInfo$a * sin(2*pi*t)
         },
         horv = {
           mu = muInfo$a * t * (1 - t)
         }
  )
  
  switch(dataInfo$basisType,
         bspline = {
           basis = create.bspline.basis(rangeval = c(0, 1), norder = 2, breaks = 0:N/N)
         },
         fourier = {
           basis = create.fourier.basis(rangeval = c(0, 1), N+1)
         }
  )
  
  process = Data2fd(argvals = t, y = values + mu, basisobj =  basis)
  
  return(process)
}

# Function generates independent data
# Input:
#   dataInfo: A list which contains all the necessary Info:
#     $n: Sample size
#     $varType: One of c("A","B"), determines the variance of the coefficients
#               of the basis functions.
#     $factor: The errors are multiplied with this factor
#     $gaussian: Boolean, determining if the coeffinients are normally
#                distibuted or t_5 distributed and scaled to have variance 1
#     $basisType: One of c("bspline", "fourier"), determines the basis functions
#     $muInfo: A list containing the Info about the expectation:
#       $type: One of c("sin", "horv") which determines the shape of the
#              expectation function
#       $a: The expectation function is multiplied with this factor
fIID <- function(dataInfo) {
  
  n = dataInfo$n
  muInfo = dataInfo$muInfo
  if(is.null(dataInfo$factor)) {
    factor = 1
  } else {factor <- dataInfo$factor}
  if(is.null(dataInfo$gaussian)) {
    gaussian = TRUE
  } else {gaussian = dataInfo$gaussian}
  
  switch(dataInfo$basisType,
         bspline = {
           basis = create.bspline.basis(rangeval=c(0, 1), nbasis=21)
         },
         fourier = {
           basis = create.fourier.basis(rangeval=c(0, 1), nbasis=21)
         }
  )
  
  # get number of coefficients
  D = basis$nbasis
  # importance of components
  switch(dataInfo$varType,
         A = { Sigma = 1/seq(1:D) },
         B = { Sigma = 1.2^(-seq(1:D)) }
  )
  
  # build coefficients matrix (initially null)
  coef = matrix(0, D, n)
  for (i in 1:n){
    if(gaussian) {
      coef[,i] = rnorm(D, 0, Sigma)
    } else {coef[,i] = rt(D, df = 5) * sqrt(3/5) * Sigma }
  }
  # generate the functional data
  fdata = fd(factor * coef, basis)
  
  if(muInfo$a != 0) {
    
    # define the mu function which is added to the sample
    N=300
    t = 0:N / N
    switch(muInfo$type,
           sin = {
             mu = muInfo$a * sin(2*pi*t)
           },
           horv = {
             mu = muInfo$a * t * (1 - t)
           }
    )
    # generate functional data object mu and add it to fdata
    mu.fd = Data2fd(argvals = t, y = mu, basisobj =  basis, lambda = 0)
    for(i in 1:n) {fdata$coefs[,i] = fdata$coefs[,i] + mu.fd$coefs}
  }
  
  return(fdata)
}


# Function generates fMA(1) data
# Input:
#   dataInfo: A list which contains all the necessary Info:
#     $n: Sample size
#     $varType: One of c("A","B"), determines the variance of the coefficients
#               of the basis functions.
#     $factor: The errors are multiplied with this factor
#     $kappa: This number determines how strong the dependency is in the fMA(1)
#     $Psi: fMA(1) operator (21x21 matrix). If NULL, it is randomly generated
#     $basisType: One of c("bspline", "fourier"), determines the basis functions
#     $muInfo: A list containing the Info about the expectation:
#       $type: One of c("sin", "horv") which determines the shape of the
#              expectation function
#       $a: The expectation function is multiplied with this factor
fMA1 <- function(dataInfo){
  
  n = dataInfo$n
  muInfo = dataInfo$muInfo
  if(is.null(dataInfo$factor)) {
    factor = 1
  } else {factor <- dataInfo$factor}
  
  switch(dataInfo$basisType,
         bspline = {
           basis = create.bspline.basis(rangeval=c(0, 1), nbasis=21)
         },
         fourier = {
           basis = create.fourier.basis(rangeval=c(0, 1), nbasis=21)
         }
  )
  
  # get number of coefficients
  D = basis$nbasis
  # importance of components
  switch(dataInfo$varType,
         A = { Sigma = 1/seq(1:D) },
         B = { Sigma = 1.2^(-seq(1:D)) }
  )
  
  if(is.null(dataInfo$Psi)) {
    Psi = matrix(0, D, D)
    for(i in 1:D){for(j in 1:D){Psi[i,j]=rnorm(1,0,Sigma[i]*Sigma[j])}}
    Psi = Psi/sqrt(eigen(Psi%*%t(Psi))$values[1])
  } else {
    Psi = dataInfo$Psi
  }
  
  Theta = dataInfo$kappa * Psi
  
  # create the fMA(1) process
  fdata = createMA(n = n, basis = basis, Theta = Theta, Sigma = Sigma)
  fdata$coefs = fdata$coefs * factor
  
  if(muInfo$a != 0) {
    # define the mu function which is added to the sample
    N = 300
    t = 0:N / N
    switch(muInfo$type,
           sin = {
             mu = muInfo$a * sin(2*pi*t)
           },
           horv = {
             mu = muInfo$a * t * (1 - t)
           }
    )
    # generate functional data object mu and add it to fdata
    mu.fd = Data2fd(argvals = t, y = mu, basisobj =  basis, lambda = 0)
    for(i in 1:n) {fdata$coefs[,i] = fdata$coefs[,i] + mu.fd$coefs}
  }
  
  return(fdata)
}

# help function for fMA1
createMA = function(n, basis, Theta, Sigma){
  
  # get number of coefficients
  D = basis$nbasis
  
  # build coefficients matrix (initially null)
  coef = matrix(0,D,n)
  
  # define the MA1
  zlag0 = matrix(0,D,n+1)
  for(i in 1:(n+1)){zlag0[,i]=rnorm(D,0,Sigma)}
  zlag1 = matrix(0,D,n)
  for (i in 2:n){zlag1[,i]=Theta %*% zlag0[,i-1]}
  
  coef=zlag1+zlag0[,2:(n+1)]
  
  # return time series
  return(fd(coef, basis = basis))
}


# Function generates fAR(1) data
# Input:
#   dataInfo: A list which contains all the necessary Info:
#     $n: Sample size
#     $varType: One of c("A","B"), determines the variance of the coefficients
#               of the basis functions.
#     $factor: The errors are multiplied with this factor
#     $gaussian: Boolean, determining if the coeffinients are normally
#                distibuted or t_5 distributed and scaled to have variance 1
#     $kappa: This number determines how strong the dependency is in the fAR(1)
#     $Psi: fAR(1) operator (21x21 matrix). If NULL, it is randomly generated
#     $basisType: One of c("bspline", "fourier"), determines the basis functions
#     $muInfo: A list containing the Info about the expectation:
#       $type: One of c("sin", "horv") which determines the shape of the
#              expectation function
#       $a: The expectation function is multiplied with this factor
fAR1 <- function(dataInfo){
  
  n = dataInfo$n
  muInfo = dataInfo$muInfo
  if(is.null(dataInfo$gaussian)) {
    gaussian = TRUE
  } else {gaussian = dataInfo$gaussian}
  
  switch(dataInfo$basisType,
         bspline = {
           basis = create.bspline.basis(rangeval=c(0, 1), nbasis=21)
         },
         fourier = {
           basis = create.fourier.basis(rangeval=c(0, 1), nbasis=21)
         }
  )
  
  # get number of coefficients
  D = basis$nbasis
  # importance of components
  switch(dataInfo$varType,
         A = { Sigma = 1/seq(1:D) },
         B = { Sigma = 1.2^(-seq(1:D)) }
  )
  
  if(is.null(dataInfo$Psi)) {
    Psi = matrix(0, D, D)
    for(i in 1:D){for(j in 1:D){Psi[i,j]=rnorm(1,0,Sigma[i]*Sigma[j])}}
    Psi = Psi/sqrt(eigen(Psi%*%t(Psi))$values[1])
  } else {
    Psi = dataInfo$Psi
  }
  
  Theta = dataInfo$kappa * Psi
  
  fdata = createAR1(n = n, basis = basis, Theta = Theta, Sigma = Sigma,
                    gaussian = gaussian)
  
  if(muInfo$a != 0) {
    # define the mu function which is added to the sample
    N = 300
    t = 0:N / N
    switch(muInfo$type,
           sin = {
             mu = muInfo$a * sin(2*pi*t)
           },
           horv = {
             mu = muInfo$a * t * (1 - t)
           }
    )
    # generate functional data object mu and add it to fdata
    mu.fd = Data2fd(argvals = t, y = mu, basisobj =  basis, lambda = 0)
    for(i in 1:n) {fdata$coefs[,i] = fdata$coefs[,i] + mu.fd$coefs}
  }
  
  return(fdata)
}

# help function for fAR1
# Input:
#   n - number of observations
#   Theta - autoregresive operator
#   Sigma - noise standard deviation
#   basis - fda basis
#   gaussian - boolean, gaussian or nongaussian data?
createAR1 = function(n, basis, Theta, Sigma, gaussian) {
  
  # get number of coefficients
  D = basis$nbasis
  
  # build coefficients matrix (initially null)
  nBurnIn <- 100
  coef = matrix(0,D,n+nBurnIn)
  
  if(gaussian) {
    coef[,1] = rnorm(D,0,Sigma)
    # recursion
    for (i in 2:(n+nBurnIn)) {
      coef[,i] = Theta %*% coef[,i-1] + rnorm(D,0,Sigma)
    }
  } else {
    coef[,1] = rt(D, df = 5) * sqrt(3/5) * Sigma
    # recursion
    for (i in 2:(n+nBurnIn)) {
      coef[,i] = Theta %*% coef[,i-1] + rt(D, df = 5) * sqrt(3/5) * Sigma
    }
  }
  
  # return time series
  return(fd(coef[,(nBurnIn+1):(nBurnIn+n)], basis = basis))
}

# Generate a sample non-gaussian data motivated by Kraus and Panaretos (2012)
# Input:
#   dataInfo: A list which contains all the necessary Info:
#     $n: Sample size
#     $factor: The errors are multiplied with this factor
#     $muInfo: A list containing the Info about the expectation:
#       $type: One of c("sin", "horv") which determines the shape of the
#              expectation function
#       $a: The expectation function is multiplied with this factor
nonGaussian <- function(dataInfo) {
  
  n = dataInfo$n
  muInfo = dataInfo$muInfo
  if(is.null(dataInfo$factor)) {
    factor = 1
  } else {factor <- dataInfo$factor}
  
  coef = matrix(nrow = 21, ncol = n)
  for(i in 1:n) {
    t.5 = rt(n = 20, df = 5)
    coef[1,i] = 0
    for (j in 2:21) {
      if(j %% 2 == 0) {
        coef[j,i] = 1/sqrt(10) * 1/(j/2)^(3/2) * t.5[j-1] * sqrt(3/5)
      } else {
        coef[j,i] = 1/sqrt(10) * 1/3^( ((j-1)/2) /2) * t.5[j-1] * sqrt(3/5)
      }
    }
  }
  basis = create.fourier.basis(rangeval = c(0,1), nbasis = 21)
  process = fd(coef = factor * coef, basisobj = basis)
  
  N = 500
  t = t = 0:(N-1) / (N-1)
  # define the mu function which is added to the sample
  switch(muInfo$type,
         sin = {
           mu = muInfo$a * sin(2*pi*t)
         },
         horv = {
           mu = muInfo$a * t * (1 - t)
         }
  )
  mu.fd = Data2fd(argvals = t, y = mu, basisobj =  basis)
  
  
  coef <- apply(process$coefs, 2, function(column){column + mu.fd$coefs})
  process = fd(coef = coef, basisobj = basis)
  
  return(process)
}


## Two sample problem specific

# Generate two dependent fMA(1) processes
# Input:
# dataInfo - list which contains necessary info about the data
#   $m, $n - sample sizes
#   $factor - number with which the errors of the second sample are multiplied
#   $basisType - one of c("bspline", "fourier")
#       $muInfo: A list containing the Info about the expectation:
#       $type: One of c("sin", "horv") which determines the shape of the
#              expectation function
#       $a: The expectation function is multiplied with this factor
#
# Output:
#   res: list containing two functional data samples:
#     $fdata1, $fdata2 - first and second sample
dependentSamples <- function(dataInfo) {
  
  m = dataInfo$m
  n = dataInfo$n
  muInfo = dataInfo$muInfo
  
  dataInfoInitSample = dataInfo
  dataInfoInitSample$m = NULL
  dataInfoInitSample$n = m+n
  dataInfoInitSample$muInfo$a = 0
  dataInfoInitSample$factor = 1
  fdata = fMA1(dataInfoInitSample)
  fdata1 = fdata[1:n]
  fdata2 = fdata[(n+1):(n+m)] * dataInfo$factor
  
  N = 300
  t = 0:N / N
  # define the mu function which is added to the sample
  switch(muInfo$type,
         sin = {
           mu = muInfo$a * sin(2*pi*t)
         },
         horv = {
           mu = muInfo$a * t * (1 - t)
         }
  )
  # generate functional data object mu and add it to fdata
  switch(dataInfo$basisType,
         bspline = {
           basis = create.bspline.basis(rangeval=c(0, 1), nbasis=21)
         },
         fourier = {
           basis = create.fourier.basis(rangeval=c(0, 1), nbasis=21)
         }
  )
  mu.fd = Data2fd(argvals = t, y = mu, basisobj =  basis, lambda = 0)
  for(i in 1:m) {fdata2$coefs[,i] = fdata2$coefs[,i] + mu.fd$coefs}
  
  res = list()
  res$fdata1 = fdata1
  res$fdata2 = fdata2
  
  return(res)
}


## Change Point Problem specific

# Use the fIID function and add a change to the sample
# Input:
#   dataInfo: A list which contains all the necessary Info:
#     $n: Sample size
#     $varType: One of c("A","B"), determines the variance of the coefficients
#               of the basis functions.
#     $s.star: number in (0,1) determining the change point
#     $factor: The errors of the sample after the change point are multiplied
#              with this factor
#     $gaussian: Boolean, determining if the coeffinients are normally
#                distibuted or t_5 distributed and scaled to have variance 1
#     $basisType: One of c("bspline", "fourier"), determines the basis functions
#     $muInfo: A list containing the Info about the expectation, which is added
#              to the sample after the change point:
#       $type: One of c("sin", "horv") which determines the shape of the
#              expectation function
#       $a: The expectation function is multiplied with this factor
fIID.cp <- function(dataInfo) {
  
  dataInfoBeforeCP = dataInfo
  dataInfoBeforeCP$muInfo$a = 0
  dataInfoBeforeCP$factor = 1
  fdata = fIID(dataInfoBeforeCP)
  
  n = dataInfo$n
  s.star = dataInfo$s.star
  muInfo = dataInfo$muInfo
  
  # add a function to the data after the cp
  N = 300
  t = 0:N / N
  switch(muInfo$type,
         sin = {
           mu = muInfo$a * sin(2*pi*t)
         },
         horv = {
           mu = muInfo$a * t * (1 - t)
         }
  )
  switch(dataInfo$basisType,
         bspline = {
           basis = create.bspline.basis(rangeval=c(0, 1), nbasis=21)
         },
         fourier = {
           basis = create.fourier.basis(rangeval=c(0, 1), nbasis=21)
         }
  )
  mu.fd = Data2fd(argvals = t, y = mu, basisobj =  basis, lambda = 0)
  for(i in (floor(n*s.star)+1):n) {
    fdata$coefs[,i] = fdata$coefs[,i]* dataInfo$factor + mu.fd$coefs
  }
  
  return(fdata)
}

# Use the fMA1 function and add a change to the sample
# Input:
#   dataInfo: A list which contains all the necessary Info:
#     $n: Sample size
#     $varType: One of c("A","B"), determines the variance of the coefficients
#               of the basis functions.
#     $s.star: number in (0,1) determining the change point
#     $factor: The errors of the sample after the change point are multiplied
#              with this factor
#     $basisType: One of c("bspline", "fourier"), determines the basis functions
#     $muInfo: A list containing the Info about the expectation, which is added
#              to the sample after the change point:
#       $type: One of c("sin", "horv") which determines the shape of the
#              expectation function
#       $a: The expectation function is multiplied with this factor
fMA1.cp <- function(dataInfo) {
  
  dataInfoBeforeCP = dataInfo
  dataInfoBeforeCP$muInfo$a = 0
  dataInfoBeforeCP$factor = 1
  fdata = fMA1(dataInfoBeforeCP)
  
  n = dataInfo$n
  s.star = dataInfo$s.star
  muInfo = dataInfo$muInfo
  
  # add a function to the data after the cp
  N = 300
  t = 0:N / N
  switch(muInfo$type,
         sin = {
           mu = muInfo$a * sin(2*pi*t)
         },
         horv = {
           mu = muInfo$a * t * (1 - t)
         }
  )
  switch(dataInfo$basisType,
         bspline = {
           basis = create.bspline.basis(rangeval=c(0,1), nbasis=21)
         },
         fourier = {
           basis = create.fourier.basis(rangeval=c(0,1), nbasis=21)
         }
  )
  mu.fd = Data2fd(argvals = t, y = mu, basisobj =  basis, lambda = 0)
  for(i in (floor(n*s.star)+1):n) {
    fdata$coefs[,i] = fdata$coefs[,i] * dataInfo$factor + mu.fd$coefs
  }
  
  return(fdata)
}


# generate a sample non-gaussian data motivated by Kraus and Panaretos (2012)
# with a change point
# Input:
#   dataInfo: A list which contains all the necessary Info:
#     $n: Sample size
#     $s.star: number in (0,1) determining the change point
#     $factor: The errors are multiplied with this factor
#     $muInfo: A list containing the Info about the expectation:
#       $type: One of c("sin", "horv") which determines the shape of the
#              expectation function
#       $a: The expectation function is multiplied with this factor
nonGaussian.cp <- function(dataInfo) {
  
  dataInfoBeforeCP = dataInfo
  dataInfoBeforeCP$muInfo$a = 0
  dataInfoBeforeCP$factor = 1
  fdata = nonGaussian(dataInfoBeforeCP)
  
  n = dataInfo$n
  s.star = dataInfo$s.star
  muInfo = dataInfo$muInfo
  
  # add a function to the data after the cp
  N = 500
  t = t = 0:(N-1) / (N-1)
  # define the mu function which is added to the sample
  switch(muInfo$type,
         sin = {
           mu = muInfo$a * sin(2*pi*t)
         },
         horv = {
           mu = muInfo$a * t * (1 - t)
         }
  )
  basis = create.fourier.basis(rangeval=c(0,1), nbasis=21)
  mu.fd = Data2fd(argvals = t, y = mu, basisobj =  basis)
  for(i in (floor(n*s.star)+1):n) {
    fdata$coefs[,i] = fdata$coefs[,i] * dataInfo$factor + mu.fd$coefs
  }
  
  return(fdata)
}
