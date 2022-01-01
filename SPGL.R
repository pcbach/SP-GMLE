rm()
library("Rmosek")
library("CVXR")
library("expm")
msk = MOSEK()
import_solver(msk)

gauss <- function(n){
  #calculate Gaussian quadrature coefficient
  #https://en.wikipedia.org/wiki/Gaussian_quadrature
  if (n == 1){
    x <- 0
    w <- 2
  }else{
    n_ <- seq(from = 1,to = n-1,by = 1)
    beta <- 0.5/sqrt(1-(2*n_)^(-2))
    T <- matrix(0,n,n)
    T[row(T)-1 == col(T)] <- beta
    T[row(T)+1 == col(T)] <- beta
    VD <- eigen(T)
    V <- VD$vectors
    D <- VD$values
    
    x <- D
    P <- 1/((1-x)*V[row(V) == 1]^2)
    ww <- 2/(n*(n+1)*P^2) 
    w <- 2*V[row(V) == 1]^2
  }
  w = w/2
  x = 1/2+1/2*x
  return(list(x,w))
}

#SIGMA: CVXR pxp variable
#optval: optimal solution value
#n: Gaussian quadrature degree
#d: dimension
#sigma_i: Sample covariance of the data
represent <- function(SIGMA,optval,n,d,sigma_i){
  
  #ids for debug access by result[[ids[["var"]]]]
  ids <- list()
  ids[["sigma"]] <- toString(id(SIGMA))
  ids[["optval"]] <- toString(id(optval))
  
  sqrtinv_Sn <- sqrtm(solve(sigma_i))
  
  #Gaussian quadrature coefficient
  xw <- gauss(n)
  t <- unlist(xw[1])
  w <- unlist(xw[2])
  
  #variable
  x <- Variable(d,d) #X
  ids[["x"]] <- toString(id(x))
  I <- Variable(d,d) #Identity matrix
  
  w_tau <- Variable(n)
  
  constraints <- list(
    optval == sum_entries(w_tau),
    x == sqrtinv_Sn %*% SIGMA %*% sqrtinv_Sn,
    I == diag(d) #Identity matrix
  )
  
  ind <- seq(from = 1, to = n,by = 1)
  for (i in ind){
    G00 <- Variable(d,d)
    G01 <- Variable(d,d)
    G02 <- Variable(d,d)
    G11 <- Variable(d,d)
    G12 <- Variable(d,d)
    G22 <- Variable(d,d)
    H00 <- Variable(d,d)
    H01 <- Variable(d,d)
    H11 <- Variable(d,d)
    
    tau <- Variable(d,d)
    ids[[paste("tau",i-1,sep = "_")]] <- toString(id(tau))
    
    G <- Variable(3*d,3*d,PSD = TRUE)
    ids[[paste("G",i-1,sep = "_")]] <- toString(id(G))
    
    H <- Variable(2*d,2*d,PSD = TRUE)
    ids[[paste("H",i-1,sep = "_")]] <- toString(id(H))
    
    constraints <- c(constraints,list(
      G == bmat(list(list(G00,G01,G02),list(t(G01),G11,G12),list(t(G02),t(G12),G22))),
      H == bmat(list(list(H00,H01),list(t(H01),H11))),
      G00 == x*(t[i]-1)^2,
      G01 + t(G01) + 3/2 * H00 == -2*(t[i]-1)*(t[i]*(x+I)-I),
      G02 + t(G02) + G11 - H00 + 3/2 * H01 + 3/2 * t(H01) == (t[i]-1)*(t[i]*(4*I+x)+tau*(t[i]-1)+x-I),
      G12 + t(G12) + 3/2 * H11 - H01 - t(H01) == -2*(t[i]-1)*(tau+I)*t[i],
      G22 - H11 == t[i]*(tau*t[i]-I),
      w_tau[i] == matrix_trace(tau)*w[i]
    ))
  }
  return(list("constraints" = constraints,"variableid" = ids))
}

