rm()
library("Rmosek")
library("CVXR")
library("mvtnorm")
library("expm")
library("psych")
#set.seed(124)
msk = MOSEK()
import_solver(msk)

gauss <- function(n){
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
g2 <- function(x,t){
  return(-2*t/(t*(x-1)+1)^3 + 2/(x^3))
}
check <- function(x,n){
  xw = gauss(n)
  
  t <- unlist(xw[1])
  w <- unlist(xw[2])
  #print(w)
  val = 0
  for (i in seq(from = 1, to = n, by = 1)){
    val = val + w[i]*g2(x,t[i])
  }
  return(val > 0)
}
bound <- function(n){
  b = 0
  e = 2
  while((e - b)>1e-13){
    m = (e+b)/2
    if(check(m,n)){
      b = m
    }else{
      e = m
    }
  }
  return((e+b)/2)
}
sampleCov <- function(data){
  N <- dim(data)[1]
  print(N)
  d <- dim(data)[2]
  total = data[1,1:d] %*% t(data[1,1:d])
  for (i in seq(from = 2, to = N,by = 1)){
    total = total + data[i,1:d] %*% t(data[i,1:d])
  }
  return(total/N)
}

represent <- function(SIGMA,optval,n,d,sigma_i){
  
  ids <- list()
  ids[["sigma"]] <- toString(id(SIGMA))
  ids[["optval"]] <- toString(id(optval))
  
  sqrtinv_Sn <- sqrtm(solve(sigma_i))
  SI_SN <- Variable(d,d)
  
  xw <- gauss(n)
  t <- unlist(xw[1])
  w <- unlist(xw[2])
  
  
  x <- Variable(d,d)
  ids[["x"]] <- toString(id(x))
  
  x_ <- Variable(d,d,PSD = TRUE)
  x_2 <- Variable(d,d,PSD = TRUE)
  I <- Variable(d,d)
  
  w_tau <- Variable(n)
  b <- bound(n)
  #cat(sprintf('%.20f',b))
  constraints <- list(
    optval == sum_entries(w_tau),
    x == sqrtinv_Sn %*% SIGMA %*% sqrtinv_Sn,
    I == diag(d)
    #x_ == I*1.5 - x,
    #x_2 == x - 0.5*I
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

likelihood <- function(data, sigmaa){
  N <- dim(data)[1]
  d <- dim(data)[2]
  total = t(data[1,1:d]) %*% solve(sigmaa) %*% data[1,1:d]
  for (i in seq(from = 2, to = N,by = 1)){
    total = total + t(data[i,1:d]) %*% solve(sigmaa) %*% data[i,1:d]
  }
  SN <- total/N
  return(d/2*log(2*pi) + 1/2*log(det(sigmaa)) + 1/2 * SN)
  
  
}

fx <- function(x){
  return(log(det(x)) + tr(solve(x)))  
} 
fx1 <- function(x){
  x <- x*diag(4)
  return(log(det(x)) + tr(solve(x)))  
}
