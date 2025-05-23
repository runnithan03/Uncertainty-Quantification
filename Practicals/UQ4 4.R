xD <- matrix(c(1,1,
               2,2,
               0,2),nrow=3,byrow=TRUE)
xD
dist(xD)

as.matrix(dist(xD))
library(pdist)
install.packages("pdist")

##############################################################################################
### Define more efficient Bayes Linear emulator for Multiple inputs in 2D ###

simple_BL_emulator_v3 <- function(xP,             # the set of emulator prediction points
                                  xD,             # the run input locations xD
                                  D,              # the run outputs D = (f(x^1),...,f(x^n))
                                  theta = 1,      # the correlation lengths (can be a vector)
                                  sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])
                                  E_f = 0,         # prior expectation of f: E(f(x)) = 0 
                                  using_pdist = 1  # if you have installed pdist package
){
  
  # store length of runs D and number of prediction points xP
  n <- length(D)
  nP <- nrow(xP)       # XXX New V3
  
  # # Rescale each input by theta. Works for different theta for each input and for same theta
  xP <- t(t(xP)/theta)     # XXX New V3: solution to Exercise 8.3, CP 3.
  xD <- t(t(xD)/theta)     # XXX New V3: solution to Exercise 8.3, CP 3.
  
  ### Define Cov structure of f(x): Cov[f(x),f(xdash)], now to act on matrix of distances ###
  # Cov_fx_fxdash <- function(x,xdash) sig^2 * exp(-sum((x-xdash)^2)/theta^2) # XXX Old V2
  Cov_fx_fxdash <- function(dist_matrix) sigma^2 * exp(-(dist_matrix)^2) # XXX New dist V3
  
  
  ### Define 5 objects needed for BL adjustment ###
  # Create E[D] vector
  E_D <- rep(E_f,n)
  
  # Create Var_D matrix:
  # Var_D <- matrix(0,nrow=n,ncol=n)                                        # XXX Old V2
  # for(i in 1:n) for(j in 1:n) Var_D[i,j] <- Cov_fx_fxdash(xD[i,],xD[j,])  # XXX Old V2
  Var_D <- Cov_fx_fxdash( as.matrix(dist(xD)) )                       # XXX New dist V3
  
  # Create E[f(x)]
  E_fx <- rep(E_f,nP)
  
  # Create Var_f(x) 
  Var_fx <- rep(sigma^2,nP)
  
  # Create Cov_fx_D row vector now using pdist() function if available, if not use dist()
  # Cov_fx_D <- matrix(0,nrow=1,ncol=n)                       # XXX Old V2
  # for(j in 1:n) Cov_fx_D[1,j] <- Cov_fx_fxdash(x,xD[j,])    # XXX Old V2
  if(using_pdist)  Cov_fx_D <- Cov_fx_fxdash( as.matrix(pdist(xP,xD)) )   # XXX New V3
  if(!using_pdist) 
    Cov_fx_D <- Cov_fx_fxdash( as.matrix(dist(rbind(xP,xD)))[1:nP,(nP+1):(nP+n)] )# XXX NewV3
  
  # find inverse of Var_D using Cholesky decomposition (check Wikipedia if interested!)
  Var_D_inv <- chol2inv(chol(Var_D))     # more efficient way to calculate inverse of Cov mat
  
  ### Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x) ###
  cov_fx_D_Var_D_inv <- Cov_fx_D %*% Var_D_inv  # Need this twice so pre-calculate here
  ED_fx   <-  E_fx + cov_fx_D_Var_D_inv %*% (D - E_D) # adj. expectation of ALL xP pts at once
  # VarD_fx <-  Var_fx - cov_fx_D_Var_D_inv %*% t(Cov_fx_D)       # XXX Old V2
  VarD_fx   <-  Var_fx - apply(cov_fx_D_Var_D_inv * Cov_fx_D,1,sum) # fast way to get diagonals 
  # and hence all nP variances (Note: does not do full np x np covariance matrix)
  
  ### return emulator expectation and variance ###
  return(cbind("ExpD_f(x)"=c(ED_fx),"VarD_f(x)"=VarD_fx))  
  
}
### End simple Bayes Linear emulator for single input in 2D ###
##############################################################################################

##############################################################################################
### Demonstrate improved efficiency 

### Redo second output from before as example to check speed ###
f_2 <- function(x) 2*sin( 2*pi*(x[,2]/2 + 0 + (x[,1]-1)^2+(x[,2]-1)^2))

### Define run locations ###
set.seed(1)
xD_w1 <- lhd_maximin(n=14)

xD <- xD_w1
xD_w2 <- matrix(
  c(0.62,0.99,
    0.98,0.28,
    0.95,0.06,
    0.04,0.05,
    0.28,0.2,
    0.04,0.4,
    0.35,0.1,
    0.22,0.36),ncol=2,byrow=TRUE
)
xD <- rbind(xD,xD_w2)   # add wave 2 runs sequentially

### Perform 22 runs of model and store as D (this would takes days for realistic example!) ###
D <- f_2(xD)

### Define 50x50 grid of prediction points xP for emulator evaluation ###
x_grid <- seq(-0.001,1.001,len=50)
xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))

### Evaluate emulator over 50x50=2500 prediction points xP and store as matrices ###

### Old Version 2 without efficiency improvements ###
system.time(em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=0.35,sigma=1,E_f=0)))  

### New Version 3 but without using pdist() ###
system.time(em_out <- simple_BL_emulator_v3(xP=xP,xD=xD,D=D,theta=0.35,sigma=1,E_f=0,using_pdist = 0))  

### New Version 3 but now using pdist() ###
system.time(em_out <- simple_BL_emulator_v3(xP=xP,xD=xD,D=D,theta=0.35,sigma=1,E_f=0))  


E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 

### Evaluate true function and store in matrix for diagnostic comparisons ###
fxP_mat <- matrix(f_2(xP),nrow=length(x_grid),ncol=length(x_grid)) 

### plot emulator expectation, variance, true function and diagnostics ###
emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
               color.palette=exp_cols,main="Emulator Adjusted Expectation E_D[f(x)]")

emul_fill_cont(cont_mat = fxP_mat,cont_levs = seq(-2,2,0.2),xD=xD, x_grid=x_grid,
               color.palette=exp_cols,plot_xD=FALSE,main="True Computer Model Function f(x)")

###############################################################################
### Demonstrate changing thetas

### Define actual 2D computer model/simulator ###
f <- function(x) sin(3*pi*x[,1]) + (1/20)*cos(2*pi*x[,2])

### Define run locations ###
set.seed(2)
xD <- lhd_maximin(n=16)

### Perform 14 runs of model and store as D ###
D <- f(xD)

### Define 50x50 grid of prediction points xP for emulator evaluation ###
x_grid <- seq(-0.001,1.001, len =50)
xP <- as.matrix(expand.grid("x1" =x_grid,"x2"=x_grid))

### Evaluate emulator over 50x50=2500 prediction points xP and store as matrices ###
### Note: now using two different theta values: theta_1=0.4 and theta_2=1.6 ###
em_out <- simple_BL_emulator_v3(xP=xP,xD=xD,D=D,theta=c(0.4,1.6),sigma=1,E_f=0,using_pdist = 1)  
E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 

### Evaluate true function and store in matrix for diagnostic comparisons ###
fxP_mat <- matrix(f(xP),nrow=length(x_grid),ncol=length(x_grid))

### plot emulator expectation, variance, true function and diagnostics ###
emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs = seq(-1.2,1.2,0.2),xD=xD,x_grid=x_grid,
               color.palette=exp_cols,main="Emulator Adjusted Expectation E_D[f(x)")

emul_fill_cont(cont_mat=fxP_mat,cont_levs=seq(-1.2,1.2,0.2),xD=xD,x_grid=x_grid,
               color.palette=exp_cols,plot_xD=FALSE,main="True Computer Model Function f(x)")


##5.1
em_out <- simple_BL_emulator_v3(xP=xP,xD=xD,D=D,theta=c(0.35,0.35),sigma=1,E_f=0,using_pdist = 1)  
E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 

fxP_mat <- matrix(f(xP),nrow=length(x_grid),ncol=length(x_grid))

emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs = seq(-1.2,1.2,0.2),xD=xD,x_grid=x_grid,
               color.palette=exp_cols,main="Emulator Adjusted Expectation E_D[f(x)")
emul_fill_cont(cont_mat=fxP_mat,cont_levs=seq(-1.2,1.2,0.2),xD=xD,x_grid=x_grid,
               color.palette=exp_cols,plot_xD=FALSE,main="True Computer Model Function f(x)")


### Evaluate emulator over 50x50=2500 prediction points xP and store as matrices ###
### Note: now using two different (small) theta values: theta_1=0.1 and theta_2=0.4 ###
em_out <- simple_BL_emulator_v3(xP=xP,xD=xD,D=D,theta=c(0.1,0.4),sigma=1,E_f=0)  
E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 

emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
               color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]")

################################################################################
### Example of Optimisation using an Implausbility Strategy


### Define 2nd 2D computer model/simulator output f(x) ###
f <- function(x) 0.5*(sin(4*pi*(x[,2]+0.5*x[,1]-0.1)) - cos(2*pi*(x[,1]-0.2))) - 3*x[,2]*(x[,2]-1)

### Create wave 1 design of n=12 runs ###
set.seed(1)
xD_w1 <- lhd_maximin(n=12)


xD <- xD_w1   # wave 1 runs

### Perform 14 runs of model and store as D ###
D <- f(xD)

### Evaluate emulator over 50x50=2500 prediction points xP and store as matrices ###
em_out <- simple_BL_emulator_v3(xP=xP,xD=xD,D=D,theta=c(0.45,0.45),sigma=1,E_f=0)
E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 

### Defined Best Run output so far = f_plus and model discrepancy ###
f_plus <- max(D)
sigma_epsilon <- 0.05

### Evaluate true function f(x) and store in matrix for (possible) diag. comparisons ###
fxP_mat <- matrix(f(xP),nrow=length(x_grid),ncol=length(x_grid)) 
### find true max (we pretend we don't know this!) ###
f_plus_true <- max(fxP_mat)

### plot emulator expectation and true function ###
emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
               color.palette=exp_cols,main="Emulator Adjusted Expectation E_D[f(x)]")
emul_fill_cont(cont_mat=fxP_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
               color.palette=exp_cols,plot_xD=FALSE,main="True Computer Model Function f(x)")


### Calculate new Implausibility Measure Over All 50x50 = 2500 input points in xP ###
Imp_mat <-  (f_plus-E_D_fx_mat) / sqrt(Var_D_fx_mat + sigma_epsilon^2) 
Imp_mat[Imp_mat<0] <- 0.0001   # replace negatives with small value for plotting purposes

### Calculate true Implausibility Measure Over All 201 input points in xP ###
Imp_true_mat <- (f_plus_true-fxP_mat) / sqrt(sigma_epsilon^2 ) 
Imp_true_mat[Imp_true_mat>40] <- 40  # limit high implausibilities for plotting purposes

imp_cols <- function(n) turbo(n,begin=0.15,end=1)
imp_levs <- c(0,seq(1,2.75,0.25),seq(3,18,2),20,30,41) 

### Plot New Optimisation Implausibility ###
emul_fill_cont_V2(cont_mat=Imp_mat,cont_levs=imp_levs,cont_levs_lines=3,
                  xD=xD,x_grid=x_grid,xD_col="purple",color.palette=imp_cols,
                  main="Max Implausibility I_M(x)")

### Plot New Optimisation Implausibility for true function with true max ###
emul_fill_cont_V2(cont_mat=Imp_true_mat,cont_levs=imp_levs,cont_levs_lines=3,xD=xD,x_grid=x_grid,
                  xD_col="purple",plot_xD=FALSE,color.palette=imp_cols,
                  main="Implausibility I(x) using True f(x)")




