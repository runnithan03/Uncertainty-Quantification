install.packages("plot3D")
install.packages("lhs")
install.packages("pdist")
install.packages("sensitivity")
install.packages("viridisLite")

library(plot3D)
library(lhs)
library(pdist)
library(sensitivity)
library(viridisLite)

simple_BL_emulator_v1 <- function(x, xD, D, theta = 1, sigma = 1, E_f = 0, nugget = 1e-8) {
  n <- length(D)
  
  Cov_fx_fxdash <- function(x, xdash) sigma^2 * exp(-(x - xdash)^2 / theta^2)
  
  E_D <- rep(E_f, n)
  
  # Construct the covariance matrix Var_D
  Var_D <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      Var_D[i, j] <- Cov_fx_fxdash(xD[i], xD[j])
    }
  }
  
  # The nugget ensures positive definiteness
  diag(Var_D) <- diag(Var_D) + nugget
  
  # Expectation and variance
  E_fx <- E_f
  Var_fx <- sigma^2
  
  # Create Cov_fx_D as a column vector (matrix form)
  Cov_fx_D <- matrix(sapply(xD, function(xD_j) Cov_fx_fxdash(x, xD_j)), nrow = 1)
  D <- matrix(D, ncol = 1)
  E_D <- matrix(E_D, ncol = 1)
  
  # Perform Bayes Linear adjustment
  ED_fx <- E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)
  VarD_fx <- Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)
  
  # Ensure variance is non-negative
  VarD_fx <- pmax(VarD_fx, 0)
  
  return(c("ExpD_f(x)" = ED_fx, "VarD_f(x)" = VarD_fx))
}

### Define simple Bayes Linear emulator for single input in 2D ###

simple_BL_emulator_v2 <- function(x,              # the emulator prediction point
                                  xD,             # the run input locations xD
                                  D,              # the run outputs D = (f(x^1),...,f(x^n))
                                  theta = 1,      # the correlation lengths
                                  sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])
                                  E_f = 0         # prior expectation of f: E(f(x)) = 0 
){
  
  # store length of runs D  
  n <- length(D)
  
  ### Define Covariance structure of f(x): Cov[f(x),f(xdash)] ###
  # Cov_fx_fxdash <- function(x,xdash) sigma^2 * exp(-(x-xdash)^2/theta^2)    # XXX Old 1D version
  Cov_fx_fxdash <- function(x,xdash) sigma^2 * exp(-sum((x-xdash)^2)/theta^2) # XXX New 2D version
  
  
  ### Define 5 objects needed for BL adjustment ###
  # Create E[D] vector
  E_D <- rep(E_f,n)
  
  # Create Var_D matrix:
  Var_D <- matrix(0,nrow=n,ncol=n)
  # for(i in 1:n) for(j in 1:n) Var_D[i,j] <- Cov_fx_fxdash(xD[i],xD[j])  # XXX Old 1D version
  for(i in 1:n) for(j in 1:n) Var_D[i,j] <- Cov_fx_fxdash(xD[i,],xD[j,])  # XXX New 2D version
  
  # Create E[f(x)]
  E_fx <- E_f
  
  # Create Var_f(x) 
  Var_fx <- sigma^2
  
  # Create Cov_fx_D row vector
  Cov_fx_D <- matrix(0,nrow=1,ncol=n)
  # for(j in 1:n) Cov_fx_D[1,j] <- Cov_fx_fxdash(x,xD[j])    # XXX Old 1D version
  for(j in 1:n) Cov_fx_D[1,j] <- Cov_fx_fxdash(x,xD[j,])    # XXX New 2D version
  
  
  ### Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x) ###
  ED_fx   <-  E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)   
  VarD_fx <-  Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)  
  
  ### return emulator expectation and variance ###
  return(c("ExpD_f(x)"=ED_fx,"VarD_f(x)"=VarD_fx))  
  
}
### End simple Bayes Linear emulator for single input in 2D ###
##############################################################################################


### Define more efficient Bayes Linear emulator for Multiple inputs in 2D ###
simple_BL_emulator_v3 <- function(xP,             # the set of emulator prediction points
                                  xD,             # the run input locations xD
                                  D,              # the run outputs D = (f(x^1),...,f(x^n))
                                  theta = 1,      # the correlation lengths (can be a vector)
                                  sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])
                                  E_f = 0,        # prior expectation of f: E(f(x)) = 0 
                                  using_pdist = 1 # if you have installed pdist package
){
  
  # store length of runs D and number of prediction points xP
  n <- length(D)
  nP <- ifelse(is.null(nrow(xP)), 1, nrow(xP))  # XXX New V3: Ensure nP is defined correctly
  
  # Ensure D and E_D are column vectors
  D <- matrix(D, ncol = 1)  # XXX New V3: Ensure proper structure for matrix multiplication
  E_D <- matrix(rep(E_f, n), ncol = 1)  # XXX New V3: Ensure E_D is a column vector
  
  # Ensure xP is a matrix
  xP <- matrix(xP, ncol = ncol(xD))  # XXX New V3: Ensure xP remains a matrix
  
  # Rescale each input by theta. Works for different theta for each input and for same theta
  xP <- t(t(xP) / theta)     # XXX New V3: solution to Exercise 8.3, CP 3.
  xD <- t(t(xD) / theta)     # XXX New V3: solution to Exercise 8.3, CP 3.
  
  ### Define Cov structure of f(x): Cov[f(x),f(xdash)], now to act on matrix of distances ###
  Cov_fx_fxdash <- function(dist_matrix) sigma^2 * exp(-(dist_matrix)^2) # XXX New dist V3
  
  ### Define 5 objects needed for BL adjustment ###
  # Create Var_D matrix:
  Var_D <- Cov_fx_fxdash(as.matrix(dist(xD)))  # XXX New dist V3
  
  # Handle potential numerical instability in Cholesky decomposition
  Var_D_inv <- tryCatch({
    chol2inv(chol(Var_D))  # more efficient way to calculate inverse of Cov mat
  }, error = function(e) {
    diag(1, n, n)  # XXX New V3: Fallback to identity matrix
  })
  
  # Ensure Var_D_inv is square
  Var_D_inv <- matrix(Var_D_inv, nrow = n, ncol = n)  # XXX New V3: Ensure correct shape
  
  # Create E[f(x)]
  E_fx <- rep(E_f, nP)  # XXX New V3: Ensure E_f replication is valid
  
  # Create Var_f(x) 
  Var_fx <- rep(sigma^2, nP)
  
  # Create Cov_fx_D row vector now using pdist() function if available, if not use dist()
  if(using_pdist) {  
    Cov_fx_D <- Cov_fx_fxdash(as.matrix(pdist(xP, xD)))  
  } else {  
    Cov_fx_D <- Cov_fx_fxdash(as.matrix(dist(rbind(xP, xD)))[1:nP, (nP+1):(nP+n), drop = FALSE])  
  }
  
  # Ensure Cov_fx_D has correct shape
  Cov_fx_D <- matrix(Cov_fx_D, nrow = nP, ncol = n)  # XXX New V3: Fix incorrect shape
  
  ### Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x) ###
  cov_fx_D_Var_D_inv <- Cov_fx_D %*% Var_D_inv  # Need this twice so pre-calculate here
  
  # Adjusted Expectation
  ED_fx <- E_fx + cov_fx_D_Var_D_inv %*% matrix(D - E_D, ncol = 1)  # XXX New V3: Fixed non-conformable multiplication
  
  # Adjusted Variance (Diagonal Only)
  VarD_fx <- Var_fx - apply(cov_fx_D_Var_D_inv * Cov_fx_D, 1, sum) # fast way to get diagonals 
  # and hence all nP variances (Note: does not do full np x np covariance matrix)
  
  ### return emulator expectation and variance ###
  return(cbind("ExpD_f(x)" = c(ED_fx), "VarD_f(x)" = VarD_fx))  
}
### End of Fixed Bayes Linear Emulator ###

braninsc <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
  
  x1bar <- 15*x1 - 5
  x2bar <- 15 * x2
  
  term1 <- x2bar - 5.1*x1bar^2/(4*pi^2) + 5*x1bar/pi - 6
  term2 <- (10 - 10/(8*pi)) * cos(x1bar)
  
  y <- (term1^2 + term2 - 44.81) / 51.95
  return(y)
}

# Model Exploration 

# Rescaled Branin Function Plot
x1 <- seq(0, 1, length.out = 100)
x2 <- seq(0, 1, length.out = 100)
grid <- expand.grid(x1 = x1, x2 = x2)
grid$z <- apply(grid, 1, braninsc)
z_matrix <- matrix(grid$z, nrow = 100, byrow = TRUE)
persp3D(x1, x2, z_matrix, theta = 45, phi = 25, col = "lightblue",
        xlab = "x1", ylab = "x2", zlab = "f(x1, x2)",
        main = "Rescaled Branin Function")

# Vary x1, keeping x2 fixed at 0.5
x1 <- seq(0, 1, length.out = 100)
x2_fixed <- 0.5
y_1D <- sapply(x1, function(x) braninsc(c(x, x2_fixed)))

# Plot 1D exploration
plot(x1, y_1D, type = "l", col = "blue", lwd = 2,
     main = "1D Exploration: Varying x1 (setting x2 = 0.5)",
     xlab = "x1", ylab = "f(x1, 0.5)")

# 2D Contour Plot
filled.contour(x1, x2, z_matrix,
               main = "2D Contour Plot of Branin Function",
               xlab = "x1", ylab = "x2")

# Emulation Preparation already done.

# 1D Emulation
xD_1D <- matrix(seq(0, 1, length.out = 10), ncol = 1)
D_1D <- apply(xD_1D, 1, function(x) braninsc(c(x, 0.5)))

xP <- matrix(seq(0, 1, length.out = 100), ncol = 1)

em_out <- t(sapply(xP, simple_BL_emulator_v1, xD = xD_1D, D = D_1D, theta = 0.2, sigma = 0.5, E_f = 0))

plot(xP, em_out[, "ExpD_f(x)"], type = "l", col = "blue", lwd = 2,
     main = "1D Emulator Output", xlab = "x1", ylab = "f(x1, 0.5)")

lines(xP, em_out[, "ExpD_f(x)"] + 2 * sqrt(em_out[, "VarD_f(x)"]), col = "red")
lines(xP, em_out[, "ExpD_f(x)"] - 2 * sqrt(em_out[, "VarD_f(x)"]), col = "red")
points(xD_1D, D_1D, col = "green", pch = 19, cex = 1.2)

true_vals <- apply(xP, 1, function(x) braninsc(c(x, 0.5)))
lines(xP, true_vals, col = "black", lwd = 2)

legend("topright", legend = c("Emulator Expectation", "Prediction Interval", "Model Evaluation", "True Function"),
       col = c("blue", "red", "green", "black"), lty = c(1, 1, NA, 1), pch = c(NA, NA, 19, NA), lwd = c(2, 2, NA, 2),
       inset = c(0.02, 0.02))

# 2D Emulation
set.seed(1)
xD_w1 <- randomLHS(14, 2)
D <- apply(xD_w1, 1, braninsc)

x_grid <- seq(0, 1, length.out = 50)
xP <- as.matrix(expand.grid(x1 = x_grid, x2 = x_grid))

em_out <- simple_BL_emulator_v3(xP = xP, xD = xD_w1, D = D, theta = c(0.3, 0.3), sigma = 1, E_f = 0)
E_D_fx_mat <- matrix(em_out[, "ExpD_f(x)"], nrow = 50)
Var_D_fx_mat <- matrix(em_out[, "VarD_f(x)"], nrow = 50)

### define filled contour plot function for emulator output ###
emul_fill_cont_V2 <- function(
    cont_mat,            # matrix of values we want contour plot of 
    cont_levs=NULL,      # contour levels (NULL: automatic selection)
    cont_levs_lines=NULL,   # contour levels for lines (NULL: automatic selection)
    nlev=20,             # approx no. of contour levels for auto select  
    plot_xD=TRUE,        # plot the design runs TRUE or FALSE
    xD=NULL,             # the design points if needed
    xD_col="green",      # colour of design runs
    x_grid,              # grid edge locations that define xP
    ...                  # extra arguments passed to filled.contour
){
  
  ### Define contour levels if necessary ###
  if(is.null(cont_levs)) cont_levs <- pretty(cont_mat,n=nlev)
  
  ### create the filled contour plot ###
  filled.contour(x_grid,x_grid,cont_mat,levels=cont_levs,xlab="x1",ylab="x2",...,  
                 plot.axes={axis(1);axis(2)                 # sets up plotting in contour box
                   if(is.null(cont_levs_lines)) contour(x_grid,x_grid,cont_mat,add=TRUE,levels=cont_levs,lwd=0.8) # plot usual contour lines 
                   if(!is.null(cont_levs_lines)) {
                     contour(x_grid,x_grid,cont_mat,add=TRUE,levels=cont_levs,lwd=0.4,labels="")   # plot thin contour lines 
                     contour(x_grid,x_grid,cont_mat,add=TRUE,levels=cont_levs_lines,lwd=2)   # plot thick contour lines
                   }
                   if(plot_xD) points(xD,pch=21,col=1,bg=xD_col,cex=1.5)})  # plot design points
}
### end define filled contour plot function for emulator output ###

exp_cols <- magma
var_cols <- function(n) hcl.colors(n, "YlOrRd", rev=TRUE)
diag_cols <- turbo

D_grid <- c(0.2,0.5,0.8)
xD <- as.matrix(expand.grid("x1"=D_grid,"x2"=D_grid)) #this expands D as a 4x4 grid 
xD

## 2D Emulation Comparison
# Evaluate true function and store in matrix for diagnostic comparisons
fxP_mat <- matrix(apply(xP, 1, function(x) braninsc(x)), 
                  nrow = length(x_grid), ncol = length(x_grid)) 

emul_fill_cont_V2(cont_mat=fxP_mat,cont_levs=seq(-2,2,0.2),xD=xD,x_grid=x_grid,
                  main="True Computer Model Function f(x)")

## 2D Emulation Expectation 
emul_fill_cont_V2(cont_mat=E_D_fx_mat,cont_levs=seq(-2,2,0.2),xD=xD,x_grid=x_grid,
                  color.palette=exp_cols, main="Emulator Adjusted Expectation (3x3 Grid)")

## 2D Emulation Variance
emul_fill_cont_V2(cont_mat=Var_D_fx_mat,cont_levs=seq(-2,2,0.2),xD=xD,x_grid=x_grid,
                  color.palette=var_cols, main="Emulator Adjusted Variance (3x3 Grid)")

## 2D Emulation Diagnostic

# Evaluate Emulator diagnostics S_D(x) and store in matrix
S_diag_mat <- (E_D_fx_mat - fxP_mat) / sqrt(Var_D_fx_mat)

emul_fill_cont_V2(cont_mat=S_diag_mat,cont_levs=seq(-3,3,0.25),xD=xD,x_grid=x_grid,
               xD_col="purple",
               color.palette=diag_cols,
               main="Emulator Diagnostics S_D[f(x)]")


## Multi-Wave History Matching
lhd_maximin <- function(nl=16){                    # nl = number of points in LHD 
  
  x_lhd <- cbind("x1"=sample(0:(nl-1)),"x2"=sample(0:(nl-1))) / nl  +  0.5/nl  # create LHD
  
  ### Maximin loop: performs swaps on 1st of two closest points with another random point
  for(i in 1:1000){
    mat <- as.matrix(dist(x_lhd)) + diag(10,nl) # creates matrix of distances between points
    # note the inflated diagonal 
    closest_runs <- which(mat==min(mat),arr.ind=TRUE)   # finds pairs of closest runs
    ind <- closest_runs[sample(nrow(closest_runs),1),1] # chooses one of close runs at random
    swap_ind <- sample(setdiff(1:nl,ind),1)       # randomly selects another run to swap with
    x_lhd2 <- x_lhd                               # creates second version of LHD
    x_lhd2[ind[1],1]   <- x_lhd[swap_ind,1] # swaps x_1 vals between 1st close run & other run
    x_lhd2[swap_ind,1] <- x_lhd[ind[1],1]   # swaps x_1 vals between 1st close run & other run
    if(min(dist(x_lhd2)) >= min(dist(x_lhd))-0.00001) {  # if min distance between points is same or better
      x_lhd <- x_lhd2                                    # we replace LHD with new LHD with the swap
      # cat("min dist =",min(dist(x_lhd)),"Iteration = ",i,"\n") # write out min dist 
    }
  }
  
  ### plot maximin LHD ###
  plot(x_lhd,xlim=c(0,1),ylim=c(0,1),pch=16,xaxs="i",yaxs="i",col="blue",
       xlab="x1",ylab="x2",cex=1.4)
  abline(h=(0:nl)/nl,col="grey60")
  abline(v=(0:nl)/nl,col="grey60")
  return(x_lhd)
}


### Define actual 2D compute model/simulator ###
### Define wave 1 run locations ###
set.seed(1)
xD_w1 <- lhd_maximin(nl=14)

### Define 50x50 grid of prediction points xP of emulator evaluation ###
x_grid <- seq(0, 1, length.out = 50)
xP <- as.matrix(expand.grid("x1" = x_grid, "x2" = x_grid))

### Defined Extra Objects for HM and Implausibility ###
z <- braninsc(c(0.8, 0.8))
sigma_e <- 0.06 
sigma_epsilon <- 0.05

### Define current run locations ###
xD <- xD_w1

### Perform 14 runs of model and store as D (this would take days for realistic example!) ###
D <- apply(xD, 1, braninsc)  

### Evaluate emulator over 50x50=2500 prediction points xP and store as matrices ###
em_out <- t(apply(xP, 1, simple_BL_emulator_v2, xD = xD, D = D, theta = 0.4, sigma = 1, E_f = 0))   
E_D_fx_mat <- matrix(em_out[, "ExpD_f(x)"], nrow = length(x_grid), ncol = length(x_grid)) 
Var_D_fx_mat <- matrix(em_out[, "VarD_f(x)"], nrow = length(x_grid), ncol = length(x_grid)) 

### Calculate Implausibility Measure Over All 50x50 = 2500 input points in xP ###
Imp_mat <- sqrt( (E_D_fx_mat - z)^2 / (Var_D_fx_mat + sigma_e^2 + sigma_epsilon^2) ) 

### Define colours and levels for implausibility plots ###
imp_cols <- function(n) turbo(n, begin = 0.15, end = 1)
imp_levs <- c(0, seq(1, 2.75, 0.25), seq(3, 18, 2), 20, 30, 41)

### plot wave 1 implausibility and wave 1 runs only ###
emul_fill_cont_V2(cont_mat = Imp_mat, cont_levs = imp_levs, cont_levs_lines = 3, xD = xD, x_grid = x_grid,
                  xD_col = "purple", color.palette = imp_cols, main = "Implausibility I(x): Wave 1")

### plot wave 1 emulator expectation ###
emul_fill_cont_V2(cont_mat = E_D_fx_mat, cont_levs = NULL, xD = xD, x_grid = x_grid,
                  color.palette = exp_cols, main = "Emulator Adjusted Expectation E_D[f(x)]")

### plot wave 1 emulator variance ###
emul_fill_cont_V2(cont_mat = Var_D_fx_mat, cont_levs = NULL, xD = xD, x_grid = x_grid,
                  color.palette = var_cols, main = "Emulator Adjusted Variance Var_D[f(x)]")

xD_w2 <- matrix(     
  c(0.05,0.3,
    0.1,0.25,
    0.1,0.1,
    0.05,0.05,
    0.5,0.9,
    0.6,0.9,
    0.7,0.8,
    0.8,0.81,
    0.92,0.9,
    0.95,0.95
    ), ncol = 2, byrow = TRUE
)

### plot current runs in purple and remaining unevaluated wave 2 runs in pink ###
emul_fill_cont_V2(cont_mat = Imp_mat, cont_levs = imp_levs, cont_levs_lines = 3, xD = rbind(xD, xD_w2),
                  x_grid = x_grid, xD_col = rep(c("purple", "pink"), c(nrow(xD), nrow(xD_w2))),
                  color.palette = imp_cols, main = "Implausibility I(x): Wave 1")

### loop over adding the wave 2 runs: add k runs ###
for (k in 0:10) {
  xD <- xD_w1
  if (k > 0) xD <- rbind(xD, xD_w2[1:k,])

  ### Perform 14 + k runs of model and store as D (this would take days for realistic example!) ###
  D <- apply(xD, 1, braninsc)

  ### Evaluate emulator over 50x50=2500 prediction points xP and store as matrices ###
  em_out <- t(apply(xP, 1, simple_BL_emulator_v2, xD = xD, D = D, theta = 0.4, sigma = 1, E_f = 0))
  E_D_fx_mat <- matrix(em_out[, "ExpD_f(x)"], nrow = length(x_grid), ncol = length(x_grid))
  Var_D_fx_mat <- matrix(em_out[, "VarD_f(x)"], nrow = length(x_grid), ncol = length(x_grid))

  ### Evaluate true function and store in matrix for diagnostic comparisons ###
  fxP_mat <- matrix(apply(xP, 1, braninsc), nrow = length(x_grid), ncol = length(x_grid))

  ### Calculate Implausibility Measure Over All 50x50 = 2500 input points in xP ###
  Imp_mat <- sqrt( (E_D_fx_mat - z)^2 / (Var_D_fx_mat + sigma_e^2 + sigma_epsilon^2) )

  ### Calculate Imp Measure for True f(x) Over All 50x50 = 2500 input points in xP ###
  Imp_true_mat <- sqrt( (fxP_mat - z)^2 / (sigma_e^2 + sigma_epsilon^2) )

  ### Define colours and levels for implausibility plots ###
  imp_cols <- function(n) turbo(n, begin = 0.15, end = 1)
  imp_levs <- c(0, seq(1, 2.75, 0.25), seq(3, 18, 2), 20, 30, 41)

  ### if k=0 plot wave 1 runs only ###
  if (k == 0) emul_fill_cont_V2(cont_mat = E_D_fx_mat, cont_levs = imp_levs, cont_levs_lines = 3, xD = xD,
                                x_grid = x_grid, xD_col = "purple", color.palette = imp_cols,
                                main = "E_D_fx: Wave 1")

  ### plot current runs in purple and remaining unevaluated wave 2 runs in pink ###
  emul_fill_cont_V2(cont_mat = E_D_fx_mat, cont_levs = imp_levs, cont_levs_lines = 3, xD = rbind(xD, xD_w2),
                    x_grid = x_grid, xD_col = rep(c("purple", "pink"), c(nrow(xD) + k, nrow(xD_w2) - k)),
                    color.palette = imp_cols, main = "Implausibility I(x): Wave 2")

  ### once last run done so k=8, plot implausibility for true function f(x) to compare ###
  if (k == nrow(xD_w2)) emul_fill_cont_V2(cont_mat = Imp_true_mat, cont_levs = imp_levs,
                                          cont_levs_lines = 3, xD = xD, x_grid = x_grid, xD_col = "purple", plot_xD = FALSE,
                                          color.palette = imp_cols, main = "Implausibility I(x) using True f(x)")
}

### Define current run locations ###
xD <- rbind(xD_w1, xD_w2)

### Perform 14 runs of model and store as D (this would take days for realistic example!) ###
D <- apply(xD, 1, braninsc)  

### Evaluate emulator over 50x50=2500 prediction points xP and store as matrices ###
em_out <- t(apply(xP, 1, simple_BL_emulator_v2, xD = xD, D = D, theta = 0.4, sigma = 1, E_f = 0))   
E_D_fx_mat <- matrix(em_out[, "ExpD_f(x)"], nrow = length(x_grid), ncol = length(x_grid)) 
Var_D_fx_mat <- matrix(em_out[, "VarD_f(x)"], nrow = length(x_grid), ncol = length(x_grid)) 

### Calculate Implausibility Measure Over All 50x50 = 2500 input points in xP ###
Imp_mat2 <- sqrt( (E_D_fx_mat - z)^2 / (Var_D_fx_mat + sigma_e^2 + sigma_epsilon^2) ) 

### Define colours and levels for implausibility plots ###
imp_cols <- function(n) turbo(n, begin = 0.15, end = 1)
imp_levs <- c(0, seq(1, 2.75, 0.25), seq(3, 18, 2), 20, 30, 41)

### plot wave 2 implausibility and wave 2 runs ###
emul_fill_cont_V2(cont_mat = Imp_mat2, cont_levs = imp_levs, cont_levs_lines = 3, xD = xD, x_grid = x_grid,
                  xD_col = "purple", color.palette = imp_cols, main = "Implausibility I(x): Wave 2")

