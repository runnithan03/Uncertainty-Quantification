install.packages("plot3D")
install.packages("lhs")
install.packages("pdist")
install.packages("sensitivity")

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

library(plot3D)

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
               color.palette = colorRampPalette(c("blue", "lightblue", "white", "red")),
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
library(lhs)
library(pdist)
set.seed(1)
xD_w1 <- randomLHS(14, 2)
D <- apply(xD_w1, 1, braninsc)

x_grid <- seq(0, 1, length.out = 50)
xP <- as.matrix(expand.grid(x1 = x_grid, x2 = x_grid))

em_out <- simple_BL_emulator_v3(xP = xP, xD = xD_w1, D = D, theta = c(0.3, 0.3), sigma = 1, E_f = 0)
E_D_fx_mat <- matrix(em_out[, "ExpD_f(x)"], nrow = 50)
Var_D_fx_mat <- matrix(em_out[, "VarD_f(x)"], nrow = 50)


##############################################################################################
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
##############################################################################################

exp_cols <- magma
var_cols <- function(n) hcl.colors(n, "YlOrRd", rev=TRUE)

D_grid <- c(0.2,0.5,0.8)
xD <- as.matrix(expand.grid("x1"=D_grid,"x2"=D_grid)) #this expands D as a 4x4 grid 
xD

# 2D Emulation Expectation 
emul_fill_cont_V2(cont_mat=E_D_fx_mat,cont_levs=seq(-2,2,0.2),xD=xD,x_grid=x_grid,
               main="Emulator Adjusted Expectation (3x3 Grid)")

filled.contour(x_grid, x_grid, E_D_fx_mat,
               color.palette = colorRampPalette(c("navy", "skyblue", "lightcyan")),
               main = "2D Emulator Expectation",
               xlab = "x1", ylab = "x2")

# Branin-Hoo Function Contour Plot Overlay with 2D Emulation Expectation 
filled.contour(x1, x2, z_matrix,
               color.palette = colorRampPalette(c("blue", "lightblue", "white", "red")),
               main = "Overlay of Branin Function and Emulator Expectation",
               xlab = "x1", ylab = "x2",
               plot.axes = {
                 axis(1); axis(2) # Add axes
                 
                 # Contours of the emulator's expectation
                 contour(x_grid, x_grid, E_D_fx_mat, 
                         add = TRUE, col = "black", lty = 2, lwd = 2)
               })

# 2D Emulation Variance
emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=seq(-2,2,0.2),xD=xD,x_grid=x_grid,
               main="Emulator Adjusted Variance E_D[f(x)]")

filled.contour(x_grid, x_grid, Var_D_fx_mat,
               color.palette = colorRampPalette(c("navy", "skyblue", "lightcyan")),
               main = "2D Emulator Variance",
               xlab = "x1", ylab = "x2")

# Branin-Hoo Function Contour Plot Overlay with 2D Emulation Variance 
filled.contour(x1, x2, z_matrix,
               color.palette = colorRampPalette(c("blue", "lightblue", "white", "red")),
               main = "Overlay of Branin Function with Emulator Variance",
               xlab = "x1", ylab = "x2",
               plot.axes = {
                 axis(1); axis(2) # Add axes
                 
                 # Contours for emulator variance
                 contour(x_grid, x_grid, Var_D_fx_mat,
                         add = TRUE, col = "green", lty = 3, lwd = 1)
               })

# Further Techniques Part 1: Implausibility
z <- -0.5
sigma_e <- 0.05
sigma_epsilon <- 0.1

Imp_mat <- sqrt((E_D_fx_mat - z)^2 / (Var_D_fx_mat + sigma_e^2 + sigma_epsilon^2))

filled.contour(x_grid, x_grid, Imp_mat,
               color.palette = colorRampPalette(c("navy", "skyblue", "lightcyan")),
               main = "Implausibility I(x)",
               xlab = "x1", ylab = "x2")

non_implausible <- which(Imp_mat < 3, arr.ind = TRUE)

# Select new points from non-implausible region
new_points <- xP[sample(non_implausible, 5), ]

xD_w2 <- rbind(xD_w1, new_points)
D_w2 <- apply(xD_w2, 1, braninsc)

# Further Techniques Part 2: Multi-Wave History Matching
library(sensitivity)
library(lhs)

# Generate New Latin Hypercube Design for Additional Wave
xD_wave2 <- randomLHS(20, 2)  # 20 new points for next wave
D_wave2 <- apply(xD_wave2, 1, braninsc)

# Combine with Existing Design
xD_combined <- rbind(xD_w1, xD_wave2)
D_combined <- c(D, D_wave2)

# Re-run Emulator with New Data Points
em_out_wave2 <- simple_BL_emulator_v3(xP = xP, xD = xD_combined, D = D_combined, theta = c(0.3, 0.3), sigma = 1)

# Plot Wave 2 Emulator Expectation
E_D_fx_wave2 <- matrix(em_out_wave2[, "ExpD_f(x)"], nrow = 50)
filled.contour(x_grid, x_grid, E_D_fx_wave2,
               color.palette = colorRampPalette(c("navy", "skyblue", "lightcyan")),
               main = "Wave 2 Emulator Expectation",
               xlab = "x1", ylab = "x2")

# # Branin-Hoo Function Contour Plot Overlay with Wave 2 Emulator Expectation 
filled.contour(x1, x2, z_matrix,
               color.palette = colorRampPalette(c("blue", "lightblue", "white", "red")),
               main = "Overlay of Branin Function and Wave 2 Emulator Expectation",
               xlab = "x1", ylab = "x2",
               plot.axes = {
                 axis(1); axis(2) # Add axes
                 
                 # Add contours of the emulator's expectation
                 contour(x_grid, x_grid, E_D_fx_wave2, 
                         add = TRUE, col = "black", lty = 2, lwd = 2)
               })
