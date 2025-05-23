### Define simple Bayes Linear emulator for single input ###
simple_BL_emulator_v1 <- function(x,              # the emulator prediction point
                                  xD,             # the run input locations xD
                                  D,              # the run outputs D = (f(x^1),...,f(x^n))
                                  theta = 1,      # the correlation lengths
                                  sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])
                                  E_f = 0         # prior expectation of f: E(f(x)) = 0 
){
  
  # store length of runs D  
  n <- length(D)
  
  ### Define Covariance structure of f(x): Cov[f(x),f(xdash)] ###
  Cov_fx_fxdash <- function(x,xdash) sigma^2 * exp(-(x-xdash)^2/theta^2) 
  
  
  ### Define 5 objects needed for BL adjustment ###
  # Create E[D] vector
  E_D <- rep(E_f,n)
  
  # Create Var_D matrix:
  Var_D <- matrix(0,nrow=n,ncol=n)
  for(i in 1:n) for(j in 1:n) Var_D[i,j] <- Cov_fx_fxdash(xD[i],xD[j])
  
  # Create E[f(x)]
  E_fx <- E_f
  
  # Create Var_f(x) 
  Var_fx <- sigma^2
  
  # Create Cov_fx_D row vector
  Cov_fx_D <- matrix(0,nrow=1,ncol=n)
  for(j in 1:n) Cov_fx_D[1,j] <- Cov_fx_fxdash(x,xD[j])
  
  
  ### Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x) ###
  ED_fx   <-  E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)   
  VarD_fx <-  Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)  
  
  ### return emulator expectation and variance ###
  return(c("ExpD_f(x)"=ED_fx,"VarD_f(x)"=VarD_fx))  
  
}
### End simple Bayes Linear emulator for single input ###
##############################################################################################


##############################################################################################
### Function to plot simple emulator output
plot_BL_emulator_V2 <- function(
    em_out,         # a nP x 2 matrix of emulator outputs
    xP,             # vector of nP inputs where emulator evaluated
    xD,             # the run input locations
    D,              # the run outputs D = (f(x^1),...,f(x^n))
    maintitle=NULL,  # title of the plot, blank if NULL
    z=NULL,         # the observed data z
    sigma_e=NULL,   # the observation errors SD s.t. Var[e] = sigma_e^2
    sigma_epsilon=NULL, # the model discrepancy SD s.t. Var[epsilon] = sigma_epsilon^2
    plot_true=FALSE   # don't plot true function unless this is TRUE
){
  
  ### plot emulator output ###
  plot(xP,em_out[,"ExpD_f(x)"],ylim=c(-1.3,1.2),ty="l",col="blue",lwd=2.5,
       xlab="Input parameter x",ylab="Output f(x)",main=maintitle) # main argument: plot title
  lines(xP,em_out[,"ExpD_f(x)"]+3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)
  lines(xP,em_out[,"ExpD_f(x)"]-3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)
  
  ### plot true function: we would not normally be able to do this! ###
  if(plot_true) lines(xP,f(xP),lwd=2,lty=1)
  
  ### Plot the runs ###
  points(xD,D,pch=21,col=1,bg="green",cex=1.5)
  
  ### Plot the Observed data plus errors due to obs and MD ###
  if(!is.null(z)){
    abline(h=z,lwd=1.4)
    abline(h=z+3*sqrt(sigma_e^2+ sigma_epsilon^2),lty=2,lwd=1.2)
    abline(h=z-3*sqrt(sigma_e^2+ sigma_epsilon^2),lty=2,lwd=1.2)
  }
  
  if(plot_true) legend('topright',legend=c("Emulator Expectation",
                                           "Emulator Prediction Interval",
                                           "True function f(x)","Model Evaluations"), 
                       lty=c(1,1,1,NA),pch=c(NA,NA,NA,16),col=c("blue","red",1,"green"),lwd=2.5,pt.cex=1.3)
  if(!is.null(z)) legend('topright',legend=c("Emulator Expectation",
                                             "Emulator Prediction Interval",
                                             "Model Evaluations","Observation z",
                                             "3 sigma interval"), 
                         lty=c(1,1,NA,1,2),pch=c(NA,NA,16,NA,NA),col=c("blue","red","green",1,1),
                         lwd=c(2.5,2.5,NA,1.4,1.2),pt.cex=1.3)
}
### Function to plot simple emulator output
##############################################################################################

library(viridisLite)

##############################################################################################
### 1D HM example wave 1

### Define actual computer model/simulator ###
f <- function(x) sin(2*pi*x)

### Evaluate emulator over 201 prediction points xP ###
xP <- seq(0.001,0.999,len=201)

### Define run locations ###
xD <- c(0,0.2,0.4,0.6,0.86,1)     # note shifted x^(5) location

### Perform 6 runs of model and store as D (this would takes days for realistic example!) ###
D <- f(xD)

### Evaluate emulator over 201 prediction points xP ###
em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))  

### Defined extra objects for HM and implausibility ###
z <- -0.8
sigma_e <- 0.01
sigma_epsilon <- 0.02

### Plot emulator output with observation errors ###
plot_BL_emulator_V2(em_out=em_out,xP=xP,xD=xD,D=D,maintitle="Wave 1 Emulator Output: 6 runs",
                    z=z,sigma_e=sigma_e,sigma_epsilon=sigma_epsilon)


### Calculate Implausibility Measure Over All 201 input points in xP ###
Implaus <- sqrt( (em_out[,"ExpD_f(x)"] - z)^2 / ( em_out[,"VarD_f(x)"] + sigma_e^2 + sigma_epsilon^2) ) 


### Plot Implausibility function ###
plot_1D_Implausibility_V1 <- function(
    Implaus,        # the vector of implausibilities to plot
    xP,             # the vector of input points where the implausibilities were evaluated
    imp_cutoff=3,   # implausibility cutoff
    maintitle=NULL  # title of the plot, blank if NULL
){
  
  imp_breaks <- c(seq(0,2.8,0.2),seq(3,9.5,0.5),seq(10,80,10))   # colour steps for imp
  
  imp_cut_index <- cut(Implaus,imp_breaks,labels=FALSE) # which imp intervals each pt lies in
  cols <- turbo(length(imp_breaks)-1,begin=0.15,end=1)  # assign reduced "turbo" colours
  
  plot(xP,Implaus,ylim=c(0,30),ty="l",lwd=2,       # plot imp with optional main title
       xlab="Input parameter x",ylab="Implausibility I(x)",main=maintitle)
  abline(h=imp_cutoff,col="green",lwd=2)    # draw horizontal implausibility cutoff
  rect_wid <- 1/(length(xP)-1)       # set up simple rectangles to colour x-axis
  for(i in 1:length(xP)) rect(xP[i]-rect_wid/2,-1,xP[i]+rect_wid/2,0,col=cols[imp_cut_index][i],border=NA)
  abline(h=0,col=1,lwd=1)     # additional black line for clarity
  
  legend('topright',legend=c("Implausibility I(x)","Implausibility Cutoff c"), 
         lty=c(1,1),col=c(1,"green"),lwd=2)   # legend
}

### Plot the implausibility ###

xP <- seq(0.001,0.999,len=50)
xD <- c(0,0.1,0.3,0.5,0.9,1)
D <- f(xD)
em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))
Implaus <- sqrt( (em_out[,"ExpD_f(x)"] - z)^2 / ( em_out[,"VarD_f(x)"] + sigma_e^2 + sigma_epsilon^2) ) 
plot_1D_Implausibility_V1(Implaus=Implaus,xP=xP)


##############################################################################################
### 1D HM example wave 2: add three points sequentially ###

for(k in c(1:2)){        # k=1: plot the emulator, k=2 plot the implausibility
  for(j in 0:3){                     # if j==0 do wave 1, if j>0 add j wave 2 points
    
    ### Define run locations ###
    xD <- c(0,0.2,0.4,0.6,0.86,1)    # the wave 1 inputs
    xD_w2 <- c(0.73,0.8,0.66)        # the wave 2 inputs
    if(j>0) xD <- c(xD,xD_w2[1:j])   # if j==0 do wave 1, if j>0 add j wave 2 points
    
    ### Perform 6 runs and store as D (this would takes days for realistic example!) ###
    D <- f(xD)
    
    ### Evaluate emulator over 201 prediction points xP ###
    em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))  
    
    ### Plot emulator output with observation errors ###
    if(k==1) plot_BL_emulator_V2(em_out=em_out,xP=xP,xD=xD,D=D,
                                 maintitle=paste("Wave",c(1,2,2,2)[j+1],"Emulator Output:",length(xD),"runs"),
                                 z=z,sigma_e=sigma_e,sigma_epsilon=sigma_epsilon)
    
    ### Calculate Implausibility Measure Over All 201 input points in xP ###
    if(k==2) Implaus <- sqrt( (em_out[,"ExpD_f(x)"] - z)^2 / 
                                ( em_out[,"VarD_f(x)"] + sigma_e^2 + sigma_epsilon^2) ) 
    
    ### Plot the implausbility ###
    if(k==2) plot_1D_Implausibility_V1(Implaus=Implaus,xP=xP,
                                       maintitle=paste("Wave",c(1,2,2,2)[j+1],"Implausibility:",length(xD),"runs"))
    
  }
}


##############################################################################################
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

### Define colour schemes for standard plots ###
exp_cols <- magma
var_cols <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
diag_cols <- turbo

##############################################################################################
### Create Type 1 LHD ###

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

### Create Type 1 LHD ###
##############################################################################################

lhd_maximin(nl=14)

set.seed(1)
lhd_maximin(nl=14)

#####################################################################################
### 2D input HM example using 2D example from before

### Define actual 2D compute model/simulator ###
f <- function(x) -sin(2*pi*x[,2]) + 0.9*sin(2*pi*(1-x[,1])*(1-x[,2]))

### Define wave 1 run locations ###
set.seed(1)
xD_w1 <-  lhd_maximin(nl=14)

### Define 50x50 grid of prediction points xP of emulator evaluation ###
x_grid <- seq(-0.001,1.001, len =50)
xP <- as.matrix(expand.grid("x1" = x_grid, "x2" = x_grid))
### Defined Extra Objects for HM and Implausibility ###
z <- -1.25
sigma_e <- 0.06 
sigma_epsilon <- 0.05

### Define current run locations ###
xD <- xD_w1

### Perform 14 runs of model and store as D (this would take days for realistic example!) ###
D <- f(xD)

### Evaluate emulator over 50x50=2500 prediction points xP and store as matrices ###
em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=0.4,sigma=1,E_f=0))   
E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 


### Calculate Implausibility Measure Over All 50x50 = 2500 input points in xP ###
Imp_mat <- sqrt( (E_D_fx_mat - z)^2 / (Var_D_fx_mat + sigma_e^2 + sigma_epsilon^2) ) 

### Define colours and levels for implausibility plots ###
imp_cols <- function(n) turbo(n,begin=0.15,end=1)
imp_levs <- c(0,seq(1,2.75,0.25),seq(3,18,2),20,30,41)

### plot wave 1 implausibility and wave 1 runs only ###
emul_fill_cont_V2(cont_mat=Imp_mat,cont_levs=imp_levs,cont_levs_lines=3,xD=xD,x_grid=x_grid,
                  xD_col="purple",color.palette=imp_cols,main="Implausibility I(x): Wave 1")

### plot wave 1 emulator expectation ###
emul_fill_cont_V2(cont_mat=E_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                  color.palette=exp_cols,main="Emulator Adjusted Expectation E_D[f(x)]")

### plot wave 1 emulator variance ###
emul_fill_cont_V2(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                  color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]")


xD_w2 <- matrix(     # the 8 point wave 2 design (chosen by hand)
  c(0.62,0.99,
    0.98,0.28,
    0.95,0.06,
    0.04,0.05,
    0.28,0.2,
    0.04,0.4,
    0.35,0.1,
    0.22,0.36),ncol=2,byrow=TRUE
)

### plot current runs in purple and remaining unevaluated wave 2 runs in pink ###
emul_fill_cont_V2(cont_mat=Imp_mat,cont_levs=imp_levs,cont_levs_lines=3,xD=rbind(xD,xD_w2),
                  x_grid=x_grid,xD_col=rep(c("purple","pink"),c(nrow(xD),nrow(xD_w2))),
                  color.palette=imp_cols,main="Implausibility I(x): Wave 1")

### loop over adding the wave 2 runs: add k runs ###
for(k in 0:8){                          # k=0: wave 1, k>0 add k wave 2 runs sequentially
  
  xD <- xD_w1                           # the 14 point wave 1 design
  if(k>0) xD <- rbind(xD,xD_w2[1:k,])   # k=0: wave 1, k>0 add wave 2 runs sequentially
  
  ### Perform 14 + k runs of model and store as D (would take days for realistic example!) ###
  D <- f(xD)
  
  ### Evaluate emulator over 50x50=2500 prediction points xP and store as matrices ###
  em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=0.4,sigma=1,E_f=0))   
  E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  ### Evaluate true function and store in matrix for diagnostic comparisons ###
  fxP_mat <- matrix(f(xP),nrow=length(x_grid),ncol=length(x_grid)) 
  
  ### Calculate Implausibility Measure Over All 50x50 = 2500 input points in xP ###
  Imp_mat <- sqrt( (E_D_fx_mat - z)^2 / (Var_D_fx_mat + sigma_e^2 + sigma_epsilon^2) ) 
  
  ### Calculate Imp Measure for True f(x) Over All 50x50 = 2500 input points in xP ###
  Imp_true_mat <- sqrt( (fxP_mat - z)^2 / (sigma_e^2 + sigma_epsilon^2) ) 
  
  ### Define colours and levels for implausibility plots ###
  imp_cols <- function(n) turbo(n,begin=0.15,end=1)
  imp_levs <- c(0,seq(1,2.75,0.25),seq(3,18,2),20,30,41)
  
  ### if k=0 plot wave 1 runs only ###
  if(k==0) emul_fill_cont_V2(cont_mat=E_D_fx_mat,cont_levs=imp_levs,cont_levs_lines=3,xD=xD,
                             x_grid=x_grid,xD_col="purple",color.palette=imp_cols,
                             main="E_D_fx: Wave 1")
  ### plot current runs in purple and remaining unevaluated wave 2 runs in pink ###
  emul_fill_cont_V2(cont_mat=E_D_fx_mat,cont_levs=imp_levs,cont_levs_lines=3,xD=rbind(xD,xD_w2),
                    x_grid=x_grid,xD_col=rep(c("purple","pink"),c(nrow(xD)+k,nrow(xD_w2)-k)),  # cover unevaluated w2 points in pink
                    color.palette=imp_cols,main="Implausibility I(x): Wave 2")
  ### once last run done so k=8, plot implausibility for true function f(x) to compare ###
  if(k==nrow(xD_w2)) emul_fill_cont_V2(cont_mat=Imp_true_mat,cont_levs=imp_levs,
                                       cont_levs_lines=3,xD=xD,x_grid=x_grid,xD_col="purple",plot_xD=FALSE,
                                       color.palette=imp_cols,main="Implausibility I(x) using True f(x)")
}
