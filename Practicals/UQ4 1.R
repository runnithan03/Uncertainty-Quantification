### Define actual computer model/simulator ###
f <- function(x) sin(2*pi*x)

#Design a set of six runs location xD evenly spread over [0,1]
### Define run locations ###
xD <- seq(0,1,0.2)
xD
### Perform 6 runs of model and store as D (this would takes days for realistic example!) ###
D <- f(xD)
D

n <- length(D)
n

##Covariance structure
theta <- 0.25                  # the correlation lengths in emulator covariance structure
sigma <- 0.5                   # the prior SD sigma sqrt(Var[f(x)])


### Define Covariance structure of weakly stationary process ###
Cov_fx_fxdash <- function(x,xdash) sigma^2 * exp(-(x-xdash)^2/theta^2) # covariance structure of f(x): Cov[f(x),f(xdash)]. 

#prior expectation of f: E(f(x)) = 0 
E_f = 0 

### Define 5 objects needed for BL adjustment ###

# 1) Create E[D] vector 
E_D <- rep(E_f,n)
E_D

# 2) 6x6 Var[D] matrix
# Create Var_D matrix:
Var_D <- matrix(0,nrow=n,ncol=n)
for(i in 1:n) for(j in 1:n) Var_D[i,j] <- Cov_fx_fxdash(xD[i],xD[j])
Var_D

# 3)
# Create E[f(x)] 
E_fx <- E_f
E_fx

# 4) Create Var_f(x)
Var_fx <- sigma^2
Var_fx

# 5) Decide on which input x point to evaluate the emulator at. 
# Then construct Cov_fx_D row vector
x <- 0.25
Cov_fx_D <- matrix(0,nrow=1,ncol=n)
for(j in 1:n) Cov_fx_D[1,j] <- Cov_fx_fxdash(x,xD[j])

### Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x) ###
# the BL adjusted expectation update equations for emulator
ED_fx  <- E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)        
# the BL adjusted variance update equations for emulator
VarD_fx <- Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)  
#solve is the inverse matrix

ED_fx
VarD_fx




##############################################################################################
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
### Define run locations and perform 6 runs of model ###
xD <- seq(0,1,0.2)
D <- f(xD)

em_out_1pt <- simple_BL_emulator_v1(x=0.25,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0)
em_out_1pt

#3 sigma prediction interval
### Prediction Interval at x=0.25 ###
low1 <- em_out_1pt["ExpD_f(x)"] - 3 * sqrt(em_out_1pt["VarD_f(x)"])
up1  <- em_out_1pt["ExpD_f(x)"] + 3 * sqrt(em_out_1pt["VarD_f(x)"])
c(low1,up1,use.names=FALSE)

### Evaluate emulator over 201 prediction points xP ###
xP <- seq(0.001,0.999,len=201)
nP <- length(xP)      # number of prediction points
em_out <- matrix(0,nrow=nP,ncol=2,dimnames=list(NULL,c("ExpD_f(x)","VarD_f(x)")))    
# set up matrix for results
for(i in 1:nP) em_out[i,]<-simple_BL_emulator_v1(x=xP[i],xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0)
head(em_out)          # look at first few rows of a matrix


em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))  
# t(): gives transpose of the matrix
head(em_out)


########### plot emulator output #############
pdf(file="fig1A_sin_emulator.pdf",height=7,width=8.5)    # uncomment if you want a pdf file
plot(xP,em_out[,"ExpD_f(x)"],ylim=c(-1.2,1.2),ty="l",col="blue",lwd=2.5,xlab="Input parameter x",ylab="Output f(x)")
lines(xP,em_out[,"ExpD_f(x)"]+3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)
lines(xP,em_out[,"ExpD_f(x)"]-3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)

### plot true function: we would not normally be able to do this! ###
lines(xP,f(xP),lwd=2,lty=1)

### plot expectation and prediction interval points at x=0.25 ###
abline(v=0.25,lty=2)
points(0.25,em_out_1pt["ExpD_f(x)"],col=1,bg="blue",pch=21)
points(c(0.25,0.25),c(low1,up1),col=1,bg="red",pch=21)

### plot the runs ###
points(xD,D,pch=21,col=1,bg="green",cex=1.5)

### add a legend to label everything ###
legend('topright',legend=c("Emulator Expectation","Emulator Prediction Interval",
                           "True function f(x)","Model Evaluations","x = 0.25 line"),
       lty=c(1,1,1,NA,2),pch=c(NA,NA,NA,16,NA),col=c("blue","red",1,"green",1),
       lwd=c(2.5,2.5,2.5,2.5,1),pt.cex=1.3)
dev.off()       # uncomment if you want to make a pdf file

##############################################################################################
### Function to plot simple emulator output
plot_BL_emulator_V1 <- function(
    em_out,         # a nP x 2 matrix of emulator outputs
    xP,             # vector of nP inputs where emulator evaluated
    xD,             # the run input locations
    D,              # the run outputs D = (f(x^1),...,f(x^n))
    maintitle=NULL  # title of the plot, blank if NULL
){
  
  ### plot emulator output ###
  plot(xP,em_out[,"ExpD_f(x)"],ylim=c(-1.2,1.2),ty="l",col="blue",lwd=2.5,
       xlab="Input parameter x",ylab="Output f(x)",main=maintitle) # main argument: plot title
  lines(xP,em_out[,"ExpD_f(x)"]+3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)
  lines(xP,em_out[,"ExpD_f(x)"]-3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)
  
  ### plot true function: we would not normally be able to do this! ###
  lines(xP,f(xP),lwd=2,lty=1)
  
  ### Plot the runs ###
  points(xD,D,pch=21,col=1,bg="green",cex=1.5)
  legend('topright',legend=c("Emulator Expectation","Emulator Prediction Interval",
                             "True function f(x)","Model Evaluations"),
         lty=c(1,1,1,NA),pch=c(NA,NA,NA,16),col=c("blue","red",1,"green"),lwd=2.5,pt.cex=1.3)
}
### Function to plot simple emulator output
##############################################################################################
plot_BL_emulator_V1(em_out,xP,xD,D,maintitle="Emulator Output")


#We add two more runs (0.7 and 0.9)
### Define run locations ###
xD <- c(seq(0,1,0.2),0.7,0.9)    # original 6 runs and two extra runs at x=0.7 and x=0.9

### Perform 8 runs of model and store as D (this would takes days for realistic example!) ###
D <- f(xD)

### Evaluate emulator over 201 prediction points xP ###
em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))  

### Plot emulator output ###
plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle="Emulator Output: 8 runs")

#We see the emulator is much more precise for low values of f(x)

### Define run locations ###
xD <- seq(0,1,0.2)    # original 6 runs

### Perform 6 runs of model and store as D (this would takes days for realistic example!) ###
D <- f(xD)

#Now we vary the prior vairance parameter sigma^2
## make sequence of sigma values to use for emulation ###
sigma_seq <- c(0.25,0.5,0.75,1,1.5)

### for loop over different sigma values in sigma_seq ###
for(i in 1:length(sigma_seq)){
  ### Evaluate emulator over 201 prediction points xP with sigma=sigma_seq[i] ###
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=sigma_seq[i],E_f=0))
  
  ### Plot emulator output in each case, note use of "paste" for plot title ###
  plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle=paste("Sigma =",sigma_seq[i]))
}


#Exercise 5.1 varying the length parameter

xD <- seq(0,1,0.2)
D <- f(xD)
theta_seq <- seq(0.0005,1,length.out =6)
for(i in 1:length(theta_seq)){ 
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=theta_seq,sigma=0.25,E_f=0))
  plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle=paste("Theta =",theta_seq[i]))
}
theta_seq  
