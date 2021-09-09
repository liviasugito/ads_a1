# installation of required packages
install.packages("reshape")
install.packages("SAVER")
install.packages("matrixcalc")
install.packages("matlib")
install.packages("plot.matrix")
install.packages("corrplot")

# Q 1.1

AV <- c(0, 20, 0, 0, 0, 0)
IV <- c(30, 45, 60, 40, 40, 40)
dur <- c(15, 20, 25, 15, 20, 25)
N <- 240
TC <- matrix(0, N, 6, byrow=FALSE)

par(mfrow=c(2,3)) 
for (col in seq(1, 6)){
  for (i in seq((AV[col]+1), N, IV[col])){
    TC[i:(i+dur[col]-1), col] <- rep(1, dur[col])
  }
  # standardise each column and plot TC
  TC[,col] <- scale(TC[,col])
  plot(TC[,col], type="l", main=sprintf("TC %i", col), xlab="", ylab="")
  
}


# Q 1.2

# correlation plot for TC
par(mfrow=c(1,1))
library(corrplot)
(CM <- cor(TC))
corrplot(CM, method="color", outline="black")

# Q 1.3
library('plot.matrix')

# setting up the slices for each SM
SM_1 <- matrix(0, 21, 21, byrow=FALSE)
SM_1[2:6,2:6] <- 1

SM_2 <- matrix(0, 21, 21, byrow=FALSE)
SM_2[2:6,15:19] <- 1

SM_3 <- matrix(0, 21, 21, byrow=FALSE)
SM_3[8:13,2:6] <- 1

SM_4 <- matrix(0, 21, 21, byrow=FALSE)
SM_4[8:13,15:19] <- 1

SM_5 <- matrix(0, 21, 21, byrow=FALSE)
SM_5[15:19,2:6] <- 1

SM_6 <- matrix(0, 21, 21, byrow=FALSE)
SM_6[15:19,15:19] <- 1


# array of SM
tmpSM <- array(c(SM_1, SM_2, SM_3, SM_4, SM_5, SM_6), dim=c(21 ,21 ,6))

# plotting SM
par(mar = c(5, 5, 5, 5), mfrow=c(2,3))
for(k in 1:6){
  plot(tmpSM[,,k], border=NA, main=sprintf("SM %i", k))
}

# reshape into 6 x 441 matrix
SM <- rbind(as.vector(tmpSM[,,1]), as.vector(tmpSM[,,2]), 
           as.vector(tmpSM[,,3]), as.vector(tmpSM[,,4]),
           as.vector(tmpSM[,,5]), as.vector(tmpSM[,,6]))

# plot the correlation plot for SM
par(mfrow=c(1,1))
cor(t(SM))
corrplot(cor(t(SM)), method="color")

# Q 1.4

Tt <- matrix(rnorm(240*6, mean=0, sd=0.5), 240, 6)
Ts <- matrix(rnorm(6*441, mean=0, sd=sqrt(0.015)), 6, 441)
Tt_Ts <- Tt%*%Ts

# plotting Tt times Ts (only the first ten)
par(mfrow=c(1,1))
corrplot(cor(Tt_Ts[1:10,1:10]), method="color", 
         title="Tt Ts", mar=c(0,0,1,0))

# Correlation Matrix for Tt and Ts
par(mar=c(0,0,1,0), mfrow=c(1,2))
corrplot(cor(Tt), method="color", 
         title="Temporal Source", mar=c(0,0,1,0))
corrplot(cor(t(Ts)), method="color", 
         title="Spatial Source", mar=c(0,0,1,0))

# Plot histogram for Tt and Ts
par(mar=c(4,4,2,2), mfrow=c(1,2))
hist(c(Tt),col="peachpuff", border="black", 
     main="Histogram for Temporal Source White Noise", 
     cex.main=1.2)
abline(v=mean(Tt), col="red")
abline(v=mean(Tt)+(1.96*sqrt(0.25)), col="blue")
abline(v=mean(Tt)-(1.96*sqrt(0.25)), col="blue")
abline(v=quantile(Tt, probs=0.025), col="green")
abline(v=quantile(Tt, probs=0.975), col="green")

hist(c(Ts),col="peachpuff", border="black", 
     main="Histogram for Spatial Source White Noise", 
     cex.main=1.2)
abline(v=mean(Ts), col="red")
abline(v=mean(Ts)+(1.96*sqrt(0.015)), col="blue")
abline(v=mean(Ts)-(1.96*sqrt(0.015)), col="blue")
abline(v=quantile(Ts, probs=0.025), col="green")
abline(v=quantile(Ts, probs=0.975), col="green")


# Q 1.5
X <- (TC + Tt) %*% (SM + Ts)

# plot 100 randomly selected X
X_sample <- X[,sample(ncol(X),size=100)]
dim(X_sample)
par(mfrow=c(1,1))
matplot(X_sample, lty=1, type="l", col=rainbow(11))

# plot the variance of 441 variables
var_X <- apply(X, 2, var)
plot(var_X)

# standardise X
scaled_X <- scale(X)
scaled_X


# Q 2.1

# declare the variables
x1 <- 21
x2 <- 21
nsrcs <- 6
V <- 441

library(matlib)
A_LSR <- abs(inv(t(TC) %*% TC) %*% t(TC) %*% scaled_X)
D_LSR <- scaled_X %*% t(A_LSR)

# plot six retrieved sources
par(mfrow=c(3,4),mar=c(2,3,2,3))
for(k in 1:6){
  plot(matrix(A_LSR[k,], nrow=21), border=NA,
       xlab="", ylab="",
       main=sprintf("Retrieved Spatial Maps %i", k))
  plot(D_LSR[,k], type="l", xlab="", ylab="",
       main=sprintf("Retrieved Time Courses %i", k))
}

# plot 3rd col of D_LSR and 30th col of standardised X
par(mfrow=c(1,2), mar=c(5,5,5,5))
plot(D_LSR[,3], scaled_X[,30], xlab="3rd col of D_LSR",
     ylab="30th col of standardised X")

# plot 4th col of D_LSR and 30th col of standardised X
plot(D_LSR[,4], scaled_X[,30], xlab="4th col of D_LSR",
     ylab="30th col of standardised X")


# Q 2.2
lambda <- seq(0, 1, 0.05)
library(SAVER)
c_TLSR <- sum(calc.maxcor(TC, D_LSR))
c_TLSR

# finding the best lambda
for(i in 1:length(lambda)) {
  A_RR <- abs(inv((t(TC) %*% TC)+ (lambda[i]*V*diag(6))) %*% 
                    t(TC) %*% scaled_X)
  D_RR <- scaled_X %*% t(A_RR)
  c_TRR <- sum(calc.maxcor(TC, D_RR))
  print(c_TRR)
  current_lambda <- lambda[i]
  if(c_TRR > c_TLSR) break
}

best_lambda <- current_lambda*V
best_lambda

# maximum absolute correlation
c_TLSR <- calc.maxcor(TC, D_LSR)
c_TRR <- calc.maxcor(TC, D_RR)

# checking whether sum of c_TRR > sum of c_TLSR
sum(c_TRR)
sum(c_TLSR)


# plot first vector of A_RR and A_LSR when lambda = 1000
new_lambda <- 1000
new_A_RR <- abs(solve((t(TC) %*% TC) + (new_lambda*V*diag(6)), 
              t(TC) %*% scaled_X))

par(mfrow=c(1,1))
plot(new_A_RR[1,], type="l", 
     col="red", ylab="Value", lty=1, lwd=2, ylim=c(0,1.6),
     main="First vector of A_RR and A_LSR", xlab="Index")
lines(A_LSR[1,], type="l",col="black",lty=1, lwd=1)
legend(300, 1.6, legend=c("A_RR", "A_LSR"),
       col=c("red", "black"),lty=1:1, cex=0.8, lwd=2)



# Q 2.3

# function to generate A_LR
generate_Alr <- function(X, rho, TC){
  step <- 1/(norm(TC %*% t(TC)) * 1.1)
  thr <- rho*N*step
  Ao <- matrix(0, nsrcs, 1)
  A <- matrix(0, nsrcs, 1)
  Alr <- matrix(0, nsrcs, x1*x2)
  
  for (k in 1:(x1*x2)) {
    A <- Ao+step*(t(TC) %*% (X[,k]-(TC%*%Ao)))
    A <- (1/(1+thr)) * (sign(A)*pmax(replicate(nsrcs, 0), abs(A)-thr))
    
    for (i in 1:10) {
      Ao <- A
      A <- Ao+step * (t(TC)%*%(X[,k]-(TC%*%Ao)))
      A <- (1/(1+thr)) * (sign(A)*pmax(replicate(nsrcs, 0), abs(A)-thr))
    }
    Alr[,k] <- A
  }
  return(Alr)
}


rho <- seq(0, 1, 0.05)
MSE_LR <- numeric(length=length(rho))

for(a in 1:length(rho)){
  MSE_per_rho <- numeric(length=10)
  # calculate MSE for this rho over 10 realisations
  for(b in 1:10) {
    # generate new X
    new_Tt <- matrix(rnorm(240*6, mean=0, sd=0.5), 240, 6)
    new_Ts <- matrix(rnorm(6*441, mean=0, sd=sqrt(0.015)), 6, 441)
    new_X <- (TC + new_Tt) %*% (SM + new_Ts)
    new_X_scaled <- as.matrix(scale(new_X))
    
    new_Alr <- matrix(0, nsrcs, x1*x2)
    new_Alr <- generate_Alr(new_X_scaled, rho[a], TC)
    Dlr <- new_X_scaled %*% t(new_Alr)
    MSE_per_rho[b] <- sum(sum((new_X_scaled - Dlr%*%new_Alr)^2))/(N*V)
  }
  MSE_LR[a] <- mean(MSE_per_rho)
}


# plot lambda vs MSE
par(mfrow=c(1,2), mar=c(5,5,5,5))
plot(rho, MSE_LR, type="l")

# zoom in plot to see where the MSE bounced back
plot(rho[10:21], MSE_LR[10:21], type="l",
     xlab="rho from 0.45 to 1", ylab="MSE_LR")

# value of rho with minimum MSE
best_rho <- rho[(which.min(MSE_LR))]
best_rho



# Q 2.4

# estimate LR params using rho = 0.6
A_LR <- abs(generate_Alr(scaled_X, best_rho, TC))
D_LR <- scaled_X %*% t(A_LR)
c_TLR <- calc.maxcor(TC, D_LR)
c_SLR <- calc.maxcor(t(SM), t(A_LR))
c_SRR <- calc.maxcor(t(SM), t(A_RR))

# checking whether sum of c_TLR > sum of c_TRR
sum(c_TLR)
sum(c_TRR)

# checking whether sum of c_SLR > sum of c_SRR
sum(c_SLR)
sum(c_SRR)

# plotting RR and LR estimates
par(mfrow=c(6,4), mar=c(2,3,2,3))
for(k in 1:6){
  plot(matrix(A_RR[k,], nrow=21), border=NA,
       xlab="", ylab="",
       main=sprintf("Retrieved Spatial Maps RR %i", k))
  plot(D_RR[,k], type="l", xlab="", ylab="",
       main=sprintf("Retrieved Time Courses RR %i", k))
  plot(matrix(A_LR[k,], nrow=21), border=NA,
       xlab="", ylab="",
       main=sprintf("Retrieved Spatial Maps LR %i", k))
  plot(D_LR[,k], type="l", xlab="", ylab="",
       main=sprintf("Retrieved Time Courses LR %i", k))
}


# Q 2.5

# find PC of TC
prcomp_TC <- prcomp(TC)
Z <- prcomp_TC$x
par(mfrow=c(3,4), mar=c(2,3,2,3))

# plotting regressor Z vs TC
for(k in 1:6){
  plot(Z[,k], type="l", xlab="", ylab="",
       main=sprintf("Z %i", k))
  plot(TC[,k], type="l", xlab="", ylab="",
       main=sprintf("TC %i", k))
}

A_PCR <- abs(generate_Alr(scaled_X, 0.001, Z))
D_PCR <- scaled_X %*% t(A_PCR)

summary(prcomp_TC)
prcomp_TC$rotation

# plot six retrieved sources
par(mfrow=c(3,4),mar=c(2,3,2,3))
for(k in 1:6){
  plot(matrix(A_PCR[k,], nrow=21), border=NA,
       xlab="", ylab="",
       main=sprintf("Retrieved Spatial Maps %i", k))
  plot(D_PCR[,k], type="l", xlab="", ylab="",
       main=sprintf("Retrieved Time Courses %i", k))
}

# the eigenvalues are sd^2
eigval <- (prcomp_TC$sdev)^2
eigval
par(mfrow=c(1,1), mar=c(5,5,5,5))
plot(x= c(1,2,3,4,5,6), y=eigval, 
     xlab="Principal Component Number", 
     type="b", ylab="Eigenvalue",)