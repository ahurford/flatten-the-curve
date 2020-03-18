rm(list=ls()) 
# Code to be embedded in the app
library(deSolve)
library(bbmle)

SIR.fit <- function(t, y, p) {
  with(as.list(c(y, p)), {
    dS <- -beta * S * I
    dI <-  beta * S * I - gamma * I - v * I
    dF <- v * I
    list(c(S = dS,I = dI, Fs = dF))
  })
}

nLL = function(beta){
  
  parms <- c(beta = beta, gamma = 1/13, v = 0.00237)
  I0 <- 0.005
  S0 <- 1 - I0
  mintime <- 0
  maxtime <- 250
  out <- ode(y = c(S = S0, I = I0, FS = 0), times = seq(mintime, maxtime, 1), SIR.fit, parms)
  df <- data.frame(out)
 
  Model.pred=NULL
  for(i in seq(1,length(Idata$I))){
    ii = which(df$time==Idata$time[i])
    Model.pred[i] = df$I[ii]
  }
  res = sum((Model.pred - Idata$I)^2)
  return(res)
}

# Fake data
Idata <- data.frame(time=c(1,5,10), I = c(.1,.12,.13))
fit = mle2(nLL, data = Idata, start=list(beta=.4))
beta.est=unname(coef(fit))
parms <- c(beta = beta.est, gamma = 1/13, v = 0.00237)
I0 <- 0.005
S0 <- 1 - I0
mintime <- 0
maxtime <- 250

out <- ode(y = c(S = S0, I = I0, FS = 0), times = seq(mintime, maxtime, 1), SIR.fit, parms)
df = data.frame(out)
plot(df$time, df$I, typ="l")
points(Idata$time, Idata$I)

# parameteric 95% CI
LowerCI = max(beta.est - 1.96*unname(stdEr(fit)),0)
UpperCI = beta.est + 1.96*unname(stdEr(fit))
parms <- c(beta = LowerCI, gamma = 1/13, v = 0.00237)
out.lower <- ode(y = c(S = S0, I = I0, FS = 0), times = seq(mintime, maxtime, 1), SIR.fit, parms)
parms <- c(beta = UpperCI, gamma = 1/13, v = 0.00237)
out.upper <-out <- ode(y = c(S = S0, I = I0, FS = 0), times = seq(mintime, maxtime, 1), SIR.fit, parms)
df.lower = data.frame(out.lower)
df.upper = data.frame(out.upper)
lines(df.lower$time, df.lower$I, lty=2)
lines(df.upper$time, df.upper$I, lty=2)
