#Hierarchical Bayesian Quantile Regression Model
#--------
#R code for: Precipitation Scaling with Temperature in the Northeast US: Variations by Weather Regime, Season, and Precipitation Intensity
#--------


# ----------------------------------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------------------------------- #

# Citation: Precipitation Scaling with Temperature in the Northeast US: Variations by Weather Regime, Season, and Precipitation Intensity
# Geophysical Research Letters, 2022
# All rights reserved: Nasser Najibi, Sudarshana Mukhopadhyay, Scott Steinschneider
# Department of Biological and Environmental Engineering, Cornell University, Ithaca 14853, NY, United States, nn289@cornell.edu

# -----------------------------------------------------------------------------------------------------------------------------------

# Version 1.0 August 2021

# Hierarchical Bayesian quantile Regression Model # 
# for Td
bayes_cc_qreg_wr_multilevel <- function() {
  for(i in 1:nsites)
  {
    for(t in 1:nevents[i])
    {
      # Level 1 #
      precip[i,t] ~ dnorm(me[i,t],pe[i,t])
      me[i,t] <- (1-2*p)/(p*(1-p))*w[i,t] + mu[i,t]
      mu[i,t] <- alpha.wrs[i,t] + beta.wrs[i,t]*td[i,t]
      pe[i,t] <- (tau.q[i]*p*(1-p))/(2*w[i,t])
      w[i,t] ~ dexp(tau.q[i])
      # Level 2 #
      beta.wrs[i,t] <- b1[i]*wr1[i,t]+b2[i]*wr2[i,t]+b3[i]*wr3[i,t]+b4[i]*wr4[i,t]      
      alpha.wrs[i,t] <- a1[i]*wr1[i,t]+a2[i]*wr2[i,t]+a3[i]*wr3[i,t]+a4[i]*wr4[i,t]
    }
    # Prior #
    lsigma[i] ~ dunif(-10,10)
    sigma[i] <- exp(lsigma[i]/2)
    tau.q[i] <- pow(sigma[i],-2)
    a1[i] ~ dnorm(0,0.01)
    a2[i] ~ dnorm(0,0.01)
    a3[i] ~ dnorm(0,0.01)
    a4[i] ~ dnorm(0,0.01)
    b1[i] ~ dnorm(mu.b1,tau.b1)
    b2[i] ~ dnorm(mu.b2,tau.b2)
    b3[i] ~ dnorm(mu.b3,tau.b3)
    b4[i] ~ dnorm(mu.b4,tau.b4)
  }
  # Prior #
  mu.b1 ~ dnorm(0,1)
  mu.b2 ~ dnorm(0,1)
  mu.b3 ~ dnorm(0,1)
  mu.b4 ~ dnorm(0,1)
  tau.b1 ~ dgamma(1,0.1)
  tau.b2 ~ dgamma(1,0.1)
  tau.b3 ~ dgamma(1,0.1)
  tau.b4 ~ dgamma(1,0.1)
  sigma.b1 <- 1/sqrt(tau.b1)
  sigma.b2 <- 1/sqrt(tau.b2)
  sigma.b3 <- 1/sqrt(tau.b3)
  sigma.b4 <- 1/sqrt(tau.b4)
}


# -----------------------------------------------------------------------------------------------------------------------------------
# Hierarchical Bayesian Quantile Regression Model # 
# for Ta
bayes_cc_qreg_wr_multilevel <- function() {
  for(i in 1:nsites)
  {
    for(t in 1:nevents[i])
    {
      # Level 1 #
      precip[i,t] ~ dnorm(me[i,t],pe[i,t])
      me[i,t] <- (1-2*p)/(p*(1-p))*w[i,t] + mu[i,t]
      mu[i,t] <- alpha.wrs[i,t] + beta.wrs[i,t]*ta[i,t]
      pe[i,t] <- (tau.q[i]*p*(1-p))/(2*w[i,t])
      w[i,t] ~ dexp(tau.q[i])
      # Level 2 #
      beta.wrs[i,t] <- b1[i]*wr1[i,t]+b2[i]*wr2[i,t]+b3[i]*wr3[i,t]+b4[i]*wr4[i,t]      
      alpha.wrs[i,t] <- a1[i]*wr1[i,t]+a2[i]*wr2[i,t]+a3[i]*wr3[i,t]+a4[i]*wr4[i,t]
    }
    # Prior #
    lsigma[i] ~ dunif(-10,10)
    sigma[i] <- exp(lsigma[i]/2)
    tau.q[i] <- pow(sigma[i],-2)
    a1[i] ~ dnorm(0,0.01)
    a2[i] ~ dnorm(0,0.01)
    a3[i] ~ dnorm(0,0.01)
    a4[i] ~ dnorm(0,0.01)
    b1[i] ~ dnorm(mu.b1,tau.b1)
    b2[i] ~ dnorm(mu.b2,tau.b2)
    b3[i] ~ dnorm(mu.b3,tau.b3)
    b4[i] ~ dnorm(mu.b4,tau.b4)
  }
  # Prior #
  mu.b1 ~ dnorm(0,1)
  mu.b2 ~ dnorm(0,1)
  mu.b3 ~ dnorm(0,1)
  mu.b4 ~ dnorm(0,1)
  tau.b1 ~ dgamma(1,1)
  tau.b2 ~ dgamma(1,1)
  tau.b3 ~ dgamma(1,1)
  tau.b4 ~ dgamma(1,1)
  sigma.b1 <- 1/sqrt(tau.b1)
  sigma.b2 <- 1/sqrt(tau.b2)
  sigma.b3 <- 1/sqrt(tau.b3)
  sigma.b4 <- 1/sqrt(tau.b4)
}


# -----------------------------------------------------------------------------------------------------------------------------------
# Bayesian Quantile Regression Model # 
# for Td
bayes_cc_qreg_wr_singlelevel <- function() {
  for(i in 1:nsites)
  {
    for(t in 1:nevents[i])
    {
      # single Level #
      precip[i,t] ~ dnorm(me[i,t],pe[i,t])
      me[i,t] <- (1-2*p)/(p*(1-p))*w[i,t] + mu[i,t]
      mu[i,t] <- alpha.wrs[i] + beta.wrs[i]*td[i,t]
      pe[i,t] <- (tau.q[i]*p*(1-p))/(2*w[i,t])
      w[i,t] ~ dexp(tau.q[i])
    }
    # Priors #
    alpha.wrs[i] ~ dnorm(0,0.01)
    beta.wrs[i] ~ dnorm(0,0.01)
    lsigma[i] ~ dunif(-10,10)
    sigma[i] <- exp(lsigma[i]/2)
    tau.q[i] <- pow(sigma[i],-2)
  }
}


# -----------------------------------------------------------------------------------------------------------------------------------
# Bayesian Quantile Regression Model # 
# for Ta
bayes_cc_qreg_wr_singlelevel <- function() {
  for(i in 1:nsites)
  {
    for(t in 1:nevents[i])
    {
      # single Level #
      precip[i,t] ~ dnorm(me[i,t],pe[i,t])
      me[i,t] <- (1-2*p)/(p*(1-p))*w[i,t] + mu[i,t]
      mu[i,t] <- alpha.wrs[i] + beta.wrs[i]*ta[i,t]
      pe[i,t] <- (tau.q[i]*p*(1-p))/(2*w[i,t])
      w[i,t] ~ dexp(tau.q[i])
    }
    # Priors #
    alpha.wrs[i] ~ dnorm(0,0.01)
    beta.wrs[i] ~ dnorm(0,0.01)
    lsigma[i] ~ dunif(-10,10)
    sigma[i] <- exp(lsigma[i]/2)
    tau.q[i] <- pow(sigma[i],-2)
  }
}

# ----------------------------------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------------------------------- #
