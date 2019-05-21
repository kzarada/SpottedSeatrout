
###################################################################################
#
#            Hierarchical Model in Jags to estimate Spawning Interval and loss
#
##################################################################################


setwd("/Users/katiezarada/Library/Mobile Documents/com~apple~CloudDocs/Seatrout/CJS_Models/Seatrout_CJS_1")

#original model in Kéry's Bayesian Population Analysis book, Chapter 7
#updated model to be hierarchical mixed effects model 
#Sex is fixed effect, Total Length is mixed

modelFilename='Seatrout_CJS_H1.txt'
cat("
model {

# hyper priors (by sex)
   for(i in 1:s){
      phi.mu.s[i] ~ dunif(-10,10)
      phi.sd.s[i] ~ dunif(0,5)
      phi.tau.s[i] = 1/(phi.sd.s[i]*phi.sd.s[i])

      p.mu.s[i] ~ dunif(-10,10)
      p.sd.s[i] ~ dunif(0,5)
      p.tau.s[i] = 1/(p.sd.s[i]*p.sd.s[i])
   }

# group priors (by length)   
for(i in 1:s){
   for (j in 1:g){
      phi.g[i,j] ~ dnorm(phi.mu.s[i],phi.tau.s[i]) # Priors for group-specific survival
      p.g[i,j] ~ dnorm(p.mu.s[i],p.tau.s[i])  # Priors for group-specific recapture

      phi.g.trans[i,j] <- exp(phi.g[i,j])/(exp(phi.g[i,j])+1) # inverse logit transformation to put on 0-1 scale
      p.g.trans[i,j] <- exp(p.g[i,j])/(exp(p.g[i,j])+1) # inverse logit transformation to put on 0-1 scale
      } #j
   } #i

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- pow(phi.g.trans[sex[i],group[i]],lag[t])
      p[i,t] <- p.g.trans[sex[i],group[i]]
      } #t
   } #i

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill=TRUE,file=modelFilename)