

# CJS_Sim Helper functions


# helper functions

# Define function to simulate a capture-history (CH) matrix
simul.cjs <- function(PHI, P, marked){
   n.occasions <- dim(PHI)[2] + 1 #gets back 6 from 5 capture events
   CH <- matrix(0, ncol = n.occasions, nrow = sum(marked)) #makes a capture matrix that is 6 capture events by number of animals 
   # Define a vector with the occasion of marking
   mark.occ <- rep(1:length(marked), marked[1:length(marked)]) #when the animals were marked
   # Fill the CH matrix
   for (i in 1:sum(marked)){
      CH[i, mark.occ[i]] <- 1       # Write an 1 at the release occasion
      if (mark.occ[i]==n.occasions) next
      for (t in (mark.occ[i]+1):n.occasions){
         # Bernoulli trial: does individual survive occasion?
         sur <- rbinom(1, 1, PHI[i,t-1])
         if (sur==0) break    # If dead, move to next individual 
         # Bernoulli trial: is individual recaptured? 
         rp <- rbinom(1, 1, P[i,t-1])
         if (rp==1) CH[i,t] <- 1
         } #t
      } #i
   return(CH)
   }

   known.state.cjs <- function(ch){
   state <- ch
   for (i in 1:dim(ch)[1]){
      n1 <- min(which(ch[i,]==1))
      n2 <- max(which(ch[i,]==1))
      state[i,n1:n2] <- 1
      state[i,n1] <- NA
      }
   state[state==0] <- NA
   return(state)
   }

   # Function to create a matrix of initial values for latent state z
cjs.init.z <- function(ch,f){
   for (i in 1:dim(ch)[1]){
      if (sum(ch[i,])==1) next
      n2 <- max(which(ch[i,]==1))
      ch[i,f[i]:n2] <- NA
      }
   for (i in 1:dim(ch)[1]){
   ch[i,1:f[i]] <- NA
   }
   return(ch)
   }

get.first <- function(x) min(which(x!=0))





   known.state.2 <- function(ch){
   state <- ch

   for (i in 1:dim(ch)[1]){
      if(sum(ch[i,]) > 1) next 
      n1 <- min(which(ch[i,]==1))
      n2 <- max(which(ch[i,]==1))
      state[i,n1:n2] <- 1
      state[i,n1] <- NA
      }
   state[state==0] <- NA
   state[state == 2] <- 1
   return(state)
   }

cjs.init.2 <- function(ch,f){
   for (i in 1:dim(ch)[1]){
      if (sum(ch[i,]) > 1) next
      n2 <- max(which(ch[i,]==1))
      ch[i,f[i]:n2] <- NA
      }
   for (i in 1:dim(ch)[1]){
   ch[i,1:f[i]] <- NA
   }
   return(ch)
   }