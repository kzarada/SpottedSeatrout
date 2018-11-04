# Code for analyzing a hierarchical bayesian CJS model for estimate 
# probability of spawning and loss from a spawning aggregation of spotted seatrout 
# K. Zarada

#load packages 
require(tidyverse)
require(data.table)
require(scales)

#set working directory 
sst = read.csv("/Users/katiezarada/Library/Mobile Documents/com~apple~CloudDocs/Seatrout/Data/sst_new_ghostremoved.csv", header= TRUE, stringsAsFactors = FALSE) #laptop

#Tidy Data 
sst$year_tag = paste(sst$year, sst$newtag)
sst = sst[-which(sst$year_tag == "2008 19")]  #removing the recaptures
sst = sst[-which(sst$year_tag == "2008 11")]  #removing the recaptures


# convert to tbl
	# make a subset for zones 2 and 3
	# a subset of just night time hours when they are spawning
	# subset to relevant columns
		sst_2007.tbl = tbl_df(sst) %>%
			  filter(zone.1 == 2| sst$zone.1 == 3) %>%
			  filter(hour > 9) %>%
			  filter(year == 2007 | year == 2008) %>%
			  rename(ID = newtag) %>%
			  rename(ImplantDate = implantdate) %>%
			  rename(RecapDate = date) %>%
			  dplyr::select(ID,TL,Sex,ImplantDate,RecapDate) %>%
			  mutate(ImplantDate=as.Date(ImplantDate, format = "%m/%d/%y")) %>%
			  mutate(RecapDate=as.Date(RecapDate, format = "%m/%d/%y")) %>%
			  arrange(ImplantDate,RecapDate) %>%
			  distinct()
	# combine ImplantDate and RecapDate into 1 column
		index1 = 1:nrow(sst_2007.tbl)
		sst_2007.tbl = rbind(sst_2007.tbl,sst_2007.tbl)
		index2 = (length(index1)+1):nrow(sst_2007.tbl)
		sst_2007.tbl$Date = c(sst_2007.tbl$ImplantDate[index1],sst_2007.tbl$RecapDate[index2])
		sst_2007.tbl = sst_2007.tbl %>%
						dplyr::select(ID,TL,Sex,Date) %>%
						 arrange(Date) %>%
						distinct()



# get unique values
	unique.ID = sort(unique(sst_2007.tbl$ID))
	unique.Dates = sort(unique(sst_2007.tbl$Date))
	lag2next = as.integer(diff(unique.Dates))

# make matrix
	CH.mat = matrix(0,nrow=length(unique.ID),ncol=length(unique.Dates)) # capture history
	rownames(CH.mat) = unique.ID
	colnames(CH.mat) = as.character(unique.Dates)

# convert to data.table for fast indexing
	sst_2007.DT = as.data.table(sst_2007.tbl)
	setkey(sst_2007.DT,ID)
# iterate across fish and recap dates
	for(i in 1:length(unique.ID))
	{
		sub = sst_2007.DT[.(unique.ID[i])]
		col.index = match(as.character(sub$Date),colnames(CH.mat))
		CH.mat[i,col.index] = 1
	}
#save.image("/Users/katiezarada/Library/Mobile Documents/com~apple~CloudDocs/Seatrout/Data/telemetry_matrix.RData")
# create grouping variable
	ID.DF = tbl_df(sst_2007.DT) %>%
			dplyr::select(ID,TL,Sex) %>%
			distinct()
	ID.DT = as.data.table(ID.DF)
	ID.DT$SexNumeric = as.numeric(as.factor(ID.DF$Sex))  #change this to switch what the grouping variable is, so Sex, TL_Group, or SexTL_Group 
	ID.DT$TLGroupNumeric = NA
	
	for(i in 1:nrow(ID.DT))
	{
		if(ID.DT$SexNumeric[i] == 1)
		{
			ID.DT$TLGroupNumeric[i] = ifelse(ID.DT$TL[i]>=500,3,ifelse(ID.DT$TL[i] >=400,2,1))
		} else {
			ID.DT$TLGroupNumeric[i] = ifelse(ID.DT$TL[i]>=500,3,ifelse(ID.DT$TL[i] >= 400 ,2,1))
		}
	}

# get vector of first marking
	source("/Users/katiezarada/Library/Mobile Documents/com~apple~CloudDocs/Seatrout/CJS_Models/CJS_Sim/CJS_Sim.helper_functions.r")
	f.vec = unname(apply(CH.mat, 1, get.first))


# run JAGS model
	require(rjags)
	require(runjags)
	# require(parallel)
	# setwd
		setwd("/Users/katiezarada/Library/Mobile Documents/com~apple~CloudDocs/Seatrout/CJS_Models/Seatrout_CJS_1")
	# Bundle data
      jags.data = list( y = CH.mat,f = f.vec, nind = dim(CH.mat)[1], n.occasions = dim(CH.mat)[2], z = known.state.cjs(CH.mat), s = length(unique(ID.DT$SexNumeric)), sex = ID.DT$SexNumeric, g = length(unique(ID.DT$TLGroupNumeric)), group = ID.DT$TLGroupNumeric, lag = lag2next)

   # Initial values
      inits = function(){list(z = cjs.init.z(CH.mat, f.vec),phi.mu.s=runif(2, 0, 1),phi.sd.s=runif(2, 0, 1),p.mu.s=runif(2, 0, 1),p.sd.s=runif(2, 0, 1) ,phi.g = matrix(runif(6, 0, 1),nrow=2,ncol=3), p.g = matrix(runif(6, 0, 1),nrow=2,ncol=3))}  
   # sink JAGS model
      source("/Users/katiezarada/Library/Mobile Documents/com~apple~CloudDocs/Seatrout/CJS_Models/Seatrout_CJS_1/Seatrout_CJS_H1.r")
   # Parameters monitored
      parameters = c("phi.mu.s","phi.sd.s","p.mu.s","p.sd.s","phi.g", "p.g")

   # MCMC settings
      ni = 50000
      nt = 3
      nb = 2000
      nc = 4

   # Call JAGS from R 
      Seatrout_CJS_H1 = run.jags(data = jags.data, 
      					inits = inits, 
      					monitor = parameters, 
      					model = "Seatrout_CJS_H1.txt", 
      					n.chains = nc, 
      					thin = nt, 
      					sample = ni, 
      					burnin = nb,
      					jags.refresh = 15) # method="parallel"

   # Summarize posteriors
      print(Seatrout_CJS_H1, digits = 3)
      plot(Seatrout_CJS_H1$mcmc)


   # diagnostics
     require(coda)
     # Gelman and Rubin Multiple Sequence Diagnostic - convergence
		gelman.diag(Seatrout_CJS_H1$mcmc)
		gelman.plot(Seatrout_CJS_H1$mcmc)

	# Geweke Diagnostic - convergence
		geweke.diag(Seatrout_CJS_H1$mcmc)

	# Heidelberg and Welch convergence diagnostic - MCMC chain length
		heidel.diag(Seatrout_CJS_H1$mcmc,eps = 0.1, pvalue = 0.05)


	# combine all chains into 1 matrix
		Seatrout_CJS_H1.MCMC = as.matrix(Seatrout_CJS_H1$mcmc)
		Seatrout_CJS_H1.params = colnames(Seatrout_CJS_H1.MCMC)
	
	#get summary on posteriors 
	info.Z <- as.data.frame(Seatrout_CJS_H1.MCMC) %>% 
	            select(contains("phi")) %>% 
	            mutate_all(funs(-log(.^365))) 
	info.Z <- apply(info.Z, 2, FUN = summary)          
	
  info <- Seatrout_CJS_H1.MCMC[,c(5,6,15:20)]
	info <- apply(info, 2, FUN = summary)
		
		
##############################################################
#
#		Violin plots for Sex_TL Bins 
#
################################################################	

data_summary <- function(x) {
   m <- mean(x)
   ymin <- m-sd(x)
   ymax <- m+sd(x)
   return(c(y=m,ymin=ymin,ymax=ymax))
}


zindex = -log((Seatrout_CJS_H1.MCMC[,c(1,2,9:14)])^365)

datZ = Seatrout_CJS_H1.MCMC[,c(1,2,9:14)]
colnames(datZ) <- c("F","M", "F_1", "M_1", "F_2", "M_2", "F_3", "M_3")
datZ = as.data.frame(datZ)
datZ$num = seq(1, 20000)

datZ = datZ %>% gather(Group, Value, 1:8)
 datZ$Value = -log(datZ$Value^365)


 ggplot(datZ, aes(x=Group, y=Value, fill = Group)) + 
  geom_violin() +  stat_summary(fun.data=data_summary) + 
  labs(y= "Instantaneous Loss", x = "Group") + ylim(0,20)+ theme_bw() +  
  theme(legend.position="none") + scale_fill_grey() + 
  scale_x_discrete(limits = c("F", "M","F_1", "F_2", "F_3","M_1", "M_2", "M_3"), labels=c("F" = "Female", "M" = "Male", "F_1" ="Small Female", "F_2" = "Med Female", "F_3"="Large Female","M_1" ="Small Male", "M_2" = "Med Male", "M_3"="Large Male" )) +
  theme(axis.text.x = element_text(face="bold", size=14),axis.text.y = element_text(face="bold", size=14), axis.title.x = element_text(face = "bold", size = 18), axis.title.y = element_text(face = "bold", size = 18))



dat = Seatrout_CJS_H1.MCMC[,c(5,6,15:20)]
colnames(dat) <- c("F","M", "F_1", "M_1", "F_2", "M_2", "F_3", "M_3")
dat = as.data.frame(dat)
dat$num = seq(1, 20000)
dat = dat %>% gather(Group, Value, 1:8)

 
 ggplot(dat, aes(x=Group, y=Value, fill = Group)) + 
  geom_violin() +  stat_summary(fun.data=data_summary) + labs(y= "Probability of Being in the Spawning Grounds", x = "Group") + theme_bw() + 
    scale_x_discrete(limits = c("F", "M","F_1", "F_2", "F_3","M_1", "M_2", "M_3"), labels=c("F" = "Female", "M" = "Male", "F_1" ="Small Female", "F_2" = "Med Female", "F_3"="Large Female","M_1" ="Small Male", "M_2" = "Med Male", "M_3"="Large Male" ))+ 
    ylim(0,0.75) + theme(legend.position="none") + scale_fill_grey() + 
      theme(axis.text.x = element_text(face="bold", size=14),axis.text.y = element_text(face="bold", size=14), axis.title.x = element_text(face = "bold", size = 18), axis.title.y = element_text(face = "bold", size = 18))



	

