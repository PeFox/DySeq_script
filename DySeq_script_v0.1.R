###############################################
### Installing & loading Packages from CRAN ###
###############################################


# install.packages("TraMineR") # for graphical analysis and research question 4
# install.packages("RColorBrewer")
# install.packages("gmodels") # must be installed!
# install.packages("MASS") # must be installed!
# install.packages("survival") # must be installed for research question 3!
# install.packages("fpc") # must be installed for research question 4!
# install.packages("cluster") # must be installed for research question 4!
# install.packages("devtools")


library(cluster)
library(fpc)
library(survival)
library(TraMineR)
library(RColorBrewer) # only used for alternative Graphics
library(devtools) # needed for DySeq-installation


####################################
### Installing DySeq from GitHub ###
####################################

# remove.packages("DySeq") # if a previous version of DySeq is installed 
# install_github("PeFox/DySeq") # make sure Devtools is installed an loaded!
library(DySeq) 
# help(DySeq)
# help(CouplesCope) # click on the help files Index for an overview of all functions and data

#############################
### Loading example data  ###
#############################

data("CouplesCope") # loading example data
View(CouplesCope) # looking into the data

# code == Couples ID
# IKCBxx == Stressrecation of female partner at the time intervall xx with {0 == no reaction; 1 == stress reaction}
# DCCBxx == dyadic copingreaction of male partner at the time intervall xx with {0 == no reaction; 1 == dyadic coping response}


mydata<-CouplesCope # mydata will be used as a working copy of the dataset


# Every line represents one observational unit (one couple)
# Sequences are spereated for every behavior
# for further analysis we will combine them using StateExpand


my.expand<-StateExpand(CouplesCope, 2:49, 50:97) # create combined sequences for further analyses



##########################
### GRAPHICAL ANALYSIS ###
##########################


# now, using the TraMineR for generating the state-distribution plot

couple.labels <-c("no reaction", "stress only", "coping only", "both reactions")  # create labels for plot
couple.seq <- seqdef(my.expand, labels = couple.labels) # create a sequence object (the way TraMineR represents sequences)
seqdplot(couple.seq)


# Alternatively a grey version.

attr(couple.seq , "cpal") <- brewer.pal(4, "Greys") # see figure 2
seqdplot(couple.seq )


# Entropy-plot: how much do states diverge within time intervals?
# Turbalance-plot how much do states diverge within a couples? both: figure 3

par (mfrow = c(1,2))
{
seqdplot(couple.seq)

Entropy <- seqstatd(couple.seq)$Entropy
plot(Entropy, main= "Entropy", col="black", xlab = "Time in 10 sec. intervall", type ="l")

Turbulance <- seqST (couple.seq)
hist (Turbulance, main = "Histrogram of turbulance")
}

par (mfrow = c(1,1)) # returning to one graphic per frame!

# If turbulance should be investigated further, one can compute a dataframe with the
# original caseIDs included:
turbID<-data.frame(mydata[,1], Turbulance)
colnames(turbID)<-c("ID", "Turbulance")
# Further analyses could be comparing Turbulances between certain types of couples
# or cases in which males were stressed with cases in which females were stressed.



##################################################################
### Research question 1:                                       ###
### Is there an association between stress communication by    ###
### one partner and dyadic coping by the other partner?        ###
##################################################################


# Step 1: Computing the frequencies of both behaviors,
# using the function NumbOccur from the DySeq-package
stress.sumscores<-NumbOccur(CouplesCope[,2:49], 1, prop=FALSE)
DC.sumscores<-NumbOccur(CouplesCope[,50:97], 1, prop=FALSE)

# Explanation:
# Argument 1: Columns 2:49 are containing the stressreactions!
# Argument 2: We are interested in the Numer of one's, because every '1' represents that this behavior was shown!
# Argument 3: prop=FALSE provides absolute frequencies, while TRUE would provide relative frequencies.
#             options would result in the same correlation, but absolute frequencies conatin more information
#             if further analysis should be conducted!


# Step 2: Correlation of both frequencies

cor.test(stress.sumscores, DC.sumscores)






##################################################################
### Research question 2:                                       ###
### Does a particular behavior by one partner evoke a specific ###
### reaction by the other? Does stress communication evoke a   ###
### dyadic coping response?                                    ###
##################################################################



# Step 1: Creating state-transition tables:
my.trans<-StateTrans(my.expand, FALSE) # Second argument FALSE means, that the second sequence will be used as dependend variable!

my.trans # inspecting the mean transitionsrates!


# Inspecting a singe case:
which(mydata==129) # getting the row number for ID 129 from the original dataset

my.trans[[41]] # inspecting couple 129 (rownumber 41)
# see table 2


# p.values are missing for single case analysis, but will be provided in the future



# Step 2: Running logit-regressions for every couple

# The LoqSeq function of the DySeq-Package can be used to simplfy this normaly cumbersome procedure

my.logseq<-LogSeq(my.trans, delta=0.5, single.case=TRUE) # Delta is a constant added to every cell to prevent cells with zero frequencies
# The higher delta, the more conservative biased the estimates (the are biased towards zero)
# A delta of 0.5 is commonly used, but if no zero frequencies are present it should be
# replaced with a zero!
# If replaced with zero and zero frequencies are present, all cases with at least one
# zero frequency cell will be excluded from analysis!
# single.case=TRUE computes p-values for every single analysis!

# Single cases can be extracted via the single.LogSeq function

single.LogSeq(my.logseq, 41) # here for couple 129 (rownumber 41)
# see table 3


# Step 3: Aggregating

# just print the output for aggregation and automatic testing of estimates for being different from zero

my.logseq # see table 4


# Plotting interaction!

plot(my.logseq)




# Further analysis: Rerun the procedure with stress as dependend variable

my.trans.stress<-StateTrans(my.expand, TRUE) # This time TRUE means: First Sequence (Stress reaction) will be used as dependend variable!

my.logseq.stress<-LogSeq(my.trans.stress, delta=0.5)

my.logseq.stress # estimates for actor- and partnereffect on stress
my.logseq # estimates for actor- and partnereffect on coping

# (see fig. 5)






########################################################################
### Research Question 3:                                             ###
### What is the typical duration of stress communication?            ###
### And does it depend on covariates, for instances men’s            ###
### dyadic coping ability rated by their partners?                   ###
########################################################################


# Calculating the last occurance of stress communication (last.stress)

last.stress<-LastOccur(mydata[,2:49],1)


# Producing hazard, survival and cumhazard for stress duration

# First a variable is needed that indicates if the event was shown (1) or not (0)
event<-c(rep(1, length(last.stress))) # it is shown for every case!
event[last.stress>=48]<-0 # exept the ones that show stress communication till time intervall 48
1-mean(event)
stress.surv<-Surv(last.stress,event) # The duration and the event variable are combined in a Surv-object

fit1 <- coxph(stress.surv~1, ties="breslow") # fittet without a covariate


# We can also plot the survival curve
par (mfrow = c(1,3)) # creats a frame for three graphics

# Survival and cumhazard are easy to produce with the survival package
plot(survfit(fit1), conf.int="none", xlab="Time", ylab="Survival Probability", xlim=c(0,48))

# Execute the following 9 lines only if Median Lifetime should be added!
x <- 45
y<-seq(0, 0.5, by=0.01)
x<-rep(x,length(y)) 
lines(x,y, lty=2)
x<-seq(0, 45, by=0.1)
y<-rep(0.5, length(x))
lines(x,y, lty=2)
x<-locator(1)
text(x$x, x$y, "ML", cex=1.2)

plot(survfit(fit1), conf.int="none", xlab="Time", ylab="cumulated hazard", fun="cumhaz")

# yet, there is no option for the uncumulated hazard,
# but the DySeq package provides one:

NonCumHaz(survfit(fit1), plot=T) # Figure 6
NonCumHaz(survfit(fit1))
survfit(fit1)$surv
survfit(fit1)$cumhaz
cbind(NonCumHaz(survfit(fit1)),survfit(fit1)$time)



# returning to one plot per frame
par (mfrow = c(1,1)) 






# Adding the covariate: 
# Dyad.Coping ability of the man rated by the woman
cent<-mydata[,98]-mean(mydata[,98])
fit2 <- coxph(stress.surv~cent, ties="breslow") # fittet without a covariate
summary(fit2)

plot(survfit(fit2)$time, NonCumHaz(survfit(fit2)), xlim=c(19,45), ylim=c(0,0.3), type="l", xlab="Time", ylab="predicted hazard")
lines(survfit(fit2)$time, NonCumHaz(survfit(fit2))*1.15^5, lty=4)
lines(survfit(fit2)$time, NonCumHaz(survfit(fit2))*1.15^-5, lty=2)

x<-locator(1)
legend(x$x, x$y, lty=c(1,4,2) )




# Hazardratio is positive, higher ratings indicate more hazard
# more hazard indicates shorter durations.
# But is very close to one (no effect) and not sign. therefore
# it is very likely that no effect exists or it is at least very small



################################################################
### Research Question 4:                                     ###
### Do couples differ from each other with respect to their  ###
### typical behavioral patterns in times of stress?          ###
################################################################


# We will use our Sequence-Object from the graphical analysis!
# That is a special type of object used within the TraMineR-Package

# This was our Sequence object from lines 43 and 44





# Now we probe for latent types of couples via OM-Distances

# first step is specifying the cost-matrix.
# in this example the Method "TRATE" is used
#
# using "CONSTANT" means that all substitutions will be weightet equaly
#
# see the TraMineR documentation for specifying custom cost-matrices!

submat <- seqsubm(couple.seq, method = "TRATE")



submat # showing the substitution matrix

seqtrate(couple.seq) # or the transition probabilities



# Now, using this cost-matrix for computation of dissimiliarities
dist.oml <- seqdist(couple.seq, method = "OM", sm = submat)



# determine the optimal number of clusters

plot(pam(dist.oml, pamk(dist.oml)$nc))


# Screeplot (Kudos to 'Ben' from stackoverflow.com for this and the last two graphics!)

wss <- (nrow(dist.oml)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")



# Dendrogramm:
clusterward1 <- agnes (dist.oml, diss=TRUE, method = "ward")
plot (clusterward1, which.plots=2)


# 2 Clustersolution [for other solution change the value behind the argument 'k']
cluster2 <- cutree (clusterward1, k=2)
cluster2fac <- factor (cluster2, labels = c("Cluster 1", "Cluster 2"))
seqdplot (couple.seq, group = cluster2fac)



# repeating the sequence analysis of Bakeman and Gottman for sperated groups
# at the moment only implented for comparsion of two groups
# if more clusters are found run the procedure several time,
# and create a subgroup vector with a one's and tow's for the two clusters
# that should be compared an zeros for the other clusters!
# and be aware of alpha cumulation!
LogSeq(my.trans, delta=0.5, subgroups=cluster2)



# coxregression with subgroups
clust.dummy<-numeric(length(cluster2))
clust.dummy[cluster2==1]<-0
clust.dummy[cluster2==2]<-1 # recoding into dummy: first group is now coded 0, second is 1!


# adding to our regression with the prolonged covariate
fit5 <- coxph(stress.surv~prolong.DC+clust.dummy)
summary(fit5)



################################################################
### Using EstFreq and EstTime to determine the number        ###
### of expected zero and low frequencies cells               ###
################################################################

# First step: Define a matrinx conatining the expected transition rates!
my.trans.table<-matrix(c(0.57, 0.13,0.05,0.05,0.05, 0.05,0.05,0.05),4,2)

# Run EstFreq: t is number of timeintervall, min.cell defines what counts as low frequencies 
# k is the number of simulations
my.cellproblems<-EstFreq(my.trans.table, t=100, min.cell=5, k=20000)

# Printing the result
my.cellproblems

# Percantage is percantage of cases with at least one Zero or low frequency
# Mean cells is the mean number of cells over all cells
# max is the maximum observed number of cells with low or zero frequencies within single cases
# min is minimum observed number of cells with low or zero frequencies within single cases (typically zero)


# plots the Results of EstFreq for several time point
# k is lowered to 5000 for demonstration purposes
my.EstTime.plot<-EstTime(my.trans.table, t=50:100, k=5000)

# printing will result in a plot of time point vs. expected number of low and zero frequencies
my.EstTime.plot

