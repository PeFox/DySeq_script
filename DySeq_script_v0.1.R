###############################################
### Installing & loading Packages from CRAN ###
###############################################

################
##  Preamble  ##
################

# This R-script provides a hands-on-tutorial for all
# analyses covered in 
# "Analyzing dyadic sequence data - research questions and implied statistical models"
# by ... , ... and ...
# 
# Please make sure to install all required packages,
# including the "Dyseq" which provides the sample data!

### Content                         Lines

## A.1 Prerequisite Steps      
#  - packages from CRAN             28-48
#  - package from Github            49-55

## A.2 Example Data
#  - loading the data
#  - details on data 

## A.3 Graphical Analysis
#  - state-distribution-plot
#  - entropy-plot
#  - Number of transitions 

## Research question 1:
#  - Pearson Correlation 
 
## Research question 2:
#  (aggregated logit models)
#  - state-transition tables
#  - multiple logit-regressions
#  - aggregating
#  - further Analysis, APIM

## Research question 3:
#  - transforming sequence-data 
#    into time-to-event-data
#  - Hazard, survival and cumhazard
#  - cox-regression

## Research question 4:
#  - OM-distances 
#  - clustering
#  - subgroup analyses


## additional DySeq functions: 
#  - needed number of time intervalls 



#############################
##  1. prerequisite steps  ##
#############################

##  installing & loading packages from CRAN

# install.packages("TraMineR")      # for graphical analysis and research question 4
# install.packages("RColorBrewer")  # for grey-shaded graphics
# install.packages("gmodels")       # must be installed!
# install.packages("MASS")          # must be installed!
# install.packages("survival")      # must be installed for research question 3!
# install.packages("fpc")           # must be installed for research question 4!
# install.packages("cluster")       # must be installed for research question 4!
# install.packages("devtools")      # must be installed for installing packages from github

# loading packages!
library(cluster)
library(fpc)
library(survival)
library(TraMineR)
library(RColorBrewer) # only used for alternative Graphics
library(devtools) # needed for DySeq-installation



##  Installing DySeq from GitHub 

# remove.packages("DySeq")          # if a previous version of DySeq is installed 
# install_github("PeFox/DySeq")     # make sure Devtools is installed an loaded!

library(DySeq)                      # loading DySeq;
                                    # for further detail use help(DySeq) and click on the help file's index!



#######################
##  2. Example data  ##
#######################

data("CouplesCope")   # loading example data
View(CouplesCope)     # Invoke data viewer for the example data

# code    ==  Couples ID
# IKCBxx  ==  Stressrecation of female partner at the time intervall xx with {0 == no reaction; 1 == stress reaction}
# DCCBxx  ==  dyadic copingreaction of male partner at the time intervall xx with {0 == no reaction; 1 == dyadic coping response}
# EDCm    ==  men’s self-assessed dyadic coping ability

mydata<-CouplesCope # mydata will be used as a working copy of the CouplesCope dataset


# Every line represents one observational unit (one couple)
# Sequences are spereated for every behavior
# for further analysis we will combine them using StateExpand


my.expand<-StateExpand(CouplesCope, 2:49, 50:97) # create combined sequences via state expand procedure
                                                 # rows == couples; columns == time intervall
                                                 # cells: {0=no SC/DC; 1=SC only; 2=DC only; 4=SC+DC}


#############################
##  3. GRAPHICAL ANALYSIS  ##
#############################

## Objects from previous sections needed: 
#  - mydata     
#  - my.expand  


# state-distribution plot (using TraMineR) 

couple.labels <-c("time interval", "SC only", "DC only", "SC+DC")  # create labels for plot
couple.seq <- seqdef(my.expand, labels = couple.labels) # create a sequence object (the way TraMineR represents sequences)
seqdplot(couple.seq)


# Alternatively a grey version (using RColorBrewer)

attr(couple.seq , "cpal") <- brewer.pal(4, "Greys") # see figure 2
seqdplot(couple.seq, cex.legend=0.8, withlegend="right")


# Figure 3 Entropy-plot and Histogramm of "Number of transitions"

par (mfrow = c(1,2)) # preparing the graphic device to show two graphics in a row
{
# Entropy-plot: how much do states diverge within time intervals?
Entropy <- seqstatd(couple.seq)$Entropy
plot(Entropy, main= "Entropy", col="black", xlab = "Time in 10 sec. intervall", type ="l")

# Histogramm of "Number of transitions": How often do couples change from one state into another
SeqL<-seqtransn(couple.seq)
hist(SeqL, main="Number of transitions", xlab="State-transitions")
}

par (mfrow = c(1,1)) # the graphic device is set back to showing only one graphik at a time. 

summary(SeqL) # provides some describtives about the "Number of transitions"

# Note:
cor(SeqL, mydata$EDCm) # the Number of transition can be investigated further.
                       # For example, the higher men’s self-assessed dyadic coping  
                       # ability the less often couples change their states!
                    



###################################################################
###  Research question 1:                                       ###
###  Is there an association between stress communication by    ###
###  one partner and dyadic coping by the other partner?        ###
###################################################################

## Objects from previous sections needed: 
#  - mydata     


# Step 1: Computing the frequencies of both behaviors,
# using the function NumbOccur from the DySeq-package
stress.sumscores<-NumbOccur(mydata[,2:49], 1, prop=FALSE)
DC.sumscores<-NumbOccur(mydata[,50:97], 1, prop=FALSE)

# Explanation:
# Argument 1: Columns 2:49 are containing the stressreactions!
# Argument 2: We are interested in the Numer of one's, because every '1' represents that this behavior was shown!
# Argument 3: prop=FALSE provides absolute frequencies, while TRUE would provide relative frequencies.
#             options would result in the same correlation, but absolute frequencies contain more information
#             if further analysis should be conducted!


# Step 2: Correlation of both frequencies
cor.test(stress.sumscores, DC.sumscores)

# plotting the correlation [Note: not shown in the article]
plot(stress.sumscores, DC.sumscores, ylab="sum of DC", xlab="sum of SC")
abline(lm(stress.sumscores~DC.sumscores))





####################################################################
###  Research question 2:                                        ###
###  Does a particular behavior by one partner evoke a specific  ###
###  reaction by the other? Does stress communication evoke a    ###
###  dyadic coping response?                                     ###
####################################################################


## Objects from previous sections needed: 
#  - mydata     
#  - my.expand  


##############################################
#  Step 1: Creating state-transition tables  # 
##############################################

my.trans<-StateTrans(my.expand, first=FALSE) # The argument "first" means:
                                             # should the variable that was entered 
                                             # first in the StateExpand() be used 
                                             # used as dependend variable?
                                             # If FALSE, the second will be used.
                                             # Therefore, DC is dependend variable
                                             # in this example. 


my.trans # calling the list containing the state-transition tables
         # shows the mean transition frequencies!


# However, if a single case should be inspected, the row number is needed. 
# For example, if the couple with ID 129 should be inspected, one could 
# use which to determine the row number of ID 129 from the dataset:
which(mydata$code==129)  
                        
my.trans[[41]] # inspecting couple 129 (rownumber 41)
               # btw: the transistion-tables are stored within a list; 
               # therefore double-brackets are needed for subscripting!

# Or, if relative frequencies are preferred: 
my.trans[[41]]/sum(my.trans[[41]]) 

# Note:
# p.values for single case analysis are missing at the moment
# but will be provided in future versions


########################################
#  Step 2: multiple logit-regressions  #
########################################

# This step can be percieved very cumbersome, 
# yet the LogSeq function from the DySeq packages 
# simplifies the process!

my.logseq<-LogSeq(my.trans, delta=0.5, single.case=TRUE) 
# Notes: 
# First Argument must be a List of dichotomous transition-tables. 
# Delta is a constant added to every cell to prevent cells with zero frequencies.
# The higher the delta the more conservative biased the estimates. 
# A value of zero should be used if no zero frequencies exist or if
# all cases with zero frequencies should be stripped from the analysis 
# (listwise deletion)
# single.case=TRUE computes p-values for every single analysis!
# Single cases can be extracted via the single.LogSeq function
# However, p-values for the single case analysis may not be 
# trustworthy at the moment (inconsistent with other software implementations).
# p-values for the overall effects are checked and are consistent with
# the results obtained by other software implementations (SPSS, Lem, etc.)

single.LogSeq(my.logseq, 41) # here for couple 129 (rownumber 41)
# see table 3



#########################
#  Step 3: Aggregating  #
#########################


# The last step is simply done by calling the output
# of the LogSeq() function
my.logseq # [Note: see table 3 in the article]


# Plotting interaction! [Note: Not shown in the articel]
plot(my.logseq)



############################
#  Further analysis: APIM  #
############################

# Rerun the procedure with stress as dependend variable!

my.trans.stress<-StateTrans(my.expand, TRUE) # This time TRUE means: First Sequence (Stress reaction) will be used as dependend variable!
my.logseq.stress<-LogSeq(my.trans.stress, delta=0.5)

my.logseq.stress # estimates for actor- and partnereffect on stress
my.logseq # estimates for actor- and partnereffect on coping

# (see fig. 4)






###############################################################
###  Research Question 3:                                   ###
###  What is the typical duration of stress communication?  ###
###  And does it depend on covariates, for instances men’s  ###
###  self-assessed dyadic coping ability?                   ###
###############################################################


# Calculating the last occurance of stress communication (last.stress)

last.stress<-LastOccur(mydata[,2:49],y=1)
# mydata[,2:49] accesses only the SC-sequences of mydata
# second Argument y specifies if the last observed 0 or 1 should be computed
# Note: other values for y are possible, but most other functions 
#       of DySeq will only support dichotomous sequences.



####################################
#  Hazard, survival and cumhazard  #
####################################

# First a variable is needed that indicates if the event was shown (1) or not (0)
event<-c(rep(1, length(last.stress))) # it is shown for every case!
event[last.stress>=48]<-0 # exept the ones that show stress communication till time intervall 48

# the mean of the new vector "event" equals the relative frequency of couples in which the 
# event occored 1 minus this frequency equals the relative frequency of censored cases
1-mean(event)

# The duration and the event variable are combined in a Surv-object
stress.surv<-Surv(last.stress,event) 

fit1 <- coxph(stress.surv~1, ties="breslow") # ~1 means: fittet without a covariate
                                             # different estimators exist to handle 
                                             # multiple events within one time-interval 
                                             # in the article "breslow" was shown


par (mfrow = c(1,3)) # Preparing R's graphic device for three graphics


# Survival and cumhazard are already contained in the object "fit1"!
plot(survfit(fit1), conf.int="none", xlab="Time", ylab="Survival Probability", xlim=c(0,48))

# Note: if you want to inspect the survival tabled against the points of time
# run:  data.frame(survfit(fit1)$surv,survfit(fit1)$time)


# Note: if Median Lifetime should be added,
# run the following 9 lines
# x <- 45
# y<-seq(0, 0.5, by=0.01)
# x<-rep(x,length(y)) 
# lines(x,y, lty=2)
# x<-seq(0, 45, by=0.1)
# y<-rep(0.5, length(x))
# lines(x,y, lty=2)
# x<-locator(1) # After this lines, click on the point within the Graphik where the ML-Label should be displayed
# text(x$x, x$y, "ML", cex=1.2)


plot(survfit(fit1), conf.int="none", xlab="Time", ylab="cumulated hazard", fun="cumhaz")

# the DySeq package provides a function to compute the non-cumulated hazard:
NonCumHaz(survfit(fit1), plot=T) # Figure 5


# set the graphics device back at displaying one plot at a time
par (mfrow = c(1,1)) 



#################################################################
##  Cox-regression: Prediction of Hazard-ratio by a covariate  ##
#################################################################



# mean-centering EDCm  (men’s self-assessed dyadic coping ability)
EDCm.cent<-scale(mydata$EDCm, TRUE, FALSE)

# Fit the coxregression with EDCm.cent used as covariate
fit2 <- coxph(stress.surv~EDCm.cent, ties="breslow") # ~cent means: the hazard-ratio should be predicted by EDCm.cent
summary(fit2)
exp.b<-unclass(summary(fit2))$coefficients[2]


# Plotting the Coxregression with simple effects [note: not shown in the article]
# run:
# plot(survfit(fit2)$time, NonCumHaz(survfit(fit2)), xlim=c(19,45), ylim=c(0,0.3), type="l", xlab="Time", ylab="predicted hazard")
# lines(survfit(fit2)$time, NonCumHaz(survfit(fit2))*exp.b^sd(EDCm.cent), lty=4)
# lines(survfit(fit2)$time, NonCumHaz(survfit(fit2))*exp.b^-sd(EDCm.cent), lty=2)
# 
# x<-locator(1) # click into the graphics device to choose the location of the legend
# legend(x$x, x$y, c("X = 0", "X = +1sd", "X = -1sd"), lty=c(1,4,2), cex=0.7 )





#####
#
#
## Stoppt here!!!
#
#
#####






##################################################################
###  Research Question 4:                                      ###
###  Do couples differ from each other with respect to their   ###
###  typical behavioral patterns in times of stress?           ###
##################################################################


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
cluster2fac <- factor (cluster2, labels = c("cluster 1; fast coper", "cluster 2; slow coper"))
seqdplot (couple.seq, group = cluster2fac)



# repeating the sequence analysis of Bakeman and Gottman for sperated groups
# at the moment only implented for comparsion of two groups
# if more clusters are found run the procedure several time,
# and create a subgroup vector with a one's and tow's for the two clusters
# that should be compared an zeros for the other clusters!
# and be aware of alpha cumulation!
LogSeq(my.trans, delta=0.5, subgroups=cluster2)

LogSeq(my.trans.stress, delta=0.5, subgroups=cluster2)


# coxregression with subgroups
clust.dummy<-numeric(length(cluster2))
clust.dummy[cluster2==1]<-0
clust.dummy[cluster2==2]<-1 # recoding into dummy: first group is now coded 0, second is 1!


# adding to our regression with the prolonged covariate
fit5 <- coxph(stress.surv~mydata[,98]*as.factor(clust.dummy))
summary(fit5)

help(coxph)

cor(mydata[,98],clust.dummy)
summary(lm(mydata[,98]~clust.dummy))







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

