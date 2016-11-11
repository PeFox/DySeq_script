##########################################################################################
### Analyzing dyadic sequence data - research questions and implied statistical models ###
##########################################################################################

################
##  Preamble  ##
################

# This R-script provides a hands-on-tutorial for all
# analyses covered in 
#
# "Analyzing dyadic sequence data - research questions and implied statistical models"
# by  --- blinded for reviewing ---
# published in --- blinded --- 
#
# 
# Please make sure to install all required packages,
# including the "Dyseq" which provides the sample data!

### Content

### Prerequisite Steps

# packages from CRAN
# package from Github
# Example Data
# loading the data
# details on data


### Graphical Analysis

# state-distribution-plot
# entropy-plot
# Number of transitions


### Research question 1

# Pearson Correlation


### Research question 2:
  
## aggregated logit models
# step 1: state-transition tables
# step 2: multiple logit-regressions
# step 3: aggregating
# step 4: APIM


## Multi-Level-Approach
# converting sequences into MLM-data-structure
# applying MLM via lme4
# Basic Markov Modell
# converting data
# obtaining the transition matrix

### Research question 3:
  
#  estimating hidden Markov model


### Research question 4:
  
#  estimating mixture Markov model
#  sequence clustering
#  OM-distances
#  clustering
#  interpret clusters


##########################
##  Prerequisite steps  ##
##########################

##  make sure the following packages are installed:

# install.packages("TraMineR")      # for graphical analysis and research question 4
# install.packages("RColorBrewer")  # for grey-shaded graphics
# install.packages("gmodels")       # must be installed!
# install.packages("MASS")          # must be installed!
# install.packages("survival")      # must be installed for research question 3!
# install.packages("fpc")           # must be installed for research question 4!
# install.packages("cluster")       # must be installed for research question 4!
# install.packages("devtools")      # must be installed for installing packages from github
# install.packages("lme4")          # must be installed for the multi-level APIM
# install.packages("lmerTest")      # must be installed for the multi-level APIM

# loading packages!
library(cluster)      # ward algorithm, opt. number of clusters
library(fpc)          # optimal number of clusters
library(survival)     # cox-regression
library(TraMineR)     # graphical analysis, state-changes, entropy, OM-procedure
library(RColorBrewer) # only used for alternative graphics
library(devtools)     # needed for DySeq-installation



##  Installing DySeq from GitHub 

# remove.packages("DySeq")          # run only if a previous version of DySeq is installed and should be updated! 
# install_github("PeFox/DySeq")     # run only once for first installation!

library(DySeq)                      # loading DySeq;
                                    # for further detail use help(DySeq) and click on the help file's index!



#######################
##  2. Example data  ##
#######################
help(CouplesCope)
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

## Objects from previous sections are needed: 
#  - mydata     
#  - my.expand  


# state-distribution plot (using TraMineR) 
citation("TraMineR")       # please cite TraMineR if you use 
citation("RColorBrewer")   # some of its functions!

# create labels for plot
couple.labels <-c("none",     # no reaction
                  "SC only",  # only stress communication
                  "DC only",  # only dyadic coping
                  "SC+DC")    # stress and dyadic

# create a stslist object (TraMineR S3-Class)
couple.seq <- seqdef(my.expand,              # the combined states 
                     labels = couple.labels) # the label

# State-Distribution plot 
seqdplot(couple.seq,
         cex.legend=0.8) # adjust size of the legend

# Alternatively a grey version (using RColorBrewer)
# And legend aligned right
attr(couple.seq , "cpal") <- brewer.pal(4, "Greys") # see figure 2
seqdplot(couple.seq, cex.legend=0.8, withlegend="right")


# Figure 3 Entropy-plot and histogramm of "Number of transitions"

# Entropy-plot: how much do states diverge within time intervals?
Entropy <- seqstatd(couple.seq)$Entropy
plot(Entropy, main= "Entropy", col="black", xlab = "Time in 10 sec. intervall", type ="l")

# Histogramm of "Number of transitions": How often do couples change from one state into another
SeqL<-seqtransn(couple.seq)
hist(SeqL, main="Number of transitions", xlab="State-transitions")


summary(SeqL) # more details on the "Number of transitions"


# Note:
cor(SeqL, mydata$EDCm) # the Number of transition can be investigated further.
                       # For example, the higher men’s self-assessed dyadic coping  
                       # ability the less often couples change their states!
                    

rm(list=c("couple.labels", "Entropy")) # these two objects are not needed in the 
                                       # subsequent analyses!




#######################################
##  4. potential research questions  ##
#######################################


###################################################################
###  Research question 1:                                       ###
###  Is there an association between a particular behavior by   ### 
###  one and the reaction by the other partner?                 ###
###################################################################

## Objects from previous sections needed: 
#  - mydata     


# Step 1: Computing the frequencies of both behaviors,
# using the function NumbOccur from the DySeq-package
SC.sumscores<-NumbOccur(mydata[,2:49], 1, prop=FALSE) 
DC.sumscores<-NumbOccur(mydata[,50:97], 1, prop=FALSE)    


# Explanation:
# Argument 1: Columns 2:49 are containing the stress reactions (SC), 50 to 97 the dyadic coping reactions (DC).
# Argument 2: We are interested in the state '1', because every '1' represents that SC / DC was shown
# Argument 3: # prop=FALSE -> absolute frequencies 
              # prop=TRUE  -> relative frequencies



# Step 2: Correlation of both frequencies
cor.test(stress.sumscores, DC.sumscores)

# plotting the correlation [Note: not shown in the article]
plot(stress.sumscores, DC.sumscores, ylab="sum of DC", xlab="sum of SC")
abline(lm(stress.sumscores~DC.sumscores))


rm(list=c("DC.sumscores", "stress.sumscores")) # these two objects are not needed in the 
                                               # subsequent analyses and therefore are removed!


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


# If a single case should be inspected, the row number is needed. 
# For example, if the couple with ID 129 should be inspected, one could 
# use which() to determine the row number of ID 129 from the dataset:
which(mydata$code==129) # gives the rownumber 41  
                        
my.trans[[41]] # inspecting couple 129 (rownumber 41)
               # btw: the transistion-tables are stored within a list; 
               # therefore double-brackets are needed for subscripting!

# Or, if relative frequencies are preferred: 
round(my.trans[[41]]/sum(my.trans[[41]]),3) 



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
# A value of zero should be used if no zero frequencies exist, or if
# all cases with zero frequencies should be stripped from the analysis 
# (listwise deletion).
# single.case=TRUE computes p-values for every single analysis!
# Single cases can be extracted via the single.LogSeq function (see l. 280)
# However, p-values for the single case analysis may not be trustworthy 
# at the moment (difer from the results of other software implementations by 
# the forth digit for this example). P-values for the overall effects are checked 
# and match the results obtained by other software implementations (SPSS, Lem)

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



# removing objects that are not necessary for the subsequent analyses
rm(list=c("my.logseq", "my.logseq.stress")) 





###############################################################
###  Research Question 3:                                   ###
###  What is the typical duration of stress communication?  ###
###  And does it depend on covariates, for instances men’s  ###
###  self-assessed dyadic coping ability?                   ###
###############################################################


## Objects from previous sections needed: 
#  - mydata     


# Calculating the last occurance of stress communication (last.stress)

last.stress<-LastOccur(mydata[,2:49],y=1)
# mydata[,2:49] accesses only the SC-sequences of mydata
# second Argument y specifies if the last observed '0' or '1' should be optained
# Note: other values for y are possible, but most other functions 
#       of the DySeq-package will only support dichotomous sequences.



####################################
#  Hazard, survival and cumhazard  #
####################################

# First a variable is needed that indicates if the event was shown (1; observed) or not (0; censored)
event<-c(rep(1, length(last.stress))) # it is shown for every case!
event[last.stress>=48]<-0 # exept the ones that show stress communication till time intervall 48

# the mean of the new vector "event" equals the relative frequency of couples in which the 
# event occored. One minus that frequency equals the relative frequency of censored cases.
1-mean(event)

# The duration and the event variable are combined in a Surv-object
stress.surv<-Surv(last.stress,event) 

fit1 <- coxph(stress.surv~1, ties="breslow") # ~1 means: fittet without a covariate.
                                             # Different estimators exist for handling 
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
##  Cox-regression: Prediction of hazard-ratio by a covariate  ##
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



# removing objects that are not necessary for the subsequent analyses
rm(list=c("event", "fit1", "fit2", "last.stress", "exp.b", "EDCm.cent")) 







##################################################################
###  Research Question 4:                                      ###
###  Do couples differ from each other with respect to their   ###
###  typical behavioral patterns in times of stress?           ###
##################################################################


## Objects from previous sections needed: 
#  - mydata     
#  - my.expand  
#  - couple.seq (see line 131)
#  - SeqL (see line 150)
#  - my.trans (see line 217)
#  - my.trans.stress(see line 293)


# Introductionary: 
summary(SeqL) # Many couples change their behaviors from time 
              # interval to time interval, but others seem to
              # remain rather stable in their behaviors (range: 5-29)




####################
##  OM-Distances  ## 
####################


# Two sequences differ to the extent to which some elements of one sequence 
# have to be changed in order to perfectly match the other sequence (cost) 
# Insertion and deletion costs 1 (by default)
# Substitution depends on the substitution-cost-matrix 
# (Substitution-cost-matrix is derived by Gabadinho's TRATE-Formula)
submat <- seqsubm(couple.seq, method = "TRATE")


# The TRATE-Formula uses the transition probabilities
# to determine the substitution-cost-matrix (table 6 in the accompanying article)
round(seqtrate(couple.seq),2) # or the transition probabilities
                              # the round() funciton is optional!

# Inspecting the substitution-cost-matrix (table 7 in the article)
round(submat,2) # the round() funciton is optional!
# state 0 == No SC/DC
# state 1 == SC
# state 2 == DC
# state 3 == SC+DC

# The substitution-cost-matrix is then used to calculate the distance-matrix
dist.oml <- seqdist(couple.seq, method = "OM", sm = submat)

# The results are stored in a distances matrix, here: "dist.oml"
# with as many rows and columns as number of observations
# cells representing the minimal cost between the associated observation units
# The dissimilarity matrix corresponds to other distance measures, 
# in an ordinary cluster analysis; 
# therefore ordinary clustering algorithms can be apllied. 



############################################ 
##  Determine optimal Number of clusters  ##
############################################


# First of all, the optimal number of clusters needs to be determined: 
plot(pam(dist.oml, pamk(dist.oml)$nc), which.plot=1)


#################################################################
##  Explanation of the two plots, generated by plot(pam(...))  ##
#######################################################################################################
#                                                                                                    ##
## First Graph: [note: not shown in the article]                                                     ##
#                                                                                                    ##
# The distances describes a number of sequences minus 1 dimensional space, and therefore it          ##
# it is impossible to plot them. However, the first plot projects the dissimilarities in a           ##
# two-dimensional space. By this, information is lost and the distances shown in the plot won't      ## 
# fit the distances in distance-matrix perfectely. Yet, this plot provides us with a first           ##
# glimpse into the cluster structure: here we see two clusters. However they are not completely      ##
# distinct, because we caan see that they overlap a little bit! Furthermore, the smaller cluster     ##
# (left) seem to be more homogenous (points are very close to each other) and the bigger one         ##
# (right) seem to be more heterogenous (more dispersion).                                            ##
#                                                                                                    ##
#                                                                                                    ##
## Second Graph: silhouette plot [note: not shown in the article]                                    ##
#                                                                                                    ##
# The silhuette (s) is the distance of one object (sequence) to its clusters centroid minus          ##
# the distance to the nearest clusters centroid devided by the maximal possible distance.            ##
# Therefore, the higher the mean s-values of a cluster, the stronger the found cluster structure.    ##
#                                                                                                    ##
# [see: Peter J. Rousseeuw, Silhouettes:                                                             ##
#       A graphical aid to the interpretation and validation of cluster analysis,                    ##
#       Journal of Computational and Applied Mathematics, Volume 20, 1987, Pages 53-65]              ##
#                                                                                                    ##
#                                                                                                    ##
# The plot contains the following information:                                                       ##
# - 2 clusters are identified                                                                        ##
# - First one includes 38 observations with s=.43                                                    ##
# - Second one includes 26 observations with s=.59                                                   ##
# - The overall mean s is .5                                                                         ##
#                                                                                                    ##
# - typical rules of thumb are: 1 >= s >.75 strong structure                                         ##
#                               0.75 >= s > 0.50 medium structure                                    ##
#                               0.50 >= s > 0.25 week structure                                      ##
#                               0.25 >= s > 0 no real structure                                      ##
#                                                                                                    ##
#######################################################################################################


# The mean silhuette value indactes an overall weak cluster structure.
# Therefore, additional methods should be used to determine the optimal 
# number of clusters. Here, for example, screeplot/dendrogramm:

# Screeplot: (Indicates 1 or 2 clusters)
wss <- (nrow(dist.oml)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")

# Dendrogramm: (Indicates 2, maybe 3, Clustersolution)
plot(agnes (dist.oml, diss=TRUE, method = "ward"), which.plots=2)



######################
##  Ward-algorithm  ##
######################

clusterward1 <- agnes (dist.oml, diss=TRUE, method = "ward") 

# 2 Clustersolution [for other solution change the value behind the argument 'k']
cluster2 <- cutree (clusterward1, k=2)
cluster2fac <- factor (cluster2, labels = c("cluster 1; fast coper", "cluster 2; slow coper"))

# [Note: The n of both clusters doesn't match exactly the n from the silhuette test! 
#        This is because the latter depends on another algorithm. Typically they provide
#        the same cluster solutions. However, in this case one sequence was assgigned differently. 
#        The first plot from line 478 shows that exactly one sequence is located in the overlapping area
#        between both clusters. Because of that, some algorithms assign it to the fist, 
#        other to the second cluster. However, changing its cluster - or deleting that
#        observation - would not lead to other results. 





#####################################################
###  Cluster interpretation and further analyses  ###
#####################################################


## separate state-distribution plots 
#  [note: figure 6 in the accompanying article]
seqdplot (couple.seq, group = cluster2fac) 

## coxregression with subgroups
fit3 <- coxph(stress.surv~mydata$EDCm+cluster2fac)
summary(fit3)

## correlation between the cluster membership 
#  and men’s self-assessed dyadic coping ability 

clust2.dummy<-as.numeric(cluster2fac) # 
clust2.dummy[clust2.dummy==1]<-0      # Dummycoding the factor as a numeric vector
clust2.dummy[clust2.dummy==2]<-1
cor.test(mydata$EDCm,clust2.dummy) # note: pearson correlation between a dichotomous
                                   #       and an interval scaled variable equals
                                   #       a point biseral correlation and is therefore
                                   #       the appropriate association measure!


## separate aggregated logit models for both clusters
#  [note: shown in table 8]

# LogSeq provides an optional argument "subgroups" that allows to compare 
# transition tables for two groups 
LogSeq(my.trans, delta=0.5, subgroups=cluster2) # Analysis for DC as dependend variable



LogSeq(my.trans.stress, delta=0.5, subgroups=cluster2) # Analysis for SC as dependend variable 
                                                       # [note: not shown in the article]


# removing surplus objects 
rm(list=c("dist.oml", "SeqL", "submat", "clust2.dummy", 
          "cluster2", "clusterward1", "fit3", "my.trans", 
          "my.trans.stress", "stress.surv", "wss", "cluster2fac")) 




##########################
##                      ##
##  Additional content  ##
##                      ##
###############################################################
###  Using EstFreq and EstTime to determine the number      ###
###  of expected zero and low frequencies cells             ###
###############################################################

# First step: Define a matrinx conatining the expected transition rates!
my.trans.table<-matrix(c(0.57, 0.13, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05),4,2)

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

# if legend position should be change:
attr(my.EstTime.plot, "pos")<-"topright" # alternatives are: bottomleft (default)
my.EstTime.plot                          #                   bottomright
                                         #                   topleft         


