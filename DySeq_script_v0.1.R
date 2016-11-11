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


### Content                                       #lines

### Prerequisite Steps                             76
# packages from CRAN
# package from Github

### Example Data                                  114
# loading the data
# details on data


### Graphical Analysis                            137
# state-distribution-plot
# entropy-plot
# Number of transitions


### Research question 1                           197
# Pearson Correlation


### Research question 2:                          233

## aggregated logit models                        246
# step 1: state-transition tables
# step 2: multiple logit-regressions
# step 3: aggregating
# step 4: APIM

## Multi-Level-Approach                           341
# converting sequences into MLM-data-structure
# applying MLM via lme4

## Basic Markov Modell                            428
# converting data
# obtaining the transition matrix


### Research question 3:                          449
#  estimating hidden Markov model


### Research question 4:                          594
  
##  estimating mixture Markov model               604

##  Sequence clustering                           547
#   OM-distances
#   clustering
#   display clusters


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
#  - couple.seq (a a stslist object from TraMineR, created in line 156)


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


##################
#  Step 4: APIM  #
##################

# Rerun the procedure with stress as dependend variable!

my.trans.SC<-StateTrans(my.expand, TRUE) # This time TRUE means: First Sequence (Stress reaction) will be used as dependend variable!
my.logseq.SC<-LogSeq(my.trans.SC, delta=0.5)

my.logseq.SC # estimates for actor- and partnereffect on stress
my.logseq # estimates for actor- and partnereffect on coping
# (see fig. 4)



# removing objects that are not necessary for the subsequent analyses
rm(list=c("my.logseq", "my.logseq.stress")) 




##########################
##########################
## Multi-Level Approach ##
##########################
##########################

# Make sure all needed packages are loaded!
#library("lme4")
#library("lmerTest")

# ML_Trans transforms dyadic sequences into multi-level data!
 ML_data<-ML_Trans(data=CouplesCope,     # The data, which should be used!
                   first=2:49,           # The sequence, which should be used as the DV
                   second=50:97)         # The sequence, which serves as the IV!



 # transitions must be recoded first into lagged actor and lagged partner effects. 
 MLAP_data<-MLAP_Trans(ML_data) 
 

 # Inspect the data:
 # View(MLAP_data)
 
 
# Adding labels:
# In example data first seq referred to females
# and second to males
names(MLAP_data)[1]<-"sex"
MLAP_data$sex<-as.factor(MLAP_data$sex)
levels(MLAP_data$sex)<-c("female", "male")
 
# Adding effect-coding
MLAP_data$Partner[MLAP_data$Partner==0]<-(-1)
MLAP_data$Actor[MLAP_data$Actor==0]<-(-1)
 
 
### Full Random Model
# Warning: Estimation can take a long time!
set.seed(1234)
glmer(DV~1+sex+Actor+Partner+Actor*Partner+
       sex*Actor+sex*Partner+sex*Actor*Partner+
       (1+sex+Actor+Partner+Actor*Partner+
       sex*Actor+sex*Partner+sex*Actor*Partner|ID),
       data=MLAP_data,
       family=binomial)
#AIC 5256.744
#BIC 5551.640
 
# The most simple MLM (Random intercept only)
set.seed(1234)
glmer(DV~1+sex+Actor+Partner+Actor*Partner+
       sex*Actor+sex*Partner+sex*Actor*Partner+
       (1|ID),
       data=MLAP_data,
       family=binomial)
# AIC 5299.386  
# BIC 5359.706
 
 
### The Model from the article (----)
 
# Random Actor und Partner Effekte
set.seed(1234)
fit<-glmer(DV~1+sex+Actor+Partner+Actor*Partner+
          sex*Actor+sex*Partner+sex*Actor*Partner+
          (1+Actor+Partner|ID),
          data=MLAP_data,
          family=binomial)
summary(fit)
AIC(fit)
BIC(fit)
# AIC 5222.063
# BIC 5315.893
 
# See Article for a close interpretation of the
# last model. 
 



##########################
##########################
## Basic Markov Model   ##
##########################
##########################

# First possibility 
round(seqtrate(couple.seq),2)

# Second possibility: Use seqHMM
# The second way, using the seqHMM-package is a little bit more complicated
# at the first glance. However, its worth to try this package too, 
# because hidden Markov and mixture Markov model follow the exact same logic.
# Please cite seqHMM if you are using it for your analysis!
citation("seqHMM")  


# starting values for transition probabilities
mytrans<-matrix(.25, 4,4)

# starting values for the initial probabilities
myinit<-c(.25,.25,.25,.25)

mybuild<-build_mm(couple.seq,
         mytrans,
         myinit)

fit<-fit_model(mybuild)

BIC(fit$model)
AIC(fit$model)
 

##################################################################
## Research Question 3: Is there an underlying dyadic process,  ##
## which might account for the observed behavior?               ##
##################################################################

# Needs following objects:
# couple.seq (created in the graphics section, using TraMineR-package)

# Following packages needed:
# seqHMM



# This question can be answered by a hidden Markov model. 

# first we have to specify starting values for the hidden chain
# by doing so, we also decide on the number of latent states!

myhtrans<-matrix(c(.50, 0, 0.5, 1), 2,2)
# In this example we created a 2*2 matrix, the second 
# state is an absorbing one, because we resticted the transition
# probabilities in a way that one cannot leave the second state!
           
my_emission<-matrix(c(.25),2,4)
# the emission matrix has one row per latent state
# and four columns for eacht observed state

myinit2<-c(.50, .50) 
# Because the hidden chain is now a 2*2 matirx, we only need two
# initial probabilities

my_hmodel<-build_hmm(couple.seq,
                       myhtrans,
                    my_emission,
                        myinit2)

fit2<-fit_model(my_hmodel)

fit2

##################################################################
###  Research Question 4:                                      ###
###  Do couples differ from each other with respect to their   ###
###  typical behavioral patterns in times of stress?           ###
##################################################################


# his question can be answered by a mixture Markov model or by the OM-procedure 

# Needs following objects:
# couple.seq (created in the graphics section, using TraMineR-package)

# Following packages needed:
# seqHMM


# First of all we have to specify starting values for the each chain
# so we need as many chains as there a latent groups. For the sake of 
# examplification we will assume two groups:

# These are transition matrices for the observed states,
# therefore they have the as many rows and columns as observed states:
mytrans1<-matrix(c(.25), 4,4)
mytrans2<-matrix(c(.25), 4,4)
# Those transition matrices must be placed into a list, before
# building the model
mymixtrans<-list(mytrans1, mytrans2)

# There are not emissions because there are no hidden states! 

# However, we need to sets of starting values. One for each chain!
myinit1<-c(.25, .25, .25, .25) 
myinit2<-c(.25, .25, .25, .25) 
# again both must be placed inside a list:
mymixinit<-list(myinit1, myinit2)


my_mmodel<-build_mmm(couple.seq, 
                    mymixtrans,
                    mymixinit)

fit3<-fit_model(my_mmodel,                   
                global_step=TRUE,     # additional arguments for the optimizier
                local_step=TRUE,      # to avoid local maximum
                control_em=list(restart = list(times=10)))

fit3




##################
## OM-procedure ##
##################

## Objects from previous sections needed: 
#  - mydata     
#  - my.expand  
#  - couple.seq (see line 131)
#  - SeqL (see line 150)
#  - my.trans (see line 217)
#  - my.trans.SC(see line 293)


# Introductionary: 
summary(SeqL) # Many couples change their behaviors from time 
              # interval to time interval, but others seem to
              # remain rather stable in their behaviors (range: 5-29)




####################
##  OM-Distances  ## 
####################


# Substitution-cost-matrix is derived by Gabadinho's TRATE-Formula
submat <- seqsubm(couple.seq, method = "TRATE")

# The substitution-cost-matrix is then used to calculate the distance-matrix
dist.oml <- seqdist(couple.seq, method = "OM", sm = submat)


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





#################################################
###  comparing clusters and further analyses  ###
#################################################


## separate state-distribution plots 
seqdplot (couple.seq, group = cluster2fac) 


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



LogSeq(my.trans.SC, delta=0.5, subgroups=cluster2) # Analysis for SC as dependend variable 
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

# First step: Define a matrix containing the expected transition rates!
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


