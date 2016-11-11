#Analyzing dyadic sequence data - research questions and implied statistical models 

Please note, a full R-script of this vignette is provided under the name DySeq_script_v0.1.R in this repository!

This Vignette provides a hands-on-tutorial for all analyses covered in 

"Analyzing dyadic sequence data - research questions and implied statistical models"
by  --- blinded for reviewing ---
published in --- still in review ---

Please make sure to install all required packages,
including the "Dyseq", which provides the sample data!


__Content__                  

<ol type="a">
  <li>Prerequisite Steps</li>
      <ul>
      <li>packages from CRAN</li>
      <li>DySeq-package from Github oranges</li>
      </ul>
      
  <li>Example Data</li>
      <ul>
      <li>Loading the Data</li>
      <li>Details on Data </li>
      </ul>
      
  <li>Graphical Analysis</li>
      <ul>
      <li>State-Distribution Plot</li>
      <li>Entropy Plot </li>
      <li>Number of Transitions</li>
      </ul>  
      
  <li>Research question 1 </li>
      <ul>
      <li>Pearson Correlation</li>
      </ul>
      
  <li>Research Question 2 </li>
      <ul>
      <li>Aggregated Logit Models</li>
          <ul>
          <li>step 1: State-Transition Tables</li>
          <li>step 2: Multiple Logit-Regressions</li>
          <li>step 3: Aggregating</li>
          <li>step 4: APIM</li>
          </ul>
      <li>Multi-Level APIM</li>
          <ul>
          <li>Data praparation</li>
          <li> Applying MLM via lme4</li>
          </ul>
      <li>pears</li>
      </ul>
      
  <li>Research question 3 </li>
  <li>Research question 4 </li>
  <li>Additional function </li>
</ol>

* (A) Prerequisite Steps (A)                
   * packages from CRAN                   
   * package from Github                  

* (B) Example Data                       
  * loading the data                     
  * details on data                      

* (C) Graphical Analysis                 
  * state-distribution-plot              
  * entropy-plot                         
  * Number of transitions                

* (D) Research question 1:
  * Pearson Correlation 
 
* (E) Research question 2:                   
  * aggregated logit models
    * step 1: state-transition tables
    * step 2: multiple logit-regressions
    * step 3: aggregating
    * step 4: APIM
  * Multi-Level-Approach
    * converting sequences into MLM-data-structure
    * applying MLM via lme4
  * Basic Markov Modell
    * converting data
    * obtaining the transition matrix

* (F) Research question 3:                 
  * estimating hidden Markov model

* (G) Research question 4:
  * estimating mixture Markov model
  * sequence clustering
    * OM-distances
    * clustering
    * interpret clusters
    
* (H) Additional functions
  * number of expected cell problems
  * find number of needed time intervals
  
---

## A Prerequisite steps

make sure the following packages are installed:

``` r
install.packages("TraMineR")      # for graphical analysis and research question 4
install.packages("RColorBrewer")  # for grey-shaded graphics
install.packages("gmodels")       # must be installed!
install.packages("MASS")          # must be installed!
install.packages("survival")      # must be installed for research question 3!
install.packages("fpc")           # must be installed for research question 4!
install.packages("cluster")       # must be installed for research question 4!
install.packages("devtools")      # must be installed for installing packages from github
install.packages("lme4")          # must be installed for the multi-level APIM
install.packages("lmerTest")      # must be installed for the multi-level APIM
install.packages("seqHMM")        # must be installed for Markov models
```

Make sure to delete previous Versions of DySeq before installing a new version from GitHub!
``` r
remove.packages("DySeq")      # remove older Version von DySeq     
install_github("PeFox/DySeq") # Install DySeq from GitHub
```

Alternatively, you can install `DySeq` from CRAN. However, the version on CRAN does not feature 
the multi-level approach until the next release!

---

## B. Example data

The data stems from a study, which was promoted as a study on close relationship and stress. 
Each row represents one dyad, each of them containing two sequences. Dyads are 64 heterosexuel
couples. The females partners were stressed using the Trier Social Stress Test 
(TSST; Kirschbaum, Pirke, Hellhammer, 1993). Directly after the stress induction, both partners joint again and the couple was left alone for eight minutes. During this period (a 'fake' waiting condition) the two partners were filmed for 8 minutes divided into 48 intervals of ten seconds length. It was coded if the female partners showed stress communication (SC) within an interval (sequence 1; Colums 50:97) and if the male partner showd dyadic coping reactions (DC; sequence 2; columns 2:49). For rurther insides about dyadic coping and/or stress communication, see Bodenmann (2015).

```r
library(DySeq)        # loading the DySeq Package
mydata<-CouplesCope   # getting the data
help(DySeq)           # get more information about the data set

# Most following approaches will need combined states!
# The function StateExpand can combine tow sequences into one!
```r
my.expand<-StateExpand(CouplesCope, 2:49, 50:97)
```

---
# C. Graphical Analysis

Following objects from the previous sections are needed:
- mydata      (the example data)
- my.expand   (the combined sequences)

Following Packages are needed:
- DySeq
- TraMineR
- RColorBrewer

```r
library(TraMineR)     # awesome package for sequence analysis in general!
library(RColorBrewer) # more colours!
citation("TraMineR") # please cite Packages if you use them!
citation("RColorBrewer")
```

### The state distribution plot 


```r
# Create labels for plot
couple.labels <-c("none",     # no reaction
                  "SC only",  # only stress communication
                  "DC only",  # only dyadic coping
                  "SC+DC")    # stress and dyadic

# Create a stslist object (TraMineR S3-Class)
couple.seq <- seqdef(my.expand,              # the combined states 
                     labels = couple.labels) # the label

# State-Distribution plot 
seqdplot(couple.seq,
         cex.legend=0.8) # adjust size of the legend

# Alternatively a grey version (using RColorBrewer)
# And legend aligned right
attr(couple.seq , "cpal") <- brewer.pal(4, "Greys") # see figure 2
seqdplot(couple.seq, cex.legend=0.8, withlegend="right")
```

### Entropy plot
```r
Entropy <- seqstatd(couple.seq)$Entropy
plot(Entropy, main= "Entropy", col="black", xlab = "Time in 10 sec. intervall", type ="l")
```

### Histogramm of transitions number
```r
SeqL<-seqtransn(couple.seq)
hist(SeqL, main="Number of transitions", xlab="State-transitions")
```

---

## D. Research Question 1: Is there an association between a particular behavior by one and the reaction by the other partner?

Following objects from the previous sections are needed:
- mydata      (the example data)

Following Packages are needed:
- DySeq

```r
# NumbOccur counts how often a certain behavior is shown within each sequence
SCstress.sumscores<-NumbOccur(x=mydata[,2:49],           # col: 2:49 represent stress communication (SC)
                                          y=1,           # 0=no SC; 1=SC was shown 
                                          prop=FALSE)    # absolute (TRUE) or relative (FALSE) frequency?
                                        
DC.sumscores<-NumbOccur(mydata[,50:97], 1, prop=FALSE)  # Same for dyadic coping (DC)

# Correlation of both frequencies
cor.test(SC.sumscores, DC.sumscores)

# plotting the correlation
plot(SC.sumscores, DC.sumscores, ylab="Number of DC", xlab="Number of SC")
abline(lm(SC.sumscores~DC.sumscores))
```

![alt text](https://github.com/PeFox/DySeq_script/blob/master/Scatter_plotDCSC.png "Scatterplot(SC,DC)")


---

## E. Research Question 2: Does the behavior of one member trigger an immediate reaction by the other?

Following objects from the previous sections are needed:
- mydata      (the example data)
- my.expand   (the combined sequences)
- couple.seq  (a a stslist object from TraMineR)

Following Packages are needed:
- DySeq
- lme4      (for MLM approach)
- lmerTest  (for MLM approach)
- TraMineR  (for the basic Markov model)
- seqHMM    (For the basic Markov model)


### Bakeman & Gottman approach: Step 1

First step is transforming the sequences into transition tables.
That can be done by DySeq's function StateTrans, which stores the transition tables in a list of the class 'state.trans'. 

```r
my.trans<-StateTrans(my.expand,   # the combined sequences
                     first=FALSE) 
```
The argument 'first' specifies which sequence should be used as a dependend variable. 
When my.expand was created, SC was defined as the first sequence and DC as the second sequence. Therefore, DC is now the dependend variable of the transition tables. 

```r
# Just print a state.trans object to get the relative frequencies accross all dyads!
my.trans

# A single case can be plottet by using the [[]] operator
my.trans[[1]] # inspects the first table!

# If the original data.frame containes a ID-variable, the following can be used:
ID<-mydata$code             # stores the ID-variable in the object ID
my.trans[[which(ID==129)]]  # inspects the dyad with ID=129

# If relative frequencies are preferred for single case analysis, just divide the transition table by its sum.
# The following shows an example for the 41th transition table, which belongs to the couple with ID 129. 
# round(x,3) rounds the frequencies to three digits
round(my.trans[[41]]/sum(my.trans[[41]],3)
```

### Bakeman & Gottman approach: Step 2 

Second step is computing a logit model for each dyad. That is a very cumbersome procedure, yet the Package DySeq provides the function LogSeq for doing so with only one function!

```r
my.logseq<-LogSeq(my.trans,         # a list containing transition plots (from step 1)
                 delta=0.5,         # adds a frequency of .5 on every cell. Needed if zero cells exist!
                 single.case=TRUE)  # if TRUE, stores single case results (and the function becomes slower)
```

If a researcher is interested in single case analysis, the function single.LogSeq() can be used to obtain these in a ready-to-interpret table! p-value test whether betas are different from zero. Note that the p-value for the intercept is not implemented yet! 

```r
# Single case analsis for transition plot 41 (aka couple 129)
single.LogSeq(my.logseq, 41)
```


### Bakeman & Gottman approach: Step 3

Next step is aggregating the results. This is automatically done by printing the results of the last step! 
No further functions are needed! Again, p-values test whether the betas are not equal zero. 

```r
my.logseq  # prints the aggrregated results!
```

Plotting the results will produce an interaction-plot. Mapping the probabilities of showing the dependend variable against the combinations of previous behavior:

```r
plot(my.logseq)
```

![alt text](https://github.com/PeFox/DySeq_script/blob/master/Interac.png "Scatterplot(SC,DC)")


### Bakeman & Gottman approach: Step 4

Rerunning the procedure a second time but switching the dependend variable. 

```r
my.trans.SC<-StateTrans(my.expand, first=TRUE)   # This time, the first sequence is the DV
my.logseq.SC<-LogSeq(my.trans.SC, delta=0.5) 

my.logseq.SC # Contains actor and partner effects for the female
```     

Our previous results from step 3 contained the effects for men, while our second results (my.logseq.SC) contain the effects for women. Thus, combining both results in a APIM as shown in the article. 

---


### Multi-level approach: Step 1

The first step is data preparation. First each transition must be recoded to be a single observation within a nested data structure. Transitions are level-1 observations, which are nested within dyads. That can be achieved by one single function from the DySeq-Package!

```r
library("lme4")      # Make sure all needed packages are loaded!
library("lmerTest")

# ML_Trans transforms dyadic sequences into multi-level data!
ML_data<-ML_Trans(data=CouplesCope,     # The data, which should be used!
                        first=2:49,     # The sequence, which should be used as the DV
                        second=50:97)   # The sequence, which serves as the IV!
```

If the data should be used to apply a multi-level APIM, transitions must be recoded first into lagged actor and lagged partner effects. The function MLAP_Trans does so. 

```r
 MLAP_data<-MLAP_Trans(ML_data) # ML_data must be the output of ML_Trans!
```

Labels should be added or else the procedure can become confusing later on. 

```r
names(MLAP_data)[1]<-"sex"
MLAP_data$sex<-as.factor(MLAP_data$sex)
levels(MLAP_data$sex)<-c("female", "male")
```    

MLAP_Trans uses dummy-coding per default. However, for the purposes of an APIM effect-coding is better, because it is easier to interpret. Furthermore, effect coding was also used in the article (-----). So for the sake of comparable results, we will stick to it for this case.  

```r
MLAP_data$Partner[MLAP_data$Partner==0]<-(-1)
MLAP_data$Actor[MLAP_data$Actor==0]<-(-1)
```

### Multi-level approach: Step 2 

The second stept is applying and testing MLM-models. There are a vast and increasing number of packages in R, wich can run multi-level modells. However, lme4 became one of the best known packages for multi-level analysis, and an increasing number of tutorials are spreading through the net. Thus, we will stick to lme4, too. Do not forget to cite lme4 and lmerTest if you use this approach!


The following shows the most complex modell, which is possible to estimate. There will be some estimation problems with this model, but it will serve as an example to explain the function's arguments. 
```r
set.seed(1234)                                             # setting SEED for replication purposes!
fit<-glmer(DV~1+sex+Actor+Partner+Actor*Partner+           # intercept, Actor, Partner and interaction effect for the referrence group
        sex*Actor+sex*Partner+sex*Actor*Partner+           # difference for the non-referrence group (here males/DC)
        (1+sex+Actor+Partner+Actor*Partner+                # Random effects for intercept, Actor, Partner and interaction effect
           sex*Actor+sex*Partner+sex*Actor*Partner|ID),    # Random effects for the differences between the DVs (SC vs. DV)
      data=MLAP_data,                                      # the actual data 
      family=binomial)                                     # provides Link-function, so logistig regression is applied!
summary(fit)                                               # Provides summary of results
AIC(fit)                                                   # Akaike information criterion (AIC; a comparative fit index)
BIC(fit)                                                   # Bayesian information criterion (BIC; a more conservative fit index)
```

Models can be compared using AIC or BIC (comparative fit-indices). Smaller values indicate better model fit. In our case, a model containting random effects for the intercept, actor, partner effects were the best fitting one. 

```r 
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
```

 For a closer interpretation of this model, inspect the article (----). 
 
 ---
 
 
### Basic Markov model

The TraMineR-package provides a function to fit a basic Markov model. The couple.seq object is needed, which was created in the 'graphic analysis' section!

```r
round(seqtrate(couple.seq),2) # the round command is optional, and rounds the transition matrix to two digits
```

The second way, using the seqHMM-package is a little bit more complicated at the first glance. However, its worth to try this package too, because hidden Markov and mixture Markov model follow the exact same logic. Please cite seqHMM if you are using it for your analysis!

```r

# First of all, starting values for the transition matrix must be specified
# each row must add up to one! In this case we assume that each transition is equally likely!
# One of the advantages of seqHMM is, that restriction can be added (e.g. setting a value to zero
# will set it zero).
mytrans<-matrix(.25, 4,4)

# Same for the initial probabilities
myinit<-c(.25,.25,.25,.25)

# Builds the acutal model
mybuild<-build_mm(couple.seq,
         mytrans,
         myinit)

# Fits the model on the data
fit<-fit_model(mybuild)

# print results
fit

# get comparative fit indices!
AIC(fit$model)   # Akaike information criterion (AIC; a comparative fit index)
BIC(fit$model)   # Bayesian information criterion (BIC; a more conservative fit index)
```

---

## Research Question 3: Is there an underlying dyadic process, which might account for the observed behavior?             

Needs following objects:
couple.seq (created in the graphics section, using TraMineR-package)

Following packages needed:
seqHMM



This question can be answered by a hidden Markov model using seqHMM. Please cite seqHMM if you are using it for your analysis! 
First we have to specify starting values for the hidden chain by doing so, we also decide on the number of latent states!

```r
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

fit2  # print the results
```

---

## Research question 4: Are there latent groups of dyads, which might account for observing different reaction patterns?

This question can be answered by a mixture Markov model or by the OM-procedure. We start with 

## Mixture Markov model: 

```r
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
```
---

## OM-procedure

Objects from previous sections needed: 
- mydata     
- my.expand  
- couple.seq 
- my.trans 
- my.trans.SC

Packages needes: 
- DySeq
- TraMineR


###  Distances   

Substitution-cost-matrix is derived by Gabadinho's TRATE-Formula
```r
submat <- seqsubm(couple.seq, method = "TRATE")
```

The substitution-cost-matrix is then used to calculate the distance-matrix
```r
dist.oml <- seqdist(couple.seq, method = "OM", sm = submat)
```



###  Determine optimal Number of clusters 

First of all, the optimal number of clusters needs to be determined: 
```r
plot(pam(dist.oml, pamk(dist.oml)$nc), which.plot=1)
```

First Graph: [note: not shown in the article]                                                     
                                                                                                   
The distances describes a number of sequences minus 1 dimensional space, and therefore it is possible to plot them. However, the first plot projects the dissimilarities in a two-dimensional space. By this, information is lost and the distances shown in the plot won't    fit the distances in distance-matrix perfectely. Yet, this plot provides us with a first limpse into the cluster structure: here we see two clusters. However they are not completely distinct, because we caan see that they overlap a little bit! Furthermore, the smaller cluster (left) seem to be more homogenous (points are very close to each other) and the bigger one (right) seem to be more heterogenous (more dispersion).                                          
                                                                                                  
Second Graph: silhouette plot [note: not shown in the article]                                   

The silhuette (s) is the distance of one object (sequence) to its clusters centroid minus          
 the distance to the nearest clusters centroid devided by the maximal possible distance.            
Therefore, the higher the mean s-values of a cluster, the stronger the found cluster structure.    
                                                                                                    
[see: Peter J. Rousseeuw, Silhouettes:                                                             
A graphical aid to the interpretation and validation of cluster analysis,                   
Journal of Computational and Applied Mathematics, Volume 20, 1987, Pages 53-65]             

The plot contains the following information:                                                       
- 2 clusters are identified                                                                       
- First one includes 38 observations with s=.43                                                   
- Second one includes 26 observations with s=.59                                                 
- The overall mean s is .5                                                                        
                                                                                                  
- typical rules of thumb are: 1 >= s >.75 strong structure                                         
                              0.75 >= s > 0.50 medium structure                                    
                              0.50 >= s > 0.25 week structure                                     
                              0.25 >= s > 0 no real structure                                    
                                                                                                   
The mean silhuette value indactes an overall weak cluster structure. Therefore, additional methods should be used to determine the optimal number of clusters. Here, for example, a screeplot and adendrogramm:

Screeplot: (Indicating 1 or 2 clusters)
```r
wss <- (nrow(dist.oml)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
```

Dendrogramm: (Indicating 2, maybe 3, Clustersolution)
```r
plot(agnes (dist.oml, diss=TRUE, method = "ward"), which.plots=2)
```

### Ward-algorithm

After determing the number of clusters, the ward-algorithm can be used for clustering. 

```r
clusterward1 <- agnes (dist.oml, diss=TRUE, method = "ward") # the algorithm 
cluster2 <- cutree (clusterward1, k=2)                       # saving the two-cluster solution

# adding labels
cluster2fac <- factor (cluster2, labels = c("cluster 1; fast coper", "cluster 2; slow coper"))
```

Note that the number of observational units of both clusters doesn't match exactly the n from the silhuette test! 
This is because the latter depends on another algorithm. Typically they provide the same cluster solutions. However, in this case one sequence was assgigned differently. The first plot from line 478 shows that exactly one sequence is located in the overlapping area
between both clusters. Because of that, some algorithms assign it to the fist, other to the second cluster. However, changing its cluster - or deleting that observation - would not lead to other results. 

### Comparing the clusters

State-distribution plot for both clusters:
```r
seqdplot (couple.seq, group = cluster2fac)
```

Correlation between the cluster membership and men’s self-assessed dyadic coping ability 
```r
clust2.dummy<-as.numeric(cluster2fac) 
clust2.dummy[clust2.dummy==1]<-0      
clust2.dummy[clust2.dummy==2]<-1
cor.test(mydata$EDCm,clust2.dummy) 
```

LogSeq from the DySeq-Package provides an optional argument "subgroups" that allows to compare transition tables for two groups. 
```r
LogSeq(my.trans, delta=0.5, subgroups=cluster2)     # Comparim aggregated logit models between clusters with DC as DV 
LogSeq(my.trans.SC, delta=0.5, subgroups=cluster2)  # Comparim aggregated logit models between clusters with SC as DV 
```


## Additional functions
  * number of expected cell problems
  * find number of needed time intervals
  
### number of expected cell problems:
1. define a matrix containing the expected transition rates
2. run EstFreq
3. print the results

```r
# First step:
my.trans.table<-matrix(c(0.57, 0.13, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05),4,2)
# Second step:
my.cellproblems<-EstFreq(my.trans.table, t=100, min.cell=5, k=20000)
# Third step
my.cellproblems
``` 

### find number of needed time intervals:
```r
my.EstTime.plot<-EstTime(my.trans.table,   # contains expected transition probabilities (from the last section)
                         t=50:100,         # limits the range of time intervals, which are testet.  
                         k=5000)           # precission of the simulation
                                           # Warning: The function is very time consuming!

# printing will result in a plot of time point vs. expected number of low and zero frequencies
my.EstTime.plot

# if legend position should be change:
attr(my.EstTime.plot, "pos")<-"topright" # alternatives are: bottomleft (default)
my.EstTime.plot
```
