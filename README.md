#Analyzing dyadic sequence data - research questions and implied statistical models 

Please note, a full R-script of this vignette is provided under the name DySeq_script_v0.1.R in this repository!

This Vignette provides a hands-on-tutorial for all analyses covered in 

"Analyzing dyadic sequence data - research questions and implied statistical models"
by  --- blinded for reviewing ---
published in --- still in review ---

Please make sure to install all required packages,
including the "Dyseq", which provides the sample data!

Content                  

A. Prerequisite Steps                
  * packages from CRAN                   
  * package from Github                  

B. Example Data                      
  * loading the data                     
  * details on data                      

C. Graphical Analysis                 
  * state-distribution-plot              
  * entropy-plot                         
  * Number of transitions                

D. Research question 1  
  * Pearson Correlation 
 
E. Research question 2:                   
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

F. Research question 3:                 
  * estimating hidden Markov model

G. Research question 4:
  * estimating mixture Markov model
  * sequence clustering
    * OM-distances
    * clustering
    * interpret clusters

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

Plotting the results will produce an interaction-plot. Mapping the probabilities of showing the dependend variable against the combinations of previous behavior.  

![alt text](https://github.com/PeFox/DySeq_script/blob/master/Interac.png "Scatterplot(SC,DC)")
