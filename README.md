#Analyzing dyadic sequence data - research questions and implied statistical models 

Please note, a full R-script of this vignette is provided under the name DySeq_script_v0.1.R in this repository!

This Vignette provides a hands-on-tutorial for all analyses covered in 

"Analyzing dyadic sequence data - research questions and implied statistical models"
by  --- blinded for reviewing ---
published in --- still in review ---

Please make sure to install all required packages,
including the "Dyseq", which provides the sample data!

Content                  

1. Prerequisite Steps                
  * packages from CRAN                   
  * package from Github                  

2. Example Data                      
  * loading the data                     
  * details on data                      

3. Graphical Analysis                 
  * state-distribution-plot              
  * entropy-plot                         
  * Number of transitions                

4. Research question 1  
  * Pearson Correlation 
 
5. Research question 2:                   
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

6. Research question 3:                 
  * estimating hidden Markov model

7. Research question 4:
  * estimating mixture Markov model
  * sequence clustering
    * OM-distances
    * clustering
    * interpret clusters

---

## Prerequisite steps

make sure the following packages are installed:

```# install.packages("TraMineR")      # for graphical analysis and research question 4
# install.packages("RColorBrewer")  # for grey-shaded graphics
# install.packages("gmodels")       # must be installed!
# install.packages("MASS")          # must be installed!
# install.packages("survival")      # must be installed for research question 3!
# install.packages("fpc")           # must be installed for research question 4!
# install.packages("cluster")       # must be installed for research question 4!
# install.packages("devtools")      # must be installed for installing packages from github
# install.packages("lme4")          # must be installed for the multi-level APIM
# install.packages("lmerTest")      # must be installed for the multi-level APIM```

