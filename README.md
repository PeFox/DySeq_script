
##Analyzing dyadic sequence data - research questions and implied statistical models 

This R-script provides a hands-on-tutorial for all
analyses covered in 
"Analyzing dyadic sequence data - research questions and implied statistical models"
by  --- blinded for reviewing ---
 
Please make sure to install all required packages,
including the "Dyseq" which provides the sample data!

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


