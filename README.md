## Table of contents
1. [About this Repository](#About-this-Repository)
2. [Meta analysis](#Meta-analysis)
3. [Mediation analysis](#Mediation-analysis)
4. [Mendelian Randomization](#Mendelian-Randomization)
4. [Overrepresentation analysis](#Overrepresentation-analysis)
5. [Associaiton to clinical phenotypes](#Associaiton-to-clinical-phenotypes)  

## About this Repository
This repository accompanies the paper __"DNA methylation signature of chronic low-grade inflammation and its role in cardio-respiratory diseases"__ 
<br/><br/>
<p align="center">
<img src="/img/FIGURE_1_ok.jpg" alt="Overview Figure" width="600"/>
<br/><br/>

> **_Abstract:_**  We performed a trans-ethnic Epigenome Wide Association study on 22,774 individuals to describe the DNA methylation signature of chronic low-grade inflammation as measured by C-Reactive protein (CRP). We found 1,511 independent loci associated with CRP. These CpG sites showed correlation structure across chromosomes, and were primarily situated in euchromatin, depleted in CpG islands and enriched in transcription factor binding sites and genomic enhancer regions. Mendelian randomisation analysis suggests altered CpG methylation is a consequence of increased CRP levels. Mediation analysis revealed obesity and smoking as important underlying factors for changed CpG methylation. Finally, we found that an activated CpG signature significantly increases the risk for cardiometabolic diseases and COPD. 
<p>
<br/>

The scripts to perform those analysis are depoisited in this Repo. Those are costum scripts, which were shared with collaborators for the paper or run on my local destop computer. Most Shell script are designed to run on a [SGE cluster](http://gridscheduler.sourceforge.net/htmlman/manuals.html). There are short example files within the folders to test scripts. 
  
## Meta analysis

## Mediation analysis
We performed a mediation analysis to better understand our CRP associated DNA methylation markers. We followed the analysis framework as suggested by [Baron Kenny](https://en.wikipedia.org/wiki/Mediation_(statistics) For example for BMI-> CRP -> DNAmeth mediation analysis, we assed the results of 4 regressions: a-path: CRP ~ BMI + covariates, b-path: DNAmeth ~ CRP + BMI + covariates, c-path: DNAmeth ~BMI + covariates (total effect) and c’-path:DNAmeth ~BMI + CRP + covariates. From this we can asses the indirect effect which is a*b. Values of this multiplication were compared to c - c', which should give the same result. Next we performed an Aroian Sobel test to assess the significance of the indirect effect. Briefly, we extracted coefficients and standard errors from a path regressions and b-path regressions and calculated a Z score like so: Zscore=(a*b)/sqrt((b^2*SEa^2)+(a^2*SEb^2))+(SEa*SEb). This was performed for every CpG separately. Z score was calculated only if the indirect effect was negative and we saw a significant association (P< 0.05) between BMI and the tested CpG. For this analysis we excluded all samples with BMI values outside of mean(BMI) +/- 4*SD(BMI).  
  
  
## Mendelian Randomization

## Overrepresentation analysis

## Associaiton to clinical phenotypes 


