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
We performed a mediation analysis to better understand our CRP associated DNA methylation markers. We followed the analysis framework as suggested by [Baron Kenny](https://en.wikipedia.org/wiki/Mediation_(statistics)). For example for BMI-> CRP -> DNAmeth mediation analysis, we assed the results of 4 regressions: a-path: CRP ~ BMI + covariates, b-path: DNAmeth ~ CRP + BMI + covariates, c-path: DNAmeth ~BMI + covariates (total effect) and c’-path:DNAmeth ~BMI + CRP + covariates. From this we can asses the indirect effect which is a*b. Values of this multiplication were compared to c - c', which should give the same result. Next we performed an Aroian Sobel test to assess the significance of the indirect effect. Briefly, we extracted coefficients and standard errors from a path regressions and b-path regressions and calculated a Z score like so: Zscore=(a*b)/sqrt((b^2*SEa^2)+(a^2*SEb^2))+(SEa*SEb). This was performed for every CpG separately. Z score was calculated only if the indirect effect was negative and we saw a significant association (P< 0.05) between BMI and the tested CpG. For this analysis we excluded all samples with BMI values outside of mean(BMI) +/- 4*SD(BMI).  
>* [script](Mediation_Analysis/1_resid_correlation_mediation_regression.R) to run 4 necessary regression analysis. This script also contains code to rsidualize your methylation data and a simple recipe to calculate a correlation matrix of your set of CpG markers.  
>* [script](Mediation_Analysis/2_mediation_analysis.R) to run the regressions, meta analysis of the individual regressions and Aroian Sobel test
<p>

#### Why residualize the data?  
  
We wanted to remove all unwanted variation from the DNA methylation values. To achieve this, we regressed out covariates known to influence the DNA methylation data from the quantile normalized DNA methylation beta values. The regression model was as follows: 
DNAmeth ~ age + sex + CD4T + NK + Bcell + Mono + Neu + Eos + batch +[..]  
For data integrated with [CPACOR pipeline](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4365767/) we added the first 10 principal components of the control probe PCs. 

## Mendelian Randomization
We performed a Mendelian Randomization analysis to better understand the reasons for the differential DNA methylation associated to CRP. We performed analysis to investigate 2 Hypothesis: 
>1. DNA methylation is causal for CRP changes  
>2. DNA methylation is a consequence of changed CRP levels  

We provided scripts to prepare the genotype data for the regression analysis necessary to run the Mendelian Randomization analysis:  
>* This QC [script](Mendelian_Randomization/prepare_dataset/QC_1_subset_to_meth_samples.sh) subsets your data set to samples required for Mendelian randomisation. This is necessary because in our experience the cohorts have DNA methylation data available just for a subset of samples.  
>* This QC [script](Mendelian_Randomization/prepare_dataset/QC_2_concat_and_make_plink.sh) combines your genetic data, which in our experience were usually stored one file per chromosome, to one file and creates a vcf file and a plink file set of your dataset  
>* This QC [script](Mendelian_Randomization/prepare_dataset/QC_3_make_kinship_matrix.sh) will calculate a kinship matrix for your dataset. You can use this in your regression analysis to protect your results being under the influence of cryptic relatedness and population stratification. This was a very important step for NFBC data, but might not be necessary in your cohort   
>* There are example files for every step of the way in this [folder](Mendelian_Randomization/example_files/). This should prevent you from getting stuck because some input file looks differnt  
<p>
  
Now that the dataset is ready we can test hypothesis 1, we needed to find instruments for DNA methylation. For this we regressed DNA methylation value of every CpG against all SNPs present in “cis” of the concordant CpG. The cis-region of every sentinel CpG was defined, as it’s chromosomal position +/- 500kb. 
Regressions were performed using [rvtests software](https://github.com/zhanxw/rvtests)  
>* [This](Mendelian_Randomization/assoc_scripts/1_assoc_CpG_cause_CRP.sh) is how we were looking for instruments for CpG methylation sites. An example file defining the genomic regions for SNP selection can be found [here](Mendelian_Randomization/example_files/)  
>* We also produced estimates for the [direct effect](Mendelian_Randomization/assoc_scripts/1a_assoc_CpG_cause_CRP_DIRECT_EFF.sh) of SNP on the CRP. This was done because we excluded SNPs with a direct effect on the CRP. That is if SNP~CRP + lnCRP gives a significant association between SNP and lnCRP if we add CRP as covariate to the model. 
>* regression outlined in the two scripts above were performed in each chohort seperately and then meta analysed. Meta analysis was done very similar in all analysis this [METAL script](Risk_Score/metal_meta.sh) can be changed to do any metanalysis job in this project.  
  
Finally we used the ratio method to determine significance in the Mendelian Randomization analysis as implemented in the [R package MendelianRandomization](https://academic.oup.com/ije/article/46/6/1734/3112150) where MRbeta=BETA CpG~SNP / BETA CRP~SNP and performed a triangulation analysis. This technically simple and was done with this [script](Mendelian_Randomization/MR_analysis/3_CpG_cause_final_version.R) the concept and motivation to do that is explained below.    

## Overrepresentation analysis

## Associaiton to clinical phenotypes 


