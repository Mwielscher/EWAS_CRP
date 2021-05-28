## Table of contents
1. [About this Repository](#About-this-Repository)
2. [Meta analysis](#Meta-analysis)
3. [Mediation analysis](#Mediation-analysis)
4. [Mendelian Randomization](#Mendelian-Randomization)
4. [Overrepresentation analysis](#Overrepresentation-analysis)
5. [Association to clinical phenotypes](#Association-to-clinical-phenotypes)  

## About this Repository
This repository accompanies the paper __"DNA methylation signature of chronic low-grade inflammation and its role in cardio-respiratory diseases"__ 
<br/><br/>
#<p align="center">
#<img src="/img/FIGURE_1_ok.jpg" alt="Overview Figure" width="600"/>
#<br/><br/>

> **_Abstract:_**  We performed a trans-ethnic Epigenome Wide Association study on 22,774 individuals to describe the DNA methylation signature of chronic low-grade inflammation as measured by C-Reactive protein (CRP). We found 1,511 independent loci associated with CRP. These CpG sites showed correlation structure across chromosomes, and were primarily situated in euchromatin, depleted in CpG islands and enriched in transcription factor binding sites and genomic enhancer regions. Mendelian randomisation analysis suggests altered CpG methylation is a consequence of increased CRP levels. Mediation analysis revealed obesity and smoking as important underlying factors for changed CpG methylation. Finally, we found that an activated CpG signature significantly increases the risk for cardiometabolic diseases and COPD. 
<p>
<br/>

The scripts to perform those analysis are depoisited in this Repo. Those are costum scripts, which were shared with collaborators for the paper or run on my local destop computer. Most Shell script are designed to run on a [SGE cluster](http://gridscheduler.sourceforge.net/htmlman/manuals.html). There are short example files within the folders to test scripts. 
  
## Meta analysis  
Meta-analysis were performed to identify DNA methylation markers associated to blood CRP levels. The analysis was restricted to autosomal markers on the Infinium Human Methylation 450K BeadChip. Further we excluded probes if they: (i) had a [SNP in last 10bp of the probe sequence](https://bioconductor.org/packages/release/data/annotation/html/IlluminaHumanMethylation450kmanifest.html) (ii) were flagged as [cross-reactive probes](https://www.tandfonline.com/doi/full/10.4161/epi.23470). This left 405,019 CpG sites for analysis in our meta-analysis. Effect sizes and standard errors of the 30 studies were combined using inverse variance weighting method as implemented in METAL software. We applied genomic control in and [genomic control](https://www.cell.com/ajhg/fulltext/S0002-9297(07)62366-0) out to control for population stratification and other unmeasured confounding factors in the analysis. Genomic control out was achieved by correcting P values for the inflation factor lambda. Genomic control procedure was also applied in all sensitivity analyses and the trans-ethnic replication analysis

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
  
Now that the dataset is ready we can **test hypothesis 1**, we needed to find instruments for DNA methylation. For this we regressed DNA methylation value of every CpG against all SNPs present in “cis” of the concordant CpG. The cis-region of every sentinel CpG was defined, as it’s chromosomal position +/- 500kb. 
Regressions were performed using [rvtests software](https://github.com/zhanxw/rvtests)  
>* This [script](Mendelian_Randomization/assoc_scripts/1_assoc_CpG_cause_CRP.sh) shows how we were looking for instruments for CpG methylation sites. An example file defining the genomic regions for SNP selection can be found [here](Mendelian_Randomization/example_files/)  
>* We produced estimates for the direct effect of SNP on the CRP with this [script](Mendelian_Randomization/assoc_scripts/1a_assoc_CpG_cause_CRP_DIRECT_EFF.sh). This was done because we excluded SNPs with a direct effect on the CRP. That is if SNP~CRP + lnCRP gives a significant association between SNP and lnCRP if we add CRP as covariate to the model. 
>* regression outlined in the two scripts above were performed in each chohort seperately and then meta analysed. Meta analysis was done very similar in all analysis this [METAL script](Risk_Score/metal_meta.sh) can be changed to do any metanalysis job in this project.  
  
Finally we used the ratio method to determine significance in the Mendelian Randomization analysis as implemented in the [R package MendelianRandomization](https://academic.oup.com/ije/article/46/6/1734/3112150)and performed a triangulation analysis. This technically simple and was done with this [script](Mendelian_Randomization/MR_analysis/3_CpG_cause_final_version.R) the concept and motivation to do that is explained at the end of this section.    

To **test hypothesis 2**, we used the latest published [GWAS on CRP](https://www.cell.com/ajhg/fulltext/S0002-9297(18)30320-3) to define a set of instruments for CRP. In contrast to testing of hypothesis 1 where we rely on large scale GWAS summary statistics for the association to our outcome (CRP). As there is currently no large scale GWAS summary statistics for Illumina 450k CpG as an outcome available, we created a CRP polygenic risk score starting 52 SNPs. To generate a beta weighted risk score we used PLINK version 1.9. Next, we regressed the risk score against every sentinel CpG under an additive model. CpG ~CRPPRS + age + sex + blood cell estimates + genetic PC [1..10].

>* [script](/Mendelian_Randomization/assoc_scripts/3_make_plink_score.sh) to create a plink risk score. Example files SNP selection intermediate files etc. can be found [here](Mendelian_Randomization/example_files/)  
>* Again we produced estimates for the [direct effect](/Mendelian_Randomization/assoc_scripts/2_assoc_CRP_cause_CpG_DIRECT_EFF.sh) of SNP on the DNA methylation. This was done because we excluded SNPs with a direct effect on the DNA methylation and rerun [plink score script](/Mendelian_Randomization/assoc_scripts/3_make_plink_score.sh) 
>* We then regressed the CRP polygenic risk score agains DNA methylation using this [script](/Mendelian_Randomization/assoc_scripts/4_run_socore_assoc.sh). 
>* we compared regression coefficients across participating cohorts using this [script}(/Mendelian_Randomization/MR_analysis/1_QC_input_data.R).  Then results were meta-analysed using a [METAL script](Risk_Score/metal_meta.sh) similar to this one. 
>* The above procedure provides us with all input data for the Mendelian Randomisation and triangulation analysis. This was done with this [R-script](/Mendelian_Randomization/MR_analysis/2_MR_triangulation_CRP_causesDNAmeth.R)

### Triangulation Analysis  
We wanted to understand if there is general trend for all CpGs to be cause or consequence of blood CRP levels. For this we applied a triangulation approach: The effect of the instrument on the outcome (observed effect) should equal the product of the effect of instrument on the exposure and the effect of the exposure on the outcome. For hypotheses 1 the observed effect is CRP ~SNP and the predicted effect is the product of the effects of DNA methylation on the SNPs and CRP on DNA methylation. For hypotheses 2 the observed effect is DNA methylation ~GRS (polygenic CRP risk score) and the predicted effect is the product of effect of CRP ~GRS and CRP ~DNA methylation. If the observed effect is mediated through the predicted effect those will correlate  
  
  

## Overrepresentation analysis
We wanted to determine the overlaps between our 1511 CRP-associated loci and publicly available lists of genomic features such as Histone Modifications, Transcription factor binding sites and chromatin states. For this we performed permutation tests. First, we assessed the number of overlaps between our 1511 CRP-associated markers and a genomic feature of interest. Then we sampled 10000 sets of 1511 markers. For every set of 1511 marker we recorded the number of overlaps to the genomic feature of interest. Those 10000 overlaps created our H0 distribution. We calculated an empirical P values separately based on either the number of entries with higher or the number of entries with lower overlap in our H0 distribution compared to our observed number of overlaps recorded within the 1511 CRP-associated loci. Additionally, we calculated a Fisher P value based on the mean of our H0 distribution and the observed overlap. We collected the mean and standard errors of for each CpG from most studies in this EWAS. 90% of standard errors were between 0.007 and 0.075. SD values in this range were binned the SD data into 0.005 intervals. For generation of 10000 random sets of markers we used the equal numbers of markers within each standard error bin as observed in the CRP-associated loci. 
We calculated overlaps to the CpG annotation as given in the Illumina Manifest file. We retrieved DNaseI-accessible sites from [Encode project](http://www.uwencode.org/proj/hotspot), specifically gapped peaks from Release 9 called with MACSv2.0.10. For histone marks (H3K9 and H3K27) were retrieved gapped peak data from Roadmap project for a collection of 127 tissues and cell lines as well as probabilities for [Roadmap 15 state chromatin model](https://www.nature.com/articles/nature14248) for selected cell lines. Encode Transcription factor binding sites were retrieved from [UCSC browser](http://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=wgEncodeRegTfbsClusteredV3). Next, we mapped Illuminas 450k probes to the GWAS catalogue hits (downloaded 20170626) of traits having 50 or more hits recorded in the catalogue. We allowed a window of 1MB to match the 450k probes to each GWAS catalogue entry and removed one or more hits in case they were overlapping. Finally, we retrieved chromosome conformation capture data (HiC) as described by [Bing Ren et al](https://www.nature.com/articles/nature11082) 
  
  
  
## Association to clinical phenotypes 

We calculated a beta weighted risk score using the coefficients from the trans-ethnic discovery analysis (similar to a polygenic risk score in GWAS).
For every participant in each study:  
CpGriskSCORE =  sum( CpGmethylation x  CoefficienttransEthnicDiscovery)  
Then we performed logistic regression model for the CpGriskSCORE against each outcome. Depending on availability of phenotypes across cohorts we combined the effect sizes and standard errors using invers variance weighted approach as implemented in METAL software. We present meta analyzed effect sizes (logODDs) and P-values from logistic regression. Odd ratios were transformed to produce adjusted relative risk estimates:  
RR=odds ratio/1- (lifetime risk) + (life time risk x odds ratio)  
We used lifetime risk estimates from current literature: [COPD 11.45%](https://erj.ersjournals.com/content/42/4/964.long), [T2D 39.9%](https://www.thelancet.com/journals/landia/article/PIIS2213-8587(14)70161-5/fulltext), [MI 24.8%](https://www.nejm.org/doi/full/10.1056/NEJMoa1804492) [CAD 40.15%](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(05)75029-7/fulltext) and [Hypertension 81%](https://jamanetwork.com/journals/jamacardiology/fullarticle/2728380).

>* [script](/Risk_Score/CRP_CpG_risk_score.R) to calculate risk score and run regressions against clinical relevant phenoptype. 
>* this is a good onlie [tool](https://clincalc.com/Stats/ConvertOR.aspx) to calculate the adjusted relative risk
