
#!/bin/bash

module load metal/2016-02-02

indir1='/rds/general/user/mwielsch/home/WORK/CRP/meta_score/data/'
outdir='/rds/general/user/mwielsch/home/WORK/CRP/meta_score/'


cat <<EOF > metal_script
SCHEME STDERR
SEPARATOR TAB 
MARKER ID
EFFECT EST
STDERR SE
PVALUELABEL P_VAL
COLUMNCOUNTING LENIENT
PROCESS ${indir1}META_framingham-pheno-risk-score.txt 
PROCESS ${indir1}META_AIRWAVE_OR_20200626_CRP_risk_scorev2.txt
PROCESS ${indir1}META_SHIP-Trend_AT_200615_CRP_risk_score.txt
PROCESS ${indir1}META_BreastCancer_KNC_June30_wlogtf_CRP_risk_score_OK.txt
PROCESS ${indir1}META_KORA_F4_BK_20200702_CRP_risk_score.txt
PROCESS ${indir1}META_YFS_PPM_22June2020_CRP_risk_score_OK.txt




OUTFILE ${outdir}CRP_causes_CpG_riskSCPRE_FINAL_ .tbl
ANALYZE HETEROGENEITY
QUIT
EOF
metal metal_script









