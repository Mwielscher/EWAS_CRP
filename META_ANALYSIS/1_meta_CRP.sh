#!/bin/bash
#PBS -l select=1:ncpus=1:mem=20GB
#PBS -l walltime=01:20:00
#PBS -N smoke_metal

module load metal/2016-02-02

cat <<EOF > metal_script

GENOMICCONTROL ON
SCHEME STDERR
SEPARATOR TAB 
MARKER probeID
EFFECT BETA 
STDERR SE
PVALUELABEL P_VAL 
PROCESS /project/eph/matthias/CRP/data/model1/AIRWAVEsub_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/ARIC_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/ARIES_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/BHS.B_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/BHS.W_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/BIOS.PAN_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/BIOS.CODAM_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/BIOS.LLS_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/BIOS.NTR_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/CARDIOGENICS_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/CHS.B_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/CHS.W_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/EGCUT.CTG_model1.txt
#PROCESS /project/eph/matthias/CRP/data/model1/EGCUT.YaO_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/EPIC.Norfolk_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/EPICOR_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/ESTHER.1a_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/ESTHER.1b_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/FHS_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/GENOA_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/KORA_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/LBC_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/LLD_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/LOLIPOP_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/NAS_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/NFBC1966_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/NFBC1986_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/ROTTERDAM_model1.txt
#PROCESS /project/eph/matthias/CRP/data/model1/SHIP_model1.txt
#PROCESS /project/eph/matthias/CRP/data/model1/YFS_model1.txt
PROCESS /project/eph/matthias/CRP/data/model1/STS_model1.txt
OUTFILE /project/eph/matthias/CRP/analysis/final_meta/COMPARE_IVW_JAN_allANCEST_het_check_1 .tbl 
ANALYZE HETEROGENEITY
QUIT

EOF

metal metal_script

