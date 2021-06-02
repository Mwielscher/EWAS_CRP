#!/usr/bin/perl


use strict;
use warnings;
use Data::Dumper;

my $dir='/project/eph/matthias/CRP/annot_hits_2/data/GWAS_cat/';
my @traits=("C-reactive_protein","3_hydroxy_1_methylpropylmercapturic_acid_levels_in_smokers","3_hydroxypropylmercapturic_acid_levels_in_smokers","Acute_lymphoblastic_leukemia__childhood_","Adiponectin_levels","Advanced_age_related_macular_degeneration","Age_related_macular_degeneration","Alzheimer_disease_and_age_of_onset","Alzheimer_s_disease","Alzheimer_s_disease__late_onset_","Amyotrophic_lateral_sclerosis","Amyotrophic_lateral_sclerosis__sporadic_","Asthma","Asthma__childhood_onset_","Atopic_dermatitis","Attention_deficit_hyperactivity_disorder","Autism_spectrum_disorder_attention_deficit_hyperactivity_disorder_bipolar_disorder_major_depressive_disorder_and_schizophrenia__combined_","Basal_cell_carcinoma","Bilirubin_levels","Bipolar_disorder","Bipolar_disorder__body_mass_index_interaction_","Bipolar_disorder_and_schizophrenia","Blood_metabolite_levels","Blood_metabolite_ratios","Blood_pressure","Body_fat_percentage","Body_mass_index","Bone_mineral_density","Breast_cancer","Breast_size","Celiac_disease","Cerebrospinal_fluid_biomarker_levels","Cholesterol_total","Chronic_lymphocytic_leukemia","Chronic_obstructive_pulmonary_disease","Cognitive_decline_rate_in_late_mild_cognitive_impairment","Cognitive_function","Cognitive_performance","Colorectal_cancer","Coronary_artery_calcification","Coronary_artery_disease","Coronary_heart_disease","Crohn_s_disease","Daytime_sleep_phenotypes","Dental_caries","Diabetic_kidney_disease","Diastolic_blood_pressure","Diisocyanate_induced_asthma","Educational_attainment","Educational_attainment__years_of_education_","Emphysema_imaging_phenotypes","Fibrinogen_levels","Glomerular_filtration_rate","Glomerular_filtration_rate__creatinine_","Glucose_homeostasis_traits","Gout","Gut_microbiota__bacterial_taxa_","Gut_microbiota__functional_units_","HDL_cholesterol","Height","Hip_circumference","Hip_circumference_adjusted_for_BMI","HIV_1_control","IgA_nephropathy","IgG_glycosylation","Inflammatory_bowel_disease","Inflammatory_skin_disease","Intraocular_pressure","Iron_status_biomarkers","Late_onset_Alzheimer_s_disease","LDL_cholesterol","Lipid_metabolism_phenotypes","Longevity","Lung_function__FEV1_FVC_","Major_depressive_disorder","Male_pattern_baldness","Mean_corpuscular_hemoglobin","Mean_corpuscular_volume","Mean_platelet_volume","Melanoma","Menarche__age_at_onset_","Metabolite_levels","Migraine","Morning_vs._evening_chronotype","Multiple_sclerosis","Myocardial_infarction","Myopia","Myopia__pathological_","Neuroticism","Night_sleep_phenotypes","Nonsyndromic_cleft_lip_with_cleft_palate","Obesity","Obesity_related_traits","Optic_cup_area","Optic_disc_area","Orofacial_clefts","Pancreatic_cancer","Parental_extreme_longevity__95_years_and_older_","Parental_longevity__combined_parental_age_at_death_","Parkinson_s_disease","Platelet_count","Polychlorinated_biphenyl_levels","Post_bronchodilator_FEV1","Post_bronchodilator_FEV1_FVC_ratio","Post_bronchodilator_FEV1_FVC_ratio_in_COPD","Post_bronchodilator_FEV1_in_COPD","PR_interval_in_Tripanosoma_cruzi_seropositivity","Primary_biliary_cirrhosis","Prostate_cancer","Psoriasis","Pulmonary_function","Pulmonary_function__smoking_interaction_","QRS_duration","QT_interval","Red_blood_cell_traits","Response_to_amphetamines","Response_to_antipsychotic_treatment","Response_to_bronchodilator_in_chronic_obstructive_pulmonary_disease__change_in_FEV1_","Response_to_paliperidone_in_schizophrenia__negative_Marder_score_","Response_to_paliperidone_in_schizophrenia__PANSS_score_","Response_to_paliperidone_in_schizophrenia__positive_Marder_score_","Resting_heart_rate","Rheumatoid_arthritis","Schizophrenia","Subjective_well_being","Sudden_cardiac_arrest","Systemic_lupus_erythematosus","Systolic_blood_pressure","Thiazide_induced_adverse_metabolic_effects_in_hypertensive_patients","Thyroid_hormone_levels","Trans_fatty_acid_levels","Triglycerides","Type_1_diabetes","Type_2_diabetes","Ulcerative_colitis","Urate_levels","Venous_thromboembolism_adjusted_for_sickle_cell_variant_rs77121243_T","Vertical_cup_disc_ratio","Visceral_adipose_tissue_adjusted_for_BMI","Waist_circumference","Waist_circumference_adjusted_for_body_mass_index","Waist_hip_ratio","Waist_to_hip_ratio_adjusted_for_body_mass_index");

#my $trait='Gut_microbiota__functional_units_';



sub make_files {
	#open(IN,$dir."HG19_gwas_cat_min50_SNPs_perTRAIT.txt") or die "cannot open GWAS cat";
	#<IN>;
	
	foreach my$trait(@traits) {
	print "now doing $trait !! \n";
	 open(OUT,">".$dir."input_data/noPRUNE_".$trait) or die "cannot open OUT";
 	open(IN,$dir."HG19_gwas_cat_min50_SNPs_perTRAIT.txt") or die "cannot open GWAS cat";
        <IN>;
	while (<IN>) {
		chomp;
		#my @line=split "\t",$_;
		#print "thats $line[8]!!!! \n"; 
		if (my ($grpline) = grep/$trait/,$_) {
		my @line=split "\t",$grpline;
		my $start = $line[1] - 500000;
		my $end = $line[1] + 500000;
		print "chr$line[0] $line[1] $line[2] $start $end $line[8] \n";	
		print OUT "chr$line[0] $line[1] $line[2] $start $end $line[8] \n";
	#	my$newline=join " ","chr$line[0]", $line[1], $line[2], $start, $end, $line[8];
	#	push @allsnp,$newline;
		}
		}
	close(OUT);

	}
	
	close(IN);
	
}
make_files();

