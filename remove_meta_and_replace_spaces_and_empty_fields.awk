# remove_meta_and_replace_spaces_and_empty_fields.awk

# takes as input a clinical file from crct of LMS samples and
# - removes the metastasis samples (5 of them, present in column 112 with their code name with T inside)
# - replaces spaces with underscores
# - replaces empty fields with NA

# example
# basedir=/work/project/bridge/workspace/l3_loanemadrid_lms
# cd $basedir/results
# pgm=$basedir/scripts/remove_meta_and_replace_spaces_and_empty_fields.awk
# awk -f $pgm Annot.ICGC.VF.tsv > Annot.ICGC.VF.wometa.ok.tsv

# input file = Annot.ICGC.VF.tsv
# Patient_Id	File_number	Centre__ConticaBase_	Birth_date	Sex	Creation_date_patient	Significant_previous_history	Significant_previous_history_det	Other_cancer	Date_of_last_contact	Lost_of_follow_up	Vital_status	Cause_of_death	Date_of_death	Primary_tumour_Id	Site_of_tumour_category	Site_of_tumour_sub_category	Site_of_tumour	Side__RESOS_	Size_of_tumour__mm_	Depth_of_tumour	Histotype_category	Histotype_sub_category	Histotype	Gradable_histotype	Date_of_original_diagnosis	Age_at_diagnosis	Diagnosis_changed	Grade_of_tumour	Radiation_area	Lymphoedema	Loco_regional_extension	Multifocal_tumour	N	M	StatusAtReferral	TTTBeforeReferral	Performance_status	Cancer_status	ImagingBeforeTTT	Surgery_of_tumour	Date_of_surgery	Time_of_surgery	Type_of_surgery_of_tumour	Surgeon	Re_excision	Re_excision_surgeon	Tumour_spillage__CTCB_	Margin_R	Isolated_limb_perfusion	DILP	SystemicTTT	ChemotherapyOfTumour	DSytTTTStart	DSystTTTEnd	IndicSystTTT	Radiotherapy_of_tumour	DRTEStart	DRTEEnd	Dose_of_radiotherapy	Hormonotherapy	DHormonoTStart	AINS	DAINSStart	OtherTTT	Clinical_trial	Locoregional_control	Complete_remission_after_treatme	Local_recurence_Id	Rank_of_local_recurrence	Date_of_local_recurrence	Size_of_tumour__mm_1	Depth_of_tumour1	Multifocal_tumour1	Radiation_area1	Performance_status1	Surgery_of_tumour1	Date_of_surgery1	Time_of_surgery1	Type_of_surgery_of_tumour1	SurgeonLR	Re_excision1	Tumour_spillage__CTCB_1	Margin_R1	SystemicTTT1	DSytTTTStart1	DSystTTTEnd1	IndicSystTTT1	Isolated_limb_perfusion1	DILP1	Radiotherapy_of_tumour1	DRTEStart1	DRTEEnd1	Dose_of_radiotherapy1	Hormonotherapy1	DHormonoTStart1	AINS1	DAINSStart1	Clinical_trial1	OtherTTT1Locoregional_control1	Complete_remission_after_treatm0	Metastasis_Id	Location_of_first_metastasis	Date_of_first_metastasis	Number_of_metastatic_sites	Location_of_further_metastasis	TTTM	Locoregional_treatment_of_metast	Sample_Id	Patient	Code_fragment	Sample_RNA	Sample_reference	Date_of_sampling	Treatment_before_sampling	Type_of_sampling	Immunohistochemistry	FISH	Result_FISH	Molecular_Biology	Result__Molecular_Biology_	Frozen_tissue	Paraffin_available	Virtual_slide	Metastasis_sample_location	sample_evnt	Chemotherapy_line_Id	Rank_of_Chemotherapy_line	Date_of_first_cycle	Performance_status_chemo	Drugs	Clinical_trial_chemo	Clinical_trial_name	Progression_chemo	Date_progression_chemo	Best_RECIST_Response_chemo	chemo_evnt
# 108783	ANGERS-000414308	Angers ICO	28/11/1960	Female	09/02/2015	No		No	01/02/2017	No	Dead	Dead of this cancer	01/02/2017	108784	Viscera	Gyneacological area	Uterus	Not applicable	27	Deep	Sarcoma	Leiomyosarcoma	Leiomyosarcoma	true	21/12/2010	50	No	NA	No	No		No	No	No	First event	No		Evidence of this cancer	Ultrasound	Yes	27/01/2011	First	Wide resection	Outside network	Yes		No	R0	No		No	No	Second			45	No		No		No	No	Yes	Yes								568758	Bone, Liver, Lung, Other	01/10/2016	4				108785	LMS28	LMS28T	LMS28R	ANGERS-000414308-11H01126	27/01/2011	No	Tumour resection	No	No		No		Yes	Yes	ICGC/ANGERS-000414308-11H01126_HES.ndpi		PRIMARY											
# 1 (49 fields)
# 1 (62 fields)
# ... 
# 1 (140 fields)
# 1 (154 fields)  *** 126 rows

# output file = Annot.ICGC.VF.wometa.ok.tsv
# Patient_Id	File_number	Centre__ConticaBase_	Birth_date	Sex	Creation_date_patient	Significant_previous_history	Significant_previous_history_det	Other_cancer	Date_of_last_contact	Lost_of_follow_up	Vital_status	Cause_of_death	Date_of_death	Primary_tumour_Id	Site_of_tumour_category	Site_of_tumour_sub_category	Site_of_tumour	Side__RESOS_	Size_of_tumour__mm_	Depth_of_tumour	Histotype_category	Histotype_sub_category	Histotype	Gradable_histotype	Date_of_original_diagnosis	Age_at_diagnosis	Diagnosis_changed	Grade_of_tumour	Radiation_area	Lymphoedema	Loco_regional_extension	Multifocal_tumour	N	M	StatusAtReferral	TTTBeforeReferral	Performance_status	Cancer_status	ImagingBeforeTTT	Surgery_of_tumour	Date_of_surgery	Time_of_surgery	Type_of_surgery_of_tumour	Surgeon	Re_excision	Re_excision_surgeon	Tumour_spillage__CTCB_	Margin_R	Isolated_limb_perfusion	DILP	SystemicTTT	ChemotherapyOfTumour	DSytTTTStart	DSystTTTEnd	IndicSystTTT	Radiotherapy_of_tumour	DRTEStart	DRTEEnd	Dose_of_radiotherapy	Hormonotherapy	DHormonoTStart	AINS	DAINSStart	OtherTTT	Clinical_trial	Locoregional_control	Complete_remission_after_treatme	Local_recurence_Id	Rank_of_local_recurrence	Date_of_local_recurrence	Size_of_tumour__mm_1	Depth_of_tumour1	Multifocal_tumour1	Radiation_area1	Performance_status1	Surgery_of_tumour1	Date_of_surgery1	Time_of_surgery1	Type_of_surgery_of_tumour1	SurgeonLR	Re_excision1	Tumour_spillage__CTCB_1	Margin_R1	SystemicTTT1	DSytTTTStart1	DSystTTTEnd1	IndicSystTTT1	Isolated_limb_perfusion1	DILP1	Radiotherapy_of_tumour1	DRTEStart1	DRTEEnd1	Dose_of_radiotherapy1	Hormonotherapy1	DHormonoTStart1	AINS1	DAINSStart1	Clinical_trial1	OtherTTT1Locoregional_control1	Complete_remission_after_treatm0	Metastasis_Id	Location_of_first_metastasis	Date_of_first_metastasis	Number_of_metastatic_sites	Location_of_further_metastasis	TTTM	Locoregional_treatment_of_metast	Sample_Id	Patient	Code_fragment	Sample_RNA	Sample_reference	Date_of_sampling	Treatment_before_sampling	Type_of_sampling	Immunohistochemistry	FISH	Result_FISH	Molecular_Biology	Result__Molecular_Biology_	Frozen_tissue	Paraffin_available	Virtual_slide	Metastasis_sample_location	sample_evnt	Chemotherapy_line_Id	Rank_of_Chemotherapy_line	Date_of_first_cycle	Performance_status_chemo	Drugs	Clinical_trial_chemo	Clinical_trial_name	Progression_chemo	Date_progression_chemo	Best_RECIST_Response_chemo	chemo_evnt
# 108783	ANGERS-000414308	Angers_ICO	28/11/1960	Female	09/02/2015	No	NA	No	01/02/2017	No	Dead	Dead_of_this_cancer	01/02/2017	108784	Viscera	Gyneacological_area	Uterus	Not_applicable	27	Deep	Sarcoma	Leiomyosarcoma	Leiomyosarcoma	true	21/12/2010	50	No	NA	No	No	NA	No	No	No	First_event	No	NA	Evidence_of_this_cancer	Ultrasound	Yes	27/01/2011	First	Wide_resection	Outside_network	Yes	NA	No	R0	No	NA	No	No	NA	NA	NA	Second	NA	NA	45	No	NA	No	NA	No	No	Yes	Yes	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	568758	Bone,_Liver,_Lung,_Other	01/10/2016	4	NA	NA	NA	108785	LMS28	LMS28T	LMS28R	ANGERS-000414308-11H01126	27/01/2011	No	Tumour_resection	No	No	NA	No	NA	Yes	Yes	ICGC/ANGERS-000414308-11H01126_HES.ndpi	NA	PRIMARY	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
# 121 (138 fields)  *** 121 rows which correspond to the 5 samples removed (LMS75T, 82T, 123MT, 129T, 130T), and ok format now
#                       however 119 samples only since LMS85T is present twice (tumor and recidive)

BEGIN{
    OFS="\t"; 
    ko["LMS82T"]=1; 
    ko["LMS75T"]=1; 
    ko["LMS123MT"]=1; 
    ko["LMS130T"]=1; 
    ko["LMS129T"]=1;
} 

NR==1{
    gsub(/\ /,"_",$0); 
    print;
} 

NR>=2{
    gsub(/\ /,"_",$0); 
    n=split($0,a,"\t"); 
    if(ko[a[112]]!=1)
    {
        s=""; 
        for(k=1; k<n; k++)
        {
            s=(s)(a[k]=="" ? "NA" : a[k])("\t");
        } 
        print (s)(a[k]=="" ? "NA" : a[k]);
    }
}
