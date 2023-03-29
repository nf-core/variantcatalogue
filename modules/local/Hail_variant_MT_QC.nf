process Hail_variant_MT_QC {
	label 'process_medium'

    conda "bioconda::hail=0.2.58"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hail:0.2.61--py37h9a982cc_1':
        'quay.io/biocontainers/hail:0.2.61--py37h9a982cc_1' }"

	input :
	path(MT_Step1_input_tsv)
	path(MT_Step2_participant_data)
	path(MT_participants_to_subset)
	path(MT_Step3_participant_data)
	file reference
        file reference_index
        file pon_prediction_table
	file artifact_prone_sites_bed
	file mitotip_predictions_table

	output :
	path('MT_post_hail_combined_sites_only.vcf.bgz')    , emit : vcf
	path('MT_post_hail_combined_sites_only.vcf.bgz.tbi'), emit : vcf_index
	path('MT_filtered_frequ_only.vcf')                  , emit : Hail_MT_frequ_only
	path('sample_annotations_MT.txt')                   , emit : sample_annotations
	path('MT_stats_pass.txt')                           , emit : stats_pass
	path('MT_stats.txt')                                , emit : stats

	script:
	"""
	mkdir -p $params.tmp_dir
	pip install gnomad
        python ${projectDir}/modules/local/Hail_variant_MT_QC.py $MT_Step1_input_tsv $MT_Step2_participant_data $MT_participants_to_subset $MT_Step3_participant_data $pon_prediction_table $artifact_prone_sites_bed $reference $refrence_index $mitotip_predictions_table $params.tmp_dir
	"""
}

