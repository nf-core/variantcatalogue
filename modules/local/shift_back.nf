process shift_back {
	tag "${meta.id}"
	label 'process_low'
	
        container = 'https://depot.galaxyproject.org/singularity/r-stringr%3A1.1.0--r3.3.1_0'

	input :
	tuple val(meta), path(MT_shifted_CollectMetrics), path(MT_CollectMetrics)
	
	output :
	tuple val(meta), path ('*per_base_coverage.tsv')  , emit : per_base_coverage
	tuple val(meta), path ('*_MT_Step1_input_tsv.tsv'), emit : Sample_MT_Step1_input_tsv

	script:
	"""
	Rscript ${projectDir}/modules/local/shift_back.R $MT_shifted_CollectMetrics ${meta.id}_sorted_chrM_Homo_sapiens_assembly38_collect_wgs_metrics_non_control_region.chrM.interval_list.tsv
	mv per_base_coverage.tsv ${meta.id}_per_base_coverage.tsv	
	"""
}
