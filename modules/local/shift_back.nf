process shift_back {
	tag "${meta.id}"
	label 'process_low'
	
    conda "bioconda::r-stringr=1.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-stringr%3A1.1.0--r3.3.1_0':
        'quay.io/biocontainers/r-stringr:1.1.0--r3.3.1_0' }"

	input :
	tuple val(meta), path(MT_shifted_CollectMetrics), path(MT_CollectMetrics)
	
	output :
	tuple val(meta), path ('*per_base_coverage.tsv')  , emit : per_base_coverage
	tuple val(meta), path ('*_MT_Step1_input_tsv.tsv'), emit : Sample_MT_Step1_input_tsv

	script:
	"""
	Rscript ${projectDir}/modules/local/shift_back.R $MT_shifted_CollectMetrics $MT_CollectMetrics
	mv per_base_coverage.tsv ${meta.id}_per_base_coverage.tsv

	echo "${meta.id}\t/Individual/MT/QC/${meta.id}_per_base_coverage.tsv\t${meta.id}" > ${meta.id}_MT_Step1_input_tsv.tsv
	
	"""
}
