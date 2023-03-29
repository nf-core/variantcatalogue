process MT_Step3_metadata_sample {
	label 'process_medium'
	
    conda "bioconda::r-stringr=1.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-stringr%3A1.1.0--r3.3.1_0':
        'quay.io/biocontainers/r-stringr:1.1.0--r3.3.1_0' }"

	input :
	tuple val(meta), path(mosdepth)
	tuple val(meta), path(haplocheck)
	
	output :
	path('*'), emit : MT_Step3_metadata_sample

	script:
	"""
	Rscript ${projectDir}/modules/local/MT_Step3_metadata_sample.R ${mosdepth} ${haplocheck} 
	mv conta_cov.tsv ${meta.id}_conta_cov.tsv
	"""
}
