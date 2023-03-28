process Picard_QualityScoreDistribution {
	tag "$meta.id"
	label 'process_low'

    conda "bioconda::picard=2.27.4 r::r-base"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.27.4--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.27.4--hdfd78af_0' }"

	input :
	tuple val(meta), path(bam), path (bai)

	output :
	tuple val(meta), path("*_qual_score_dist.*"), emit: report

	script :
	"""
	picard "-Xmx2G" QualityScoreDistribution \
		I=${bam} \
		O=${meta.id}_qual_score_dist.txt \
		CHART= ${meta.id}_qual_score_dist.pdf
	"""
}
