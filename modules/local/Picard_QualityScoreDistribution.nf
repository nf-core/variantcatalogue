process Picard_QualityScoreDistribution {
	tag "$meta.id"
	label 'process_low'

    conda "bioconda::picard=3.0.0 r::r-base"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.0.0--hdfd78af_1' :
        'quay.io/biocontainers/picard:3.0.0--hdfd78af_1' }"

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
