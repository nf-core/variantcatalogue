process Picard_CollectAlignmentSummaryMetrics {
	ag "$meta.id"
	label 'process_medium'

    conda "bioconda::gatk4=4.3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.3.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0' }"
	
	input :
	tuple val (meta), path (bam), path (bai)

	output :
	tuple val(meta), path("*"), emit: report

	script :
	"""
	gatk CollectAlignmentSummaryMetrics \
		--java-options "-Xmx2000M" \
		-I ${bam} \
		-O ${meta.id}_Picard_Alignment
	"""
}
