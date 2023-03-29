process GATK4_VARIANTFILTRATION {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  fasta
    path  fai
    path  dict

    output:
    tuple val(meta), path("*.vcf.gz")                        , emit: vcf
    tuple val(meta), path("*.tbi")                           , emit: tbi
    tuple val(meta), path('*_MT_Step2_participant_data.tsv') , emit: sample_MT_Step2_participant_data
    tuple val(meta), path('*_list.txt')                      , emit: Sample_list
    path "versions.yml"		                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK VariantFiltration] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M" VariantFiltration \\
        --variant $vcf \\
        --output ${prefix}.vcf.gz \\
        --reference $fasta \\
        --tmp-dir . \\
        --mask-name "GATK_artifact" \\
	--mask ${params.blacklist_sites_hg38_MT_file} \\
        $args

		echo "${meta.id}\t${meta.id}\tIndividual/MT/Sample/${meta.id}_filtered_sites.vcf.gz" > ${meta.id}_MT_Step2_participant_data.tsv
		echo "${meta.id}" > ${meta.id}_list.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
