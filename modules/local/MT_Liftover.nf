process MT_Liftover {
        tag "${meta.id}"
	label 'process_low'

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"
	
	input :
        tuple val(meta), path(MT_call_variants_shifted), path(MT_call_variants_shifted_index)
	tuple val(meta2), path(ref_genome_MT_file)
	tuple val(meta3), path(ref_genome_MT_file_dict)
	tuple val(meta4), path(bwa_index_ref_genome)
	file ShiftBack_chain_MT_file

        output :
        tuple val(meta), path ('*_lifted_over.vcf')      , emit : lifted_vcf
	tuple val(meta), path ('*_rejected_variants.vcf'), emit : rejected_vcf

        script :
        """
	gatk LiftoverVcf \
		I=${MT_call_variants_shifted} \
		O=${meta.id}_lifted_over.vcf \
		CHAIN=${ShiftBack_chain_MT_file} \
		REJECT=${meta.id}_rejected_variants.vcf \
		R=${ref_genome_MT_file}
	"""
}
