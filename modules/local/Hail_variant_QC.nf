process Hail_variant_QC {
	label 'process_medium'

	conda "bioconda::hail=0.2.58"

	input :
	tuple val(meta) , path(vcf_sample_filtered)
	tuple val(meta2), path(sample_sex_file)

	output :
	tuple val(meta), path('*.html')                             , emit : graph
	tuple val(meta), path('SNV_indel_QC_report.txt')            , emit : SNV_QC_report

	tuple val(meta), path('SNV_filtered_with_geno.vcf.bgz*')    , emit : SNV_filtered_with_geno
        tuple val(meta), path('SNV_filtered_frequ_only.vcf.bgz')    , emit : vcf_SNV_filtered_frequ_only
        tuple val(meta), path('SNV_filtered_frequ_only.vcf.bgz.tbi'), emit : index_SNV_filtered_frequ_only

	script:
	"""
	mkdir -p $params.tmp_dir
        python ${projectDir}/modules/Hail_variant_QC.py $vcf_sample_filtered $sample_sex_file $params.tmp_dir $params.assembly
	"""
}
