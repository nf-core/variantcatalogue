process Hail_sample_QC {
	label 'process_medium'

	conda "bioconda::hail=0.2.58"
    
	input :
	tuple val (meta), path (SNV_vcf)

	output :
	tuple val(meta), path('*.html')	                  , emit : graph               , optional : true
	tuple val(meta), path('*filtered_samples.vcf.bgz'), emit : vcf_sample_filtered
	tuple val(meta), path('*filtered_samples_sex.tsv'), emit : filtered_sample_sex

	script:
	"""
	mkdir -p $params.tmp_dir
	python ${projectDir}/modules/Hail_sample_QC.py $SNV_vcf $params.tmp_dir $params.assembly
	"""
}
