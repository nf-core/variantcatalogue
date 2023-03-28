process list_vcfs_txt {
	label 'process_low'
        
	input :
        file individual_vcf
	val var_type

        output :
        tuple val('combined'), path('*.txt')

        script:
	if(var_type == "MT") {	
		"""
                find $params.outdir_ind/${assembly}/${batch}/${run}/${var_type}/Sample/ -name "*_filtered_sites.vcf.gz" > MT_vcfs.txt
		"""
	} else if (var_type == "SNV") {
                """
                find $params.outdir_ind/${assembly}/${batch}/${run}/${var_type}/Sample/ -name "*.g.vcf.gz" > ${var_type}_vcfs.txt
                """
	} else if (var_type == "SV") {
                """
                find $params.outdir_ind/${assembly}/${batch}/${run}/${var_type}/Sample/paragraph/ -name "*.vcf.gz" > ${var_type}_vcfs.txt
		"""
        } else {
        	"""
                find $params.outdir_ind/${assembly}/${batch}/${run}/${var_type}/Sample/ -name "*.vcf.gz" > ${var_type}_vcfs.txt
        	"""
	}

}

