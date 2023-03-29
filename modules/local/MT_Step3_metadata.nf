process MT_Step3_metadata {
	label 'process_medium'
        
	input :
        path(MT_Step3_metadata_sample)

        output :
        path('*.tsv'), emit : tsv

        script:
	"""
	echo "entity:participant_id\ts\tcontamination\twgs_mean_coverage\tmt_mean_coverage" > header_MT_metadata_step3
	cat header_MT_metadata_step3 $MT_Step3_metadata_sample > MT_Step3_participant_data.tsv
	"""
}

