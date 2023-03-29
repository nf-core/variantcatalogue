process MT_Step2_participant_data {
	label 'process_medium'

	input :
        tuple val(meta), path(Sample_MT_Step2_participant_data)
	tuple val(meta), path(Sample_list)

        output :
        tuple val(meta), path('*.tsv'), emit : MT_Step2_participant_data_tsv
	tuple val(meta), path('*.txt'), emit : participants_to_subset_txt

        script:
	"""
	echo "entity:participant_id\ts\tvcf_path" > header_MT_Step2_participant_data
	cat header_MT_Step2_participant_data $Sample_MT_Step2_participant_data > MT_Step2_participant_data.tsv
	
	echo "participant" > header_participants_to_subset
	cat header_participants_to_subset $Sample_list > MT_participants_to_subset.txt
	"""
}

