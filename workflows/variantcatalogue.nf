/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowVariantcatalogue.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { INPUT_CHECK } from '../subworkflows/local/input_check'

//
// For the mapping subworkflow
//

include { Picard_CollectAlignmentSummaryMetrics } from '../modules/local/Picard_CollectAlignmentSummaryMetrics'
include { Picard_QualityScoreDistribution       } from '../modules/local/Picard_QualityScoreDistribution'


//
// For the SNV/indel subworkflow
//

include { DEEPVARIANT     } from '../modules/local/deepvariant/main'
include { Hail_sample_QC  } from '../modules/local/Hail_sample_QC'
include { Hail_variant_QC } from '../modules/local/Hail_variant_QC'


include { GATK4_PRINTREADS } from '../modules/local/gatk4/printreads/main'
include { shift_back       } from '../modules/local/shift_back' 
include { PICARD_COLLECTHSMETRICS as PICARD_COLLECTHSMETRICS_MT; PICARD_COLLECTHSMETRICS as PICARD_COLLECTHSMETRICS_MT_shifted } from '../modules/local/collecthsmetrics/main'
include { MT_Liftover      } from '../modules/local/MT_Liftover'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

//
// For the mapping subworkflow
//

include { BWA_INDEX                } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM                  } from '../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_INDEX           } from '../modules/nf-core/samtools/index/main'
include { TRIMMOMATIC              } from '../modules/nf-core/trimmomatic/main'

include { MOSDEPTH                 } from '../modules/nf-core/mosdepth/main'
include { PICARD_COLLECTWGSMETRICS } from '../modules/nf-core/picard/collectwgsmetrics/main'
include { FASTQC                   } from '../modules/nf-core/fastqc/main'
include { MULTIQC                  } from '../modules/nf-core/multiqc/main'

//
// For the SNV/indel subworkflow
//

include { GLNEXUS           } from '../modules/nf-core/glnexus/main'
include { TABIX_TABIX       } from '../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_NORM     } from '../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_ANNOTATE } from '../modules/nf-core/bcftools/annotate/main'
include { ENSEMBLVEP_VEP    } from '../modules/nf-core/ensemblvep/vep/main'

//
// For the MT subworkflow
//

include { BWA_INDEX as BWA_INDEX_MT; BWA_INDEX as BWA_INDEX_MT_shifted       } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM as BWA_MEM_MT; BWA_MEM as BWA_MEM_MT_shifted               } from '../modules/nf-core/bwa/mem/main'
include { GATK4_SAMTOFASTQ                                                   } from '../modules/nf-core/gatk4/samtofastq/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MT; SAMTOOLS_INDEX as SAMTOOLS_INDEX_MT_shifted                         } from '../modules/nf-core/samtools/index/main'
include { GATK4_MARKDUPLICATES as GATK4_MARKDUPLICATES_MT; GATK4_MARKDUPLICATES as GATK4_MARKDUPLICATES_MT_shifted } from '../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_MUTECT2; GATK4_MUTECT2 as GATK4_MUTECT2_shifted               } from '../modules/nf-core/gatk4/mutect2/main'
include { TABIX_TABIX as TABIX_TABIX_MT                                       } from '../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_MT                                   } from '../modules/nf-core/bcftools/norm/main'

include { GATK4_VARIANTFILTRATION } from '../modules/nf-core/gatk4/variantfiltration/main'  
include { GATK4_FILTERMUTECTCALLS } from '../modules/nf-core/gatk4/filtermutectcalls/main' 
include { GATK4_LEFTALIGNANDTRIMVARIANTS } from '../modules/nf-core/gatk4/leftalignandtrimvariants/main'
include { GATK4_MERGEVCFS } from '../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEMUTECTSTATS } from '../modules/nf-core/gatk4/mergemutectstats/main'  

include { HAPLOCHECK } from '../modules/nf-core/haplocheck/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow VARIANTCATALOGUE {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)


    reference = [
        [ id:'reference', single_end:true ],
        file(params.reference)
	]

    reference_index = [
        [ id:'reference', single_end:true ],
        file(params.reference_index)
        ]

    //
    // Subworkflow : Mapping
    //

    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    TRIMMOMATIC    (INPUT_CHECK.out.reads )
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions.first())

    BWA_INDEX      ( reference )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions.first())

    BWA_MEM        (TRIMMOMATIC.out.trimmed_reads, BWA_INDEX.out.index, true )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    SAMTOOLS_INDEX (BWA_MEM.out.bam ) 
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_bam = BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai)       // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_bam.map { meta, bam, bai ->
                        return [meta, bam, bai, []]
            }
            .set { ch_mosdepth_in }

    MOSDEPTH ( ch_mosdepth_in, [[:],[]]) 
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())

    PICARD_COLLECTWGSMETRICS ( BWA_MEM.out.bam
                                        .join(SAMTOOLS_INDEX.out.bai),
				reference,
				reference_index,
				[]
			)
    ch_versions = ch_versions.mix(PICARD_COLLECTWGSMETRICS.out.versions.first())

    Picard_CollectAlignmentSummaryMetrics ( BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai))
    Picard_QualityScoreDistribution       ( BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai))


    //
    // Subworkflow : SNV/indel
    //

    ch_bam = BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai)       // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_bam.map { meta, bam, bai ->
                        return [meta, bam, bai, file(params.test_bed)]
            }
            .set { ch_deepvar_in }


    DEEPVARIANT       (ch_deepvar_in, file(params.reference), file(params.reference_index))
    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions.first())

    ch_case_info = Channel.from([id:'SNV'])
 
        DEEPVARIANT.out.gvcf
            .collect{it[1]}
            .toList()
            .collect()
            .set { ch_file_list }

        ch_case_info
            .combine(ch_file_list)
            .set { ch_gvcfs }

    GLNEXUS           (ch_gvcfs)
    ch_versions = ch_versions.mix(GLNEXUS.out.versions.first())

    TABIX_TABIX       (GLNEXUS.out.bcf)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    BCFTOOLS_NORM     (GLNEXUS.out.bcf.join(TABIX_TABIX.out.csi), file(params.reference))
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())

    BCFTOOLS_ANNOTATE (BCFTOOLS_NORM.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())

    Hail_sample_QC    (BCFTOOLS_ANNOTATE.out.vcf)
    Hail_variant_QC   (Hail_sample_QC.out.vcf_sample_filtered, Hail_sample_QC.out.filtered_sample_sex)
//    ENSEMBLVEP_VEP    (Hail_variant_QC.out.vcf_SNV_filtered_frequ_only, params.genome, params.species, params.cache_version, params.cache_path, file(params.reference), params.vep_extra_files )
//    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions.first())


    //
    // Subworkflow : MT
    //

    reference_MT = [
        [ id:'reference_MT', single_end:true ],
        file(params.reference_MT)
        ]

    reference_MT_index = [
        [ id:'reference_MT_index', single_end:true ],
        file(params.reference_MT_index)
        ]

    reference_MT_dict = [
        [ id:'reference_MT_dict', single_end:true ],
        file(params.reference_MT_dict)
        ]

    reference_MT_shifted = [
        [ id:'reference_MT_shifted', single_end:true ],
        file(params.reference_MT_shifted)
        ]

    reference_MT_shifted_index = [
        [ id:'reference_MT_shifted_index', single_end:true ],
        file(params.reference_MT_shifted_index)
        ]

    reference_MT_shifted_dict = [
        [ id:'reference_MT_shifted_dict', single_end:true ],
        file(params.reference_MT_shifted_dict)
        ]

    BWA_INDEX_MT                       (reference_MT)
    BWA_INDEX_MT_shifted               (reference_MT_shifted)
    GATK4_PRINTREADS                   (BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai), reference_MT, file(params.reference_MT_index), file(params.reference_MT_dict))
    ch_versions = ch_versions.mix(GATK4_PRINTREADS.out.versions.first())

    GATK4_SAMTOFASTQ                   (GATK4_PRINTREADS.out.bam)
    ch_versions = ch_versions.mix(GATK4_SAMTOFASTQ.out.versions.first())

    BWA_MEM_MT                         (GATK4_SAMTOFASTQ.out.fastq, BWA_INDEX_MT.out.index, true)
    BWA_MEM_MT_shifted                 (GATK4_SAMTOFASTQ.out.fastq, BWA_INDEX_MT_shifted.out.index, true)
    SAMTOOLS_INDEX_MT                  (BWA_MEM_MT.out.bam)
    SAMTOOLS_INDEX_MT_shifted          (BWA_MEM_MT_shifted.out.bam)
    GATK4_MARKDUPLICATES_MT            (BWA_MEM_MT.out.bam, [], [])
    GATK4_MARKDUPLICATES_MT_shifted    (BWA_MEM_MT_shifted.out.bam, [], [])
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES_MT.out.versions.first())

    PICARD_COLLECTHSMETRICS_MT         (GATK4_MARKDUPLICATES_MT.out.bam.join(GATK4_MARKDUPLICATES_MT.out.bai), reference_MT, reference_MT_index, file(params.non_control_region_interval_list))
    PICARD_COLLECTHSMETRICS_MT_shifted (GATK4_MARKDUPLICATES_MT_shifted.out.bam.join(GATK4_MARKDUPLICATES_MT_shifted.out.bai), reference_MT_shifted, reference_MT_shifted_index, file(params.control_region_shifted_reference_interval_list))
    ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS_MT.out.versions.first())

    ch_bam_MT = GATK4_MARKDUPLICATES_MT.out.bam.join(GATK4_MARKDUPLICATES_MT.out.bai)       // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_bam_MT.map { meta, bam, bai ->
                        return [meta, bam, bai, []]
            }
            .set { ch_mutect2_in }

    GATK4_MUTECT2(ch_mutect2_in, file(params.reference_MT), file(params.reference_MT_index), file(params.reference_MT_dict), [], [], [], [])

//    ch_bam_MT_shifted = GATK4_MARKDUPLICATES_MT_shifted.out.bam.join(GATK4_MARKDUPLICATES_MT_shifted.out.bai)       // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
//    ch_bam_MT_shifted.map { meta, bam, bai ->
//                        return [meta, bam, bai, []]
 //           }
//            .set { ch_mutect2_shifted_in }

//    GATK4_MUTECT2_shifted(ch_mutect2_shifted_in, reference_MT_shifted, reference_MT_index_shifted, reference_MT_dict_shifted, [], [], [], [])

//    MT_Liftover(GATK4_MUTECT2_shifted.out.vcf.join(GATK4_MUTECT2_shifted.out.tbi), reference_MT, reference_dict, BWA_INDEX_MT.out.index, file(params.ShiftBack_chain))
//    GATK4_MERGEVCFS(MT_Liftover.out.vcf.join(GATK4_MUTECT2.out.vcf), reference_dict)
//    TABIX_TABIX_MT (GATK4_MERGEVCFS.out.vcf)
//    BCFTOOLS_NORM_MT(GATK4_MERGEVCFS.out.vcf.join(TABIX_TABIX_MT.out.tbi), file(params.reference_MT))
//    TABIX_TABIX_MT2 (BCFTOOLS_NORM_MT.out.vcf)
//    GATK4_MERGEMUTECTSTATS(GATK4_MUTECT2.out.stats.join(GATK4_MUTECT2_shifted.out.stat))


//    ch_vcf_MT = BCFTOOLS_NORM_MT.out.vcf
//		.join(TABIX_TABIX_MT2.out.tbi)
//		.join(GATK4_MERGEMUTECTSTATS.out.stat)
//    ch_vcf_MT.map { meta, vcf, tbi, stat ->
//                        return [meta, vcf, tbi, stat, [], [], [], []]
 //           }
//            .set { ch_vcf_MT_in }


//    GATK4_FILTERMUTECTCALLS(ch_vcf_MT_in, file(params.reference_MT), file(params.reference_MT_index), file(params.reference_MT_dict) )

//    ch_vcf_MT_filtered = GATK4_FILTERMUTECTCALLS.out.vcf
//              .join(GATK4_FILTERMUTECTCALLS.out.tbi)
//    ch_vcf_MT_filtered.map { meta, vcf, tbi ->
//                        return [meta, vcf, tbi, []]
 //           }
//            .set { ch_vcf_MT_filtered_in }
//    GATK4_LEFTALIGNANDTRIMVARIANTS(ch_vcf_MT_filtered_in, file(params.reference_MT), file(params.reference_MT_index), file(params.reference_MT_dict))









  CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowVariantcatalogue.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowVariantcatalogue.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTWGSMETRICS.out.metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(Picard_CollectAlignmentSummaryMetrics.out.report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(Picard_QualityScoreDistribution.out.report.collect{it[1]}.ifEmpty([]))


    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
