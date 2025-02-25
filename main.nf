#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    INPUT FILES
======================================================================================== */
params.input = null

if (!params.input) {
    error "Input not specified. Use --input to specify the input."
}

input_files = file(params.input)


/* ========================================================================================
    OUTPUT DIRECTORY
======================================================================================== */
params.outdir = false

if(params.outdir){
    outdir = params.outdir
} else {
    outdir = '.'
}


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.genome            = 'Mus_musculus_GRCm39'
params.single_end        = false  // default mode is auto-detect.
params.seq_method        = 'WGBS'


/* ========================================================================================
    PACKAGE PARAMETERS
======================================================================================== */
params.fastqc_args                        = ''
params.fastq_screen_args                  = ''
params.trim_galore_args                   = ''
params.bismark_args                       = ''
params.filter_non_conversion_args         = ''
params.deduplicate_bismark_args           = ''
params.bismark_methylation_extractor_args = ''
params.coverage2cytosine_args             = ''
params.bismark2summary_args               = ''
params.bismark2report_args                = ''
params.multiqc_args                       = ''

fastqc_args                        = params.fastqc_args
fastq_screen_args                  = params.fastq_screen_args
trim_galore_args                   = params.trim_galore_args
bismark_args                       = params.bismark_args
filter_non_conversion_args         = params.filter_non_conversion_args
deduplicate_bismark_args           = params.deduplicate_bismark_args
bismark_methylation_extractor_args = params.bismark_methylation_extractor_args
coverage2cytosine_args             = params.coverage2cytosine_args
bismark2summary_args               = params.bismark2summary_args
bismark2report_args                = params.bismark2report_args
multiqc_args                       = params.multiqc_args


/* ========================================================================================
    SKIP STEPS
======================================================================================== */
params.skip_qc                   = false
params.skip_fastq_screen         = false
params.skip_deduplication        = false
params.add_filter_non_conversion = false


/* ========================================================================================
    SEQUENCING METHOD PARAMETERS
======================================================================================== */

// Trim Galore parameters
params.clip_r1 = false
clip_r1        = params.clip_r1

params.clip_r2 = false
clip_r2        = params.clip_r2

// Bismark methylation extractor parameters
params.ignore_r2 = false
ignore_r2        = params.ignore_r2


// WGBS
if (params.seq_method == 'WGBS'){

    fastq_screen_args += " --bisulfite "
    ignore_r2          = 2 // Bismark methylation extractor 

}

// PBAT
if (params.seq_method == 'PBAT'){

    fastq_screen_args += " --bisulfite "
    bismark_args      += " --pbat "
    clip_r1            = 9       // Trim Galore
    clip_r2            = clip_r1 // Trim Galore
    trim_galore_args  += " --clip_R1 ${clip_r1} "

}

// RRBS
if (params.seq_method == 'RRBS'){

    fastq_screen_args += " --bisulfite "
    trim_galore_args  += " --rrbs "

}

// NOMe-Seq
if (params.seq_method == 'NOMe-Seq'){

    fastq_screen_args      += " --bisulfite "
    coverage2cytosine_args += " --nome "
    ignore_r2               = 2 // Bismark methylation extractor

}

// scBS-Seq
if (params.seq_method == 'scBS-Seq'){

    fastq_screen_args += " --bisulfite "
    bismark_args      += " --non_directional "
    clip_r1            = 6       // Trim Galore
    clip_r2            = clip_r1 // Trim Galore
    trim_galore_args  += " --clip_R1 ${clip_r1} "

}

// scNMT-Seq
if (params.seq_method == 'scNMT-Seq'){

    fastq_screen_args                  += " --bisulfite "
    bismark_args                       += " --non_directional "
    bismark_methylation_extractor_args += " --CX "
    clip_r1                             = 6       // Trim Galore
    clip_r2                             = clip_r1 // Trim Galore
    trim_galore_args                   += " --clip_R1 ${clip_r1} "

}


/* ========================================================================================
    FASTQ SCREEN PARAMETERS
======================================================================================== */
params.fastq_screen_conf = "/cluster/work/nme/software/config/fastq_screen.conf" // FastQ Screen config file directory


/* ========================================================================================
    BISMARK PARAMETERS
======================================================================================== */

// --unmapped
// Write all reads that could not be aligned to a file in the output directory.
params.unmapped  = false
if (params.unmapped){
		bismark_args += " --unmapped "
	}

// --ambiguous
// Write all reads which produce more than one valid alignment with the same number of
// lowest mismatches or other reads that fail to align uniquely to a file in the output directory.
params.ambiguous = false 
if (params.ambiguous){
        bismark_args += " --ambiguous "
    }

// --local
// In this mode, it is not required that the entire read aligns from one end to the other. Rather, some
// characters may be omitted (“soft-clipped”) from the ends in order to achieve the greatest possible
// alignment score.
params.local = false 
if (params.local){
        bismark_args += " --local "
    }
        

/* ========================================================================================
    GENOMES
======================================================================================== */
params.custom_genome_file = '' // Option to add a directory for a custom genome file

include { getGenome } from './modules/genomes.mod.nf' params(custom_genome_file: params.custom_genome_file)
genome = getGenome(params.genome)


/* ========================================================================================
    FILES CHANNEL
======================================================================================== */
include { makeFilesChannel; getFileBaseNames } from './modules/files.mod.nf'
file_ch = makeFilesChannel(input_files)


/* ========================================================================================
    WORKFLOW
======================================================================================== */
include { FASTQC }                               from './modules/fastqc.mod.nf'                         
include { FASTQC as FASTQC2 }                    from './modules/fastqc.mod.nf' 
include { FASTQ_SCREEN }                         from './modules/fastq_screen.mod.nf'                  params(fastq_screen_conf: params.fastq_screen_conf)
include { TRIM_GALORE }                          from './modules/trim_galore.mod.nf'                   params(clip_r2: clip_r2)
include { BISMARK }                              from './modules/bismark.mod.nf'                       params(genome: genome)
include { BISMARK_FILTER_NON_CONVERSION }        from './modules/bismark_filter_non_conversion.mod.nf'
include { BISMARK_DEDUPLICATION }                from './modules/bismark_deduplication.mod.nf'
include { BISMARK_METHYLATION_EXTRACTOR }        from './modules/bismark_methylation_extractor.mod.nf' params(ignore_r2: ignore_r2)
include { COVERAGE2CYTOSINE }                    from './modules/coverage2cytosine.mod.nf'             params(genome: genome)
include { BISMARK2REPORT }                       from './modules/bismark2report.mod.nf'
include { BISMARK2SUMMARY }                      from './modules/bismark2summary.mod.nf'
include { MULTIQC }                              from './modules/multiqc.mod.nf' 

workflow {

    main:
        // QC conditional
        if (!params.skip_qc){ 
            FASTQC                          (file_ch, outdir, fastqc_args)
            
            // FastQ Screen conditional
            if (!params.skip_fastq_screen){ 
            FASTQ_SCREEN                    (file_ch, outdir, fastq_screen_args)
            }

            TRIM_GALORE                     (file_ch, outdir, trim_galore_args)
            FASTQC2                         (TRIM_GALORE.out.reads, outdir, fastqc_args)
            BISMARK                         (TRIM_GALORE.out.reads, outdir, bismark_args)
        } else {
            BISMARK                         (file_ch, outdir, bismark_args)
        }

        // Deduplication conditional
        if (!params.skip_deduplication){ 
        
                // Filter non-conversion conditional
                if (params.add_filter_non_conversion){ 
                    BISMARK_FILTER_NON_CONVERSION   (BISMARK.out.bam, outdir, filter_non_conversion_args)
                    BISMARK_DEDUPLICATION           (BISMARK_FILTER_NON_CONVERSION.out.bam, outdir, deduplicate_bismark_args)
                    BISMARK_METHYLATION_EXTRACTOR   (BISMARK_DEDUPLICATION.out.bam, outdir, bismark_methylation_extractor_args)
                } else {
                    BISMARK_DEDUPLICATION           (BISMARK.out.bam, outdir, deduplicate_bismark_args)
                    BISMARK_METHYLATION_EXTRACTOR   (BISMARK_DEDUPLICATION.out.bam, outdir, bismark_methylation_extractor_args)
                }

        } else {

                // Filter non-conversion conditional
                if (params.add_filter_non_conversion){ 
                    BISMARK_FILTER_NON_CONVERSION   (BISMARK.out.bam, outdir, filter_non_conversion_args)
                    BISMARK_METHYLATION_EXTRACTOR   (BISMARK_FILTER_NON_CONVERSION.out.bam, outdir, bismark_methylation_extractor_args)
                } else {
                    BISMARK_METHYLATION_EXTRACTOR   (BISMARK.out.bam, outdir, bismark_methylation_extractor_args)
                }
                
        }

        // scNMT-Seq conditional
        if (params.seq_method == 'scNMT-Seq'){ 
            COVERAGE2CYTOSINE             (BISMARK_METHYLATION_EXTRACTOR.out.coverage, outdir, coverage2cytosine_args)
        }


        /* ========================================================================================
            Reports
        ======================================================================================== */

        // Merging channels for MultiQC
        if (!params.skip_qc){

            multiqc_ch = FASTQC.out.report.mix(
                            TRIM_GALORE.out.report.ifEmpty([]),
                            FASTQC2.out.report.ifEmpty([])
                            ).collect()

            if (!params.skip_fastq_screen){
                multiqc_ch = multiqc_ch.mix(
                            FASTQ_SCREEN.out.report.ifEmpty([])
                            ).collect()
            }

            multiqc_ch = multiqc_ch.mix(
                            BISMARK.out.report.ifEmpty([])
                            ).collect()

        } else {

            multiqc_ch = BISMARK.out.report.ifEmpty([]).collect()

        }

        if (!params.skip_deduplication){

            if (params.add_filter_non_conversion){ 
                
                // with deduplication & with filter non-conversion
                multiqc_ch = multiqc_ch.mix(
                            BISMARK_FILTER_NON_CONVERSION.out.report.ifEmpty([]),
                            BISMARK_DEDUPLICATION.out.report.ifEmpty([]),
                            BISMARK_METHYLATION_EXTRACTOR.out.report.ifEmpty([]),
                            BISMARK_METHYLATION_EXTRACTOR.out.mbias.ifEmpty([])
                            ).collect() 

                bismark_report_ch = BISMARK.out.bam.mix(
                            BISMARK.out.report.ifEmpty([]),
                            BISMARK_FILTER_NON_CONVERSION.out.report.ifEmpty([]),
                            BISMARK_DEDUPLICATION.out.report.ifEmpty([]),
                            BISMARK_METHYLATION_EXTRACTOR.out.report.ifEmpty([]),
                            BISMARK_METHYLATION_EXTRACTOR.out.mbias.ifEmpty([])
                            ).collect()                            

            } else { 

                // with deduplication & without filter non-conversion
                multiqc_ch = multiqc_ch.mix(
                            BISMARK_DEDUPLICATION.out.report.ifEmpty([]),
                            BISMARK_METHYLATION_EXTRACTOR.out.report.ifEmpty([]),
                            BISMARK_METHYLATION_EXTRACTOR.out.mbias.ifEmpty([])
                            ).collect()

                bismark_report_ch = BISMARK.out.bam.mix(
                            BISMARK.out.report.ifEmpty([]),
                            BISMARK_DEDUPLICATION.out.report.ifEmpty([]),
                            BISMARK_METHYLATION_EXTRACTOR.out.report.ifEmpty([]),
                            BISMARK_METHYLATION_EXTRACTOR.out.mbias.ifEmpty([])
                            ).collect()   

            }   

        } else { 
            if (params.add_filter_non_conversion){ 

                // without deduplication & with filter non-conversion
                multiqc_ch = multiqc_ch.mix(
                        BISMARK_FILTER_NON_CONVERSION.out.report.ifEmpty([]),
                        BISMARK_METHYLATION_EXTRACTOR.out.report.ifEmpty([]),
                        BISMARK_METHYLATION_EXTRACTOR.out.mbias.ifEmpty([])
                        ).collect()

                bismark_report_ch = BISMARK.out.bam.mix(
                        BISMARK.out.report,
                        BISMARK_FILTER_NON_CONVERSION.out.report.ifEmpty([]),
                        BISMARK_METHYLATION_EXTRACTOR.out.report.ifEmpty([]),
                        BISMARK_METHYLATION_EXTRACTOR.out.mbias.ifEmpty([])
                        ).collect()     


            } else {

                // without deduplication & without filter non-conversion
                multiqc_ch = multiqc_ch.mix(
                            BISMARK_METHYLATION_EXTRACTOR.out.report.ifEmpty([]),
                            BISMARK_METHYLATION_EXTRACTOR.out.mbias.ifEmpty([])
                            ).collect()

                bismark_report_ch = BISMARK.out.bam.mix(
                            BISMARK.out.report,
                            BISMARK_METHYLATION_EXTRACTOR.out.report.ifEmpty([]),
                            BISMARK_METHYLATION_EXTRACTOR.out.mbias.ifEmpty([])
                            ).collect()                  
                
            }
        }

        // with coverage2cytosine
        if (params.seq_method == 'scNMT-Seq'){

            multiqc_ch = multiqc_ch.mix(
                        COVERAGE2CYTOSINE.out.report.ifEmpty([])
                        ).collect()

            bismark_report_ch = bismark_report_ch.mix(
                        COVERAGE2CYTOSINE.out.report.ifEmpty([])
                        ).collect()             

        }

        /* ========================================================================================
        ======================================================================================== */

        BISMARK2REPORT      (bismark_report_ch, outdir, bismark2report_args)
        BISMARK2SUMMARY     (bismark_report_ch, outdir, bismark2summary_args)
        MULTIQC             (multiqc_ch, outdir, multiqc_args)
}
