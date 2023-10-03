#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    INPUT FILES
======================================================================================== */
params.input = null
input_files  = file(params.input)


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
    SKIP STEPS
======================================================================================== */
params.skip_qc                   = false
params.skip_fastq_screen         = false
params.skip_deduplication        = false
params.add_filter_non_conversion = false


/* ========================================================================================
    SEQUENCING TYPE PARAMETERS
======================================================================================== */
params.seqtype     = 'WGBS'
params.scNMT       = false
params.pbat        = false
params.rrbs        = false
params.singlecell  = false
params.nonCG       = false
params.nome        = false
params.pbat_trim   = false

scNMT       = params.scNMT     
pbat        = params.pbat      
rrbs        = params.rrbs      
singlecell  = params.singlecell
nonCG       = params.nonCG     
nome        = params.nome      
pbat_trim   = params.pbat_trim 

// PBAT (with Trim Galore --clip_R1 9 --clip_R2 9)
if (params.seqtype == 'PBAT'){
    pbat = true
    pbat_trim = '9'
}

// PBAT (without Trim Galore --clip_R1 9 --clip_R2 9)
else if (params.seqtype == 'PBAT without default TrimGalore clip'){
    pbat = true
}

// RRBS
else if (params.seqtype == 'RRBS'){
    rrbs = true
}

// NOME-seq
else if (params.seqtype == 'NOME-seq'){
    nome = true
}

// scBS-seq
else if (params.seqtype == 'scBS-seq'){
    singlecell = '6'
}

// scNMT-seq
else if (params.seqtype == 'scNMT-seq'){
    singlecell = '6'
    nome       = true
    nonCG      = true
    scNMT      = true
}


/* ========================================================================================
    HELPER PARAMETERS
======================================================================================== */
params.fastq_screen_conf = "/cluster/work/nme/software/config/fastq_screen.conf" // FastQ Screen config file directory
params.genome            = ''
params.verbose           = false
params.single_end        = false  // default mode is auto-detect. NOTE: params are handed over automatically
params.help              = false
params.list_genomes      = false


/* ========================================================================================
    WORKFLOW STEPS PARAMETERS
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


/* ========================================================================================
    MESSAGES
======================================================================================== */
if (params.list_genomes){
    println ("[WORKLFOW] List genomes selected")
}

if (params.verbose){
    println ("[WORKFLOW] FASTQC ARGS: "                            + params.fastqc_args)
    println ("[WORKFLOW] FASTQ SCREEN ARGS ARE: "                  + params.fastq_screen_args)
    println ("[WORKFLOW] TRIM GALORE ARGS: "                       + params.trim_galore_args)
    println ("[WORKFLOW] BISMARK ARGS ARE: "                       + params.bismark_args)
    println ("[WORKFLOW] BISMARK FILTER NON-CONVERSION ARGS ARE: " + params.filter_non_conversion_args)
    println ("[WORKFLOW] BISMARK DEDUPLICATION ARGS ARE: "         + params.deduplicate_bismark_args)
    println ("[WORKFLOW] BISMARK METHYLATION EXTRACTOR ARGS ARE: " + params.bismark_methylation_extractor_args)
    println ("[WORKFLOW] BISMARK COVERAGE2CYTOSINE ARGS ARE: "     + params.coverage2cytosine_args)
    println ("[WORKFLOW] BISMARK2SUMMARY ARGS ARE: "               + params.bismark2summary_args)
    println ("[WORKFLOW] BISMARK2REPORT ARGS ARE: "                + params.bismark2report_args)
    println ("[WORKFLOW] MULTIQC ARGS: "                           + params.multiqc_args)
}


/* ========================================================================================
    GENOMES
======================================================================================== */
include { getGenome; listGenomes } from './modules/genomes.mod.nf'

if (params.list_genomes){
    listGenomes()
}

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
include { FASTQ_SCREEN }                         from './modules/fastq_screen.mod.nf'                  params(fastq_screen_conf: params.fastq_screen_conf, bisulfite: '--bisulfite')
include { TRIM_GALORE }                          from './modules/trim_galore.mod.nf'                   params(pbat: pbat_trim, singlecell: singlecell, rrbs: rrbs)
include { BISMARK }                              from './modules/bismark.mod.nf'                       params(genome: genome, pbat: pbat, singlecell: singlecell)
include { BISMARK_FILTER_NON_CONVERSION }        from './modules/bismark_filter_non_conversion.mod.nf'
include { BISMARK_DEDUPLICATION }                from './modules/bismark_deduplication.mod.nf'
include { BISMARK_METHYLATION_EXTRACTOR }        from './modules/bismark_methylation_extractor.mod.nf' params(pbat: pbat, singlecell: singlecell, rrbs: rrbs, nonCG: nonCG)
include { COVERAGE2CYTOSINE }                    from './modules/coverage2cytosine.mod.nf'             params(genome: genome, nome: nome)
include { BISMARK2REPORT }                       from './modules/bismark2report.mod.nf'
include { BISMARK2SUMMARY }                      from './modules/bismark2summary.mod.nf'
include { MULTIQC }                              from './modules/multiqc.mod.nf' 

workflow {

    main:
        // QC conditional
        if (!params.skip_qc){ 
            FASTQC                          (file_ch, outdir, params.fastqc_args, params.verbose)
            TRIM_GALORE                     (file_ch, outdir, params.trim_galore_args, params.verbose)
            // FastQ Screen conditional
            if (!params.skip_fastq_screen){ 
            FASTQ_SCREEN                    (TRIM_GALORE.out.reads, outdir, params.fastq_screen_args, params.verbose)
            }
            FASTQC2                         (TRIM_GALORE.out.reads, outdir, params.fastqc_args, params.verbose)
            BISMARK                         (TRIM_GALORE.out.reads, outdir, params.bismark_args, params.verbose)
        } else {
            BISMARK                         (file_ch, outdir, params.bismark_args, params.verbose)
        }

        // Deduplication conditional
        if (!params.skip_deduplication){ 
        
                // Filter non-conversion conditional
                if (params.add_filter_non_conversion){ 
                    BISMARK_FILTER_NON_CONVERSION   (BISMARK.out.bam, outdir, params.filter_non_conversion_args, params.verbose)
                    BISMARK_DEDUPLICATION           (BISMARK_FILTER_NON_CONVERSION.out.bam, outdir, params.deduplicate_bismark_args, params.verbose)
                    BISMARK_METHYLATION_EXTRACTOR   (BISMARK_DEDUPLICATION.out.bam, outdir, params.bismark_methylation_extractor_args, params.verbose)
                } else {
                    BISMARK_DEDUPLICATION           (BISMARK.out.bam, outdir, params.deduplicate_bismark_args, params.verbose)
                    BISMARK_METHYLATION_EXTRACTOR   (BISMARK_DEDUPLICATION.out.bam, outdir, params.bismark_methylation_extractor_args, params.verbose)
                }
        } else {
                // Filter non-conversion conditional
                if (params.add_filter_non_conversion){ 
                    BISMARK_FILTER_NON_CONVERSION   (BISMARK.out.bam, outdir, params.filter_non_conversion_args, params.verbose)
                    BISMARK_METHYLATION_EXTRACTOR   (BISMARK_FILTER_NON_CONVERSION.out.bam, outdir, params.bismark_methylation_extractor_args, params.verbose)
                } else {
                    BISMARK_METHYLATION_EXTRACTOR   (BISMARK.out.bam, outdir, params.bismark_methylation_extractor_args, params.verbose)
                }
        }

        // scNME conditional
        if (scNMT){ 
            COVERAGE2CYTOSINE             (BISMARK_METHYLATION_EXTRACTOR.out.coverage, outdir, params.coverage2cytosine_args, params.verbose)
        }


        /* ========================================================================================
            Reports
        ======================================================================================== */

        // Merging channels for MultiQC
        if (!params.skip_qc){

            multiqc_ch = FASTQC.out.report.mix(
                         TRIM_GALORE.out.report,
                         FASTQC2.out.report.ifEmpty([])
                         ).collect()

            if (!params.skip_fastq_screen){
                multiqc_ch = multiqc_ch.mix(
                            FASTQ_SCREEN.out.report.ifEmpty([])
                            ).collect()
            }

        } else {

            multiqc_ch = BISMARK.out.report.ifEmpty([])

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
                            BISMARK.out.report,
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
                            BISMARK.out.report,
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
        if (scNMT){

            multiqc_ch = multiqc_ch.mix(
                         COVERAGE2CYTOSINE.out.report.ifEmpty([])
                         ).collect()

            bismark_report_ch = bismark_report_ch.mix(
                         COVERAGE2CYTOSINE.out.report.ifEmpty([])
                         ).collect()             

        }

        /* ========================================================================================
        ======================================================================================== */

        BISMARK2REPORT      (bismark_report_ch, outdir, params.bismark2report_args, params.verbose)
        BISMARK2SUMMARY     (bismark_report_ch, outdir, params.bismark2summary_args, params.verbose)
        MULTIQC             (multiqc_ch, outdir, params.multiqc_args, params.verbose)
}

workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Jobname     : ${workflow.runName}
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
    .stripIndent()

    sendMail(to: "${workflow.userName}@ethz.ch", subject: 'Minimal pipeline execution report', body: msg)
}
