#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.bam_output = true // Setting if the bam file should be published

/* ========================================================================================
    PROCESSES
======================================================================================== */
process BISMARK_DEDUPLICATION {
    
    label 'bismarkDeduplication'
    tag "$bam" // Adds name to job submission
    
    input:
        tuple val(name), path(bam)
        val(outputdir)
        val(deduplicate_bismark_args)
        
    output:
        path "*report.txt",             emit: report
        tuple val(name), path("*bam"), emit: bam

        publishDir "$outputdir/aligned/logs",              mode: "link", overwrite: true, pattern: "*report.txt"
        publishDir "$outputdir/aligned/bam/deduplicated",  mode: "link", overwrite: true, pattern: "*bam", enabled: params.bam_output

    script:
        """
        module load bismark

        deduplicate_bismark --bam ${deduplicate_bismark_args} ${bam}
        """
}