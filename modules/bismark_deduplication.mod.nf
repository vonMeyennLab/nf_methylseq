#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
	DEFAULT PARAMETERS
======================================================================================== */
params.verbose    = true
params.bam_output = true // Setting if the bam file should be published


/* ========================================================================================
	PROCESSES
======================================================================================== */
process BISMARK_DEDUPLICATION {
	
	label 'bismarkDeduplication'
	tag "$bam" // Adds name to job submission
  	
    input:
	    tuple val(name), path(bam)
		val (outputdir)
		val (deduplicate_bismark_args)
		val (verbose)

	output:
		path "*report.txt",             emit: report
		tuple val(name), path ("*bam"), emit: bam

		publishDir "$outputdir/aligned/logs",              mode: "link", overwrite: true, pattern: "*report.txt"
		publishDir "$outputdir/aligned/bam/deduplicated",  mode: "link", overwrite: true, pattern: "*bam", enabled: params.bam_output

    script:
	
		/* ==========
			Verbose
		========== */
		if (verbose){
			println ("[MODULE] BISMARK DEDUPLICATION ARGS: " + deduplicate_bismark_args)
		}

		/* ==========
			Arguments
		========== */
		deduplicate_bismark_args += " --bam "

		"""
		module load bismark

		deduplicate_bismark ${deduplicate_bismark_args} ${bam}
		"""
}
