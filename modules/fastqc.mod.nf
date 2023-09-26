#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.verbose = true


/* ========================================================================================
    PROCESSES
======================================================================================== */
process FASTQC {

	tag "$name" // Adds name to job submission

	input:
	  	tuple val(name), path(reads)
		val(outputdir)
		val(fastqc_args)
		val(verbose)

	output:
	  	tuple val(name), path ("*fastqc*"), emit: all
		path "*.zip", 						emit: report
		
		publishDir "$outputdir/qc/fastqc", mode: "link", overwrite: true

	script:

		/* ==========
			Verbose
		========== */
		if (verbose){
			println ("[MODULE] FASTQC ARGS: " + fastqc_args)
		}

		"""
		module load fastqc

		fastqc $fastqc_args -q -t 2 ${reads}
		"""
}
