#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.verbose = true


/* ========================================================================================
    PROCESSES
======================================================================================== */
process BISMARK2REPORT {

	input:
		file(file)
		val(outputdir)
		val(bismark2report_args)
		val(verbose)

	output:
		path "*html", emit: html
		
		publishDir "$outputdir/aligned/logs", mode: "link", overwrite: true

	script:

		/* ==========
			Verbose
		========== */
		if (verbose){
			println ("[MODULE] BISMARK2REPORT ARGS: " + bismark2report_args)
		}


		"""
		module load bismark

		bismark2report
		"""
}
