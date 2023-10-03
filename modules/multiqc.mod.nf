#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.verbose = true
params.prefix  = "" 


/* ========================================================================================
    PROCESSES
======================================================================================== */
process MULTIQC {

	label 'multiQC'

	input:
		path(file)
		val(outputdir)
		val(multiqc_args)
		val(verbose)

	output:
		path "*html", emit: html
		
		publishDir "$outputdir/qc", mode: "link", overwrite: true

	script:

		/* ==========
			Verbose
		========== */
		if (verbose){
			println ("[MODULE] MULTIQC ARGS: " + multiqc_args)
		}


		"""
		module load multiqc

		multiqc $multiqc_args -x work --filename ${params.prefix}multiqc_report.html .
		"""
}
