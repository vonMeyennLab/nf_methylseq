#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.verbose  = true


/* ========================================================================================
    PROCESSES
======================================================================================== */
process SRADOWNLOADER {

	tag "$geoaccession" // Adds name to job submission
    label 'sradownloader'

	input:
	  	file(sra_metadata)
        val(geoaccession)
		val(outputdir)
		val(sradownloader_args)
		val(verbose)

	output:
	  	path "*gz", emit: fastq
		publishDir "$outputdir/fastq", mode: "link", overwrite: true

	script:

		/* ==========
			Verbose
		========== */
		if (verbose){
			println ("[MODULE] SRADOWNLOADER ARGS: " + sradownloader_args)
		}

		"""
		module load sradownloader

		sradownloader ${sradownloader_args} --noena ${sra_metadata} 
		"""
}
