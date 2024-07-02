#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    PROCESSES
======================================================================================== */
process BISMARK2REPORT {

	container 'docker://josousa/bismark:0.24.2'

	input:
		file(file)
		val(outputdir)
		val(bismark2report_args)

	output:
		path "*html", emit: html
		
		publishDir "$outputdir/aligned/logs", mode: "link", overwrite: true

	script:
		"""
		bismark2report
		"""
}
