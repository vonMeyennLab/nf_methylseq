#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    PROCESSES
======================================================================================== */
process BISMARK2REPORT {

	input:
		file(file)
		val(outputdir)
		val(bismark2report_args)

	output:
		path "*html", emit: html
		
		publishDir "$outputdir/aligned/logs", mode: "link", overwrite: true

	script:
		"""
		export MODULEPATH=/cluster/work/nme/software/modules:$MODULEPATH

		module load bismark

		bismark2report
		"""
}
