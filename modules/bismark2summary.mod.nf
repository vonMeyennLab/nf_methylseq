#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    PROCESSES
======================================================================================== */
process BISMARK2SUMMARY {

	input:
		file(file)
		val(outputdir)
		val(bismark2summary_args)

	output:
		path "*html", emit: html
		path "*txt",  emit: report

		publishDir "$outputdir/aligned/logs", mode: "link", overwrite: true

	script:
		/* ==========	
		We need to replace single quotes in the arguments so that they are not getting passed in as a single string
		This is only a temporary workaround until Paolo has fixed the Nextflow bug.
		https://github.com/nextflow-io/nextflow/issues/1519
		========== */
		bismark2summary_args = bismark2summary_args.replaceAll(/'/,"")
		
		"""
		export MODULEPATH=/cluster/work/nme/software/modules:$MODULEPATH

		module load bismark

		bismark2summary
		"""

}
