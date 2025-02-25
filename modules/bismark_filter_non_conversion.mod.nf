#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    PROCESSES
======================================================================================== */
process BISMARK_FILTER_NON_CONVERSION {
	
	label 'bismarkFilterNonConversion'
	tag "$bam" // Adds name to job submission

	container 'docker://josousa/bismark:0.24.2'

    input:
		tuple val(name), path(bam)
		val (outputdir)
		val (filter_non_conversion_args)

	output:
		path "*.txt",                                      emit: report
		tuple val(name), path ("*nonCG_filtered.bam"),     emit: bam
		tuple val(name), path ("*nonCG_removed_seqs.bam"), emit: nonCG_removed_seqs_bam

		publishDir "$outputdir/aligned/logs", mode: "link", overwrite: true, pattern: "*report.txt"
		publishDir "$outputdir/aligned/bam",  mode: "link", overwrite: true, pattern: "*nonCG_filtered.bam"
		publishDir "$outputdir/aligned/bam",  mode: "link", overwrite: true, pattern: "*nonCG_removed_seqs.bam"

    script:
		"""
		filter_non_conversion ${filter_non_conversion_args} ${bam}
		"""
}
