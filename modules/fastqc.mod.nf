#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    PROCESSES
======================================================================================== */
process FASTQC {

	label 'fastqc'
	tag "$name" // Adds name to job submission

	container 'docker://staphb/fastqc:0.12.1'

	input:
	  	tuple val(name), path(reads)
		val(outputdir)
		val(fastqc_args)

	output:
	  	tuple val(name), path ("*fastqc*"), emit: all
		path "*.zip", 						emit: report
		
		publishDir "$outputdir/qc/fastqc", mode: "link", overwrite: true

	script:
		"""
		ls
		echo $SINGULARITY_BIND
		echo $GENOMES
		echo $WORK
		ls $GENOMES
		ls /cluster

		fastqc ${fastqc_args} --threads ${task.cpus} ${reads}
		"""
}
