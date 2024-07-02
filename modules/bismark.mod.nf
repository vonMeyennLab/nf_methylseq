#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.bam_output = true // Setting if the bam file should be published


/* ========================================================================================
    PROCESSES
======================================================================================== */
process BISMARK {

	label 'bismark'
	tag "$name" // Adds name to job submission
	
	container 'docker://josousa/bismark:0.24.2'

    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (bismark_args)

	output:
	    tuple val(name), path("*bam"), emit: bam
		path "*report.txt",  	       emit: report

		publishDir "$outputdir/aligned/bam", 	 mode: "link", overwrite: true, pattern: "*bam", enabled: params.bam_output
		publishDir "$outputdir/aligned/logs",    mode: "link", overwrite: true, pattern: "*report.txt"
		publishDir "$outputdir/unaligned/fastq", mode: "link", overwrite: true, pattern: "*.fq.gz"

    script:
		/* ==========
			Paired-end
		========== */
		readString = ""
		if (reads instanceof List) {
			readString = "-1 " + reads[0] + " -2 " + reads[1]
		}
		else {
			readString = reads
		}

		/* ==========
			Index
		========== */
		index = "--genome " + params.genome["bismark"]


		/* ==========
			Basename
		========== */
		// HISAT2
		if (bismark_args =~ /-hisat/){ 
			bismark_name = name + "_" + params.genome["name"] + "_bismark_ht2"
		}
		// Bowtie 2
		else {
			bismark_name = name + "_" + params.genome["name"] + "_bismark_bt2"
		}

		"""
		bismark --parallel 1 -p ${task.cpus} --basename ${bismark_name} ${index} ${bismark_args} ${readString}
		"""
}
