#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.verbose       = true
params.bam_output    = true // Setting if the bam file should be published

params.pbat 	     = false
params.unmapped      = false
params.ambiguous     = false
params.singlecell    = false
params.read_identity = ''


/* ========================================================================================
    PROCESSES
======================================================================================== */
process BISMARK {
	
	label 'bismark'
	tag "$name" // Adds name to job submission
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (bismark_args)
		val (verbose)

	output:
	    tuple val(name), path("*bam"), emit: bam
		path "*report.txt",  	       emit: report
		
		tuple val(name), path("*unmapped_reads_1.fq.gz"),  optional: true, emit: unmapped1
		tuple val(name), path("*unmapped_reads_2.fq.gz"),  optional: true, emit: unmapped2
        tuple val(name), path("*ambiguous_reads_1.fq.gz"), optional: true, emit: ambiguous1
		tuple val(name), path("*ambiguous_reads_2.fq.gz"), optional: true, emit: ambiguous2

		publishDir "$outputdir/aligned/bam", 	 mode: "link", overwrite: true, pattern: "*bam", enabled: params.bam_output
		publishDir "$outputdir/aligned/logs",    mode: "link", overwrite: true, pattern: "*report.txt"
		publishDir "$outputdir/unaligned/fastq", mode: "link", overwrite: true, pattern: "*.fq.gz"

    script:

		/* ==========
			Verbose
		========== */
		if (verbose){
			println ("[MODULE] BISMARK ARGS: " + bismark_args)
		}


		/* ==========
			Single-cell
		========== */
		if (params.singlecell){
			bismark_args += " --non_directional "
		}


		/* ==========
			PBAT
		========== */
		if (params.pbat){
			bismark_args += " --pbat "
		}


		/* ==========
			File names
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
			Unmapped
		========== */
		if (params.unmapped){
			bismark_args += " --unmapped "
		}
		

		/* ==========
			Ambiguous
		========== */
		if (params.ambiguous){
			bismark_args += " --ambiguous "
		}


		/* ==========
			Basename
		========== */
		unmapped_name  = ''
		ambiguous_name = ''
		if (params.read_identity == "1" || params.read_identity == "2"){


			// Unmapped and ambiguous reads
			if (params.read_identity == "1"){
				unmapped_name  = name + "_unmapped_R1"
				ambiguous_name = name + "_ambiguous_R1"
			}
			else {
				unmapped_name  = name + "_unmapped_R2"
				ambiguous_name = name + "_ambiguous_R2"
			}

			// HISAT2
			if (bismark_args =~ /-hisat/){
				bismark_name = unmapped_name  + "_" + params.genome["name"] + "_bismark_ht2"
				bismark_name = ambiguous_name + "_" + params.genome["name"] + "_bismark_ht2"
			}
			// Bowtie 2
			else {
				bismark_name = unmapped_name  + "_" + params.genome["name"] + "_bismark_bt2"
				bismark_name = ambiguous_name + "_" + params.genome["name"] + "_bismark_bt2"
			}


		}
		else {

			// HISAT2
			if (bismark_args =~ /-hisat/){ 
				bismark_name = name + "_" + params.genome["name"] + "_bismark_ht2"
			}
			// Bowtie 2
			else {
				bismark_name = name + "_" + params.genome["name"] + "_bismark_bt2"
			}
		}	
		

		"""
		module load bismark

		bismark --parallel 1 -p ${task.cpus} --basename $bismark_name $index $bismark_args $readString
		"""
}
