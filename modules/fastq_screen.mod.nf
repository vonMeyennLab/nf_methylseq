#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.verbose           = true
params.fastq_screen_conf = "/cluster/work/nme/software/config/fastq_screen.conf"
params.bisulfite         = false
params.single_end        = false


/* ========================================================================================
    PROCESSES
======================================================================================== */
process FASTQ_SCREEN {

	label 'fastqScreen'
	tag "$name" // Adds name to job submission
	
	input:
		tuple val(name), path(reads)
		val(outputdir)
		val(fastq_screen_args)
		val(verbose)

	output:
	  	path "*html", emit: html
	  	path "*txt",  emit: report
		  
		publishDir "$outputdir/qc/fastq_screen", mode: "link", overwrite: true

	script:

		/* ==========
			Verbose
		========== */		
		if (verbose){
			println ("[MODULE] FASTQ SCREEN ARGS: "+ fastq_screen_args)
		}


		/* ==========
			File names
		========== */	
		if (params.single_end){
			// TODO: Add single-end parameter
		}
		else{
			// for paired-end files we only use Read 1 (as Read 2 tends to show the exact same thing)
			if (reads instanceof List) {
				reads = reads[0]
			}
		}
		

		/* ==========
			Bisulfite
		========== */		
		if (params.bisulfite){
			fastq_screen_args += " --bisulfite "
		}


		"""
		module load fastq_screen

		fastq_screen --conf $params.fastq_screen_conf $fastq_screen_args $reads
		"""
}
