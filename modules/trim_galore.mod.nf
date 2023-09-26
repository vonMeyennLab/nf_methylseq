#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.verbose    = true

params.singlecell = false
params.rrbs		  = false
params.pbat		  = false


/* ========================================================================================
    PROCESSES
======================================================================================== */
process TRIM_GALORE {

	label 'trimGalore'
	tag "$name" // Adds name to job submission

	input:
		tuple val(name), path(reads)
		val(outputdir)
		val(trim_galore_args)
		val(verbose)

	output:
		tuple val(name), path("*fq.gz"), 			 emit: reads
		path "*trimming_report.txt", optional: true, emit: report
		
		publishDir "$outputdir/unaligned/fastq", mode: "link", overwrite: true, pattern: "*fq.gz"
		publishDir "$outputdir/unaligned/logs",  mode: "link", overwrite: true, pattern: "*trimming_report.txt"

	script:

		/* ==========
			Verbose
		========== */	
		if (verbose){
			println ("[MODULE] TRIM GALORE ARGS: " + trim_galore_args)
		}


		/* ==========
			Paired-end
		========== */
		pairedString = ""
		if (reads instanceof List) {
			pairedString = "--paired"
		}


		/* ==========
			Single-cell
		========== */
		if (params.singlecell == '6'){
			trim_galore_args = trim_galore_args + " --clip_r1 $params.singlecell "
			if (pairedString == "--paired"){
				trim_galore_args = trim_galore_args + " --clip_r2 $params.singlecell "
			}
		}

		/* ==========
			RRBS
		========== */
		if (params.rrbs){
			trim_galore_args = trim_galore_args + " --rrbs "
		}


		/* ==========
			PBAT
		========== */
		if (params.pbat == '9'){
			trim_galore_args = trim_galore_args + " --clip_r1 $params.pbat "
			if (pairedString == "--paired"){
				trim_galore_args = trim_galore_args + " --clip_r2 $params.pbat "
			}
		}


		"""
		module load trimgalore

		trim_galore --cores ${task.cpus} $trim_galore_args ${pairedString} ${reads}
		"""
}
