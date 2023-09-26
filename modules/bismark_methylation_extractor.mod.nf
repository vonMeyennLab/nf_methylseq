#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.verbose    = true

params.singlecell = false
params.rrbs       = false
params.pbat       = false
params.nonCG      = false


/* ========================================================================================
    PROCESSES
======================================================================================== */
process BISMARK_METHYLATION_EXTRACTOR {
	
	label 'BismarkMethylationExtractor'
	tag "$bam" // Adds name to job submission

    input:
	    tuple val(name), path(bam)
		val (outputdir)
		val (bismark_methylation_extractor_args)
		val (verbose)

	output:
	    tuple val(name), path("CpG*.txt.gz"),   emit: context_files_CG
		path "CH*.txt.gz",                      emit: context_files_nonCG
		path "*report.txt",                     emit: report
		path "*M-bias.txt",                     emit: mbias
		path "*cov.gz",                         emit: coverage
		
		publishDir "$outputdir/aligned/methylation_calls",    mode: "link", overwrite: true, pattern: "CpG*.txt.gz"
		publishDir "$outputdir/aligned/methylation_calls",    mode: "link", overwrite: true, pattern: "CH*.txt.gz"
		publishDir "$outputdir/aligned/methylation_coverage", mode: "link", overwrite: true, pattern: "*cov.gz"
		publishDir "$outputdir/aligned/logs", 				  mode: "link", overwrite: true, pattern: "*.txt"
    
	script:

		/* ==========
			Verbose
		========== */
		if (verbose){
			println ("[MODULE] BISMARK METHYLATION EXTRACTOR ARGS: " + bismark_methylation_extractor_args)
		}


		/* ==========
			Arguments
		========== */
		bismark_methylation_extractor_args += " --gzip "
		

		/* ==========
			Single-cell
		========== */
		if (params.singlecell){
			println ("FLAG SINGLE CELL SPECIFIED: PROCESSING ACCORDINGLY")
		}


		/* ==========
			Non-CpG methylation
		========== */
		if (params.nonCG){
			if (verbose){
				println ("FLAG nonCG specified: adding flag --CX ")
			}
			bismark_methylation_extractor_args +=  " --CX "
		}


		/* ==========
			Paired-end
		========== */
		isPE = isPairedEnd(bam)
		if (isPE){
			// not perform any ignoring behaviour for RRBS or single-cell libraries
			if (!params.rrbs && !params.singlecell && !params.pbat){
				bismark_methylation_extractor_args +=  " --ignore_r2 2 "
			}
		}
		else {
			println("File seems to be single-end")
		}

		
		"""
		module load bismark

		bismark_methylation_extractor --bedGraph --buffer 10G -parallel ${task.cpus} ${bismark_methylation_extractor_args} ${bam}
		"""

}


/* ========================================================================================
    FUNCTIONS
======================================================================================== */
def isPairedEnd(bamfile) {

	// Transforms the nextflow.processor.TaskPath object to String
	bamfile = bamfile.toString()
	if (params.verbose){
		println ("Processing file: " + bamfile)
	}
	
	if (bamfile =~ /_pe/){
		if (params.verbose){
			println ("File is paired-end")
		}
		return true
	}
	else{
	 	if (params.verbose){
			 println ("File is single-end")
		 }
		return false
	}
}