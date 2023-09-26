#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.verbose    = true

params.singlecell = false
params.rrbs       = false
params.verbose    = false
params.pbat       = false
params.nome       = false


/* ========================================================================================
    PROCESSES
======================================================================================== */
process COVERAGE2CYTOSINE {

	label 'coverage2Cytosine'
	tag "$coverage_file" // Adds name to job submission
	
	input:
		path(coverage_file)
		val(outputdir)
		val(coverage2cytosine_args)
		val(verbose)

	output:
		path "*{report.txt.gz,report.txt}", emit: report
		path "*{.cov.gz,.cov}",             emit: coverage

		publishDir "$outputdir/aligned/logs", 				  mode: "link", overwrite: true, pattern: "*report.txt*"
		publishDir "$outputdir/aligned/methylation_coverage", mode: "link", overwrite: true, pattern: "*.cov*"

	script:

		/* ==========
			Genome
		========== */
		genome = params.genome["bismark"]


		/* ==========
			Verbose
		========== */
		if (verbose){
			println ("[MODULE] BISMARK COVERAGE2CYTOSINE ARGS: " + coverage2cytosine_args)
			println ("Bismark Genome is: " + genome)
		}


		/* ==========
			Basename
		========== */
		// Removing the file extension from the input file name
		// (https://www.nextflow.io/docs/latest/script.html#removing-part-of-a-string)
		outfile_basename = coverage_file.toString()  // Important to convert nextflow.processor.TaskPath object to String first
		outfile_basename = (outfile_basename - ~/.bismark.cov.gz$/)
		outfile_basename = (outfile_basename - ~/.cov.gz$/)
		outfile_basename = (outfile_basename - ~/.cov$/)


		/* ==========
			Arguments
		========== */
		coverage2cytosine_args = coverage2cytosine_args + " --gzip "


		/* ==========
			NOMe-Seq
		========== */
		if (params.nome){
			if (verbose){
				println ("NOMe-seq outfile basename: $outfile_basename")
			}
			coverage2cytosine_args = coverage2cytosine_args + " --nome"
		}

		"""
		module load bismark

		coverage2cytosine --genome $genome $coverage2cytosine_args --output ${outfile_basename} $coverage_file
		"""
}
