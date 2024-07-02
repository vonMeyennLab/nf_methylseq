#!/usr/bin/env nextflow
nextflow.enable.dsl=2


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
			Basename
		========== */
		// Removing the file extension from the input file name
		// (https://www.nextflow.io/docs/latest/script.html#removing-part-of-a-string)
		outfile_basename = coverage_file.toString()  // Important to convert nextflow.processor.TaskPath object to String first
		outfile_basename = (outfile_basename - ~/.bismark.cov.gz$/)
		outfile_basename = (outfile_basename - ~/.cov.gz$/)
		outfile_basename = (outfile_basename - ~/.cov$/)

		"""
		module load bismark

		coverage2cytosine --gzip --genome ${genome} ${coverage2cytosine_args} --output ${outfile_basename} ${coverage_file}
		"""
}
