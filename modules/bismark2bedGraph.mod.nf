#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.verbose	   = true

params.dirty_harry = false


/* ========================================================================================
    PROCESSES
======================================================================================== */
process BISMARK2BEDGRAPH {
	
	label 'bismark2bedGraph'
	tag "$name" // Adds name to job submission
			
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (bismark2bedGraph_args)
		val (verbose)

	output:
	    path "*cov.gz",        emit: coverage
		path "*bedGraph.gz",   emit: bedGraph

		publishDir "$outputdir/aligned/methylation_coverage", mode: "link", overwrite: true, pattern: "*cov.gz"
		publishDir "$outputdir/aligned/methylation_bedgraph", mode: "link", overwrite: true, pattern: "*bedGraph.gz"

    script:

		/* ==========
			Verbose
		========== */
		if (verbose){
			println ("[MODULE] BISMARK2BEDGRAPH ARGS: " + bismark2bedGraph_args)
		}


		/* ==========
			Output name
		========== */
		if (params.dirty_harry){
			output_name = name + "_DH.bedGraph.gz"
		} 
		else {
			output_name = name + ".bedGraph.gz"
		}


		"""
		module load bismark

		bismark2bedGraph --buffer 15G -o $output_name $bismark2bedGraph_args $reads
		"""
}
