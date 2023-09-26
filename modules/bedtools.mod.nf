#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
	DEFAULT PARAMETERS
======================================================================================== */
params.verbose = true
params.cutntag = false


/* ========================================================================================
	PROCESSES
======================================================================================== */
process BEDTOOLS_GENOMECOV {

	tag "$bam" // Adds name to job submission

	input:
		path(bam)
		val(outputdir)
		val(bedtools_genomecov_args)
		val(verbose)

	output:
		path("*bedgraph"), emit: bedgraph
		publishDir "$outputdir/aligned/bedgraph", mode: "link", overwrite: true

    script:

		/* ==========
			Verbose
		========== */
		if (verbose){
			println ("[MODULE] BEDTOOLS GENOMECOV ARGS: " + bedtools_genomecov_args)
		}

		/* ==========
			Parameters for the CUT&Tag pipeline
		========== */
        if (params.cutntag) {
			bedtools_genomecov_args += " -bga " 
			/* ==========
			Report depth in BedGraph format, as above (-bg).
			However with this option, regions with zero
			coverage are also reported. This allows one to
			quickly extract all regions of a genome with 0
			coverage by applying: "grep -w 0$" to the output.
			========== */
		}

		"""
		module load bedtools2

		bedtools genomecov ${bedtools_genomecov_args} -ibam ${bam} > ${bam}.bedgraph

		rename .bam.bedgraph .bedgraph *
    	"""
}
