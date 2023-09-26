#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.verbose  = true
params.metadata = false


/* ========================================================================================
    PROCESSES
======================================================================================== */
process GEOFETCH {

	tag "$geoaccession" // Adds name to job submission

	input:
        val(geoaccession)
		val(outputdir)
		val(geofetch_args)
		val(verbose)

	output:
	  	path "*/*soft",     emit: soft
        path "*/*_SRA.csv", emit: sra_metadata
        path "*/*/*_PEP*",  emit: pep
		publishDir "$outputdir/metadata", mode: "link", overwrite: true

	script:

		/* ==========
			Verbose
		========== */
		if (verbose){
			println ("[MODULE] GEOFETCH ARGS: " + geofetch_args)
		}

        if (params.metadata){
            geofetch_args += " --just-metadata "
		}

		"""
		module load geofetch

		geofetch -i ${geoaccession} ${geofetch_args}
		"""
}
