nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.verbose = true

params.dual    = false


/* ========================================================================================
    PROCESSES
======================================================================================== */
process UMIBAM {
    
	tag "$bam" // Adds name to job submission

	input:
	    tuple val(name), path(bam)
		val (outputdir)
		val (umibam_args)
		val (verbose)

	output:
		path "*report.txt", 			emit: report
		tuple val(name), path ("*bam"), emit: bam

		publishDir "$outputdir/aligned/logs", mode: "link", overwrite: true, pattern: "*report.txt"
		publishDir "$outputdir/aligned/bam",  mode: "link", overwrite: true, pattern: "*bam"

    script:

		/* ==========
			Verbose
		========== */
		if (verbose){
			println ("[MODULE] UMIBAM ARGS: " + umibam_args)
		}
		

		/* ==========
			Dual
		========== */	
		if (params.dual){
			umibam_args += " --double_umi "	
		}
		else {
			umibam_args += " --umi "	
		}
		

		"""
		module load umibam
		
		umibam $umibam_args $bam
		
		rename UMI_d d *
		"""
}