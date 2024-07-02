#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.ignore_r2 = false


/* ========================================================================================
    PROCESSES
======================================================================================== */
process BISMARK_METHYLATION_EXTRACTOR {
	
	label 'BismarkMethylationExtractor'
	tag "$bam" // Adds name to job submission

	container 'docker://josousa/bismark:0.24.2'

    input:
	    tuple val(name), path(bam)
		val (outputdir)
		val (bismark_methylation_extractor_args)

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
			Paired-end
		========== */
		isPE = isPairedEnd(bam)
		
		if (isPE && params.ignore_r2){
			bismark_methylation_extractor_args += " --ignore_r2 ${params.ignore_r2} "
		}
		
		"""
		bismark_methylation_extractor --gzip --bedGraph --buffer 10G --parallel ${task.cpus} ${bismark_methylation_extractor_args} ${bam}
		"""

}


/* ========================================================================================
    FUNCTIONS
======================================================================================== */
def isPairedEnd(bamfile) {

	bamfile = bamfile.toString()
	if (bamfile =~ /_pe/){
		return true
	}
	else{
		return false
	}
	
}