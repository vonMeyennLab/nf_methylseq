#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.verbose = true


/* ========================================================================================
    PROCESSES
======================================================================================== */
process MACS_CALLPEAK {	

	input:
		path(files)
		val(outputdir)
		val(macs_callpeak_args)
		val(verbose)

	output:
		path "*_peaks.xls",        emit: peaks_xls,				optional: true
		path "*_peaks.narrowPeak", emit: peaks_narrowPeak,		optional: true
		path "*_summits.bed",      emit: summits_bed,			optional: true
		path "*_model.r",          emit: model_r,				optional: true
		path "*_peaks.gappedPeak", emit: peaks_gappedPeak,		optional: true
		path "*_peaks.broadPeak",  emit: peaks_broadPeakPeak,	optional: true
		path "*_pileup.bdg",       emit: pileup_bdg,			optional: true
		path "*_lambda.bdg",       emit: lambda_bdg,			optional: true

		publishDir "$outputdir/aligned/macs", mode: "link", overwrite: true

    script:

		/* ==========
			Verbose
		========== */
		if (verbose){
			println ("[MODULE] MACS CALLPEAK ARGS: " + macs_callpeak_args)
		}


		/* ==========
			File names and output name
		========== */
		if (files instanceof List) {
			files_command = "-t " + files[0] + " -c " + files[1]
			output_name   = files[0].toString()[0..<files[0].toString().lastIndexOf('.')]
		} else {
			files_command = "-t " + files
			output_name   = files.toString()[0..<files.toString().lastIndexOf('.')]
		}


		/* ==========
			Effective genome size
		========== */
        if (params.genome["name"] == 'GRCh37' || params.genome["name"] == 'GRCh38') {
			gsize = 'hs'
		}
        if (params.genome["name"] == 'GRCm38' || params.genome["name"] == 'GRCm39') {
			gsize = 'mm'
		}
        if (params.genome["name"] == 'WBcel235') {
			gsize = 'ce'
		}
        if (params.genome["name"] == 'BDGP6') {
			gsize = 'dm'
		}

		"""
		module load macs

		macs3 callpeak ${files_command} -g ${gsize} -n ${output_name} $macs_callpeak_args
    	"""
}