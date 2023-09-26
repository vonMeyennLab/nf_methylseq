#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.verbose = true


/* ========================================================================================
    PROCESSES
======================================================================================== */
process FEATURECOUNTS {	

    label "featureCounts"
	tag "$bam" // Adds name to job submission

	input:
		path(bam)
        val(single_end)
		val(outputdir)
		val(featurecounts_args)
		val(verbose)

	output:
        path("*featureCounts.txt")        , emit: counts
        path("*featureCounts.txt.summary"), emit: summary
        
	 	publishDir "$outputdir/aligned/counts", mode: "link", overwrite: true

	script:

  		/* ==========
			Verbose
		========== */  
		if (verbose){
			println ("[MODULE] FEATURECOUNTS ARGS: " + featurecounts_args)
		}


		/* ==========
			Arguments
		========== */
            // -B  Only count read pairs that have both ends aligned.
            // -C  Do not count read pairs that have their two ends mapping
            //     to different chromosomes or mapping to same chromosome
            //     but on different strands.
        featurecounts_args = featurecounts_args + " -B -C "


		/* ==========
			Strandedness
		========== */
        if (params.strand == 'forward') {
            strandedness = 1
        } else if (params.strand == 'reverse') {
            strandedness = 2
        } else if (params.strand == 'unstranded' || params.strand == 'smartseq2') {
            strandedness = 0
        }


		/* ==========
			Paired-end or single-end
		========== */
        paired_end = single_end ? '' : '-p'


		/* ==========
			Basename
		========== */
        basename = bam.toString() - ".bam"


		/* ==========
			Annotation file
		========== */
        annotation = params.genome["gtf"]
        
        
		"""
		module load subread

		featureCounts \\
            $featurecounts_args \\
            $paired_end \\
            -T $task.cpus \\
            -a $annotation \\
            -s $strandedness \\
            -o ${basename}.featureCounts.txt \\
            ${bam}
		"""
}




process FEATURECOUNTS_MERGE_COUNTS {

    input:
        path('counts/*')
        val(outputdir)
        val(verbose)

    output:
        path "*.txt", emit: merged_counts

        publishDir "$outputdir/aligned/counts", mode: "link", overwrite: true

    script:
    """
    mkdir -p tmp/counts

    cut -f 1 `ls ./counts/* | head -n 1` | grep -v "^#" > ids.tsv

    for fileid in `ls ./counts/*featureCounts.txt`; do
        samplename=`basename \$fileid | sed s/\\.featureCounts.txt\$//g`
        echo \$samplename > tmp/counts/\$samplename.featureCounts.txt
        grep -v "^#" \${fileid} | cut -f 7 | tail -n+2 >> tmp/counts/\$samplename.featureCounts.txt
    done
    
    paste ids.tsv tmp/counts/* > gene_counts_merged.txt
    """
}
