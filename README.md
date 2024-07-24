# Bisulfite Sequencing Pipeline

<img width="30%" src="https://raw.githubusercontent.com/nextflow-io/trademark/master/nextflow-logo-bg-light.png" />
<img width="30%" src="https://tower.nf/assets/nf-tower-black.svg" />

A Nextflow pipeline to align and quantify Methylation (Bisulfite) sequencing data.

>The pipeline was created to run on the [ETH Euler cluster](https://scicomp.ethz.ch/wiki/Euler) and it relies on the server's genome files. Thus, the pipeline needs to be adapted before running it in a different HPC cluster.

## Pipeline steps
1. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
2. [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
3. [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
4. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
5. [Bismark](https://felixkrueger.github.io/Bismark/)
6. [Bismark filter non-conversion](https://felixkrueger.github.io/Bismark/bismark/filter_nonconverted_reads/) _[Optional]_
7. [Bismark deduplication](https://felixkrueger.github.io/Bismark/bismark/deduplication/)
8. [Bismark methylation extractor](https://felixkrueger.github.io/Bismark/bismark/methylation_extraction/)
9. [coverage2cytosine](https://felixkrueger.github.io/Bismark/bismark/methylation_extraction/) _[Optional]_
10. [Bismark2report](https://felixkrueger.github.io/Bismark/bismark/processing_report/)
11. [Bismark2summary](https://felixkrueger.github.io/Bismark/bismark/summary_report/)
12. [MultiQC](https://multiqc.info/)

## Required parameters
Path to the folder where the FASTQ files are located. 
```bash
--input /cluster/work/nme/data/josousa/project/fastq/*fastq.gz
```

Output directory where the files will be saved.
```bash
--outdir /cluster/work/nme/data/josousa/project
```

## Genomes
- Reference genome used for alignment.

    `--genome`

    Available genomes:
    ``` bash
        GRCm39 # Default
        GRCm38
        GRCh38
        GRCh37 
        panTro6
        CHIMP2.1.4
        BDGP6
        susScr11
        Rnor_6.0
        R64-1-1
        TAIR10
        WBcel235
        E_coli_K_12_DH10B
        E_coli_K_12_MG1655
        Vectors
        Lambda
        PhiX
        Mitochondria
    ```

- Option to use a custom genome for alignment by providing an absolute path to a custom genome file.

    ``` bash
    --custom_genome_file '/cluster/work/nme/data/josousa/project/genome/CHM13.genome'
    ```

    Example of a genome file:
    ``` bash
    name           GRCm39
    species        Mouse
    bismark        /cluster/work/nme/genomes/Mus_musculus/Ensembl/GRCm39/Sequence/BismarkIndex/          
    ```

## FastQ Screen optional parameters
- Option to provide a custom FastQ Screen config file.
    ``` bash
    --fastq_screen_conf '/cluster/work/nme/software/config/fastq_screen.conf' # Default
    ```


## Bismark optional parameters
- Option to set the alignment mode to _local_.

    `--local`
    > In this mode, it is not required that the entire read aligns from one end to the other. Rather, some characters may be omitted (“soft-clipped”) from the ends in order to achieve the greatest possible alignment score.

- Option to write all reads that could not be aligned to a file in the output directory.
    
    `--unmapped`

- Option to write all reads which produce more than one valid alignment with the same number of lowest mismatches or other reads that fail to align uniquely to a file in the output directory.

    `--ambiguous`


## Skipping and adding options
- Option to skip FastQC, TrimGalore, and FastQ Screen. The first step of the pipeline will be the Bismark alignment. 
`--skip_qc`

- Option to skip FastQ Screen. 
`--skip_fastq_screen`

- Option to skip Bismark deduplication. 
`--skip_deduplication`

- Option to add Bismark filter non-conversion before deduplication (if selected) and before Bismark methylation extractor. 
`--add_filter_non_conversion`


## Extra arguments
- Option to add extra arguments to [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
`--fastqc_args`

- Option to add extra arguments to [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/).
`--fastq_screen_args`

- Option to add extra arguments to [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).
`--trim_galore_args`

- Option to add extra arguments to [Bismark](https://felixkrueger.github.io/Bismark/).
`--bismark_arg`

- Option to add extra arguments to [Bismark filter non-conversion](https://felixkrueger.github.io/Bismark/bismark/filter_nonconverted_reads/).
`--filter_non_conversion_args`

- Option to add extra arguments to [Bismark deduplication](https://felixkrueger.github.io/Bismark/bismark/deduplication/).
`--deduplicate_bismark_args`

- Option to add extra arguments to [Bismark methylation extractor](https://felixkrueger.github.io/Bismark/bismark/methylation_extraction/).
`--bismark_methylation_extractor_args`

- Option to add extra arguments to [Bismark coverage2cytosine](https://felixkrueger.github.io/Bismark/bismark/methylation_extraction/).
`--coverage2cytosine_args`

- Option to add extra arguments to [Bismark2summary](https://felixkrueger.github.io/Bismark/bismark/summary_report/).
`--bismark2summary_args`

- Option to add extra arguments to [Bismark2report](https://felixkrueger.github.io/Bismark/bismark/processing_report/).
`--bismark2report_args`

- Option to add extra arguments to [MultiQC](https://multiqc.info/).
`--multiqc_args`

## Acknowledgements
This pipeline was adapted from the Nextflow pipelines created by the [Babraham Institute Bioinformatics Group](https://github.com/s-andrews/nextflow_pipelines) and from the [nf-core](https://nf-co.re/) pipelines. We thank all the contributors for both projects. We also thank the [Nextflow community](https://nextflow.slack.com/join) and the [nf-core community](https://nf-co.re/join) for all the help and support.
