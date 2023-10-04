# Bisulfite Sequencing Pipeline

<img width="40%" src="https://raw.githubusercontent.com/nextflow-io/trademark/master/nextflow2014_no-bg.png" /></br>

A Nextflow pipeline to align and quantify Methylation (Bisulfite) sequencing data.


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
`--input` _Path to the folder where the FASTQ files are located. For example: /cluster/work/nme/data/josousa/project/fastq/*fastq.gz_</br>

`--outdir` _output directory where the files will be saved._

## Genomes
- Reference genome used for alignment.</br>
`--genome` _[Default: GRCm39]_
    > Available genomes:</br>
    GRCm39</br>
    GRCm38</br>
    GRCh38</br>
    GRCh37</br> 
    panTro6</br>
    CHIMP2.1.4</br>
    BDGP6</br>
    susScr11</br>
    Rnor_6.0</br>
    R64-1-1</br>
    TAIR10</br>
    WBcel235</br>
    E_coli_K_12_DH10B</br>
    E_coli_K_12_MG1655</br>
    Vectors</br>
    Lambda</br>
    PhiX</br>
    Mitochondria

- Option to use a custom genome for alignment by providing an absolute path to a custom genome file. </br>
`--custom_genome_file`

> |                |                                                                              |
> |----------------|------------------------------------------------------------------------------|
> | name           | GRCm39                                                                       |
> | species        | Mouse                                                                        |
> | fasta          | ${GENOMES}/Mus_musculus/Ensembl/GRCm39/Sequence/WholeGenomeFasta/            |
> | bismark        | ${GENOMES}/Mus_musculus/Ensembl/GRCm39/Sequence/BismarkIndex/                |
> | bowtie         | ${GENOMES}/Mus_musculus/Ensembl/GRCm39/Sequence/BowtieIndex/genome           |
> | bowtie2        | ${GENOMES}/Mus_musculus/Ensembl/GRCm39/Sequence/Bowtie2Index/genome          |
> | star           | ${GENOMES}/Mus_musculus/Ensembl/GRCm39/Sequence/STARIndex/genome             |
> | bwa            | ${GENOMES}/Mus_musculus/Ensembl/GRCm39/Sequence/BWAIndex/genome              |
> | hisat2         | ${GENOMES}/Mus_musculus/Ensembl/GRCm39/Sequence/Hisat2Index/genome           |
> | hisat2_splices | ${GENOMES}/Mus_musculus/Ensembl/GRCm39/Sequence/Hisat2Index/splice_sites.txt |
> | gtf            | ${GENOMES}/Mus_musculus/Ensembl/GRCm39/Annotation/Genes/genes.gtf            |



## FastQ Screen optional parameters
- Option to provide a custom FastQ Screen config file.</br>
`--fastq_screen_conf` _[default: /cluster/work/nme/software/config/fastq_screen.conf]_</br>


## Bismark optional parameters
- Option to set the minimum insert size for valid paired-end alignments.</br>
`--minins` _[Defaul: 0]_</br>

- Option to set the maximum insert size for valid paired-end alignments.</br>
`--maxins` _[Defaul: 500]_</br>

- Option to set the alignment mode to _local_.</br>
`--local`</br>
    > In this mode, it is not required that the entire read aligns from one end to the other. Rather, some characters may be omitted (“soft-clipped”) from the ends in order to achieve the greatest possible alignment score.

- Option to write all reads that could not be aligned to a file in the output directory.</br>
`--unmapped`</br>

- Option to write all reads which produce more than one valid alignment with the same number of lowest mismatches or other reads that fail to align uniquely to a file in the output directory.</br>
`--ambiguous`


## Skipping and adding options
- Option to skip FastQC, TrimGalore, and FastQ Screen. The first step of the pipeline will be the Bismark alignment. </br>
`--skip_qc`

- Option to skip FastQ Screen. </br>
`--skip_fastq_screen`</br>

- Option to skip Bismark deduplication. </br>
`--skip_deduplication`

- Option to add Bismark filter non-conversion before deduplication (if selected) and before Bismark methylation extractor. </br>
`--add_filter_non_conversion`


## Extra arguments
- Option to add extra arguments to the package [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) .</br>
`--fastqc_args`</br>

- Option to add extra arguments to the package [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/).</br>
`--fastq_screen_args`</br>

- Option to add extra arguments to the package [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).</br>
`--trim_galore_args`</br>

- Option to add extra arguments to the package [Bismark](https://felixkrueger.github.io/Bismark/).</br>
`--bismark_arg`</br>

- Option to add extra arguments to the funcion [Bismark filter non-conversion](https://felixkrueger.github.io/Bismark/bismark/filter_nonconverted_reads/).</br>
`--filter_non_conversion_args`</br>

- Option to add extra arguments to the funcion [Bismark deduplication](https://felixkrueger.github.io/Bismark/bismark/deduplication/).</br>
`--deduplicate_bismark_args`</br>

- Option to add extra arguments to the funcion [Bismark methylation extractor](https://felixkrueger.github.io/Bismark/bismark/methylation_extraction/).</br>
`--bismark_methylation_extractor_args`</br>

- Option to add extra arguments to the funcion [Bismark coverage2cytosine](https://felixkrueger.github.io/Bismark/bismark/methylation_extraction/).</br>
`--coverage2cytosine_args`</br>

- Option to add extra arguments to the funcion [Bismark2summary](https://felixkrueger.github.io/Bismark/bismark/summary_report/).</br>
`--bismark2summary_args`</br>

- Option to add extra arguments to the funcion [Bismark2report](https://felixkrueger.github.io/Bismark/bismark/processing_report/).</br>
`--bismark2report_args`</br>

- Option to add extra arguments to the package [MultiQC](https://multiqc.info/).</br>
`--multiqc_args`</br>