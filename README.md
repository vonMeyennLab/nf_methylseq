# Bisulfite Sequencing Pipeline

<img width="40%" src="https://raw.githubusercontent.com/nextflow-io/trademark/master/nextflow2014_no-bg.png" /></br>

A Nextflow pipeline to align and quantify Methylation (Bisulfite) sequencing data.


### Pipeline steps
1. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
2. [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
3. [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
4. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
5. [Bismark](https://felixkrueger.github.io/Bismark/)
6. [Bismark filter non-conversion](https://felixkrueger.github.io/Bismark/bismark/filter_nonconverted_reads/) _[Optional]_
7. [Bismark deduplication](https://felixkrueger.github.io/Bismark/bismark/deduplication/)
8. [Bismark methylation extractor](https://felixkrueger.github.io/Bismark/bismark/methylation_extraction/)
9. [Bismark2report](https://felixkrueger.github.io/Bismark/bismark/processing_report/)
10. [Bismark2summary](https://felixkrueger.github.io/Bismark/bismark/summary_report/)
11. [MultiQC](https://multiqc.info/)

### Required parameters
`--input` _Path to the folder where the FASTQ files are located. For example: /cluster/work/nme/data/josousa/project/fastq/*fastq.gz_</br>

`--outdir` _output directory where the files will be saved._


### Optional parameters
- Option to provide a custom FastQ Screen config file.</br>
`--fastq_screen_conf` _[default: /cluster/work/nme/software/config/fastq_screen.conf]_</br>

- Option to pass the option --bisulfite to FastQ Screen. This will result with FastQ Screen having the arguments: fastq_screen_args='--bisulfite'.</br>
`--bisulfite`</br>

- Option to define sequencing data as **single-cell**. This will result with Trim Galore having the arguments: trim_galore_args='--clip_R1 6 --clip_R2 6'.</br>
`--singlecell`</br>

- Option to define sequencing data as **PBAT**. This will result with Trim Galore having the arguments: trim_galore_args='--clip_R1 9 --clip_R2 9'.</br>
`--pbat`</br>
    > If this option is selected, the optional parameter `--bisulfite` will automatically be set to `true`.

- Option to define sequencing data as **RRBS**. This will result with Trim Galore having the arguments: trim_galore_args='--rrbs'.</br>
`--rrbs`</br>


### Skipping options
- Option to skip FastQ Screen. </br>
`--skip_fastq_screen`</br>

- Option to skip Trim Galore. </br>
`--skip_trim_galore`


### Extra arguments
- Option to add extra arguments to the package [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</br>
`--fastqc_args`</br>

- Option to add extra arguments to the package [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/).</br>
`--fastq_screen_args`</br>

- Option to add extra arguments to the package [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).</br>
`--trim_galore_args`

- Option to add extra arguments to the package [MultiQC](https://multiqc.info/).</br>
`--multiqc_args`

