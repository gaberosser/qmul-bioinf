# Marino Lab processing pipelines for Next Gen Seq data
_Author: Gabriel Rosser_
_Last updated: 09/03/2018_

## Preamble
All `commands` can be run at a Linux terminal, providing the appropriate software has been installed and - in the case of Apocrita - `module load` has been called.

I've created some Python classes to automate running most of pipelines locally or on Apocrita. On Apocrita, jobs are always submitted in an array. The Python scripts are purely a convenience: it is always possible to generate a script and associated parameters file.

I organise files as follows. 

- Raw data (fastq files) are located in `<data_dir>/<experiment_type>/<project_id>/<lane_id>`.
- Where processing is lane-specific, files are located in the lane subdirectory, e.g. `<data_dir>/<experiment_type>/<project_id>/<lane_id>/trim_galore/`.
- Where processing combines multiple lanes, the files are located in a subdirectory within the project directory, e.g. `<data_dir>/<experiment_type>/<project_id>/human/star_alignment/`.
- Where processing is reference-specific, files are located in a subdirectory naming the reference, e.g. `<data_dir>/<experiment_type>/<project_id>/<lane_id>/mouse/bwa_alignment` (lane-specific) or `<data_dir>/<experiment_type>/<project_id>/mouse/bwa_alignment` (combined lanes).

## Downloading data
The Wellcome Trust Centre for Human Genetics (WTCHG) always serves data over ftp. The ftp server requires a username and password, which are supplied in the URL they send. To download the entire directory from a URL, run 

```bash
wget -nH -np -r ftp://<username>:<password>@bsg-ftp.well.ox.ac.uk/<run_id>
```

The entire URL can just be copied from the email. This will create a directory called `run_id`.

## FastQC

Always run `fastqc` on the raw fastq files first. If desirable, the reports can be aggregated into a single summary HTML file using `multiqc`.

On Apocrita: `module load fastqc`. MultiQC should be installed in a Python virtual environment.

```bash
fastqc --noextract -t <num_threads> -o <output_dir> <file.fastq.gz>
```
To automate on Apocrita, the following will run `fastqc` on  all fastq files in the current directory: 
```bash
python /path/to/pyscript/apocrita/fastqc.py
```

## Small RNA-Seq
### Typical experimental parameters
- 51 bp reads
- Single end
- 30mi reads per sample
- Split over a large number of lanes (6 - 8)?
- WTCHG usually run their own trim - align - feature count pipeline, which is in a subdirectory `miRNA`. I don't trust it, but keep it for reference.
### Process
This has not been used enough to be fully automated. First load modules if on Apocrita:
```bash
module load trimgalore
module load bwa
module load samtools
```

sRNA-Seq reads always need trimming to remove adapters:
```bash
trimgalore -o <output_dir> file.fastq.gz
```

Now we can align using an **ungapped** aligner like `bwa`. WE should be fairly strict about disallowing gaps and mismatches.

Before we can run `bwa`, we need to build an indexed reference genome. Start with a reference file (in fasta format) - I have been using Ensembl GrCh38 for all human work:
```bash
bwa index /path/to/reference.fa
```
This dumps some extra files in the same directory as the `.fa` file. Now we can use this to align the raw sequencing data:
```bash
bwa aln -n 1 -o 0 -e 0 -k 1 -t 2 /path/to/reference.fa reads.fastq.gz > reads.sai
bwa samse /path/to/reference.fa reads.fastq.gz reads.sai reads.fastq.gz | samtools view -hb - > reads.bam
samtools sort -@ <num_threads> reads.bam > reads.sorted.bam
rm reads.sai reads.bam
```

This is implemented in a Python script:
```bash
python /path/to/pyscript/apocrita/bwa_se.py -i /path/to/reference.fa -p <num_threads> --srna -o <output_dir>
```
Running this script from within a directory will submit an alignment job for every fastq and fastq.gz file in that directory.

At present, I run this separately for each lane. This means we end up with `n` reads.sorted.bam files (where `n` is the number of lanes) for every sample. We need to merge these at the end. For WTCHG,  the common ID linking each sample is the final number (e.g. `WTCHG_xxxxxx_101.bam`, where `xxxxxx` is specific to the lane and 101 is the sample ID). I'm assuming that you have a subdirectory `/human/bwa_alignment` in the project directory, a subdirectory `trim_galore/human/bwa_alignment` in each of the lane subdirectories and that the lane IDs all start with the number 1:

```bash
arr=(101 102 ... 201) # all of the sample IDs
for i in "${arr[@]}"; do
	echo $i
	samtools merge human/bwa_alignment/$i.sorted.bam 1*/trim_galore/human/bwa_alignment/*_$i.sorted.bam
done
```

This generates sorted bam files (by coordinate), as the output name suggests.

Since BWA is light on logging, we might want to generate our own reports about the alignment success:

```bash
samtools flagstat <alignment.bam>
```

Now we need to run `featureCounts` using a `miRBase` annotation file to count against known miRNA transcripts. The annotation file can be obtained in `.gff3` format from the [site](http://www.mirbase.org/ftp.shtml). I found I had to convert this to `.gtf` first, which can be achieved using the R package `rtracklayer`:
```r
library(rtracklayer)
gff <- import("path/to/file.gff3")
export(gff, "path/to/file.gtf", "gtf")
```
The chromosome names in this `gtf` file _must_ match those in the original reference fasta sequence (and the version must be the same). If you are using an Ensembl reference, you may need to rename chromosomes in the `gtf` file to remove the prefix `chr`:
```bash
sed 's/chr\([1-9XY]*\)/\1/' annotation.chr.gtf > annotation.gtf
```
Finally, we run `featureCounts`:
```bash
featureCounts -a path/to/annotation.gtf -o featurecounts_output -t miRNA -T <num_threads> -g ID sample1.bam sample2.bam ... sampleN.bam
```
Note that this can accept a list of bam files separated by a space. If run in this way, all the results are included in a single output file. A summary file, named `featurecounts_output.summary` gives an overview of the number of reads assigned to features.

**NB** The default operation mode in `featureCounts` assumes that the reads are *unstranded*, which they probably are in single read mode(?). However, if using this for paired end data an additional parameter `-s 0/1/2` allows different configurations.