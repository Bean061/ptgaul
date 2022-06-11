# chloroplast Genome Assembly Using Long reads data (cpgaul)

This pipeline is used for plastid genome assembly based on long reads data, including both Nanopore and PacBio. It will help assemble the cpgenome with many repeat regions and reduce the assembly path number from short reads data. Our data is not published yet [Zhou et al. (unpublished)]. We introduced this pipeline in [BAGGs workshop]().

## Prerequisites and Software/dependencies

1. [minimap2](https://github.com/lh3/minimap2).
```
check if minimap2 is installed successfully by typing "minimap2 -h" in terminal.
```

2. [seqkit](https://bioinf.shenwei.me/seqkit/) for Step2.

```
check if seqkit is installed successfully by typing "seqkit -h" in terminal.
```
3. [assembly-stats](https://github.com/sanger-pathogens/assembly-stats) for Step2.

```
check if assembly-stats is installed successfully by typing "assembly-stats" in terminal.
```
4. [seqtk](https://github.com/lh3/seqtk) for Step2.
 
```
check if seqtk is installed successfully by typing "seqtk" in terminal.
```
5. [flye](https://github.com/fenderglass/Flye) for Step2.

```
check if flye is installed successfully by typing "flye -h" in terminal.
```
6. [python3](https://www.python.org/downloads/) and its dependencies for Step3.

```
check if python3 is installed successfully by typing "python3 -h" in terminal.
```

## Environment
Examples can be run on Mac and Linux.

## Quick run
The basic things that you need to run cpGAUL are 1) a cpgenome from a closely related cpgenome data and 2) longread data (works for all fasta, fastq, and fq.gz files).

  
  ```
  Example: bash cpGAUL1_0_3.sh -r [PATH]/[reference_genome]/ -l [PATH]/[long_read_data]
  ```
  
## Output


#### EXAMPLE
##### The command for the example data.
  ```bash
  bash cpGAUL_1.0.3.sh -r Beta.fasta -l SRR1980665.1 -t 8 -f 3000 
  ```
##### This script is used to find putative paralogs in 353 enrichment data. It requires two input parameters. The "-r", "-l" are required arguments.
  
  To check all parameters in cpGAUL using:
  ```bash
  bash cpGUAL.sh -h
  ```
  
##### Parameters in Detail
```bash
Usage: this script is used for ONT cp genome assembly.
cpGAUL.sh [options]
Options:
-r <MANDATORY: contigs or scaffolds in fasta format>
-l <MANDATORY: long reads in fasta/fastq/fq.gz format>
-t <number of threads, default:1>
-f <filter the long reads less than this number; default: 3000>
-h <help manu>
(base)
```

## repeat runs
```
for n in `seq 3000 1000 15000`; do
	bash cpGAUL_1.0.3.sh -r Beta.fasta -l SRR1980665.1 -t 8 -f $n
done
```
## Citation

coming soon.
