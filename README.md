# PlasTid Genome Assembly Using Long reads data (ptGAUL)
```
===========================================================================
                       _____           _        _         _    _
    ___      _       /  ___  \       / _ \     | |       | |  | |
   / _ \    | |     / /     \ \     / / \ \    | |       | |  | |
  / / \ \ __| |__  | |       \_|    / / \ \    | |       | |  | |
  ||   |||__   __| | |             / / _ \ \   | |       | |  | |
  | \_/ /   | |    | |      ___    /  ___  \   | |       | |  | |
  |  __/    | |_   | |     |__ |  / /     \ \  \ \       / /  | |        _
  | |       |   |   \ \ ___ / /   / /     \ \   \ \ ___ / /   | | _____ | |
  |_|       |__/     \ _____ /   /_/       \_\   \ _____ /    | _________ |

===========================================================================
```
This pipeline is used for plastid genome assembly based on long reads data, including both Nanopore and PacBio. It will help assemble the cpgenome with many repeat regions and reduce the assembly path number from short reads data. Our data is not published yet [Zhou et al. (unpublished)]. We introduced this pipeline in [BAGGs workshop](https://tarheels.live/baags/).

## Prerequisites and Software/dependencies

1. [minimap2](https://github.com/lh3/minimap2) or use [conda](https://anaconda.org/bioconda/minimap2) to install.
```
check if minimap2 is installed successfully by typing "minimap2 -h" in terminal.
```

2. [seqkit](https://bioinf.shenwei.me/seqkit/) or use [conda](https://anaconda.org/bioconda/seqkit) to install.for Step2.

```
check if seqkit is installed successfully by typing "seqkit -h" in terminal.
```
3. [assembly-stats](https://github.com/sanger-pathogens/assembly-stats) or use [conda](https://anaconda.org/bioconda/assembly-stats) to install.for Step2.

```
check if assembly-stats is installed successfully by typing "assembly-stats" in terminal.
```
4. [seqtk](https://github.com/lh3/seqtk) or use [conda](https://anaconda.org/bioconda/seqtk) to install. for Step2.
 
```
check if seqtk is installed successfully by typing "seqtk" in terminal.
```
5. [flye](https://github.com/fenderglass/Flye) or use [conda](https://anaconda.org/bioconda/flye) to install. for Step2.

```
check if flye is installed successfully by typing "flye -h" in terminal.
```
6. [python3](https://www.python.org/downloads/) and its dependencies for Step3.

```
check if python3 is installed successfully by typing "python3 -h" in terminal.
```

#### Optional Software for polishing step
7. [ropebwt2](https://github.com/lh3/ropebwt2) or use [conda](https://anaconda.org/bioconda/ropebwt2) to install.
```
check if ropebwt2 is installed successfully by typing "ropebwt2 -h" in terminal.
```
8. [msbwt](https://github.com/holtjma/msbwt) or use [conda](https://anaconda.org/kbchoi/msbwt) to install.
```
check if msbwt is installed successfully by typing "msbwt -h" in terminal.
```
9. [fmlrc](https://github.com/holtjma/fmlrc) or use [conda](https://anaconda.org/bioconda/fmlrc) to install.
Use the fmlrc instead of fmlrc2.
```
check if fmlrc is installed successfully by typing "fmlrc -h" in terminal.
```

10. put combine_gfa.py and ptGAUL.sh in the same directory.

## Environment
Examples can be run on Mac and Linux.

![](ptGAUL_image.png)

## Quick run
The basic things that you need to run ptGAUL are 1) a cpgenome from a closely related cpgenome data and 2) longread data (works for all fasta, fastq, and fq.gz files).

  
  ```
  Example: bash ptGAUL.sh -r [PATH]/[reference_genome]/ -l [PATH]/[long_read_data]
  ```
  
## Output


#### EXAMPLE
##### The command for the example data.
  ```bash
  bash ptGAUL.sh -r Beta.fasta -l SRR1980665.1 -t 8 -f 3000 
  ```
##### This script is used to find putative paralogs in 353 enrichment data. It requires two input parameters. The "-r", "-l" are required arguments.
  
  To check all parameters in ptGAUL using:
  ```bash
  bash ptGAUL.sh -h
  ```
  
##### Parameters in Detail
```bash
Usage: this script is used for plastome assembly using long read data.
ptGAUL.sh [options]
Options:
-r <MANDATORY: contigs or scaffolds in fasta format>
-l <MANDATORY: long reads in fasta/fastq/fq.gz format>
-t <number of threads, default:1>
-f <filter the long reads less than this number; default: 3000>
-h <help manu>
(base)
```

## (Optional) If the edge number does not equal 3
You should manually check the assembled data using [BANDAGE](https://rrwick.github.io/Bandage/). Then, manually run the python script again to get the assembly results including two paths.
```
python3 ./combine_gfa.py -e ./PATH_OF_EDGES_FILE/edges.fa -d ./PATH_OF_SORTED_DEPTH_FILE/sorted_depth
```

## (Optional) Final assembly polish using long reads data
install racon using conda.
```minimap2 -x ava-ont -t $n $asm $nanopore > ${racon_outdir}/map.paf
```racon -t $n $nanopore_fq ${racon_outdir}/map.paf $asm > ${racon_outdir}/asm.racon.fasta

## (Optional) Final assembly polish using short reads data
Highly recommend to use fmlrc for polishing step.
files illumina_* are the fq.gz file of illumina reads. Change the output path directory "/PATH/msbwt".

```
gunzip -c $illumina_1r1 $illumina_1r2 $illumina_2r1 $illumina_2r2 | awk 'NR % 4 == 2' | sort | tr NT TN | ropebwt2 -LR | tr NT TN | msbwt convert /PATH/msbwt
```
Once you finished msbwt run. $N means thread number. $assembled_cp is assembled plastome from ptGAUL. Change the output path of "/PATH/fmlrc/corrected.fasta"

```
fmlrc -p $N /PATH/msbwt/comp_msbwt.npy $assembled_cp /PATH/fmlrc/corrected_cp.fasta
```

## Citation

(coming soon) Zhou et al., Plastid Genome Assembly Using Long-read data (ptGAUL).

If you are using fmlrc, please cite Wang, Jeremy R. and Holt, James and McMillan, Leonard and Jones, Corbin D. FMLRC: Hybrid long read error correction using an FM-index. BMC Bioinformatics, 2018. 19 (1) 50.
