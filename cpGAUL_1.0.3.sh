#!/bin/bash

# cpGAUL represents chloroplast genome assembly using long reads data.
# We prefer using the long reads corrected by fmlrc (https://github.com/holtjma/fmlrc).
# This pipeline is used for chloroplast genome assembly using long reads data.
# It can help resolve the fragmented contigs results from getorganelles.
# It is more suitable for long reads data with large N50.
# Please contact wenbin.evolution@gmail.com



######## get prepared ########
#### using conda install the following modules ####
# 1. minimap2
# 2. seqkit
# 3. assembly-stats
# 4. seqtk
# 5. flye

function filter_fasta () {
	awk 'BEGIN{keep=0;}
	NR==FNR{remove[$1]=1}
	NR!=FNR{
	if(substr($1,1,1)==">"){if(substr($1,2) in remove){keep=1}else{keep=0}}
	if(keep==1){print}
	}' $1 $2
}

function filter_fastq () {
	awk 'BEGIN{keep=0;}
	NR==FNR{remove[$1]=1}
	NR!=FNR{
	if(FNR%4==1){if(substr($1,2) in remove){keep=1}else{keep=0}}
	if(keep==1){print}
	}' $1 $2
}


function filter_fastq_gz (){
	awk 'BEGIN{keep=0;}
	NR==FNR{remove[$1]=1}
	NR!=FNR{
	if(FNR%4==1){if(substr($1,2) in remove){keep=1}else{keep=0}}
	if(keep==1){print}
	}' $1 <(zcat $2)
}

function usage () {
  echo "Usage: this script is used for ONT cp genome assembly."
  echo "cp_genome_ONT.sh [options]"
  echo "Options: "
  echo "-r <MANDATORY: contigs or scaffolds in fasta format>"
  echo "-l <MANDATORY: long reads in fasta/fastq/fq.gz format>"
  echo "-t <number of threads, default:1>" 
  echo "-f <filter the long reads less than this number; default: 3000>"
}

NUM_THREADS=1
FRL=3000

if [ $# -lt 1 ]; then
	usage
fi

while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -t|--threads)
            NUM_THREADS="$2"
            shift
            ;;
        -l|--longreads)
            LR="$2"
            shift
            ;;
        -r|--reference)
            REF="$2"
            shift
            ;;
        -f|--filtered)
            FRL="$2"
            shift
            ;;
	-h|--help|-u|--usage)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option $1"
            exit 1        # unknown option
            ;;
    esac
    shift
done


if [ ! -s $REF ];then
error_exit "reference $REF does not exist or size zero"
fi

if [ ! -s $LR ];then
error_exit "Long reads of ONT or PACBIO $LR does not exist or size zero"
fi

#echo $LR
#echo $NUM_THREADS
#echo $REF

if [ -e result_$FRL ];then
	rm -rf result_$FRL
fi

mkdir result_$FRL

module add minimap2
### Step1 minimap to generate a table including the matching reads heads.
if [ ! -s result_$FRL/filter_name ];then
	#minimap2 -cx map-ont $REF $LR | awk '{print $1}' > result_$FRL/filter_name 2>error.txt
	# more filters
	minimap2 -cx map-ont $REF $LR > result_$FRL/filter.paf
	awk '{print $1, $10, $11, $10/$11}' result_$FRL/filter.paf > result_$FRL/filter_1.paf
	awk '{if (($4>=0.7) && ($3 >=1000)) {print}}' result_$FRL/filter_1.paf > result_$FRL/filter_2.paf
	awk '{print $1}' result_$FRL/filter_2.paf > result_$FRL/filter_name

	echo "########## Step1 result ###########"
	echo $LR
	echo "Step1 has been finished!! The filter_name includes all names of chloroplast long reads match to the reference."
fi

module add seqkit
### Step2 keep the matching reads. 
### $1 is bash script for the filter.
if [ ! -s result_$FRL/filter_reads ];then
	if [[ $LR = *.fasta || $LR = *.fa ]];then		
		filter_fasta result_$FRL/filter_name $LR > result_$FRL/filter_reads

	elif [[ $LR = *.fastq || $LR = *.fq ]];then
		filter_fastq result_$FRL/filter_name $LR > result_$FRL/filter_reads.fq
		seqkit fq2fa result_$FRL/filter_reads.fq > result_$FRL/filter_reads
		rm result_$FRL/filter_reads.fq
	elif [[ $LR = *.fq.gz || $LR = *.fastq.gz ]];then
		filter_fastq_gz result_$FRL/filter_name $LR > result_$FRL/filter_reads.fq
		seqkit fq2fa result_$FRL/filter_reads.fq > result_$FRL/filter_reads
		rm result_$FRL/filter_reads.fq
	fi
fi

echo "########## Step2 result ###########"
echo "filter_reads"
echo "Step2 has been finished!! The filter_reads includes all chloroplast long reads match to the reference."

# module add seqkit
#if [ ! -s filter_reads_final ]; then
### Step3 filter reads > $FRL bp, and keep a fraction of reads
### make sure loaded seqkit, seqtk, and assembly-stats module.


if [ ! -s result_$FRL/new_filter_gt$FRL.fa ];then
	seqkit seq -g -m $FRL result_$FRL/filter_reads > result_$FRL/new_filter_gt$FRL.fa 2>>error.txt
	echo "########## Step3 result ###########"
	echo "new_filter_gt$FRL.fa"
	echo "Step3 has been finished!! The new_filter_gt$FRL.fa keeps all chloroplast long reads greater than $FRL bp."

fi

module add seqtk
module add flye
module add anaconda/2021.11
source activate grasstool
if [ ! -s result_$FRL/total_length ];then
	assembly-stats result_$FRL/new_filter_gt$FRL.fa | awk 'FNR == 2 {print $3}' | sed 's/,//' > result_$FRL/total_length 
	num=$(awk '$1/160000>15 && $1/160000<50{print $1/160000}' result_$FRL/total_length)
	echo "########## Step4 result ###########"
	echo $num
	if [[ $num = "" ]];then
		echo "Step4 has been finished!! The coverage of chloroplast genome is greater than 50x coverage"
		echo "It needs to reduce reads to 50x coverage"
		frac=$(awk '{print 160000*50/$1}' result_$FRL/total_length)
		echo "The fraction number is" $frac
		if [ $frac != 0 ];then
			seqtk sample -s100 result_$FRL/new_filter_gt$FRL.fa $frac > result_$FRL/filter_reads_final.fa 
			echo "########## Step5 result ###########"
			echo "filter_reads_final.fa"
			echo "Step5 has been finished!! The filter_reads_final.fa includes about 50x coverage of long reads."
			echo "It is ready for flye assembly."
			flye --nano-raw result_$FRL/filter_reads_final.fa --genome-size 0.16m --out-dir ./result_$FRL/flye_cpONT --threads $NUM_THREADS 2>>error.txt

		fi
	else
		echo "Step4 has been finished!! The coverage of chloroplast genome is between 15-50x coverage"
		echo "It is ready for flye assembly"

		flye --nano-raw result_$FRL/new_filter_gt$FRL.fa --genome-size 0.16m --out-dir ./result_$FRL/flye_cpONT --threads $NUM_THREADS 2>>error.txt

	fi
fi

if [ -s ./result_$FRL/flye_cpONT/assembly.fasta ];then
	echo "########## Step6 result ###########"
	echo "Assembly is done. ./flye_cpONT/assembly.fasta."
	echo "Please check the fragments and length in ./result_$FRL/flye_cpONT/flye.log"
	echo "You can use blast the aseembly sequence itself to check the IR region, and use either long reads and short reads to check the evenness of the coverages."
fi
