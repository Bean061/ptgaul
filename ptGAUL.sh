#!/bin/bash
set -e

# ptGAUL represents plastid genome assembly using long reads data.
# We prefer using the long reads corrected by fmlrc (https://github.com/holtjma/fmlrc).
# This pipeline is used for plastome assembly using long reads data.
# It can help resolve the fragmented contigs results from getorganelles.
# It is more suitable for long reads data with large N50.
# Please contact Wenbin (Bean) Zhou. Email: wenbin.evolution@gmail.com

start=$(date +%s)

######## get prepared ########
#### using conda install the following modules ####
# 1. minimap2
# 2. seqkit
# 3. assembly-stats
# 4. seqtk
# 5. flye 2.7

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
  echo "Usage: ptGAUL.sh -r (REFERENCE FILE) -l (LONG READ FILE)"
  echo ""
  echo "                 [-t threads int] [-g genome size int]"
  echo "                 [-c coverage int] [-f filter threshold int]"
  echo "                 [-o output directory string]"
  echo ""
  echo "this pipeline is used for plastome assembly using long read data."
  echo ""
  echo "optional arguments: "
  echo "-h, --help            <show this help message and exit>"
  echo "-r, --reference       <MANDATORY: reference contigs or scaffolds in fasta format>"
  echo "-l, --longreads       <MANDATORY: raw long reads in fasta/fastq/fq.gz format>"
  echo "-t, --threads         <number of threads, default:1>" 
  echo "-g, --genomesize      <expected genome size of plastome (bp), default:160000>"
  echo "-c, --coverage        <a rough coverage of data used for plastome assembly, default:50>"
  echo "-f, --filtered        <the raw long reads will be filtered if the lengths are less than this number (bp); default: 3000>"
  echo "-o, --outputdir       <output directory of results, defult is current directory>"
  echo ""
  echo ""
  echo "                     _____           _        _         _    _"
  echo "  ___      _       /  ___  \       / _ \     | |       | |  | |"
  echo " / _ \    | |     / /     \ \     / / \ \    | |       | |  | |"
  echo "/ / \ \ __| |__  | |       \_|    / / \ \    | |       | |  | |"
  echo "||   |||__   __| | |             / / _ \ \   | |       | |  | |"
  echo "| \_/ /   | |    | |      ___    /  ___  \   | |       | |  | |"
  echo "|  __/    | |_   | |     |__ |  / /     \ \  \ \       / /  | |        _"
  echo "| |       |   |   \ \ ___ / /   / /     \ \   \ \ ___ / /   | | _____ | |"
  echo "|_|       |__/     \ _____ /   /_/       \_\   \ _____ /    | _________ |"
  echo ""
  echo ""
  echo ""
}


NUM_THREADS=1
FRL=3000
GEM=160000
COV=50
OPD=./

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
        -c|--coverage)
            COV="$2"
            shift
            ;;
        -g|--genomesize)
            GEM="$2"
            shift
            ;;
        -o|--outputdir)
            OPD="$2"
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


echo ""
echo ""
echo ""
echo "                     _____           _        _         _    _"
echo "  ___      _       /  ___  \       / _ \     | |       | |  | |"
echo " / _ \    | |     / /     \ \     / / \ \    | |       | |  | |"
echo "/ / \ \ __| |__  | |       \_|    / / \ \    | |       | |  | |"
echo "||   |||__   __| | |             / / _ \ \   | |       | |  | |"
echo "| \_/ /   | |    | |      ___    /  ___  \   | |       | |  | |"
echo "|  __/    | |_   | |     |__ |  / /     \ \  \ \       / /  | |        _"
echo "| |       |   |   \ \ ___ / /   / /     \ \   \ \ ___ / /   | | _____ | |"
echo "|_|       |__/     \ _____ /   /_/       \_\   \ _____ /    | _________ |"
echo ""
echo ""
echo ""


if [ ! -s $REF ];then
error_exit "reference $REF does not exist or size is 0 bp"
fi

if [ ! -s $LR ];then
error_exit "Long reads of ONT or PACBIO $LR does not exist or size is 0 bp"
fi

#echo $LR
#echo $NUM_THREADS
#echo $REF
#echo $GEM
#echo $COV

#if [ -e result_$FRL ];then
#	rm -rf result_$FRL
#fi

mkdir -p $OPD/result_$FRL/ptGAUL_final_assembly

asm_cov=$(( $COV - 20 )) 

#module add minimap2
### Step1 minimap to generate a table including the matching reads heads.
if [ ! -s $OPD/result_$FRL/filter_name ];then
	#minimap2 -cx map-ont $REF $LR | awk '{print $1}' > $OPD/result_$FRL/filter_name 2>error.txt
	# more filters
	minimap2 -cx map-ont $REF $LR > $OPD/result_$FRL/filter.paf 2>>error.txt
	awk '{print $1, $10, $11, $10/$11}' $OPD/result_$FRL/filter.paf > $OPD/result_$FRL/filter_1.paf
	awk '{if (($4>=0.7) && ($3 >=1000)) {print}}' $OPD/result_$FRL/filter_1.paf > $OPD/result_$FRL/filter_2.paf
	awk '{print $1}' $OPD/result_$FRL/filter_2.paf > $OPD/result_$FRL/filter_name

fi

echo ""
echo "####################################"
echo "########## Step1 results ###########"
echo "####################################"
echo ">> $LR is the sequence input of long read data."
echo ">> $REF is the reference data of closely related species."
echo ""
echo 'Step1 is finished!! The filter_name includes all names of plastid long reads that match to the reference.'

checkpoint1=$(date +%s)
echo "Step1 took $(($checkpoint1 - $start)) sec."

#module add seqkit
### Step2 keep the matching reads. 
### $1 is bash script for the filter.
if [ ! -s $OPD/result_$FRL/filter_reads ];then
	if [[ $LR = *.fasta || $LR = *.fa ]];then		
		filter_fasta $OPD/result_$FRL/filter_name $LR > $OPD/result_$FRL/filter_reads

	elif [[ $LR = *.fastq || $LR = *.fq ]];then
		filter_fastq $OPD/result_$FRL/filter_name $LR > $OPD/result_$FRL/filter_reads.fq
		seqkit fq2fa $OPD/result_$FRL/filter_reads.fq > $OPD/result_$FRL/filter_reads 2>>error.txt
		rm $OPD/result_$FRL/filter_reads.fq

	elif [[ $LR = *.fq.gz || $LR = *.fastq.gz ]];then
		filter_fastq_gz $OPD/result_$FRL/filter_name $LR > $OPD/result_$FRL/filter_reads.fq
		seqkit fq2fa $OPD/result_$FRL/filter_reads.fq > $OPD/result_$FRL/filter_reads 2>>error.txt
		rm $OPD/result_$FRL/filter_reads.fq
	fi
fi

echo ""
echo "####################################"
echo "########## Step2 results ###########"
echo "####################################"
echo ""
echo 'Step2 is finished!! filter_reads is your filtered sequence reads with all plastid sequences matching to the reference.'
checkpoint2=$(date +%s)
echo "Step2 took $(($checkpoint2 - $checkpoint1)) sec."


if [ ! -s $OPD/result_$FRL/new_filter_gt$FRL.fa ];then
	seqkit seq -g -m $FRL $OPD/result_$FRL/filter_reads \
	> $OPD/result_$FRL/new_filter_gt$FRL.fa 2>>error.txt

	echo ""
	echo "####################################"
	echo "########## Step3 results ###########"
	echo "####################################"
	echo ""
	echo 'Step3 is finished!! new_filter_gt'"$FRL"'.fa stores all plastid long reads greater than '"$FRL"' bp.'

fi
checkpoint3=$(date +%s)
echo "Step3 took $(($checkpoint3 - $checkpoint2)) sec."


if [ ! -s $OPD/result_$FRL/total_length_of_plastid_reads ];then
	assembly-stats $OPD/result_$FRL/new_filter_gt$FRL.fa 2>>error.txt \
	| awk 'FNR == 2 {print $3}' \
	| sed 's/,//' > $OPD/result_$FRL/total_length_of_plastid_reads

	num=$(awk '$1/'"$GEM"'>10 && $1/'"$GEM"'<'"$COV"'{print $1/'"$GEM"'}' $OPD/result_$FRL/total_length_of_plastid_reads)
	real_cov=$(awk '{print $1/'"$GEM"'}' $OPD/result_$FRL/total_length_of_plastid_reads)
	echo ""
	echo "####################################"
	echo "########## Step4 results ###########"
	echo "####################################"
	
	if [[ $num = "" ]];then
		echo ""
		echo 'Step4 is finished!! It has '"$real_cov"'x coverage of plastome data.'
		echo 'We need to reduce reads to '"$COV"'x coverage.'
		frac=$(awk '{print '"$GEM"'*'"$COV"'/$1}' $OPD/result_$FRL/total_length_of_plastid_reads)
		echo 'The fraction number is '"$frac"'.'
		

		if [ $frac != 0 ];then
			seqtk sample -s100 $OPD/result_$FRL/new_filter_gt$FRL.fa $frac \
			> $OPD/result_$FRL/filter_reads_final.fa 
			echo ""
			echo "####################################"
			echo "########## Step5 results ###########"
			echo "####################################"
			echo ""
			echo 'Step5 is finished!! The filter_reads_final.fa includes about '"$COV"'x coverage of long reads.'
			echo "It is ready for flye assembly."
			
			flye --nano-raw $OPD/result_$FRL/filter_reads_final.fa --genome-size 0.16m \
			--out-dir $OPD/result_$FRL/flye_cpONT --threads $NUM_THREADS --asm-coverage $asm_cov 2>>error.txt
			
		fi
	else
		echo ""
		echo 'It has '"$num"'x coverage of plastome data.'
		echo 'Step4 is finished!! The coverage of plastome is between 10-'"$COV"'x coverage.'
		echo "It is ready for flye assembly."

		flye --nano-raw $OPD/result_$FRL/new_filter_gt$FRL.fa --genome-size 0.16m \
		--out-dir $OPD/result_$FRL/flye_cpONT --threads $NUM_THREADS --asm-coverage $asm_cov 2>>error.txt

	fi
fi

checkpoint6=$(date +%s)
echo "Step4, 5, and 6 take $(($checkpoint6 - $checkpoint3)) sec."

if [ -s $OPD/result_$FRL/flye_cpONT/assembly.fasta ];then
	echo ""
	echo "####################################"
	echo "########## Step6 results ###########"
	echo "####################################"
    echo ""
    echo "Assembly is done. All results are in ./flye_cpONT/"
	echo "Please use the Bandage to visualize the assembly result."
	echo 'If the edges are three, it is a well-assembled result! Congratulations!!'
	echo "Otherwise, you need to manually check the graph results."
	echo "You can use blast the assembly sequence itself to check the IR region, and use either long reads or short reads to check the evenness of the coverages."
fi


if [ ! -s $OPD/result_$FRL/edges.fa ];then
	awk '/^S/{print ">"$2"\n"$3}' $OPD/result_$FRL/flye_cpONT/assembly_graph.gfa \
	| fold > $OPD/result_$FRL/edges.fa
	awk -F ":" '/^S/{print substr($1,3,7)"\t"$3}' $OPD/result_$FRL/flye_cpONT/assembly_graph.gfa \
	| sort -k 2 -n -r > $OPD/result_$FRL/sorted_depth
fi	

number_edge=$(grep ">" $OPD/result_$FRL/edges.fa | wc -l)

echo ""
echo "####################################"
echo "########## Step7 results ###########"
echo "####################################"
if [[ $number_edge -eq 3 ]]; then
	echo ""
	echo '==============================================='
	echo '-------------- Congratulations ----------------'
	echo '==============================================='
	echo ""
	echo "The result_$FRL/flye_cpONT/assembly_graph.gfa file has 3 edges for a plastome. Check two genotypes of plastome: ./path1.fasta and ./path2.fasta."
	combine_gfa.py -e $OPD/result_$FRL/edges.fa -d $OPD/result_$FRL/sorted_depth -o $OPD/result_$FRL/ptGAUL_final_assembly/ 2>>error.txt
	end=$(date +%s)
	echo "It took $(($end - $start)) sec in total."

elif [[ $number_edge -eq 1 ]]; then
	echo ""
	echo '==============================================='
	echo '---------------- Attentions -------------------'
	echo '==============================================='
	echo ""
	echo "The gfa file only has 1 edge. Please manually check the result_$FRL/flye_cpONT/assembly_graph.gfa file in Bandage."
	echo "Confirm the assembly result is a circle with a reasonable length (~160k bp)."
	echo "Congratulations!!!"
	ln $OPD/result_$FRL/flye_cpONT/assembly.fasta $OPD/result_$FRL/ptGAUL_final_assembly/final_assembly.fasta
	echo "result_$FRL/ptGAUL_final_assembly/final_assembly.fasta is the final plastome. You can check the details in the flye output."

else
	echo "edge number is $number_edge" '!'
	echo ""
	echo '==============================================='
	echo '-------- Oops! Detect a weird result ----------'
	echo '==============================================='
	echo ""
	echo "The gfa file has more than/less than three edges. Please manually check the result_$FRL/flye_cpONT/assembly_graph.gfa file in Bandage."
	echo "After confirming the assembly result, keep three edges in the edges.fa and three lines in sorted_depth files."
	echo "Then run the python script as follows manually:"
	echo "combine_gfa.py -e result_$FRL/edges.fa -d result_$FRL/sorted_depth -o result_$FRL/ptGAUL_final_assembly"

fi


end=$(date +%s)
echo ''
echo '==============================================='
echo '------------- total running time --------------'
echo '==============================================='
echo ""
echo "It took $(($end - $start)) sec in total."
