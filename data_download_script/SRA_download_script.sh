#Download Long read data or Illumina data
module add sratoolkit/3.0.0

fasterq-dump --split-files SRR14883332
#Download reference genome
#Long read data
esearch -db nucleotide -query "NC_035584.1" | efetch -format fasta > NC_035584.1.fasta

# Replace SRA accession numbers (SRR14883332) and reference number (NC_035584.1) with the following numbers.
# Long read data for ptGAUL validation:

# Arctostaphylos glauca (SRR14883332) (Ref: NC_035584.1/NC_042713.1/NC_047438.1)
# Lepidium sativum (SRR18079816) (Ref: NC_047178.1)
# Chaetoceros muellerii (SRR13238610) (Ref: MW004650.1)
# Potentilla micrantha (ERR338629) (Ref: NC_015206.1)
# Durio zibethinus (SRR11547303) (Ref: MT321069)
# Beta vulgaris (SRR1980665) (Ref: KR230391.1)
# Eleocharis dulcis - available by asking authors (Ref: NC_047447.1)
# Eucalyptus pauciflora (SRR7153095) (Ref: MZ670598.1/HM347959.1/NC_014570.1/AY780259.1/ NC_039597.1)
# Leucanthemum vulgare (SRR10948618) (Ref: NC_047460.1)
# Oryza glaberrima (SRR8989349) (Ref: NC_024175.1)
# Cenchrus americanus (SRR8989348) (Ref: NC_024171.1)
# Digitaria exilis (SRR8989347) (Ref: NC_024176.1)
# Podococcus acaulis (SRR8989346) (Ref: NC_027276.1)
# Raphia textilis (SRR8989345) (Ref: NC_020365.1)
# Phytelephas aequatorialis (SRR8989344) (Ref: NC_029957.1)
# Picea glauca (SRR21038753) (Ref: NC_021456.1).


# Illumina data of Juncus for GetOrganelle:

# J. inflexus (SRR14298745)
# J. effusus (SRR13309655 (Lu et al. 2021) and SRR14298746)
# J. roamerianus (SRR21976092)
# J. validus (SRR21976091)

# example of using GetOrganelle for J. effusus
fasterq-dump --split-files SRR13309655
# load GetOrganelle
get_organelle_from_reads.py -1 SRR13309655.1_1.fastq -2 SRR13309655.1_2.fastq -o SSR55 -R 15 -k 21,39,45,65,85,105,115,125 -F embplant_pt -t 20
get_organelle_from_reads.py -1 SRR13309655.1_1.fastq -2 SRR13309655.1_2.fastq -o SSR55_w075 -R 15 -k 21,39,45,65,85,105,115,125 -F embplant_pt -t 20 -w 0.75 
