#!/bin/bash

#SBATCH --job-name=1Aformat
##SBATCH --array=0-59%32
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH --partition=fast
#SBATCH --export=ALL

#lancement script pour chromosome 1A:
#sbatch --array=0-59%32 /home/palasser/bin/formatting_clariTE_outputs_RENAN_V2.sh chr1A

#for chr in "1A" "1B" "1D" "2A" "2B" "2D" "3A" "3B" "3D" "4A" "4B" "4D" "5A" "5B" "5D" "6A" "6B" "6D" "7A" "7B" "7D"; do ls /home/napapon/Renan/Fasta/chr$chr*.fa |wc -l; done
#            60   71   50   80   82   67   76   87   64   76   67   53   72   71   52   63   72   50   78   76   64

module load gcc/4.8.4 triannotTools/1.2 gdecTools/1.1 bioperl/1.7.0_rc5 GenomeTools

##### DEFINITION VARIABLES
INPUT='/home/napapon/Renan/ClariTE'
OUTPUT='/home/palasser/results/RepeatMasker_clariTE/RENAN_v2/gff'
mkdir $OUTPUT

files=($(find /home/napapon/Renan/Fasta/$1*.fa))
i=$(basename ${files[$SLURM_ARRAY_TASK_ID]})

##### embl to gff
/storage/groups/gdec/bin/scripts/embl2gff.clarite.pl $INPUT/${i}.out_without_Unspecified_anno.embl > $OUTPUT/${i}.annoTE.gff

#chrom=$(cut -d'#' -f1 <(echo $i))
gffsource=$(cut -d'.' -f1 <(echo $i))
sed 's/'$gffsource'/'$1'/' $OUTPUT/${i}.annoTE.gff |sed 's/EMBL2GFF/clariTE/' |grep -v 'EMBL'$'\t''region' \
|gt gff3 -sort -tidy -retainids 1> $OUTPUT/${i}.annoTE_gtsort.gff3  2> $OUTPUT/${i}.annoTE_gt.log

# recalcul des coordonnees: passage coordonnees relatives aux chunks a celles relatives au chrom
start=$(cut -d'#' -f2 <(echo $i))
endchrom=$(grep $1 /storage/databanks/bio/triannot/databaseTemp/TraesRenanV2pseudoV1.fasta.fai |cut -f2)

gawk -v chr=$1 -v end=$endchrom 'BEGIN{FS="\t";OFS="\t"} { if ($0~"##sequence-region") $0="##sequence-region\t"chr"\t1\t"end; print }' $OUTPUT/${i}.annoTE_gtsort.gff3 \
|gawk -v strt=$start 'BEGIN{FS="\t";OFS="\t"} { if ($0!~"^#") $4=$4+strt ; if ($0!~"^#") $5=$5+strt; print }' > $OUTPUT/${i}.annoTE.gff3

gt gff3validator $OUTPUT/${i}.annoTE.gff3

#rm $OUTPUT/${i}.out_anno_addID.gff3 $OUTPUT/${i}.out_anno_addID_gtsort.gff3

#verif
#cat slurm* |grep -c 'input is valid GFF3' 
#si chiffre egal au chiffre ligne 13 indiquant le nb de fichier fasta chunk, tout s'est bien passe
#rm slurm* *gt*

## merge gff3
# for chr in "1A" "1B" "1D" "2A" "2B" "2D" "3A" "3B" "3D" "4A" "4B" "4D" "5A" "5B" "5D" "6A" "6B" "6D" "7A" "7B" "7D";
# do
#     gt gff3 -sort -tidy -retainids chr${chr}*.fa.annoTE.gff3 1> TraesRenanV2pseudoV1_clariTE_${chr}.gff3 2> TraesRenanV2pseudoV1_clariTE_${chr}_gt.log
# done

## merge gff3 whole genome
# gt gff3 -sort -tidy -retainids TraesRenanV2pseudoV1_clariTE_*.gff3 1> TraesRenanV2pseudoV1_clariTE.gff3 2> TraesRenanV2pseudoV1_clariTE_gt.log
# gt gff3validator TraesRenanV2pseudoV1_clariTE.gff3

## genome fasta masking
# grep -v '#' TraesRenanV2pseudoV1_clariTE.gff3 |gawk -F'\t' '{print $1"\t"$4-1"\t"$5"\t"$9}' > TraesRenanV2pseudoV1_clariTE.bed
# sbatch -p smp -c 1 --mem=128G --wrap="ml bedtools; bedtools maskfasta -fi /home/celmonat/Renan_v13_v2.pseudo.v1.fa -bed /home/palasser/results/RepeatMasker_clariTE/RENAN_v2/gff/TraesRenanV2pseudoV1_clariTE.bed -fo /home/palasser/results/RepeatMasker_clariTE/RENAN_v2/gff/TraesRenanV2pseudoV1_clariTE_masked.fasta"

# ml exonerate
# fastacomposition /home/palasser/results/RepeatMasker_clariTE/RENAN_v2/gff/TraesRenanV2pseudoV1_clariTE_masked.fasta