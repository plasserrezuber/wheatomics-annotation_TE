#!/bin/bash

#SBATCH --job-name=annoTE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 1
#SBATCH --array=0-20
#SBATCH --mem=8G
#SBATCH --partition=normal

module load bedtools/2.27.1 exonerate samtools GenomeTools bioperl gdecTools

####################################################################################
##### DEFINITION VARIABLES

DATA='/home/palasser/data/RENAN_v2_pseudo'
OUTPUT='/home/palasser/results/RepeatMasker_clariTE/RENANV2_pseudoV2'

####################################################################################
chromosomes=("1A" "1B" "1D" "2A" "2B" "2D" "3A" "3B" "3D" "4A" "4B" "4D" "5A" "5B" "5D" "6A" "6B" "6D" "7A" "7B" "7D")
chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}
chrom='chr'$chr
mkdir -p $OUTPUT/${chrom}
cd $OUTPUT/${chrom}

###### Chunks
## create one fasta per chromosome
samtools faidx $DATA/Renan_v13_v2.pseudo.v2.fa $chrom > $OUTPUT/${chrom}/${chrom}.fasta
##other methods
##gawk -v seq='chr'$chrom 'BEGIN { RS=">" } { if ($0 ~ seq) print RS $0 }' $DATA/Renan_v13_v2.pseudo.v2.fa > $OUTPUT/${chrom}/${chrom}.fasta
##echo $chrom > chrom.txt; seqtk subseq $DATA/Renan_v13_v2.pseudo.v2.fa chrom.txt > $OUTPUT/${chrom}/${chrom}.fasta

samtools faidx $OUTPUT/$chrom/${chrom}.fasta
bedtools makewindows -g $OUTPUT/${chrom}/${chrom}.fasta.fai -w 10000000 > $OUTPUT/${chrom}/${chrom}.windows.bed

bedtools getfasta -bed $OUTPUT/${chrom}/${chrom}.windows.bed -fi $OUTPUT/${chrom}/${chrom}.fasta > $OUTPUT/${chrom}/${chrom}.windows.fasta
fastaexplode -f $OUTPUT/${chrom}/${chrom}.windows.fasta -d $OUTPUT/${chrom}

files=($(find $OUTPUT/${chrom}/${chrom}*.fa))
k=$((${#files[@]}-1))  #recupere le nb d'item dans la liste 'files'


###### RepeatMasker
RMjid=`sbatch --wait --parsable -J ${chr}_RM --array=0-${k}%20 /home/palasser/bin/repeatMasker_RENAN_V2.sh ${chrom}`

echo "$RMjid : RepeatMasker on ${chrom}"

if [ $(find ${chrom}:*.fa.out.xm |wc -l) = $(($k+1)) ];
    then echo ${chrom}": all RepeatMasker output files there";
    else echo ${chrom}": missing "$(($(find $OUTPUT/${chrom}/${chrom}:*.fa |wc -l) - $(find $OUTPUT/${chrom}/${chrom}:*.fa.out.xm |wc -l) ))" RepeatMasker output files";
fi

# # ## rescue:
# # sbatch -J 6A_RMrescue -p normal -c 16 --mem=8G --wrap="module load gcc/4.8.4 triannotTools/1.2 gdecTools/1.1 RepeatMasker; cd /home/palasser/results/RepeatMasker_clariTE/RENANV2_pseudoV2/chr6A; RepeatMasker -e crossmatch -lib /storage/groups/gdec/triannot/tools/ClariTE/clariTE_1.0/databank/CLARIwheat.fasta -xsmall -nolow -xm -pa 16 -q /home/palasser/results/RepeatMasker_clariTE/RENANV2_pseudoV2/chr6A/chr6A:80000000-90000000.fa"
# # sbatch -J 7B_RMrescue -p normal -c 16 --mem=8G --wrap="module load gcc/4.8.4 triannotTools/1.2 gdecTools/1.1 RepeatMasker; cd /home/palasser/results/RepeatMasker_clariTE/RENANV2_pseudoV2/chr7B; RepeatMasker -e crossmatch -lib /storage/groups/gdec/triannot/tools/ClariTE/clariTE_1.0/databank/CLARIwheat.fasta -xsmall -nolow -xm -pa 16 -q /home/palasser/results/RepeatMasker_clariTE/RENANV2_pseudoV2/chr7B/chr7B:440000000-450000000.fa"


###### clariTE
## option `--wait` pour que le script principal attende la fin du job array clariTE avant cmd suivante
# CLARITEjid=`sbatch --wait --parsable -J ${chr}_clariTE --dependency=afterok:$RMjid --array=0-$k%20 /home/palasser/bin/clariTE_RENAN_V2.sh ${chrom}`
CLARITEjid=`sbatch --wait --parsable -J ${chr}_clariTE --array=0-$k%20 /home/palasser/bin/clariTE_RENAN_V2.sh ${chrom}`

echo "$CLARITEjid : clariTE on ${chrom}"

if [ $(find ${chrom}:*.fa.out_anno.embl |wc -l) = $(($k+1)) ];
    then echo ${chrom}": all clariTE output files there";
    else echo ${chrom}": missing "$(($(find $OUTPUT/${chrom}/${chrom}:*.fa |wc -l) - $(find $OUTPUT/${chrom}/${chrom}:*.fa.out_anno.embl |wc -l) ))" clariTE output files";
fi

###### embl2gff: takes all chunks embl and produces one merged .gff
## le script embl2gff.generic.pl numerote les TE par chrom (et non par chunk) en repartant de zero pour chaque chrom a condition
## de fournir la liste ordonnee des fichiers pour l'ensemble des chunks par chromosome
embl_files=$(\ls -1 $OUTPUT/${chrom}/${chrom}*.embl |sort -t ':' -k2,2n |tr -s '\n' ' ')
/storage/groups/gdec/bin/scripts/embl2gff.generic.pl -renan -note -l 10 $embl_files > $OUTPUT/TaeRenan_refseq_v2.0_${chr}_clariTE.gff

###### gff to gff3: recalcul des coordonnees (transfo de relatives aux chunks vers relatives au chrom) + mise en forme gff3 par GenomeTools sort
grep -v $'\t''region' $OUTPUT/TaeRenan_refseq_v2.0_${chr}_clariTE.gff |sed 's/embl2gff.generic.pl-2.1/clariTE/' \
|gawk -v FS='\t' -v OFS='\t' '/^chr/ {match($1, /:[0-9]+-/); $4=$4+substr($0, RSTART+1, RLENGTH-2); $5=$5+substr($0, RSTART+1, RLENGTH-2); print}' \
|gt gff3 -sort -tidy -retainids 1> $OUTPUT/TaeRenan_refseq_v2.0_${chr}_clariTE.gff3  2> $OUTPUT/TaeRenan_refseq_v2.0_${chr}_clariTE_gt.log

###### friendly gff3
endchrom=$(grep $chrom $DATA/TaeRenan_refseq_v2.0.fa.fai |cut -f2)
grep -v -P '##sequence-region *'$chrom':[1-9]' $OUTPUT/TaeRenan_refseq_v2.0_${chr}_clariTE.gff3 \
|gawk -v chr=$chrom -v end=$endchrom 'BEGIN{FS="\t";OFS="\t"} { if ($0~"##sequence-region") $0="##sequence-region\t"chr"\t1\t"end; print }' \
|sed -E 's/'$chrom':[0-9]*-[0-9]*/'$chrom'/' |sed -E 's/Compo:.* (Family)/\1/' | sed -E 's/Post:.* (Status)/\1/' > $OUTPUT/TaeRenan_refseq_v2.0_${chr}_clariTE_friendly.gff3

echo "Check clariTE_gff3 validity for ${chrom}"
gt gff3validator $OUTPUT/TaeRenan_refseq_v2.0_${chr}_clariTE_friendly.gff3


rm $OUTPUT/${chrom}/${chrom}.fasta $OUTPUT/${chrom}/${chrom}.fasta.fai $OUTPUT/${chrom}/${chrom}.windows.fasta $OUTPUT/${chrom}/${chrom}.windows.bed

############################################################################################################################################################################
###### MAKE MERGED GFF FOR WHOLE GENOME
OUTPUT='/home/palasser/results/RepeatMasker_clariTE/RENANV2_pseudoV2'
# gt gff3 -sort -tidy -retainids $OUTPUT/TaeRenan_refseq_v2.0_*_clariTE_friendly.gff3 1>  $OUTPUT/TaeRenan_refseq_v2.0_clariTE.gff3 2> $OUTPUT/TaeRenan_refseq_v2.0_clariTE_gt.log

###### gff3 to bed
grep '\srepeat_region' $OUTPUT/TaeRenan_refseq_v2.0_clariTE.gff3 |gawk -v OFS='\t' '{print $1,$4-1,$5,$9}' > $OUTPUT/TaeRenan_refseq_v2.0_clariTE.bed
grep '\srepeat_region' $OUTPUT/unanchored/TaeRenan_refseq_v2.0_unanchored_scaffolds_clariTE.gff3 |gawk -v OFS='\t' '{print $1,$4-1,$5,$9}' > $OUTPUT/unanchored/TaeRenan_refseq_v2.0_unanchored_scaffolds_clariTE.bed

###### masking before Triannot
# ml bedtools
# bedtools maskfasta -fi $DATA/TaeRenan_refseq_v2.0.fa -bed $OUTPUT/TaeRenan_refseq_v2.0_clariTE.bed -fo $OUTPUT/TaeRenan_refseq_v2.0_TEmasked.fasta

###### generate TEs fasta file
# ml bedtools
bedtools getfasta -fi $DATA/TaeRenan_refseq_v2.0.fa -bed $OUTPUT/TaeRenan_refseq_v2.0_clariTE.bed -fo $OUTPUT/TaeRenan_refseq_v2.0_clariTE.fasta
bedtools getfasta -fi $DATA/TaeRenan_refseq_v2.0_unanchored_scaffolds.fa -bed $OUTPUT/unanchored/TaeRenan_refseq_v2.0_unanchored_scaffolds_clariTE.bed -fo $OUTPUT/unanchored/TaeRenan_refseq_v2.0_unanchored_scaffolds_clariTE.fasta

# ml exonerate
## genome entier
# fastacomposition $OUTPUT/TaeRenan_refseq_v2.0_TEmasked.fasta
## chrom par chrom
# fastacomposition -s --fasta $OUTPUT/TaeRenan_refseq_v2.0_TEmasked.fasta

###### TOTAL GENOME LENGTH repeat_region:
## bedtools merge: Merges overlapping BED/GFF/VCF entries into a single interval, the input file (-i) file must be sorted by chrom, then start.
bedtools merge -i <(grep '\srepeat_region\s' $OUTPUT/TaeRenan_refseq_v2.0_clariTE.gff3 |cut -f1,4,5 |sort -k1,1 -k2,2n) |gawk '{ sum+=$3-$2; print sum }' |tail -n1
### 12395112757

###### TOTAL GENOME LENGTH match_part:
bedtools merge -i <(grep '\smatch_part\s' $OUTPUT/TaeRenan_refseq_v2.0_clariTE.gff3 |cut -f1,4,5 |sort -k1,1 -k2,2n) |gawk '{ sum+=$3-$2; print sum }' |tail -n1
### 11967394667 (diff avec repeat_region=427 665 657 s'explique par le fait que les match_part ne se suivent pas a la base pres (intervalle entre deux par ex), alors que repeat_region englobe tout)
gawk '{sum+=$3; print sum}' Length_TE_Renan_SUPERFamLevel_Renan.tsv |tail -n1
### 11967447100 (diff 52433 bp)

###### TOTAL GENOME LENGTH ALL FEATURES (match_part overlapping repeat_region):
bedtools merge -i <(grep '^chr' $OUTPUT/TaeRenan_refseq_v2.0_clariTE.gff3 |cut -f1,4,5 |sort -k1,1 -k2,2n) |gawk '{ sum+=$3-$2; print sum }' |tail -n1
### 12395112982 (diff avec repeat_region=225)

/home/palasser/bin/length_Family_TE.py $OUTPUT/TaeRenan_refseq_v2.0_clariTE.gff3 |sort -k1,1 -k4,4n > $OUTPUT/TaeRenan_refseq_v2.0_clariTE_For_Length_TE_calculation.gff3

for chr in "1A" "1B" "1D" "2A" "2B" "2D" "3A" "3B" "3D" "4A" "4B" "4D" "5A" "5B" "5D" "6A" "6B" "6D" "7A" "7B" "7D";
do
    chrom='chr'$chr
    # grep ^$chrom $OUTPUT/TaeRenan_refseq_v2.0_clariTE_For_Length_TE_calculation.gff3 |gawk -v OFS='\t' '{print $NF,$5-$4}' |sed -E 's/\.[0-9]*//' \
    # |gawk -v grp=$chrom -v OFS='\t' '{a[$1]+=$2}END{for(i in a) print grp,i,a[i]}'  |sort -k1,1 >> Length_TE_Renan_FamLevel_Renan.tsv

    grep ^$chrom $OUTPUT/TaeRenan_refseq_v2.0_clariTE_For_Length_TE_calculation.gff3 |gawk -v OFS='\t' '{print $NF,$5-$4}' |sed -E 's/_.*\t/\t/' \
    |gawk -v grp=$chrom -v OFS='\t' '{a[$1]+=$2}END{for(i in a) print grp,i,a[i]}'  |sort -k1,1 >> Length_TE_Renan_SUPERFamLevel_Renan.tsv

    grep ^$chrom $OUTPUT/TaeRenan_refseq_v2.0_clariTE.gff3 |grep '\srepeat_region\s' |cut -f9 |sed -E 's/ID=.*Family:(.*)[_|w].* Matching.*/\1/' \
    |sort |uniq -c |gawk -v grp=$chrom -v OFS='\t' '{print grp,$2,$1}' >> Nb_TE_Renan_SUPERFamLevel_Renan.tsv
done
