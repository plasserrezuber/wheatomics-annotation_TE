#!/bin/bash
#SBATCH --job-name=clariUn
#SBATCH --array=0-337%20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH --partition=fast
#SBATCH --export=ALL

module load gcc/4.8.4 triannotTools/1.2 gdecTools/1.1 bioperl/1.7.0_rc5 GenomeTools

##### DEFINITION VARIABLES
INPUT='/home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fasta'
DATABANK='/storage/groups/gdec/triannot/tools/ClariTE/clariTE_1.0'
OUTPUT='/home/palasser/results/RepeatMasker_clariTE/RENANV2_pseudoV2/unanchored'

cd $OUTPUT

files=($(find $INPUT/Tae*.fa))
i=${files[$SLURM_ARRAY_TASK_ID]}
prefix=$(basename $i)

# mise en forme
sed -i "s/#Unspecified/#Unknown/" $prefix.out.xm

# ##### clariTE
$DATABANK/bin/clariTE.pl -dir $OUTPUT/ \
-LTR $DATABANK/databank/CLARIwheat.LTR_position \
-classi $DATABANK/databank/CLARIwheat.classification \
-fasta $i \
$prefix.out.xm


###### embl2gff
## le script embl2gff.generic.pl numerote les TE
/storage/groups/gdec/bin/scripts/embl2gff.generic.pl -renan -note -l 10 ${prefix}.out_anno.embl > ${prefix}.out_anno.gff

###### mise en forme gff3 par GenomeTools sort
sed 's/.fa.out_anno.embl.1//' ${prefix}.out_anno.gff |sed 's/embl2gff.generic.pl-2.1/clariTE/' \
|gt gff3 -sort -tidy -retainids 1> ${prefix}_clariTE.gff3  2> ${prefix}_clariTE_gt.log

###### friendly gff3
sed -E 's/Compo:.* (Family)/\1/' ${prefix}_clariTE.gff3 |sed -E 's/Post:.* (Status)/\1/' > ${prefix}_clariTE_friendly.gff3

###### echo "Check clariTE_gff3 validity"
gt gff3validator ${prefix}_clariTE_friendly.gff3

###### merge gff3
# gt gff3 -tidy -retainids $(\ls -1 Tae*.fa_clariTE_friendly.gff3 |sort -V |tr '\n' ' ') \
# 1> TaeRenan_refseq_v2.0_unanchored_scaffolds_clariTE.gff3 2> TaeRenan_refseq_v2.0_unanchored_scaffolds_clariTE_gt.log

###### verif
# grep $'\t''region'$'\t' TaeRenan_refseq_v2.0_unanchored_scaffolds_clariTE.gff3 |wc -l
# gt gff3validator TaeRenan_refseq_v2.0_unanchored_scaffolds_clariTE.gff3