#!/usr/bin/python3

##SYNOPSIS:
#/home/palasser/bin/length_Family_TE_Renan.py TaeRenan_refseq_v2.0_clariTE.gff3
#Permet d'obtenir fichier pour calcul longueur de chaque famille de TE

#Permet de restituter l'information de la famille provenant des features "repeat_region" pour chaque "match_part" associe
#exemple : chr1A:0-10000000        clariTE repeat_region   17498   19831   .       -       1       ID=TraesRe_chr1A_repeat_region_000001;Note=Compo:RLG_famc1.4 100.00 Family:RLG_famc1.4 Matching_repeat:ctgI_rep_0052 Post:RLG_famc1.4 9066bp 6722..9066 Status:fragmented
#chr1A:0-10000000        clariTE match_part      17498   19831   .       -       1       ID=TraesRe_chr1A_repeat_region_000001_1;Parent=TraesRe_chr1A_repeat_region_000001
#devient une ligne: chr1A:0-10000000        clariTE match_part      17498   19831   .       -       1       ID=TraesRe_chr1A_repeat_region_000001_1;Parent=TraesRe_chr1A_repeat_region_000001 RLG_famc1.4

import sys
from collections import defaultdict

match_part=defaultdict(list)
repeat_region=defaultdict(list)

with open(sys.argv[1],"r") as f1:
	for line in f1:
		if line.startswith("chr") :
			lig_final=line.rstrip()
			li=line.rstrip().split()

			if li[2] == "match_part" :
				id_match=li[-1].split("Parent=")[1]
				#dico cle=ID of Parent repeat_region, value=list of match_part lines
				match_part[id_match].append(lig_final)

			if li[2] == "repeat_region" :
				id_repeat=li[8].split(";Note=")[0][3:]
				family_TE=li[8].split(";Note=")[1][7:]
				#dico cle=ID of repeat_region, value=list of family
				repeat_region[id_repeat].append(family_TE)
				
for key in repeat_region :
	if key in match_part:
		for elem in match_part[key]:
			print(elem, repeat_region[key][0])
