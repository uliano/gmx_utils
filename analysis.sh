# to select a subset of the trajectory make an index with a group
# containing the atoms you want to keep 
# gmx trjconv -s MD.tpr -f MD.xtc -o small.xtc -n the_index_you_just_made.ndx
# gmx convert-tpr -s MD.tpr -o small.tpr -n the_index_you_just_made.ndx

# these values are for graphics sampled every 100 ps
SKIP_ENERGY=50    
DT=100

BEGIN_EQUILIBRIUM=50000

SS_SIZE=500

MAX_EIGENVAL=10 


#   0 System              :  7599 atoms
#   1 Protein             :  7537 atoms
#   2 Protein-H           :  3687 atoms
#   3 C-alpha             :   462 atoms
#   4 Backbone            :  1388 atoms
#   5 MainChain           :  1851 atoms
#   6 MainChain+Cb        :  2290 atoms
#   7 MainChain+H         :  2304 atoms
#   8 SideChain           :  5233 atoms
#   9 SideChain-H         :  1836 atoms
#  10 Prot-Masses         :  7537 atoms
#  11 non-Protein         :    62 atoms
#  12 Other               :    62 atoms
#  13 LIG                 :    62 atoms

echo $'1 | 13\n name\n r 15-331\n name 15 Transmembrane\n r 1002-1161\n name 16 Lysozyme\n r 16-20\n name 17 Interface_Nter\n r 71 82 115 158 161 165 169 172 | r 179-181 | r 183-181\n name 18 Interface\n 17 | 18\n name 19 Interface_All\n r 138\n name 20 I138_3x46\n r258\n name 21 V258_6x36\n r 294\n name 22 E294_7x35\n r 311\n name 23 Y311_7x53\n r 295\n name 24 Y295_7x36\n r 101\n name 25 N101_2x60\n r 105\n name 26 S105_EL2\n r 121\n name 27 E121_3x29\n r 98\n name 28 Y98_2x57\n r 151\n name 29 Appunti_ciaraia1\n r154\n name 30 Appunti_ciaraia2\n 20 | 21 |  22 | 23 | 24 | 25 | 26 | 27 |  28 | 29 | 30\n name 31 Interesting_residues\n 15 & 8\n name 32 Transmembrane_sidechains\n 15 & 3 \n name 33 Transmembrane_c-alpha\n 13 | 32\n name 34 Transmembrane_sidechains_LIG\n q' | gmx make_ndx -f small.tpr -o index.ndx

#   0 System              :  7599 atoms
#   1 Protein             :  7537 atoms
#   2 Protein-H           :  3687 atoms
#   3 C-alpha             :   462 atoms
#   4 Backbone            :  1388 atoms
#   5 MainChain           :  1851 atoms
#   6 MainChain+Cb        :  2290 atoms
#   7 MainChain+H         :  2304 atoms
#   8 SideChain           :  5233 atoms
#   9 SideChain-H         :  1836 atoms
#  10 Prot-Masses         :  7537 atoms
#  11 non-Protein         :    62 atoms
#  12 Other               :    62 atoms
#  13 LIG                 :    62 atoms
#  14 Protein_LIG         :  7599 atoms
#  15 Transmembrane       :  4967 atoms
#  16 Lysozyme            :  2570 atoms
#  17 Interface_Nter      :    76 atoms
#  18 Interface           :   184 atoms
#  19 Interface_All       :   260 atoms
#  20 I138_3x46           :    19 atoms
#  21 V258_6x36           :    16 atoms
#  22 E294_7x35           :    15 atoms
#  23 Y311_7x53           :    21 atoms
#  24 Y295_7x36           :    21 atoms
#  25 N101_2x60           :    14 atoms
#  26 S105_EL2            :    11 atoms
#  27 E121_3x29           :    15 atoms
#  28 Y98_2x57            :    21 atoms
#  29 Appunti_ciaraia1    :    19 atoms
#  30 Appunti_ciaraia2    :     7 atoms
#  31 Interesting_residues:   179 atoms
#  32 Transmembrane_sidechains:  3460 atoms
#  33 Transmembrane_c-alpha:   302 atoms

# temperature 
echo "15 0" | gmx energy -f ener.edr  -o temperature.xvg -skip $SKIP_ENERGY

# min periodic distance 
echo "1" | gmx mindist -s small.tpr -f small.xtc -pi -od min_periodic_distance.xvg -dt $DT

# gyration
echo "1" | gmx gyrate -s small.tpr -f small.xtc -o gyrate.xvg -dt $DT

# rmsf
echo "3" | gmx rmsf -s small.tpr -f small.xtc  -o rmsf.xvg -res

# rmsd
echo "3 3" | gmx rms -s small.tpr -f small.xtc -o rmsd.xvg -dt $DT

# secondary structure
echo "1" |gmx do_dssp -s small.tpr -f small.xtc -o ss.xpm -sc scount.xvg -dt $DT 
# QUI BISOGNA PRODURRE L'EPS CORRETTO ++ OPPURE importare l'xpm da qualche parte

# sasa
echo "1" | gmx sasa -s small.tpr -f small.xtc -o sasa.xvg -dt $DT

# hbonds

echo "15 15" | gmx hbond -f short.xtc -dt 100 -s small.tpr -n index.ndx -num hbnum.xvg

echo "15 13" | gmx hbond -f short.xtc -dt 100 -s small.tpr -n index.ndx -num hbnum_lig.xvg

do_hbonds -f short.xtc -dt 100 -s small.tpr -o hbonds_tm -sel 15 -t 5 -n index.ndx -min-res-dist 5

# after running do_hbonds on both directory the following command should be issued once
# do_missing_dist -dirs /Users/uliano/Desktop/MDproc/S1PR1_Lys+FP /Users/uliano/Desktop/MDproc/S1PR1_Lys+S1P -hdir hbonds_tm

do_hbonds -f short.xtc -dt 100 -s small.tpr -o hbonds_lig -sel 15 -sel2 13 -t 5 -n index.ndx

# missing bonds should be computed by hand forming groups 
# (see S1PR1_Lys+S1P/saltbr_hbond_lig.ndx) and then running distance

gmx distance -s first.pdb -f short.xtc -n saltbr_hbond_lig.ndx -oav hbonds_lig/hbond-LIG1515-SER105.csv -select 'com of group 15 plus com of group 20'

# then they can be converted from .xvg to .csv using xvg2csv from within the bond directory

# also hbonds between ligand and protein could be splitted in 2 files (swapping donor and acceptor) and should be merged 

merge2csv -f1 hbond-SER105-LIG1515.csv -f2 hbond-LIG1515-SER105.csv
merge2csv -f1 hbond-ASN101-LIG1515.csv -f2 hbond-LIG1515-ASN101.csv
merge2csv -f1 hbond-ARG120-LIG1515.csv -f2 hbond-LIG1515-ARG120.csv

# salt bridges 

# for protein-protein bridges use the tool from VMD and save to salbr directory
# missing bridge should be processed by hand forming groups (see S1PR1_Lys+FP/saltbr.ndx)

gmx distance -s first.pdb -f short.xtc -n saltbr.ndx -oav saltbr/ASP18-ARG27.xvg -select 'com of group 14 plus com of group 15'
gmx distance -s first.pdb -f short.xtc -n saltbr.ndx -oav saltbr/ASP279-LYS200.xvg -select 'com of group 16 plus com of group 17'
gmx distance -s first.pdb -f short.xtc -n saltbr.ndx -oav saltbr/ASP288-ARG292.xvg -select 'com of group 18 plus com of group 19'
gmx distance -s first.pdb -f short.xtc -n saltbr.ndx -oav saltbr/ASP288-LYS285.xvg -select 'com of group 18 plus com of group 20'
gmx distance -s first.pdb -f short.xtc -n saltbr.ndx -oav saltbr/GLU294-LYS34.xvg -select 'com of group 21 plus com of group 22'
gmx distance -s first.pdb -f short.xtc -n saltbr.ndx -oav saltbr/GLU317-ARG247.xvg -select 'com of group 23 plus com of group 24'
gmx distance -s first.pdb -f short.xtc -n saltbr.ndx -oav saltbr/GLU317-ARG78.xvg -select 'com of group 23 plus com of group 25'

# for ligand-protein it is all done by hand forming groups (see S1PR1_Lys+S1P/saltbr_hbond_lig.ndx)

gmx distance -s first.pdb -f short.xtc -n saltbr_hbond_lig.ndx -oav saltbr_lig/saltbr-LIG1515-ARG120.xvg -select 'com of group 16 plus com of group 14'
gmx distance -s first.pdb -f short.xtc -n saltbr_hbond_lig.ndx -oav saltbr_lig/saltbr-LIG1515-LYS34.xvg -select 'com of group 17 plus com of group 14'
gmx distance -s first.pdb -f short.xtc -n saltbr_hbond_lig.ndx -oav saltbr_lig/saltbr-GLU121-LIG1515.xvg -select 'com of group 18 plus com of group 15'
gmx distance -s first.pdb -f short.xtc -n saltbr_hbond_lig.ndx -oav saltbr_lig/saltbr-GLU294-LIG1515.xvg -select 'com of group 19 plus com of group 15'


# 3x46 switch
echo "20 21" |gmx mindist -s small.tpr -f small.xtc  -n  -od mindist_i138_v258.xvg -on ncontact_i138_v258.xvg -group -dt $DT

echo "20 23" |gmx mindist -s small.tpr -f small.xtc  -n  -od mindist_i138_y311.xvg -on ncontact_i138_y311.xvg -group -dt $DT

# cluster

echo "15 15" | gmx cluster -f short.xtc -dt 100 -s small.tpr -n index.ndx -cutoff 0.2 -method gromos -o rmsd-clust-0.2.xpm -g cluster-0.2.log -dist rmsd-distribution-0.2.xvg -minstruct 20 -tr rmsd-clust-transitions-0.2.xpm -ntr rmsd-clust-transitions-0.2.xvg -clid rmsd-clust-id-over-time-0.2.xvg

# essential dynamics

echo "33 33" | gmx covar -f small.xtc -dt 10 -s small.tpr -n index.ndx -l -ascii -av -o -v

echo "33 33" |gmx anaeig -f small.xtc -s small.tpr -v eigenvec.trr -n index.ndx -eig eigenval.xvg  -extr extreme.pdb -first 1 -last 1 -nframes 50 

