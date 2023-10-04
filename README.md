# PlastomeBurstAndAlign
Extract and align coding regions, introns and intergenic spacers across a set of plastomes

#### Installation
##### MAFFT
```
# Debian:
apt install mafft
```
##### Other dependencies
```
pip install coloredlogs
```

#### Exemplary usage
```
MYSCRIPT=~/git/PlastomeBurstAndAlign/PlastomeRegionBurstAndAlign.py

folder_CDS=./output_CDS
mkdir -p $folder_CDS
mkdir -p $folder_CDS/01_unalign
mkdir -p $folder_CDS/02_aligned
mkdir -p $folder_CDS/02_aligned/fasta
mkdir -p $folder_CDS/02_aligned/nexus
python $MYSCRIPT -i . -o $folder_CDS -s cds 1>${folder_CDS}/${folder_CDS}.log 2>&1
mv $folder_CDS/*.unalign.fasta $folder_CDS/01_unalign
mv $folder_CDS/*.aligned.fasta $folder_CDS/02_aligned/fasta
mv $folder_CDS/*.aligned.nexus $folder_CDS/02_aligned/nexus

folder_INT=./output_INT
mkdir -p $folder_INT
mkdir -p $folder_INT/01_unalign
mkdir -p $folder_INT/02_aligned
mkdir -p $folder_INT/02_aligned/fasta
mkdir -p $folder_INT/02_aligned/nexus
python $MYSCRIPT -i . -o $folder_INT -s int 1>${folder_INT}/${folder_INT}.log 2>&1
mv $folder_INT/*.unalign.fasta $folder_INT/01_unalign
mv $folder_INT/*.aligned.fasta $folder_INT/02_aligned/fasta
mv $folder_INT/*.aligned.nexus $folder_INT/02_aligned/nexus

folder_IGS=./output_IGS
mkdir -p $folder_IGS
mkdir -p $folder_IGS/01_unalign
mkdir -p $folder_IGS/02_aligned
mkdir -p $folder_IGS/02_aligned/fasta
mkdir -p $folder_IGS/02_aligned/nexus
python $MYSCRIPT -i . -o $folder_IGS -s igs 1>${folder_IGS}/${folder_IGS}.log 2>&1
mv $folder_IGS/*.unalign.fasta $folder_IGS/01_unalign
mv $folder_IGS/*.aligned.fasta $folder_IGS/02_aligned/fasta
mv $folder_IGS/*.aligned.nexus $folder_IGS/02_aligned/nexus

```

#### Exemplary exclude lists

###### Exclude-list for genes (CDS)
```
--excllist ['rpl12']
```

###### Exclude-list for introns (INT)
```
--excllist ['rpl12']
```

###### Exclude-list for intergenic spacers (IGS)
```
--excllist ['rpl12_rpl2']
--excllist ['trnI_CAU_ycf2', 'rpl23_trnI_CAU', 'trnM_CAU_ycf2', 'rpl23_trnM_CAU', 'ndhI_ndhH', 'clpP_rpl20', 'rpl2_trnI_CAU', 'ycf4_cemA', 'ycf3_trnS_GGA', 'trnV_GAC_trnI_GAU', 'trnI_CAU_trnI_CAU', 'trnH_GUG_trnK_UUU', 'trnH_GUG_rpl2', 'psaI_ycf4', 'psaA_ycf3', 'trnM_CAU_rps14', 'trnl_CAU_ycf2', 'trnI_CAU_trnL_CAA', 'trnI_CAU_rpl12', 'trnG_UCC_trnM_CAU', 'trnG_UCC_atpA', 'rrn16_rps7', 'rps19_trnI_CAU', 'rps19_rpl23', 'rps15_ndhF', 'rpl23_trnl_CAU', 'rpl12_rpl2', 'psbZ_trnG_UCC', 'psbT_pbf1', 'psbH_rpoA', 'psbH_petB', 'petB_petD', 'pbf1_psbH', 'matK_rps16']
```

#### Setting up benchmark
```
apt install ncbi-entrez-direct  # on Debian
cd benchmark
while read accn; do
  efetch -db nucleotide -id $accn -format gb > ${accn}.gb;
done <Accession_numbers_Asteraceae_YangEtAl2022.txt  # doi: 10.3389/fpls.2022.808156
```
