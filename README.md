# PlastomeBurstAndAlign
Extract and align coding regions, introns and intergenic spacers across a set of plastomes

#### Exemplary usage
```
MYSCRIPT=~/git/PlastomeBurstAndAlign/PlastomeRegionBurstAndAlign.py

folder_CDS=./output_CDS
mkdir -p $folder_CDS
python $MYSCRIPT -i . -o $folder_CDS -s cds

folder_INT=./output_INT
mkdir -p $folder_INT
python $MYSCRIPT -i . -o $folder_INT -s int

folder_IGS=./output_IGS
mkdir -p $folder_IGS
python $MYSCRIPT -i . -o $folder_IGS -s igs
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
