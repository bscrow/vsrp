# vsrp

This is a bash pipeline that converts NetAffx annotation and sample files
(TSV format, from either NA29 or NA32), into sorted, Gunzipped VCF files 
for each chromosome - the required input file format for the Michigan
Imputation Server.

The proposed pipeline builds on top of the data processing pipeline as 
described on [Michigan Imputation Server Docs](https://imputationserver.readthedocs.io/en/latest/prepare-your-data/)

Steps taken by the merge command:
---------------------------------
1. Remove probesets with unknwon rsID, chromosome or position
2. Compare probesets from NA29 and NA32: 
    - Probesets with the same ID but different rsID are removed
    - For probesets with the same ID and rsID, chromosome and position from 
    NA32 are used. If NA29 and NA32 alleles differ, this program checks if 
    it is due to a strand reversal. If so, NA29 alleles are updated to 
    match NA32 alleles. If not, both sets of alleles are kept. 
3. Write the merged annotation into a TSV file with the following columns:
    1. Probe Set ID
    2. rsID
    3. Chromosome
    4. Chromosomal Position
    5. NA29 A Allele
    6. NA29 B Allele
    7. NA32 A Allele 
    8. NA32 B Allele

Steps taken by the run command:
-------------------------------
1. Convert sample and annotation files from TSV format to PLINK formats
2. Perform quality control* and convert to single chromosome VCF files using
    the pipeline recommended by Michigan Imputation Server.
3. Sort and gunzip VCF files using BCFTools

Requirements:
-------------
- Python 3.5 or newer
- plink and perl installed and executable added to path
- HRC-1000G-check-bim.pl (v4.3.0) & the tab delimited HRC reference,
    or wget and about 700Mb of disk space
- BCFTools installed and executable added to path

Usage:
------

#### Merge:
```
python NetAffx2Imputation.py merge [-h] [-o O] anno29 anno32

positional arguments:
  anno29      CSV annotation file from NA29
  anno32      CSV annotation file from NA32

optional arguments:
  -h, --help  show this help message and exit
  -o O        Name of the output merged annotation file. Default is 'merged_annotation.tsv'

```

#### Run: 
```
python NetAffx2Imputation.py run [-h] [--hrc HRC] [-t T] [-o O] [-r] anno samples group

positional arguments:
  anno        TSV annotation file
  samples     Directory containing sample TSV files
  group       File specifying the phenotype group of each sample

optional arguments:
  -h, --help  show this help message and exit
  --hrc HRC   Directory containing HRC-1000G-check-bim-v4.3.0.zip & the tab delimited HRC reference. If not included, the required files will be
              downloaded into the current directory in hrc/.
  -t T        Allele difference threshold value. Default is 0.2.
  -o O        Output directory of vcf.gz files
  -r          Remove temporary files
```

Sample Run:
-----------
`python NetAffx2Imputation.py merge Axiom_GW_Hu_SNP.r2.na29.annot.csv Axiom_GW_Hu_SNP.na32.annot.csv -o merged_anno.tsv`

`python NetAffx2Imputation.py run <annotation_tsv_file> <control_sample_dir> <case_annotation_tsv_file> 
--hrc ./hrc/`

Execution Time:
---------------

| No. of Samples | Test 1 /s| Test 2 /s| Test 3 /s|
| :-------------:|:--------:|:--------:|:--------:|
| 40             |4m25.254s |4m42.964s |4m29.313s |
| 100            |6m45.031s |6m53.969s |6m38.285s |
| 200            |12m5.040s |12m52.375s|12m39.313s|


Results are from the bash time command (real time)

*Quality Control checks:
------------------------
1. Checks that the number of variants in the PLINK files match; variants 
are coded as ACTG bases.

2. Match the position and rsID of SNPs in the PLINK files to the 
Haplotype Reference Consortium reference panel (v1.1 HRC.r1-1.GRCh37.wgs.mac5.sites.tab):
    - For a SNP, if its rsID exists in the reference but its position is 
    mismatched, the position of the SNP will be corrected if the rsID is 
    unique in the list.  SNPs with other mismatch cases are removed. 
    - Otherwise, if a SNP has a valid position in the reference but its 
    rsID cannot be found, its rsID will be updated based on the reference.
    - SNPs with other cases of mismatches/duplicated SNPs are removed. 

3. Based on the reference, SNPs identified as indels are excluded.

4. Match SNP alleles with that of the reference.
    - SNPs with swapped ref/alt alleles are corrected 
    - SNPs that are found to be mapped to the wrong strand are corrected
    - SNPs with differing alleles are removed
    - Palindromic SNPs with MAF > 0.4 are removed (Can be disabled)
    - SNPs exhibiting MAF that differs from the reference MAF by > 0.2 are 
    removed (Can be disabled / Threshold can be adjusted by specifying -t)