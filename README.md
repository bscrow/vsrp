# vsrp

This is a bash pipeline that converts NetAffx annotation and sample files
(TSV format, from NA29 and/or NA32), into sorted, Gunzipped VCF files for
each chromosome - the required input file format for the Michigan
Imputation Server.

The proposed pipeline builds on top of the data processing pipeline as 
described on [Michigan Imputation Server Docs](https://imputationserver.readthedocs.io/en/latest/prepare-your-data/)

Steps:
------
1. Merge annotation files. If there is a conflict between the two files,
    the more up-to-date NA32 annotation will be kept
2. Convert sample and annotation files from TSV format to PLINK formats
3. Perform quality control* and convert to single chromosome VCF files using
    the pipeline recommended by Michigan Imputation Server.
4. Sort and gunzip VCF files using BCFTools

Requirements:
-------------
- Python 3.5 or newer
- plink and perl installed and executable added to path
- HRC-1000G-check-bim.pl (v4.3.0) & the tab delimited HRC reference,
    or wget and about 700Mb of disk space
- BCFTools installed and executable added to path

Usage:
------

```
python NetAffx2Imputation.py [-h] [--dir29 DIR29] [--anno29 ANNO29] [--hrc HRC] [-o O] [-r] dir32 anno32

This is a bash pipeline that converts NetAffx annotation and sample files (TSV format, from NA29 and/or NA32), into 
sorted, Gunzipped VCF files for each chromosome - the required input file format for the Michigan Imputation Server.

positional arguments:
  dir32            Directory containing sample TSV files from NA32
  anno32           TSV annotation file from NA32

optional arguments:
  -h, --help       show this help message and exit
  --dir29 DIR29    Directory containing sample TSV files from NA29
  --anno29 ANNO29  TSV annotation file for NA29
  --hrc HRC        Directory containing HRC-1000G-check-bim-v4.3.0.zip & the tab delimited HRC reference. If not 
                   included, the required files will be downloaded into the current directory in hrc/.
  -o O             Output directory of vcf.gz files
  -r               Remove temporary files
```

Sample Run:
-----------
`python NetAffx2Imputation.py <na32_sample_dir> <na32_annotation_tsv_file> 
--dir29 <na29_sample_dir> --anno29 <na32_annotation_tsv_file> -o study_vcf`

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
    removed (Can be disabled / Threshold can be adjusted)